package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.StringLineReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.StorageLevels;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkPipelineProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.AddContextDataToReadSpark;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.JoinStrategy;
import org.broadinstitute.hellbender.engine.spark.SparkReadShard;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSink;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.ApplyBQSRUniqueArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.transforms.ApplyBQSRSparkFn;
import org.broadinstitute.hellbender.tools.spark.transforms.BaseRecalibratorSparkFn;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.broadinstitute.hellbender.tools.spark.bwa.NativeBwaSparkEngine;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import scala.Tuple2;

import java.io.IOException;
import java.io.Serializable;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;


import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;

@CommandLineProgramProperties(
        summary = "Takes unaligned fastq reads and runs Burrows-Wheeler Aligner, Mark Duplicates, Base Quality Score Recalibration, and Haplotype Caller. The final result is a VCF file.",
        oneLineSummary = "Takes unaligned fastq reads and runs BWA, MD, BQSR, and HC. The final result is a VCF file.",
        usageExample = "ReadsPipelineSpark --inputFastq1 input.fastq.gz [--inputFastq2 input2.fastq.gz] --bwamemIndexImage bwaIndexImage -R referenceURL --knownSites variants.vcf -O file:///tmp/output.vcf",
        programGroup = SparkPipelineProgramGroup.class
)

/**
 * ReadsPipelineSpark is our standard pipeline that takes unaligned fastq reads
 * and runs BWA, MarkDuplicates, BQSR, and HaplotypeCaller.
 */
@DocumentedFeature
@BetaFeature
public class ReadsPipelineSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReference() { return true; }

    @Argument(doc = "the known variants", shortName = "knownSites", fullName = "knownSites", optional = false)
    protected List<String> baseRecalibrationKnownVariants;

    @Argument(doc = "input [gzipped] fastq file prefix", fullName = "inputFastq1")
    protected String inputFastq1 = null;

    @Argument(doc = "another input [gzipped] fastq file prefix", fullName = "inputFastq2", optional = true)
    protected String inputFastq2 = null;

    @Argument(doc = "dump bwa results in BAM (currently broken)", fullName = "bwaResultBAM", optional = true)
    protected String bwaResultBAM = null;

    @Argument(doc = "dump bwa results in SAM", fullName = "bwaResultSAM", optional = true)
    protected String bwaResultSAM = null;

    @Argument(doc = "Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’.", fullName = "readGroupHeaderLine", optional = true)
    protected String readGroupHeaderLine = null;

    @Override
    public SAMFileHeader getHeaderForReads() {
        if(readsHeader == null) {
            if(readGroupHeaderLine!=null) {
                readsHeader = new SAMTextHeaderCodec().decode(new StringLineReader(readGroupHeaderLine), null);
            }
        }
        return readsHeader;
    }

    public void setHeaderForReads(String headerString) {
        readsHeader = new SAMTextHeaderCodec().decode(new StringLineReader(headerString), null);
    }

    protected SAMFileHeader readsHeader = null;

    @Argument(doc = "Single file to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    protected String output;

    @Argument(doc = "the bwa mem index image file name that you've distributed to each executor",
              fullName = "bwamemIndexImage")
    private String indexImageFile;

    @ArgumentCollection
    public final ShardingArgumentCollection shardingArgs = new ShardingArgumentCollection();

		public static class ShardingArgumentCollection implements Serializable {
        private static final long serialVersionUID = 1L;

        @Argument(fullName="readShardSizeForHaplotypeCaller", shortName="readShardSizeForHaplotypeCaller", doc = "Maximum size of each read shard, in bases. For good performance, this should be much larger than the maximum assembly region size.", optional = true)
        public int readShardSize = HaplotypeCaller.DEFAULT_READSHARD_SIZE;

        @Argument(fullName="readShardPaddingForHaplotypeCaller", shortName="readShardPaddingForHaplotypeCaller", doc = "Each read shard has this many bases of extra context on each side. Read shards must have as much or more padding than assembly regions.", optional = true)
        public int readShardPadding = HaplotypeCaller.DEFAULT_READSHARD_PADDING;

        @Argument(fullName = "minAssemblyRegionSize", shortName = "minAssemblyRegionSize", doc = "Minimum size of an assembly region", optional = true)
        public int minAssemblyRegionSize = HaplotypeCaller.DEFAULT_MIN_ASSEMBLY_REGION_SIZE;

        @Argument(fullName = "maxAssemblyRegionSize", shortName = "maxAssemblyRegionSize", doc = "Maximum size of an assembly region", optional = true)
        public int maxAssemblyRegionSize = HaplotypeCaller.DEFAULT_MAX_ASSEMBLY_REGION_SIZE;

        @Argument(fullName = "assemblyRegionPadding", shortName = "assemblyRegionPadding", doc = "Number of additional bases of context to include around each assembly region", optional = true)
        public int  assemblyRegionPadding = HaplotypeCaller.DEFAULT_ASSEMBLY_REGION_PADDING;

        @Advanced
        @Argument(fullName = "activeProbabilityThreshold", shortName = "activeProbabilityThreshold", doc="Minimum probability for a locus to be considered active.", optional = true)
        public double activeProbThreshold = HaplotypeCaller.DEFAULT_ACTIVE_PROB_THRESHOLD;

        @Advanced
        @Argument(fullName = "maxProbPropagationDistance", shortName = "maxProbPropagationDistance", doc="Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions", optional = true)
        public int maxProbPropagationDistance = HaplotypeCaller.DEFAULT_MAX_PROB_PROPAGATION_DISTANCE;

    }

		@ArgumentCollection
    public HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

    @Argument(doc = "the join strategy for reference bases and known variants", shortName = "joinStrategy", fullName = "joinStrategy", optional = true)
    private JoinStrategy joinStrategy = JoinStrategy.BROADCAST;

    @Argument(shortName = "DS", fullName ="duplicates_scoring_strategy", doc = "The scoring strategy for choosing the non-duplicate among candidates.")
    public MarkDuplicatesScoringStrategy duplicatesScoringStrategy = MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES;

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private final RecalibrationArgumentCollection bqsrArgs = new RecalibrationArgumentCollection();

    @Argument(fullName="readShardSize", shortName="readShardSize", doc = "Maximum size of each read shard, in bases. Only applies when using the OVERLAPS_PARTITIONER join strategy.", optional = true)
    public int readShardSize = 10000;

    @Argument(fullName="readShardPadding", shortName="readShardPadding", doc = "Each read shard has this many bases of extra context on each side. Only applies when using the OVERLAPS_PARTITIONER join strategy.", optional = true)
    public int readShardPadding = 1000;

    /**
     * command-line arguments to fine tune the apply BQSR step.
     */
    @ArgumentCollection
    public ApplyBQSRUniqueArgumentCollection applyBqsrArgs = new ApplyBQSRUniqueArgumentCollection();

    @Override
    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return BaseRecalibrationEngine.BQSR_REFERENCE_WINDOW_FUNCTION;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        // read fastq inputs
        List<Tuple2<String, String> > fastqRecordList = new ArrayList<Tuple2<String, String> >();
        Set<Integer> fileSet = new HashSet<Integer>();
        try {
            FileSystem fs = FileSystem.get(new URI(inputFastq1), new Configuration(true));
            FileStatus[] inputFastq1Files = fs.globStatus(new Path(new URI(inputFastq1+".part*")));
            if(inputFastq1Files.length==0) {
                // TODO: try to do the split
                throw new UserException("Cannot find file: "+inputFastq1+".part*");
            }
            for(FileStatus file:inputFastq1Files)
            {
                String pathString = file.getPath().toString();
                fileSet.add(Integer.parseInt(pathString.substring(pathString.lastIndexOf(".part")+5)));
            }
        } catch (IOException e) {
            throw new UserException("Cannot open file: "+e.getMessage());
        } catch (URISyntaxException e) {
            throw new UserException("URI syntax error: "+e.getMessage());
        }

        // stop if there is a gap
        for(int i = 0;; ++i)
        {
            if(fileSet.contains(new Integer(i)))
            {
                fastqRecordList.add(
                    new Tuple2<>(inputFastq1+".part"+i, inputFastq2!=null ? inputFastq2+".part"+i : null)
                );
            } else {
                numReducers = i;
                break;
            }
        }

        JavaRDD<Tuple2<String, String> > fastqRecords = ctx.parallelize(fastqRecordList, fastqRecordList.size());

        // bwa
        final NativeBwaSparkEngine bwaEngine = new NativeBwaSparkEngine(indexImageFile, readGroupHeaderLine);
        setHeaderForReads(bwaEngine.getHeaderString());
        bwaEngine.setHeaderObject(getHeaderForReads());

        final WellformedReadFilter wellformedReadFilter = new WellformedReadFilter(getHeaderForReads());
        final JavaRDD<GATKRead> rawReads = bwaEngine.align(fastqRecords);

        rawReads.persist(org.apache.spark.api.java.StorageLevels.MEMORY_AND_DISK);

        final JavaRDD<GATKRead> initialReads = rawReads.filter(read -> wellformedReadFilter.test(read));
        initialReads.persist(org.apache.spark.api.java.StorageLevels.MEMORY_AND_DISK);

        if(bwaResultBAM != null || bwaResultSAM != null)
        {
            System.err.println("Get "+rawReads.count()+" reads");
            System.err.println("Get "+initialReads.count()+" wellformed reads");

            if(bwaResultBAM != null) {
                String outputFile = bwaResultBAM;
                try {
                    ReadsSparkSink.writeReads(ctx, outputFile,
                            hasReference() ? referenceArguments.getReferenceFile().getAbsolutePath() : null,
                            rawReads, getHeaderForReads(), shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE,
                            getRecommendedNumReducers());
                } catch (IOException e) {
                    throw new UserException.CouldNotCreateOutputFile(outputFile,"writing failed", e);
                }
            }

            if(bwaResultSAM != null) {
                String outputFile = bwaResultSAM;
                BufferedWriter out = null;
                try {
                    out = new BufferedWriter(new FileWriter(bwaResultSAM, false));
                    for(GATKRead read: rawReads.collect()) {
                        out.write(read.getSAMString());
                    }
                    out.close();
                } catch (IOException e)
                {
                    System.err.println("Error: " + e.getMessage());
                }
            }
        }

        rawReads.unpersist();

        if (joinStrategy == JoinStrategy.BROADCAST && ! getReference().isCompatibleWithSparkBroadcast()){
            throw new UserException.Require2BitReferenceForBroadcast();
        }

        // mark duplicates
        final JavaRDD<GATKRead> markedReadsWithOD = MarkDuplicatesSpark.mark(initialReads, getHeaderForReads(), duplicatesScoringStrategy, new OpticalDuplicateFinder(), getRecommendedNumReducers());
        final JavaRDD<GATKRead> markedReads = MarkDuplicatesSpark.cleanupTemporaryAttributes(markedReadsWithOD);

        // The markedReads have already had the WellformedReadFilter applied to them, which
        // is all the filtering that MarkDupes and ApplyBQSR want. BQSR itself wants additional
        // filtering performed, so we do that here.
        //NOTE: this doesn't honor enabled/disabled commandline filters
        final ReadFilter bqsrReadFilter = ReadFilter.fromList(BaseRecalibrator.getBQSRSpecificReadFilterList(), getHeaderForReads());

        JavaRDD<GATKRead> markedFilteredReadsForBQSR = markedReads.filter(read -> bqsrReadFilter.test(read));

        if (joinStrategy.equals(JoinStrategy.OVERLAPS_PARTITIONER)) {
            // the overlaps partitioner requires that reads are coordinate-sorted
            final SAMFileHeader readsHeader = getHeaderForReads().clone();
            readsHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            markedFilteredReadsForBQSR = SparkUtils.coordinateSortReads(markedFilteredReadsForBQSR, readsHeader, numReducers);
        }

        final VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        final JavaRDD<GATKVariant> bqsrKnownVariants = variantsSparkSource.getParallelVariants(baseRecalibrationKnownVariants, getIntervals());

        final JavaPairRDD<GATKRead, ReadContextData> rddReadContext = AddContextDataToReadSpark.add(ctx, markedFilteredReadsForBQSR, getReference(), bqsrKnownVariants, baseRecalibrationKnownVariants, joinStrategy, getHeaderForReads().getSequenceDictionary(), readShardSize, readShardPadding);
        final RecalibrationReport bqsrReport = BaseRecalibratorSparkFn.apply(rddReadContext, getHeaderForReads(), getReferenceSequenceDictionary(), bqsrArgs);

        final Broadcast<RecalibrationReport> reportBroadcast = ctx.broadcast(bqsrReport);
        final JavaRDD<GATKRead> calibratedReads = ApplyBQSRSparkFn.apply(markedReads, reportBroadcast, getHeaderForReads(), applyBqsrArgs.toApplyBQSRArgumentCollection(bqsrArgs.PRESERVE_QSCORES_LESS_THAN));

        final List<SimpleInterval> intervals = hasIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(getHeaderForReads().getSequenceDictionary());
        final ReadFilter haplotypeCallerFilter = (new GATKReadFilterPluginDescriptor(HaplotypeCallerEngine.makeStandardHCReadFilters())).getMergedReadFilter(getHeaderForReads());
        final JavaRDD<GATKRead> filteredReads = calibratedReads.filter(read -> haplotypeCallerFilter.test(read));

        // hyplotype caller
        final JavaRDD<VariantContext> variants = callVariantsWithHaplotypeCaller(getAuthHolder(), ctx, filteredReads, getHeaderForReads(), getReference(), intervals, hcArgs, shardingArgs);
        if (hcArgs.emitReferenceConfidence == ReferenceConfidenceMode.GVCF) {
            // VariantsSparkSink/Hadoop-BAM VCFOutputFormat do not support writing GVCF, see https://github.com/broadinstitute/gatk/issues/2738
            writeVariants(variants);
        } else {
            final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgs, getHeaderForReads(), new ReferenceMultiSourceAdapter(getReference(), getAuthHolder()));
            variants.cache(); // without caching, computations are run twice as a side effect of finding partition boundaries for sorting
            try {
                VariantsSparkSink.writeVariants(ctx, output, variants, hcEngine.makeVCFHeader(getHeaderForReads().getSequenceDictionary(), new HashSet<>()));
            } catch (IOException e) {
                throw new UserException.CouldNotCreateOutputFile(output, "writing failed", e);
            }
        }
    }

    @Override
    public int getRecommendedNumReducers() {
      if (numReducers != 0) {
          return numReducers;
      }
      return 32;  // TODO: make it work
    }

		/**
     * Call Variants using HaplotypeCaller on Spark and return an RDD of  {@link VariantContext}
     *
     * This may be called from any spark pipeline in order to call variants from an RDD of GATKRead
     *
     * @param authHolder authorization needed for the reading the reference
     * @param ctx the spark context
     * @param reads the reads variants should be called from
     * @param header the header that goes with the reads
     * @param reference the reference to use when calling
     * @param intervals the intervals to restrict calling to
     * @param hcArgs haplotype caller arguments
     * @param shardingArgs arguments to control how the assembly regions are sharded
     * @return an RDD of Variants
     */
    public static JavaRDD<VariantContext> callVariantsWithHaplotypeCaller(
            final AuthHolder authHolder,
            final JavaSparkContext ctx,
            final JavaRDD<GATKRead> reads,
            final SAMFileHeader header,
            final ReferenceMultiSource reference,
            final List<SimpleInterval> intervals,
            final HaplotypeCallerArgumentCollection hcArgs,
            final ShardingArgumentCollection shardingArgs) {
        Utils.validateArg(hcArgs.dbsnp.dbsnp == null, "HaplotypeCallerSpark does not yet support -D or --dbsnp arguments" );
        Utils.validateArg(hcArgs.comps.isEmpty(), "HaplotypeCallerSpark does not yet support -comp or --comp arguments" );
        Utils.validateArg(hcArgs.bamOutputPath == null, "HaplotypeCallerSpark does not yet support -bamout or --bamOutput");
        if ( !reference.isCompatibleWithSparkBroadcast()){
            throw new UserException.Require2BitReferenceForBroadcast();
        }

        final Broadcast<ReferenceMultiSource> referenceBroadcast = ctx.broadcast(reference);
        final Broadcast<HaplotypeCallerArgumentCollection> hcArgsBroadcast = ctx.broadcast(hcArgs);
        final OverlapDetector<ShardBoundary> overlaps = getShardBoundaryOverlapDetector(header, intervals, shardingArgs.readShardSize, shardingArgs.readShardPadding);
        final Broadcast<OverlapDetector<ShardBoundary>> shardBoundariesBroadcast = ctx.broadcast(overlaps);

        final JavaRDD<Shard<GATKRead>> readShards = createReadShards(shardBoundariesBroadcast, reads);

        final JavaRDD<Tuple2<AssemblyRegion, SimpleInterval>> assemblyRegions = readShards
                .mapPartitions(shardsToAssemblyRegions(authHolder, referenceBroadcast, hcArgsBroadcast, shardingArgs, header));

        return assemblyRegions.mapPartitions(callVariantsFromAssemblyRegions(authHolder, header, referenceBroadcast, hcArgsBroadcast));
    }

    /**
     * Call variants from Tuples of AssemblyRegion and Simple Interval
     * The interval should be the non-padded shard boundary for the shard that the corresponding AssemblyRegion was
     * created in, it's used to eliminate redundant variant calls at the edge of shard boundaries.
     */
    private static FlatMapFunction<Iterator<Tuple2<AssemblyRegion, SimpleInterval>>, VariantContext> callVariantsFromAssemblyRegions(
            final AuthHolder authHolder,
            final SAMFileHeader header,
            final Broadcast<ReferenceMultiSource> referenceBroadcast,
            final Broadcast<HaplotypeCallerArgumentCollection> hcArgsBroadcast) {
        return regionAndIntervals -> {
            //HaplotypeCallerEngine isn't serializable but is expensive to instantiate, so construct and reuse one for every partition
            final ReferenceMultiSourceAdapter referenceReader = new ReferenceMultiSourceAdapter(referenceBroadcast.getValue(), authHolder);
            final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgsBroadcast.value(), header, referenceReader);
            return iteratorToStream(regionAndIntervals).flatMap(regionToVariants(hcEngine)).iterator();
        };
    }

    private static <T> Stream<T> iteratorToStream(Iterator<T> iterator) {
        Iterable<T> regionsIterable = () -> iterator;
        return StreamSupport.stream(regionsIterable.spliterator(), false);
    }

    private static Function<Tuple2<AssemblyRegion, SimpleInterval>, Stream<? extends VariantContext>> regionToVariants(HaplotypeCallerEngine hcEngine) {
        return regionAndInterval -> {
            final List<VariantContext> variantContexts = hcEngine.callRegion(regionAndInterval._1(), new FeatureContext());
            final SimpleInterval shardBoundary = regionAndInterval._2();
            return variantContexts.stream()
                .filter(vc -> shardBoundary.contains(new SimpleInterval(vc.getContig(), vc.getStart(), vc.getStart())));
        };
    }

    /**
     * WriteVariants, this is currently going to be horribly slow and explosive on a full size file since it performs a collect.
     *
     * This will be replaced by a parallel writer similar to what's done with {@link org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink}
     */
    private void writeVariants(JavaRDD<VariantContext> variants) {
        final List<VariantContext> collectedVariants = variants.collect();
        final SAMSequenceDictionary referenceDictionary = getReferenceSequenceDictionary();

        final List<VariantContext> sortedVariants = collectedVariants.stream()
            .sorted((o1, o2) -> IntervalUtils.compareLocatables(o1, o2, referenceDictionary))
            .collect(Collectors.toList());

        final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgs, getHeaderForReads(), new ReferenceMultiSourceAdapter(getReference(), getAuthHolder()));
        try(final VariantContextWriter writer = hcEngine.makeVCFWriter(output, getBestAvailableSequenceDictionary())) {
            hcEngine.writeHeader(writer, getHeaderForReads().getSequenceDictionary(), new HashSet<>());
            sortedVariants.forEach(writer::add);
        }
    }

    /**
     * Create an RDD of {@link Shard} from an RDD of {@link GATKRead}
     * @param shardBoundariesBroadcast  broadcast of an {@link OverlapDetector} loaded with the intervals that should be used for creating ReadShards
     * @param reads Rdd of {@link GATKRead}
     * @return a Rdd of reads grouped into potentially overlapping shards
     */
    private static JavaRDD<Shard<GATKRead>> createReadShards(final Broadcast<OverlapDetector<ShardBoundary>> shardBoundariesBroadcast, final JavaRDD<GATKRead> reads) {
        final JavaPairRDD<ShardBoundary, GATKRead> paired = reads.flatMapToPair(read -> {
            final Collection<ShardBoundary> overlappingShards = shardBoundariesBroadcast.value().getOverlaps(read);
            return overlappingShards.stream().map(key -> new Tuple2<>(key, read)).iterator();
        });
        final JavaPairRDD<ShardBoundary, Iterable<GATKRead>> shardsWithReads = paired.groupByKey();
        return shardsWithReads.map(shard -> new SparkReadShard(shard._1(), shard._2()));
    }

    /**
     * @return an {@link OverlapDetector} loaded with {@link ShardBoundary}
     * based on the -L intervals
     */
    private static OverlapDetector<ShardBoundary> getShardBoundaryOverlapDetector(final SAMFileHeader header, final List<SimpleInterval> intervals, final int readShardSize, final int readShardPadding) {
        final OverlapDetector<ShardBoundary> shardBoundaryOverlapDetector = new OverlapDetector<>(0, 0);
        intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, readShardSize, readShardPadding, header.getSequenceDictionary()).stream())
                .forEach(boundary -> shardBoundaryOverlapDetector.addLhs(boundary, boundary.getPaddedInterval()));
        return shardBoundaryOverlapDetector;
    }

    /**
     * @return and RDD of {@link Tuple2<AssemblyRegion, SimpleInterval>} which pairs each AssemblyRegion with the
     * interval it was generated in
     */
    private static FlatMapFunction<Iterator<Shard<GATKRead>>, Tuple2<AssemblyRegion, SimpleInterval>> shardsToAssemblyRegions(
            final AuthHolder authHolder,
            final Broadcast<ReferenceMultiSource> reference,
            final Broadcast<HaplotypeCallerArgumentCollection> hcArgsBroadcast,
            final ShardingArgumentCollection assemblyArgs,
            final SAMFileHeader header) {
        return shards -> {
            final ReferenceMultiSource referenceMultiSource = reference.value();
            final ReferenceMultiSourceAdapter referenceSource = new ReferenceMultiSourceAdapter(referenceMultiSource, authHolder);
            final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgsBroadcast.value(), header, referenceSource);

            return iteratorToStream(shards).flatMap(shardToRegion(assemblyArgs, header, referenceSource, hcEngine)).iterator();
        };
    }

    private static Function<Shard<GATKRead>, Stream<? extends Tuple2<AssemblyRegion, SimpleInterval>>> shardToRegion(
            ShardingArgumentCollection assemblyArgs,
            SAMFileHeader header,
            ReferenceMultiSourceAdapter referenceSource,
            HaplotypeCallerEngine evaluator) {
        return shard -> {
            final ReferenceContext refContext = new ReferenceContext(referenceSource, shard.getPaddedInterval());

            //TODO load features as a side input
            final FeatureContext features = new FeatureContext();

            final Iterable<AssemblyRegion> assemblyRegions = AssemblyRegion.createFromReadShard(
                    shard, header, refContext, features, evaluator,
                    assemblyArgs.minAssemblyRegionSize, assemblyArgs.maxAssemblyRegionSize,
                    assemblyArgs.assemblyRegionPadding, assemblyArgs.activeProbThreshold,
                    assemblyArgs.maxProbPropagationDistance);

            return StreamSupport.stream(assemblyRegions.spliterator(), false)
                    .map(a -> new Tuple2<>(a, shard.getInterval()));
        };
    }

    /**
     * Adapter to allow a 2bit reference to be used in HaplotypeCallerEngine.
     * This is not intended as a general purpose adapter, it only enables the operations needed in {@link HaplotypeCallerEngine}
     * This should not be used outside of this class except for testing purposes.
     */
    @VisibleForTesting
    public static final class ReferenceMultiSourceAdapter implements ReferenceSequenceFile, ReferenceDataSource, Serializable{
        private static final long serialVersionUID = 1L;

        private final ReferenceMultiSource source;
        private final AuthHolder auth;
        private final SAMSequenceDictionary sequenceDictionary;

        public ReferenceMultiSourceAdapter(final ReferenceMultiSource source, final AuthHolder auth) {
            this.source = source;
            this.auth = auth;
            sequenceDictionary = source.getReferenceSequenceDictionary(null);
        }

        @Override
        public ReferenceSequence queryAndPrefetch(final String contig, final long start, final long stop) {
           return getSubsequenceAt(contig, start, stop);
        }

        @Override
        public SAMSequenceDictionary getSequenceDictionary() {
            return source.getReferenceSequenceDictionary(null);
        }

        @Override
        public ReferenceSequence nextSequence() {
            throw new UnsupportedOperationException("nextSequence is not implemented");
        }

        @Override
        public void reset() {
            throw new UnsupportedOperationException("reset is not implemented");
        }

        @Override
        public boolean isIndexed() {
            return true;
        }

        @Override
        public ReferenceSequence getSequence(final String contig) {
            throw new UnsupportedOperationException("getSequence is not supported");
        }

        @Override
        public ReferenceSequence getSubsequenceAt(final String contig, final long start, final long stop) {
            try {
                final ReferenceBases bases = source.getReferenceBases(auth.asPipelineOptionsDeprecated(), new SimpleInterval(contig, (int) start, (int) stop));
                return new ReferenceSequence(contig, sequenceDictionary.getSequenceIndex(contig), bases.getBases());
            } catch (final IOException e) {
                throw new GATKException(String.format("Failed to load reference bases for %s:%d-%d", contig, start, stop));
            }
        }

        @Override
        public void close() {
            // doesn't do anything because you can't close a two-bit file
        }

        @Override
        public Iterator<Byte> iterator() {
            throw new UnsupportedOperationException("iterator is not supported");
        }
    }
}
