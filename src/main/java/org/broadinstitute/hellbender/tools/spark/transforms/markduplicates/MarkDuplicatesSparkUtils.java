package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import com.google.common.collect.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.*;
import scala.Tuple2;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Utility classes and functions for Mark Duplicates.
 */
public class MarkDuplicatesSparkUtils {

    // Used to set an attribute on the GATKRead marking this read as an optical duplicate.
    public static final String OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME = "OD";

    /**
     * (0) filter: remove unpaired reads and reads with an unmapped mate.
     * (1) keyReadsByName: label each read with its read group and read name.
     * (2) GroupByKey: group together reads with the same group and name.
     * (3) keyPairedEndsWithAlignmentInfo:
     *   (a) Sort each group of reads (see GATKOrder below).
     *   (b) Pair consecutive reads into PairedEnds. In most cases there will only be two reads
     *       with the same name. TODO: explain why there might be more.
     *   (c) Label each read with alignment information: Library, reference index,
     *       stranded unclipped start and reverse strand.
     *   (d) Leftover reads are emitted, unmodified, as an unpaired end.
     * (4) GroupByKey: Group PairedEnds that share alignment information. These pairs
     *     are duplicates of each other.
     * (5) markDuplicatePairs:
     *   (a) For each group created by (4), sort the pairs by score and mark all but the
     *       highest scoring as duplicates.
     *   (b) Determine which duplicates are optical duplicates and increase the overall count.
     */

//convert String to hashcode
    static long myHashCodeForString (String a){
	if ( a == null )
                return 0;
	long result = 0;
	for (int i = 0; i < a.length(); i++){
		result = result * 31 *31 + a.charAt(i);
	}
	return result;
    }

// Convert byte array key into hashcode
    static long myHashCode (byte[] a){
	if ( a == null )
		return 0;
	
	long result = 0;
	for (byte e : a)result = 31 * 31 * result + e;
	return result;
    }

// Personal version of Readskey.keyForRead
    static long keyForRead2 (final SAMFileHeader header, GATKRead read){
	String rg = read.getReadGroup();
	String name = read.getName();
	return myHashCodeForString(rg) + myHashCodeForString(name) * 31 * 31;
    }

    static byte[] key2(PairedEnds pr, final SAMFileHeader header, short library) {
    GATKRead first = pr.first();
    GATKRead second = pr.second();
    int r1 = ReadUtils.getReferenceIndex(first,header);
    int s1 = ReadUtils.getStrandedUnclippedStart(first);
    int r2 = ReadUtils.getReferenceIndex(second,header);
    int s2 = ReadUtils.getStrandedUnclippedStart(second);
    if(first.isReverseStrand()){s1 = ~s1 + 1;}
    if(second.isReverseStrand()){s2 = ~s2 + 1;}
    return new byte[] {
        (byte)(library >>> 8),
        (byte)(library),
        (byte)(r1 >>> 24),
        (byte)(r1 >>> 16),
        (byte)(r1 >>> 8),
        (byte)(r1),
        (byte)(s1 >>> 24),
        (byte)(s1 >>> 16),
        (byte)(s1 >>> 8),
        (byte)(s1),
        (byte)(r2 >>> 24),
        (byte)(r2 >>> 16),
        (byte)(r2 >>> 8),
        (byte)(r2),
        (byte)(s2 >>> 24),
        (byte)(s2 >>> 16),
        (byte)(s2 >>> 8),
        (byte)(s2)};
  }
  static byte[] keyForFragment2(GATKRead first, final SAMFileHeader header, short library) {
    int r1 = ReadUtils.getReferenceIndex(first,header);
    int s1 = ReadUtils.getStrandedUnclippedStart(first);
    if(first.isReverseStrand()){s1 = ~s1 + 1;}
    return new byte[] {
        (byte)(library >>> 8),
        (byte)(library),
        (byte)(r1 >>> 24),
        (byte)(r1 >>> 16),
        (byte)(r1 >>> 8),
        (byte)(r1),
        (byte)(s1 >>> 24),
        (byte)(s1 >>> 16),
        (byte)(s1 >>> 8),
        (byte)(s1)};
  }
  static String getLibraryName(final SAMFileHeader header, String readGroupId) {
      if (readGroupId != null) {
            final SAMReadGroupRecord rg = header.getReadGroup(readGroupId);
            if (rg != null) {
                final String libraryName = rg.getLibrary();
                if (null != libraryName) return libraryName;
            }
        }

        return "Unknown Library";
    }

  static short getLibraryId(final SAMFileHeader header, Map<String, Short> libraryIds, 
				Short nextLibraryId, String readGroupId){
        final String library = getLibraryName(header, readGroupId);
        Short libraryId = libraryIds.get(library);

        if (libraryId == null) {
            libraryId = nextLibraryId++;
            libraryIds.put(library, libraryId);
        }

        return libraryId;
  }
	
    static JavaRDD<GATKRead> transformReads(final SAMFileHeader header, final MarkDuplicatesScoringStrategy scoringStrategy, final OpticalDuplicateFinder finder, final JavaRDD<GATKRead> reads, final int numReducers) {
	JavaPairRDD<Long, Iterable<GATKRead>> keyedReads;
        if (SAMFileHeader.SortOrder.queryname.equals(header.getSortOrder())) {
            // reads are already sorted by name, so perform grouping within the partition (no shuffle)
            keyedReads = spanReadsByKey(header, reads);
        } else {
            // sort by group and name (incurs a shuffle)
            JavaPairRDD<Long, GATKRead> keyReadPairs = reads.mapToPair(read -> new Tuple2<>(keyForRead2(header, read), read));
            keyedReads = keyReadPairs.groupByKey(numReducers);
        }

	Map<String, Short> libraryIds = new LinkedHashMap<>(); // from library string to library id
	short nextLibraryId = 1;	
	JavaPairRDD<Long, Iterable<Tuple2<Long, Integer>>> keyedValues = keyedReads.flatMapToPair(keyedRead -> {
            List<Tuple2<Long, Tuple2<Long, Integer>>> out = Lists.newArrayList();
            // Write each read out as a pair with only the first slot filled
            for (GATKRead read : keyedRead._2()) {
                read.setIsDuplicate(false);
		int score = scoringStrategy.score(read);
		if (!ReadUtils.readHasMappedMate(read)){score = ~score;} // Effective way to distinguish fragment. largest:-1
                out.add(new Tuple2<>((myHashCode(keyForFragment2(read, header, getLibraryId(header, libraryIds, nextLibraryId,
			read.getReadGroup()))) & 0x3fffffffffffffffL) | 0x8000000000000000L, // Fragment key top two bits 10
			new Tuple2<>(keyedRead._1(),score)));
            }
            // Write each paired read with a mapped mate as a pair
            final List<GATKRead> sorted = Lists.newArrayList(Iterables.filter(keyedRead._2(), read -> ReadUtils.readHasMappedMate(read)));
            sorted.sort(new GATKOrder(header));
            int score = -1;
	    PairedEnds pair = null;
            //Records are sorted, we iterate over them and pair them up.
            for (final GATKRead record : sorted) {
                if (score == -1) {                                //first in pair
		    pair = PairedEnds.of(record);
		    score = scoringStrategy.score(record);
                } else {                                           //second in pair
                    pair.and(record);
		    score += scoringStrategy.score(record);
                    out.add(new Tuple2<>((myHashCode(key2(pair, header, getLibraryId(header, libraryIds, nextLibraryId, 
					record.getReadGroup()))) & 0x3fffffffffffffffL), new Tuple2<>(keyedRead._1(),score))); // Pair key for pair top two bits 00
                    score = -1;                                   //back to first
                }
            }
            if (score != -1) {                                    //left over read
                out.add(new Tuple2<>((myHashCode(key2(pair, header, getLibraryId(header, libraryIds, nextLibraryId, 
				pair.first().getReadGroup()))) & 0x3fffffffffffffffL) | 0x4fffffffffffffffL ,
				 new Tuple2<>(keyedRead._1(),score))); // Pair key for single read top two bits 01
            }
            return out.iterator();
        }).groupByKey(numReducers);
	// If duplicate for each RG|name
	JavaPairRDD<Long, Boolean> keyedFlags = markPairedEnds2(keyedValues);
	
	return keyedReads.join(keyedFlags).flatMap(keyedReadFlag -> {
		List<GATKRead> out = Lists.newArrayList();
		Tuple2<Iterable<GATKRead>, Boolean> tmp = keyedReadFlag._2();
		if (!tmp._2){
			for (GATKRead read : tmp._1){
				read.setIsDuplicate(true);
			}
		}
		for (Iterator<GATKRead> it = tmp._1.iterator();it.hasNext();){
			out.add(it.next());
		}
		return out.iterator();
	});
    }
    
    static JavaPairRDD<Long, Boolean> markPairedEnds2(final JavaPairRDD<Long, Iterable<Tuple2<Long, Integer>>> keyedValues) {
        return keyedValues.flatMapToPair(keyedValue -> {
            Iterable<Tuple2<Long, Integer>> keyedScores = keyedValue._2();
            //final ImmutableListMultimap<Boolean, PairedEnds> paired = Multimaps.index(pairedEnds, pair -> pair.second() != null);

            // Each key corresponds to either fragments or paired ends, not a mixture of both.

            if (keyedValue._1() < 0) { // fragments
                return handleFragments2(keyedScores).iterator();
            }

            List<Tuple2<Long, Boolean>> out = Lists.newArrayList();

            // As in Picard, unpaired ends left alone.
            if((keyedValue._1() & 0x4000000000000000L) != 0){
		for (Iterator<Tuple2<Long, Integer>> it = keyedScores.iterator(); it.hasNext();){
			out.add(new Tuple2<Long, Boolean>(it.next()._1, false));
		}
		return out.iterator();
	    }
	    
	    int count=0;
	    int bestScore=-1;
	    int bestScoreId=0;
	    for (Iterator<Tuple2<Long, Integer>> it = keyedScores.iterator(); it.hasNext();){
		Tuple2<Long, Integer> tmp = it.next();
		if (tmp._2 >= bestScore) {
			bestScore = tmp._2;
			bestScoreId = count;
		}
		count++;
	    }
	    int index=0;
	    for (Iterator<Tuple2<Long, Integer>> it = keyedScores.iterator(); it.hasNext();) {
		boolean flag = false;
		if(index == bestScoreId) {flag = true;}
		out.add(new Tuple2<Long, Boolean>(it.next()._1, flag));
		index++;
	    }
            return out.iterator();
        });
    }
    
    private static List<Tuple2<Long, Boolean>> handleFragments2(Iterable<Tuple2<Long, Integer>> keyedScores) {
        
	boolean ifPairFlag = false;
	for(Iterator<Tuple2<Long, Integer>> it = keyedScores.iterator(); it.hasNext();){
		if(it.next()._2 >= 0){
			ifPairFlag = true;
			break;
		}
	}
	List<Tuple2<Long, Boolean>> flags = Lists.newArrayList();	
        // Note the we emit only fragments from this mapper.
        if (ifPairFlag == false) {
            // There are no paired reads, mark all but the highest scoring fragment as duplicate.
        	int id = 0;
		int bestScore = 0;
		int bestScorePlace = 0;
		for (Iterator<Tuple2<Long, Integer>> it = keyedScores.iterator(); it.hasNext();){
			Tuple2<Long, Integer> tmp = it.next();
			if (tmp._2 <= bestScore){
				bestScore = tmp._2;
				bestScorePlace = id;
			}
			id++;
		}
		id = 0;
		for (Iterator<Tuple2<Long, Integer>> it = keyedScores.iterator(); it.hasNext();){
			boolean flag = false;
			if (id == bestScorePlace) {flag = true;}
			flags.add(new Tuple2<Long, Boolean>(it.next()._1, flag));
			id++;
		}
	
	} else {
            // There are paired ends so we mark all fragments as duplicates.
        	for (Iterator<Tuple2<Long, Integer>> it = keyedScores.iterator(); it.hasNext();){
			Tuple2<Long, Integer> tmp = it.next();
			if(tmp._2 < 0){
				flags.add(new Tuple2<Long, Boolean>(tmp._1, false));
			}
		}
	}
        return flags;
    }
    static JavaPairRDD<Long, Iterable<GATKRead>> spanReadsByKey(final SAMFileHeader header, final JavaRDD<GATKRead> reads) {
        JavaPairRDD<String, GATKRead> nameReadPairs = reads.mapToPair(read -> new Tuple2<>(read.getName(), read));
        return spanByKey(nameReadPairs).flatMapToPair(namedRead -> {
            // for each name, separate reads by key (group name)
            List<Tuple2<Long, Iterable<GATKRead>>> out = Lists.newArrayList();
            ListMultimap<Long, GATKRead> multi = LinkedListMultimap.create();
            for (GATKRead read : namedRead._2()) {
                multi.put(keyForRead2(header, read), read);
            }
            for (Long key : multi.keySet()) {
                // list from Multimap is not serializable by Kryo, so put in a new array list
                out.add(new Tuple2<>(key, Lists.newArrayList(multi.get(key))));
            }
            return out.iterator();
        });
    }

    /**
     * Like <code>groupByKey</code>, but assumes that values are already sorted by key, so no shuffle is needed,
     * which is much faster.
     * @param rdd the input RDD
     * @param <K> type of keys
     * @param <V> type of values
     * @return an RDD where each the values for each key are grouped into an iterable collection
     */
    static <K, V> JavaPairRDD<K, Iterable<V>> spanByKey(JavaPairRDD<K, V> rdd) {
        return rdd.mapPartitionsToPair(iter -> spanningIterator(iter));
    }

    /**
     * An iterator that groups values having the same key into iterable collections.
     * @param iterator an iterator over key-value pairs
     * @param <K> type of keys
     * @param <V> type of values
     * @return an iterator over pairs of keys and grouped values
     */
    static <K, V> Iterator<Tuple2<K, Iterable<V>>> spanningIterator(Iterator<Tuple2<K, V>> iterator) {
        final PeekingIterator<Tuple2<K, V>> iter = Iterators.peekingIterator(iterator);
        return new AbstractIterator<Tuple2<K, Iterable<V>>>() {
            @Override
            protected Tuple2<K, Iterable<V>> computeNext() {
                K key = null;
                List<V> group = Lists.newArrayList();
                while (iter.hasNext()) {
                    if (key == null) {
                        Tuple2<K, V> next = iter.next();
                        key = next._1();
                        V value = next._2();
                        group.add(value);
                        continue;
                    }
                    K nextKey = iter.peek()._1(); // don't advance...
                    if (nextKey.equals(key)) {
                        group.add(iter.next()._2()); // .. unless the keys match
                    } else {
                        return new Tuple2<>(key, group);
                    }
                }
                if (key != null) {
                    return new Tuple2<>(key, group);
                }
                return endOfData();
            }
        };
    }
    private static int countOpticalDuplicates(OpticalDuplicateFinder finder, List<PairedEnds> scored) {
        final boolean[] opticalDuplicateFlags = finder.findOpticalDuplicates(scored);
        int numOpticalDuplicates = 0;
        for (final boolean b : opticalDuplicateFlags) {
            if (b) {
                numOpticalDuplicates++;
            }
        }
        return numOpticalDuplicates;
    }
    static JavaPairRDD<String, DuplicationMetrics> generateMetrics(final SAMFileHeader header, final JavaRDD<GATKRead> reads) {
        return reads.filter(read -> !read.isSecondaryAlignment() && !read.isSupplementaryAlignment())
                .mapToPair(read -> {
                    final String library = LibraryIdGenerator.getLibraryName(header, read.getReadGroup());
                    DuplicationMetrics metrics = new DuplicationMetrics();
                    metrics.LIBRARY = library;
                    if (read.isUnmapped()) {
                        ++metrics.UNMAPPED_READS;
                    } else if (!read.isPaired() || read.mateIsUnmapped()) {
                        ++metrics.UNPAIRED_READS_EXAMINED;
                    } else {
                        ++metrics.READ_PAIRS_EXAMINED;
                    }

                    if (read.isDuplicate()) {
                        if (!read.isPaired() || read.mateIsUnmapped()) {
                            ++metrics.UNPAIRED_READ_DUPLICATES;
                        } else {
                            ++metrics.READ_PAIR_DUPLICATES;
                        }
                    }
                    if (read.hasAttribute(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME)) {
                        metrics.READ_PAIR_OPTICAL_DUPLICATES +=
                                read.getAttributeAsInteger(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME);
                    }
                    return new Tuple2<>(library, metrics);
                })
                .foldByKey(new DuplicationMetrics(), (metricsSum, m) -> {
                    if (metricsSum.LIBRARY == null) {
                        metricsSum.LIBRARY = m.LIBRARY;
                    }
                    // This should never happen, as we grouped by key using library as the key.
                    if (!metricsSum.LIBRARY.equals(m.LIBRARY)) {
                        throw new GATKException("Two different libraries encountered while summing metrics: " + metricsSum.LIBRARY
                                + " and " + m.LIBRARY);
                    }
                    metricsSum.UNMAPPED_READS += m.UNMAPPED_READS;
                    metricsSum.UNPAIRED_READS_EXAMINED += m.UNPAIRED_READS_EXAMINED;
                    metricsSum.READ_PAIRS_EXAMINED += m.READ_PAIRS_EXAMINED;
                    metricsSum.UNPAIRED_READ_DUPLICATES += m.UNPAIRED_READ_DUPLICATES;
                    metricsSum.READ_PAIR_DUPLICATES += m.READ_PAIR_DUPLICATES;
                    metricsSum.READ_PAIR_OPTICAL_DUPLICATES += m.READ_PAIR_OPTICAL_DUPLICATES;
                    return metricsSum;
                })
                .mapValues(metrics -> {
                    DuplicationMetrics copy = metrics.copy();
                    // Divide these by 2 because they are counted for each read
                    // when they should be counted by pair.
                    copy.READ_PAIRS_EXAMINED = metrics.READ_PAIRS_EXAMINED / 2;
                    copy.READ_PAIR_DUPLICATES = metrics.READ_PAIR_DUPLICATES / 2;

                    copy.calculateDerivedMetrics();
                    if (copy.ESTIMATED_LIBRARY_SIZE == null) {
                        copy.ESTIMATED_LIBRARY_SIZE = 0L;
                    }
                    return copy;
                });
    }
    /**
     * Saves the metrics to a file.
     * Note: the SamFileHeader is needed in order to include libraries that didn't have any duplicates.
     * @param result metrics object, potentially pre-initialized with headers,
     */
    public static void saveMetricsRDD(final MetricsFile<DuplicationMetrics, Double> result, final SAMFileHeader header, final JavaPairRDD<String, DuplicationMetrics> metricsRDD, final String metricsOutputPath, AuthHolder authHolder) {
        final LibraryIdGenerator libraryIdGenerator = new LibraryIdGenerator(header);

        final Map<String, DuplicationMetrics> nonEmptyMetricsByLibrary = metricsRDD.collectAsMap();           //Unknown Library
        final Map<String, DuplicationMetrics> emptyMapByLibrary = libraryIdGenerator.getMetricsByLibraryMap();//with null

        final List<String> sortedListOfLibraryNames = new ArrayList<>(Sets.union(emptyMapByLibrary.keySet(), nonEmptyMetricsByLibrary.keySet()));
        sortedListOfLibraryNames.sort(Utils.COMPARE_STRINGS_NULLS_FIRST);
        for (final String library : sortedListOfLibraryNames){
            //if a non-empty exists, take it, otherwise take from the the empties. This is done to include libraries with zero data in them.
            //But not all libraries are listed in the header (esp in testing data) so we union empty and non-empty
            final DuplicationMetrics metricsToAdd = nonEmptyMetricsByLibrary.containsKey(library) ? nonEmptyMetricsByLibrary.get(library) : emptyMapByLibrary.get(library);
            metricsToAdd.calculateDerivedMetrics();
            result.addMetric(metricsToAdd);
        }

        if (nonEmptyMetricsByLibrary.size() == 1) {
            result.setHistogram(nonEmptyMetricsByLibrary.values().iterator().next().calculateRoiHistogram());
        }

        MetricsUtils.saveMetrics(result, metricsOutputPath);
    }
    /**
     * GATKRead comparator that compares based on mapping position followed by SAM flags.
     */
    final static class GATKOrder implements Comparator<GATKRead>, Serializable {
        private static final long serialVersionUID = 1l;
        private final SAMFileHeader header;
        // TODO: Unify with other comparators in the codebase

        public GATKOrder(final SAMFileHeader header) {
            this.header = header;
        }

        @Override
        public int compare(final GATKRead lhs, final GATKRead rhs) {
            if (rhs == lhs) return 0; //shortcut

            final int res1 = Integer.compare(ReadUtils.getReferenceIndex(lhs, header), ReadUtils.getReferenceIndex(rhs, header));
            if (res1 != 0) return res1;

            final int res2 = Long.compare(lhs.getStart(), rhs.getStart());
            if (res2 != 0) return res2;

            final int res3 = Boolean.compare(lhs.isDuplicate(), rhs.isDuplicate());
            if (res3 != 0) return res3;

            final int res4 = Boolean.compare(lhs.failsVendorQualityCheck(), rhs.failsVendorQualityCheck());
            if (res4 != 0) return res4;

            final int res5 = Boolean.compare(lhs.isPaired(), rhs.isPaired());
            if (res5 != 0) return res5;

            final int res6 = Boolean.compare(lhs.isProperlyPaired(), rhs.isProperlyPaired());
            if (res6 != 0) return res6;

            //Note: negate the result because we want first-of-pair to be before second
            //ie, want 'second' to be sorted after first, so want to return -1 for (true, false)
            final int res7 = -Boolean.compare(lhs.isFirstOfPair(), rhs.isFirstOfPair());
            if (res7 != 0) return res7;

            final int res8 = Boolean.compare(lhs.isSecondaryAlignment(), rhs.isSecondaryAlignment());
            if (res8 != 0) return res8;

            final int res9 = Boolean.compare(lhs.isSupplementaryAlignment(), rhs.isSupplementaryAlignment());
            if (res9 != 0) return res9;

            final int res10 = Integer.compare(lhs.getMappingQuality(), rhs.getMappingQuality());
            if (res10 != 0) return res10;

            final int res11 = Integer.compare(ReadUtils.getMateReferenceIndex(lhs, header), ReadUtils.getMateReferenceIndex(rhs, header));
            if (res11 != 0) return res11;

            final int res12 = Long.compare(lhs.getMateStart(), rhs.getMateStart());
            return res12;
        }
    }
}
