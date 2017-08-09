package org.broadinstitute.hellbender.tools.spark.bwa;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import scala.Tuple2;
import scala.Tuple3;

import java.io.Serializable;
import java.util.*;


/**
 * The NativeBwaSparkEngine provides a simple interface for transforming a
 * JavaRDD<Tuple3<String fastqFile1, String fastqFile2, Long start_idx> in
 * which the reads are unaligned, into a JavaRDD<GATKRead> of aligned reads,
 * and does so using a native binary.
 *
 * Use it like this:
 *     Make one, call the align method for each of your input RDDs in
 *     a pipeline that runs some action, close it when the content of result
 *     RDD is no longer needed.
 *
 * The reason that the pipeline must culminate in some action, is because Spark
 * implements this class as a lazy transform and nothing will happen otherwise.
 */
public final class NativeBwaSparkEngine implements Serializable, AutoCloseable {
    private static final long serialVersionUID = 0L;
    private final String indexFileName;
    private final String readGroupHeaderLine;
    private String headerString = null;
    private SAMFileHeader headerObject = null;

    public NativeBwaSparkEngine(final String indexFileName,
                                final String readGroupHeaderLine) {
        Utils.nonNull(indexFileName);
        Utils.nonNull(readGroupHeaderLine);
        this.indexFileName = indexFileName;
        this.readGroupHeaderLine = readGroupHeaderLine;

        // needed for generating header
        System.loadLibrary("gatkbwa");
        initNativeAlign();
        closeNativeAlign();
    }

    /* When this class initializes, it will generate header using native bwa.
     * The generated header text string will be stored in headerString field.
     * To actually run align, the text header string must be converted to 
     * SAMFileHeader externally and set with setHeaderObject.
     */

    public String getHeaderString() {
        return headerString;
    }

    public void setHeaderObject(SAMFileHeader header) {
        headerObject = header;
        headerString = null;
    }

    public JavaRDD<GATKRead> align(JavaRDD<Tuple3<String, String, Long> > fastqRecords) {
        Utils.nonNull(headerObject);
        return fastqRecords.flatMap((Tuple3<String, String, Long> e) -> {
            /* Only this part is executed on worker nodes;
             * repeated call of initNativeAlign() will be ignored. */
            System.loadLibrary("gatkbwa");
            initNativeAlign();
            List<SAMRecord> samRecordList =  doNativeAlign(e._1(), e._2(), e._3().longValue());
            return samRecordList.iterator();
        }).map(e -> new SAMRecordToGATKReadAdapter(e));
    }

    @Override
    public void close() {
        closeNativeAlign();
    }

    private native List<SAMRecord> doNativeAlign(String fastqFile1, String fastq2File, long start_idx);

    private native void initNativeAlign();
    private native void closeNativeAlign();
}
