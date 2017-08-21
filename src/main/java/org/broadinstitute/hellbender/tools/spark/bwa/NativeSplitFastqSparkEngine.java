package org.broadinstitute.hellbender.tools.spark.bwa;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;

import scala.Tuple2;
import scala.Tuple3;

import java.io.Serializable;
import java.util.*;


public final class NativeSplitFastqSparkEngine implements Serializable {
    private static final long serialVersionUID = 0L;

    public static JavaRDD<Tuple3<String, String, Long> > split(JavaRDD<Tuple2<String[], long[]> > fastqSplitParams) {
        return fastqSplitParams.flatMap(e -> doSplit(e).iterator());
    }

    public static List<Tuple3<String, String, Long> > doSplit(Tuple2<String[], long[]> fastqSplitParamsTuple) {
        String inputFastq1 = fastqSplitParamsTuple._1()[0];
        String inputSplittedFastq1 = fastqSplitParamsTuple._1()[1];
        String inputFastq2 = fastqSplitParamsTuple._1()[2];
        String inputSplittedFastq2 = fastqSplitParamsTuple._1()[3];
        long fastqSplitSize = fastqSplitParamsTuple._2()[0];
        long fastqSplitReplication = fastqSplitParamsTuple._2()[1];
        long fastqSplitCompressionLevel = fastqSplitParamsTuple._2()[2];
        System.loadLibrary("splitfastq");
        int[] nativeSplitResult = doNativeSplit(
            inputFastq1, inputSplittedFastq1,
            inputFastq2, inputSplittedFastq2,
            fastqSplitSize, fastqSplitReplication, fastqSplitCompressionLevel
        );
        int numSplits = nativeSplitResult[0];
        int fastqSplitSizeInSeq = (int)(fastqSplitSize/nativeSplitResult[1]);
        List<Tuple3<String, String, Long> > fastqSplits = new ArrayList<Tuple3<String, String, Long> >(numSplits);
        for(int i = 0; i<numSplits; ++i) {
            fastqSplits.add(new Tuple3<String, String, Long>(
                inputSplittedFastq1+".part"+i,
                inputSplittedFastq2==null ? null : inputSplittedFastq2+".part"+i,
                new Long((long)i*(fastqSplitSizeInSeq))
            ));
        }
        return fastqSplits;
    }

    private static native int[] doNativeSplit(
        String inputFastq1, String inputSplittedFastq1,
        String inputFastq2, String inputSplittedFastq2,
        long fastqSplitSize, long fastqSplitReplication, long fastqSplitCompressionLevel
    );
}
