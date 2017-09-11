#include<string>

#include<jni.h>

#include"org_broadinstitute_hellbender_tools_spark_bwa_NativeSplitFastqSparkEngine.h"

using namespace std;

int SplitFASTQ(const int kVerboseFlag, const size_t kBatchSize, int* seq_length, const string& kInputFastq1, const string& kOutputFastq1, const string& kInputFastq2, const string& kOutputFastq2, const int kHDFSBufferSize = 0, const short kHDFSReplication = 0, const size_t kHDFSBlockSize = 0, const int8_t kCompressionLevel = 1);

/*
 * Class:     org_broadinstitute_hellbender_tools_spark_bwa_NativeSplitFastqSparkEngine
 * Method:    doNativeSplit
 * Signature: (Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;JJJ)[I
 */
JNIEXPORT jintArray JNICALL Java_org_broadinstitute_hellbender_tools_spark_bwa_NativeSplitFastqSparkEngine_doNativeSplit
  (JNIEnv* env, jclass obj, jstring input_fastq1, jstring input_splitted_fastq1, jstring input_fastq2, jstring input_splitted_fastq2, jlong fastq_split_size, jlong fastq_split_replication, jlong fastq_split_compression_level)
{
    const char* input_fastq1_chars = env->GetStringUTFChars(input_fastq1, nullptr);
    const char* input_splitted_fastq1_chars = env->GetStringUTFChars(input_splitted_fastq1, nullptr);
    const string kInputFastq1(input_fastq1_chars);
    const string kInputSplittedFastq1(input_splitted_fastq1_chars);
    env->ReleaseStringUTFChars(input_fastq1, input_fastq1_chars);
    env->ReleaseStringUTFChars(input_splitted_fastq1, input_splitted_fastq1_chars);

    const char* input_fastq2_chars = input_fastq2==nullptr ? nullptr : env->GetStringUTFChars(input_fastq2, nullptr);
    const char* input_splitted_fastq2_chars = input_fastq2==nullptr ? nullptr : env->GetStringUTFChars(input_splitted_fastq2, nullptr);
    const string kInputFastq2(input_fastq2!=nullptr ? input_fastq2_chars : "");
    const string kInputSplittedFastq2(input_fastq2!=nullptr ? input_splitted_fastq2_chars : "");
    if(input_fastq2!=nullptr)
    {
        env->ReleaseStringUTFChars(input_fastq2, input_fastq2_chars);
        env->ReleaseStringUTFChars(input_splitted_fastq2, input_splitted_fastq2_chars);
    }

    const int kHDFSReplication = int(fastq_split_replication);
    const size_t kBatchSize = size_t(fastq_split_size);
    const int8_t kCompressionLevel = int8_t(fastq_split_compression_level);

    jintArray result = env->NewIntArray(2);
    int result_array[2];
    result_array[0] = SplitFASTQ(true, kBatchSize, result_array+1, kInputFastq1, kInputSplittedFastq1, kInputFastq2, kInputSplittedFastq2, 0, kHDFSReplication, 0, kCompressionLevel);
    env->SetIntArrayRegion(result, 0, 2, result_array);
    return result;
}

