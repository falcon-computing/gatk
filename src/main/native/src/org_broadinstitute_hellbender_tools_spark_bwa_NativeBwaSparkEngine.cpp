#include<string.h>

#include<iostream>
#include<mutex>
#include<string>
#include<thread>
#include<vector>

#include<bwa.h>
#include<hdfs.h>
#include<kseq.h>
#include<zlib.h>

#include"bwa_wrapper.h"
#include"org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine.h"

#ifndef  __cplusplus
#error "This file requires a C++ compiler."
#endif //__cplusplus

using namespace std;

/* These global variables are shared on the same node and must be initialized
 * before calling doNativeAlign by calling initNativeAlign. When no longer
 * needed, call closeNativeAlign to release native memory.
 *
 * Note that JVM does not manage objects created in native code here (mainly
 * BAMRecords); holding these objects requires a lot of memory but releasing
 * them too soon will crash the JVM.
 *
 * Calling initNativeAlign and closeNativeAlign multiple times is safe; aux is
 * protected with a mutex lock and repeated calls will be ignored.
 */

ktp_aux_t* aux = nullptr;
mutex aux_mutex;
hdfsFS hdfs_fs1 = nullptr;
mutex hdfs_fs1_mutex;
hdfsFS hdfs_fs2 = nullptr;
mutex hdfs_fs2_mutex;
size_t total_output_count = 0;
mutex total_output_count_mutex;
size_t total_input_count = 0;
mutex total_input_count_mutex;

int hdfs_buffer_size = 0;
short hdfs_replication = 0;
size_t hdfs_block_size = 0;
int bwa_chunk_size = 10000000;

static void ReadHDFS(hdfsFS hdfs_fs, hdfsFile hdfs_file, int write_fd)
{
    const size_t kBufferSize = 1<<16;
    char* buffer = new char[kBufferSize];
    tSize read_bytes = 0;
    for(;;)
    {
        read_bytes = hdfsRead(hdfs_fs, hdfs_file, buffer, kBufferSize);
        if(read_bytes < 0)
        {
            if(errno==EINTR)
            {
                clog<<"libgatkbwa:WARNING Retry read from HDFS: "<<strerror(errno)<<endl;
                continue;
            }
            clog<<"libgatkbwa:ERROR Cannot read from HDFS: "<<strerror(errno)<<endl;
            exit(2);
        }
        read_bytes = write(write_fd, buffer, read_bytes);
        if(read_bytes==0)
        {
            break;
        }
        else if(read_bytes < 0)
        {
            clog<<"libgatkbwa:ERROR Cannot write to pipe: "<<strerror(errno)<<endl;
            exit(2);
        }
    }
    if(hdfsCloseFile(hdfs_fs, hdfs_file) < 0)
    {
        clog<<"libgatkbwa:ERROR Cannot close HDFS file: "<<strerror(errno)<<endl;
        exit(2);
    }
    close(write_fd);
    delete[] buffer;
}

extern "C" {

/*
 * Class:     org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine
 * Method:    doNativeAlign
 * Signature: (Ljava/lang/String;Ljava/lang/String;)Ljava/util/List;
 */
JNIEXPORT jobject JNICALL Java_org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_doNativeAlign
    (JNIEnv * env, jobject obj, jstring fastq_file1, jstring fastq_file2, jlong start_idx)
{
    clog<<"libgatkbwa:INFO Start native calculation task @ "<<start_idx<<"\n";

    // java.util.ArrayList
    jclass java_util_ArrayList;
    jmethodID java_util_ArrayList_;
    jmethodID java_util_ArrayList_add;

    java_util_ArrayList      = env->FindClass("java/util/ArrayList");
    java_util_ArrayList_     = env->GetMethodID(java_util_ArrayList, "<init>", "(I)V");
    java_util_ArrayList_add  = env->GetMethodID(java_util_ArrayList, "add", "(Ljava/lang/Object;)Z");

    // htsjdk.samtools.BAMRecord
    jclass htsjdk_samtools_BAMRecord;
    jmethodID htsjdk_samtools_BAMRecord_;

    htsjdk_samtools_BAMRecord   = env->FindClass("htsjdk/samtools/BAMRecord");
    htsjdk_samtools_BAMRecord_  = env->GetMethodID(htsjdk_samtools_BAMRecord, "<init>", "(Lhtsjdk/samtools/SAMFileHeader;IISSIIIIIII[B)V");

    // org.broadinstitute.hellbender.tools.spark.bwa.NativeBwaSparkEngine
    jclass org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine;
    jfieldID org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_headerObject_;

    org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine = 
        env->FindClass("org/broadinstitute/hellbender/tools/spark/bwa/NativeBwaSparkEngine");
    if(nullptr==org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine)
    {
        clog<<"libgatkbwa:ERROR Cannot find class org.broadinstitute.hellbender.tools.spark.bwa.NativeBwaSparkEngine"<<std::endl;
    }
    org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_headerObject_ =
        env->GetFieldID(org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine,
            "headerObject", "Lhtsjdk/samtools/SAMFileHeader;");
    if(nullptr==org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_headerObject_)
    {
        clog<<"libgatkbwa:ERROR Cannot find field headerObject"<<std::endl;
    }
    jobject org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_headerObject = 
        (jobject)env->GetObjectField(obj, 
            org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_headerObject_);

    // get C++ filename from Java String
    const char* fastq_file1_chars = env->GetStringUTFChars(fastq_file1, nullptr);
    string fastq_file1_string(fastq_file1_chars);
    env->ReleaseStringUTFChars(fastq_file1, fastq_file1_chars);
    string fastq_file2_string;
    const bool kIsPaired = fastq_file2!=nullptr;
    if(kIsPaired)
    {
        const char* fastq_file2_chars = env->GetStringUTFChars(fastq_file2, nullptr);
        fastq_file2_string = string(fastq_file2_chars);
        env->ReleaseStringUTFChars(fastq_file2, fastq_file2_chars);
    }

    // open input file
    int read_count = 0;
    gzFile fastq_gzfile1 = nullptr;
    gzFile fastq_gzfile2 = nullptr;
    kseq_t* fastq_kseq1 = nullptr;
    kseq_t* fastq_kseq2 = nullptr;
    string hdfs1_nodename;
    string hdfs2_nodename;
    thread* hdfs_read_thread1 = nullptr;
    thread* hdfs_read_thread2 = nullptr;
    const string kHDFSProto = "hdfs://";
    const bool kInput1IsHDFS = fastq_file1_string.substr(0, kHDFSProto.size()) == kHDFSProto;
    const bool kInput2IsHDFS = fastq_file2_string.substr(0, kHDFSProto.size()) == kHDFSProto;

    size_t local_output_count = 0;

    if(kInput1IsHDFS)
    {
        // connect if not already connected
        hdfs1_nodename = fastq_file1_string.substr(0, fastq_file1_string.find('/', kHDFSProto.size()));
        hdfs_fs1_mutex.lock();
        if(hdfs_fs1==nullptr)
        {
            clog<<"libgatkbwa:INFO Connect to "<<hdfs1_nodename<<endl;
            hdfsBuilder* hdfs_builder = hdfsNewBuilder();
            hdfsBuilderSetNameNode(hdfs_builder, hdfs1_nodename.c_str());
            hdfs_fs1 = hdfsBuilderConnect(hdfs_builder);
            if(hdfs_fs1==nullptr)
            {
                clog<<"libgatkbwa:ERROR Cannot connect to "<<hdfs1_nodename<<endl;
            }
        }
        hdfs_fs1_mutex.unlock();

        // start read thread
        hdfsFile hdfs_file = hdfsOpenFile(hdfs_fs1, fastq_file1_string.c_str(), O_RDONLY, hdfs_buffer_size, hdfs_replication, hdfs_block_size);
        if(hdfs_file==nullptr)
        {
            clog<<"libgatkbwa:ERROR Cannot open file \""<<fastq_file1_string<<"\" on HDFS"<<endl;
        }
        int pipe_fd[2];
        pipe(pipe_fd);
        hdfs_read_thread1 = new thread(ReadHDFS, hdfs_fs1, hdfs_file, pipe_fd[1]);
        fastq_gzfile1 = gzdopen(pipe_fd[0], "r");
        if(fastq_gzfile1==nullptr)
        {
            clog<<"libgatkbwa:ERROR Cannot open fastq 1 from pipe\n";
        }
    }
    else
    {
        fastq_gzfile1 = gzopen(fastq_file1_string.c_str(), "r");
        if(fastq_gzfile1==nullptr)
        {
            clog<<"libgatkbwa:ERROR Cannot open fastq 1: "<<fastq_file1_string<<endl;
        }
    }
    clog<<"libgatkbwa:INFO Read fastq 1 from "<<fastq_file1_string<<endl;
    fastq_kseq1 = kseq_init(fastq_gzfile1);
    if(nullptr == fastq_gzfile1)
    {
        clog<<"libgatkbwa:ERROR Cannot read fastq 1\n";
    }

    if(kIsPaired)
    {
        if(kInput2IsHDFS)
        {
            hdfs2_nodename = fastq_file2_string.substr(0, fastq_file2_string.find('/', kHDFSProto.size()));
            hdfs_fs2_mutex.lock();
            if(hdfs_fs2==nullptr)
            {
                if(hdfs1_nodename==hdfs2_nodename)
                {
                    hdfs_fs2 = hdfs_fs1;
                }
                else
                {
                    clog<<"libgatkbwa:INFO Connect to "<<hdfs2_nodename<<endl;
                    hdfsBuilder* hdfs_builder = hdfsNewBuilder();
                    hdfsBuilderSetNameNode(hdfs_builder, hdfs2_nodename.c_str());
                    hdfs_fs2 = hdfsBuilderConnect(hdfs_builder);
                    if(hdfs_fs2==nullptr)
                    {
                        clog<<"libgatkbwa:ERROR Cannot connect to "<<hdfs2_nodename<<endl;
                    }
                }
            }
            hdfs_fs2_mutex.unlock();

            // start read thread
            hdfsFile hdfs_file = hdfsOpenFile(hdfs_fs2, fastq_file1_string.c_str(), O_RDONLY, hdfs_buffer_size, hdfs_replication, hdfs_block_size);
            if(hdfs_file==nullptr)
            {
                clog<<"libgatkbwa:ERROR Cannot open file \""<<fastq_file2_string<<"\" on HDFS"<<endl;
            }
            int pipe_fd[2];
            pipe(pipe_fd);
            hdfs_read_thread2 = new thread(ReadHDFS, hdfs_fs2, hdfs_file, pipe_fd[1]);
            fastq_gzfile2 = gzdopen(pipe_fd[0], "r");
            if(fastq_gzfile2==nullptr)
            {
                clog<<"libgatkbwa:ERROR Cannot open fastq 1 from pipe\n";
            }
        }
        else
        {
            fastq_gzfile2 = gzopen(fastq_file2_string.c_str(), "r");
            if(fastq_gzfile2==nullptr)
            {
                clog<<"libgatkbwa:ERROR Cannot open fastq 2: "<<fastq_file2_string<<endl;
            }
        }
        clog<<"libgatkbwa:INFO Read fastq 2 from "<<fastq_file2_string<<endl;
        fastq_kseq2 = kseq_init(fastq_gzfile2);
        if(nullptr == fastq_gzfile2)
        {
            clog<<"libgatkbwa:ERROR Cannot read fastq 2\n";
        }
    }

    jobject bam_record_list = env->NewObject(java_util_ArrayList, java_util_ArrayList_, 1<<16); // TODO: better guess on init capacity

    bool all_right = true;
    for(;;)
    {
        bseq1_t* seqs = bseq_read(bwa_chunk_size, &read_count, fastq_kseq1, fastq_kseq2);
        if(seqs==nullptr)
        {
            break;
        }
        clog<<"libgatkbwa:INFO Get "<<read_count<<" sequences\n";
        total_input_count_mutex.lock();
        total_input_count += read_count;
        total_input_count_mutex.unlock();

        /* Most likely comments in fastq files is not compatible with SAM file,
         * so do not copy comments from fastq. */
        if(!aux->copy_comment)
        {
            for(int i = 0; i < read_count; i++)
            {
                free(seqs[i].comment);
                seqs[i].comment = 0;
            }
        }

        mem_alnreg_v* alnreg = new mem_alnreg_v[read_count];

        for (int i = 0; i < read_count; i++)
        {
            mem_chain_v chains = seq2chain(aux, &seqs[i]);
            kv_init(alnreg[i]);
            for (int j = 0; j < int(chains.n); j++)
            {
                mem_chain2aln(
                    aux->opt, 
                    aux->idx->bns, 
                    aux->idx->pac,
                    seqs[i].l_seq,
                    (uint8_t*)seqs[i].seq,
                    &chains.a[j],
                    &alnreg[i]);
                free(chains.a[j].seeds);
            }
            free(chains.a);

            // Post-process each chain before output
            alnreg[i].n = mem_sort_dedup_patch(
                aux->opt,
                aux->idx->bns,
                aux->idx->pac,
                (uint8_t*)seqs[i].seq,
                alnreg[i].n,
                alnreg[i].a);

            for (int j = 0; j < int(alnreg[i].n); j++)
            {
                mem_alnreg_t *p = &alnreg[i].a[j];
                if (p->rid >= 0 && aux->idx->bns->anns[p->rid].is_alt)
                    p->is_alt = 1;
            }
        }

        if(kIsPaired)
        {
            mem_pestat_t pes[4];
            mem_pestat(aux->opt, aux->idx->bns->l_pac, read_count, alnreg, pes);

#ifdef USE_HTSLIB
            for (int i =0; i< read_count/2; i++) {
                seqs[i<<1].bams = bams_init();
                seqs[1+(i<<1)].bams = bams_init();
                mem_sam_pe(
                    aux->opt,
                    aux->idx->bns,
                    aux->idx->pac,
                    pes,
                    (start_idx>>1)+i,
                    &seqs[i<<1],
                    &alnreg[i<<1],
                    aux->h);
            }
#else//ifndef USE_HTSLIB
            for (int i = 0; i < read_count/2; i++) {
                mem_sam_pe(
                    aux->opt,
                    aux->idx->bns,
                    aux->idx->pac,
                    pes,
                    (start_idx>>1)+i,
                    &seqs[i<<1],
                    &alnreg[i<<1]);
            }
#endif//ifndef USE_HTSLIB
        }
        else
        {
            for (int i=0; i<read_count; i++)
            {
#ifdef USE_HTSLIB
                seqs[i].bams = bams_init();
#endif//ifndef USE_HTSLIB
                mem_mark_primary_se(
                    aux->opt,
                    alnreg[i].n,
                    alnreg[i].a,
                    start_idx+i
                    );
                mem_reg2sam(
                    aux->opt,
                    aux->idx->bns,
                    aux->idx->pac,
                    &seqs[i],
                    &alnreg[i],
                    0,
                    0,
                    aux->h
                    );
            }
        }

        start_idx += read_count;

        freeAligns(alnreg, read_count);
        
        // Free fields in seq except sam
        for (int i = 0; i < read_count; i++) {
            free(seqs[i].name);
            free(seqs[i].comment);
            free(seqs[i].seq);
            free(seqs[i].qual);
        }
        
        for(int i = 0; i < read_count; ++i)
        {
            for(int j = 0; j < seqs[i].bams->l && all_right; ++j)
            {
                const bam1_core_t& bam = seqs[i].bams->bams[j]->core;
                jbyteArray rest_of_bam = env->NewByteArray(seqs[i].bams->bams[j]->l_data-bam.l_extranul);
                env->SetByteArrayRegion((jbyteArray)rest_of_bam, 0, bam.l_qname-bam.l_extranul, (jbyte*)(seqs[i].bams->bams[j]->data));
                env->SetByteArrayRegion((jbyteArray)rest_of_bam, bam.l_qname-bam.l_extranul, seqs[i].bams->bams[j]->l_data-bam.l_qname, (jbyte*)(seqs[i].bams->bams[j]->data+bam.l_qname));
                jobject bam_record = env->NewObject(htsjdk_samtools_BAMRecord, htsjdk_samtools_BAMRecord_,
                    org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_headerObject,
                    jint(bam.tid),
                    jint(bam.pos+1),
                    jshort(bam.l_qname-bam.l_extranul),
                    jshort(bam.qual),
                    jint(bam.bin),
                    jint(bam.n_cigar),
                    jint(bam.flag),
                    jint(bam.l_qseq),
                    jint(bam.mtid),
                    jint(bam.mpos+1),
                    jint(bam.isize),
                    rest_of_bam
                    );
                if(env->ExceptionCheck())
                {
                    clog<<"libgatkbwa:ERROR Cannot create BAMRecord\n";
                    env->ExceptionDescribe();
                    all_right = false;
                }
                if(env->IsSameObject(bam_record, nullptr))
                {
                    clog<<"libgatkbwa:ERROR Object bam_record (type:BAMRecord) is null\n";
                    all_right = false;
                }
                if(env->IsSameObject(bam_record_list, nullptr))
                {
                    clog<<"libgatkbwa:ERROR Object bam_record_list (type:ArrayList<BAMRecord>) is null\n";
                    all_right = false;
                }
                env->CallBooleanMethod(bam_record_list, java_util_ArrayList_add, bam_record);
                local_output_count++;
                if(env->ExceptionCheck())
                {
                    clog<<"libgatkbwa:ERROR Failed to add bam_record (type:BAMRecord) to bam_record_list (type:ArrayList<BAMRecord>)\n";
                    env->ExceptionDescribe();
                    all_right = false;
                }
                env->DeleteLocalRef(rest_of_bam);
                env->DeleteLocalRef(bam_record);
            }
            bams_destroy(seqs[i].bams);
        }

        free(seqs);
    }

    total_output_count_mutex.lock();
    total_output_count += local_output_count;
    clog<<"libgatkbwa:INFO Total records processed: "<<total_output_count<<"\n";
    total_output_count_mutex.unlock();

    if(hdfs_read_thread1!=nullptr)
    {
        hdfs_read_thread1->join();
        delete hdfs_read_thread1;
        hdfs_read_thread1 = nullptr;
    }
    if(hdfs_read_thread2!=nullptr)
    {
        hdfs_read_thread2->join();
        delete hdfs_read_thread2;
        hdfs_read_thread2 = nullptr;
    }

    if(fastq_gzfile1!=nullptr)
    {
        gzclose(fastq_gzfile1);
    }
    if(fastq_gzfile2!=nullptr)
    {
        gzclose(fastq_gzfile2);
    }
    
    kseq_destroy(fastq_kseq1);
    kseq_destroy(fastq_kseq2);

    if(all_right)
    {
        clog<<"libgatkbwa:INFO Native calculation task finished\n";
    }
    else
    {
        clog<<"libgatkbwa:INFO Native calculation went wrong; go back to Java\n";
    }
    return bam_record_list;
}

/*
 * Class:     org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine
 * Method:    initNativeAlign
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_initNativeAlign(JNIEnv* env, jobject obj)
{
    aux_mutex.lock();
    if(aux!=nullptr)
    {
        clog<<"libgatkbwa:INFO Native resources already initialized\n";
        aux_mutex.unlock();
        return;
    }
    clog<<"libgatkbwa:INFO Initialize native resources\n";

    // org.broadinstitute.hellbender.tools.spark.bwa.NativeBwaSparkEngine
    const char* index_file_name;
    const char* read_group_header_line;
    jclass org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine;
    jfieldID org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_indexFileName_;
    jfieldID org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_readGroupHeaderLine_;
    jfieldID org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_headerString_;
    jstring org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_indexFileName;
    jstring org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_readGroupHeaderLine;

    org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine = 
        env->FindClass("org/broadinstitute/hellbender/tools/spark/bwa/NativeBwaSparkEngine");
    if(nullptr==org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine)
    {
        clog<<"libgatkbwa:ERROR Cannot find class org.broadinstitute.hellbender.tools.spark.bwa.NativeBwaSparkEngine"<<std::endl;
    }
    org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_indexFileName_ =
        env->GetFieldID(org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine,
            "indexFileName", "Ljava/lang/String;");
    if(nullptr==org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_indexFileName_)
    {
        clog<<"libgatkbwa:ERROR Cannot find field indexFileName"<<std::endl;
    }
    org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_readGroupHeaderLine_ =
        env->GetFieldID(org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine,
            "readGroupHeaderLine", "Ljava/lang/String;");
    if(nullptr==org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_readGroupHeaderLine_)
    {
        clog<<"libgatkbwa:ERROR Cannot find field readGroupHeaderLine"<<std::endl;
    }
    org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_headerString_ =
        env->GetFieldID(org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine,
            "headerString", "Ljava/lang/String;");
    if(nullptr==org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_headerString_)
    {
        clog<<"libgatkbwa:ERROR Cannot find field headerString"<<std::endl;
    }
    org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_indexFileName = 
        (jstring)env->GetObjectField(obj, 
            org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_indexFileName_);
    org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_readGroupHeaderLine = 
        (jstring)env->GetObjectField(obj, 
            org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_readGroupHeaderLine_);

    // bwa parameters
    index_file_name = env->GetStringUTFChars(
        org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_indexFileName,
        nullptr);
    clog<<"libgatkbwa:INFO Using index file: "<<index_file_name<<std::endl;
    read_group_header_line = env->GetStringUTFChars(
        org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_readGroupHeaderLine,
        nullptr);
    clog<<"libgatkbwa:INFO Using read group header: "<<read_group_header_line<<std::endl;

    aux = new ktp_aux_t;
    memset(aux, 0, sizeof(ktp_aux_t));

    vector<const char*> bwa_args;
    bwa_args.push_back("mem");
    bwa_args.push_back("-Ma");
    bwa_args.push_back("-R");
    bwa_args.push_back(read_group_header_line);
    bwa_args.push_back(index_file_name);
    // placeholders, prevent preprocess from complaining
    bwa_args.push_back("");

    const int header_string_max_len = 1<<20;
    char* header_string = new char[header_string_max_len]();

    // Catch header string from preprocess stdout.
    int out_pipe[2];
    int saved_stdout_fd;
    saved_stdout_fd = dup(STDOUT_FILENO);
    if(pipe(out_pipe) != 0)
    {
        clog<<"libgatkbwa:ERROR Cannot make pipe for redirecting header output\n";
        aux_mutex.unlock();
        return;
    }
    dup2(out_pipe[1], STDOUT_FILENO);
    close(out_pipe[1]);

    pre_process(bwa_args.size(), (char**)&bwa_args[0], aux, true);
    fflush(stdout);

    read(out_pipe[0], header_string, header_string_max_len);
    close(out_pipe[0]);
    dup2(saved_stdout_fd, STDOUT_FILENO);
    jstring header_jstring = env->NewStringUTF(header_string);
    env->SetObjectField(obj, org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_headerString_, header_jstring);
    delete[] header_string;
    env->ReleaseStringUTFChars(
        org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_indexFileName,
        index_file_name);
    env->ReleaseStringUTFChars(
        org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_readGroupHeaderLine,
        read_group_header_line);
    clog<<"libgatkbwa:INFO Native resources initialized\n";
    aux_mutex.unlock();
    total_input_count_mutex.lock();
    total_input_count = 0;
    total_input_count_mutex.unlock();
    total_output_count_mutex.lock();
    total_output_count = 0;
    total_output_count_mutex.unlock();
}

/*
 * Class:     org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine
 * Method:    closeNativeAlign
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine_closeNativeAlign(JNIEnv* env, jobject obj)
{
    aux_mutex.lock();
    if(aux == nullptr)
    {
        clog<<"libgatkbwa:WARNING No native resource found\n";
        aux_mutex.unlock();
        return;
    }
    clog<<"libgatkbwa:INFO Release native resources\n";

    if(hdfs_fs1!=nullptr)
    {
        hdfsDisconnect(hdfs_fs1);
    }
    if(hdfs_fs2!=nullptr && hdfs_fs1!=hdfs_fs2)
    {
        hdfsDisconnect(hdfs_fs2);
    }

    free(aux->opt);
    bwa_idx_destroy(aux->idx);
    kseq_destroy(aux->ks);
    if(aux->ks2!=nullptr)
    {
        kseq_destroy(aux->ks2);
    }

    delete aux;
    aux = nullptr;

    clog<<"libgatkbwa:INFO Native resources released\n";
    aux_mutex.unlock();
}

}// extern "C"

// These two functions are copied from bwa-flow code and require C++ linkage.
mem_chain_v seq2chain(
    ktp_aux_t *aux,
    bseq1_t *seqs
)
{
    int i;
    mem_chain_v chn;
    for (i = 0; i < seqs->l_seq; ++i) // convert to 2-bit encoding if we have not done so
    {
        seqs->seq[i] = seqs->seq[i] < 4? seqs->seq[i] : nst_nt4_table[(int)seqs->seq[i]];
    }
    chn = mem_chain(aux->opt, aux->idx->bwt, aux->idx->bns, seqs->l_seq, (uint8_t*)seqs->seq, 0);
    // the 0 should be reconsidered
    chn.n = mem_chain_flt(aux->opt, chn.n, chn.a);
    mem_flt_chained_seeds(aux->opt, aux->idx->bns, aux->idx->pac, seqs->l_seq, (uint8_t*)seqs->seq, chn.n, chn.a);
    return chn;
}

void freeAligns(mem_alnreg_v* alnreg, int batch_num)
{
    for (int i = 0; i < batch_num; i++)
    {
        free(alnreg[i].a);
    }
    delete[] alnreg;
}

