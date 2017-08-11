#include<iostream>
#include<string>

#include<cxxopts.hpp>

using namespace std;

void SplitFASTQ(const int kVerboseFlag, const size_t kBatchSize, const string& kInputFastq1, const string& kOutputFastq1, const string& kInputFastq2, const string& kOutputFastq2, const int kHDFSBufferSize = 0, const short kHDFSReplication = 0, const size_t kHDFSBlockSize = 0);

int main(int argc, char* argv[])
{
    cxxopts::Options options(argv[0], "Align fastq file(s) to a fixed block size");
    options.add_options()
        ("v,verbose", "Enable verbose logging")
        ("h,help", "Print this help message")
        ("s,batch-size", "Number of reads per batch", cxxopts::value<size_t>()->default_value("100000"))
        ("i,input-fastq", "Input fastq file name", cxxopts::value<string>())
        ("o,output-fastq", "Output fastq file prefix", cxxopts::value<string>())
        ("I,paired-input-fastq", "Paired input fastq file name", cxxopts::value<string>())
        ("O,paired-output-fastq", "Paired output fastq file prefix", cxxopts::value<string>())
        ("hdfs-buffer-size", "HDFS buffer size", cxxopts::value<int>()->default_value("0"))
        ("hdfs-replication", "HDFS replication", cxxopts::value<short>()->default_value("3"))
        ("hdfs-block-size", "HDFS block size", cxxopts::value<size_t>()->default_value("0"))
        ;

    options.parse(argc, argv);

    // validate options
    if(options.count("help")>0)
    {
        clog<<options.help()<<endl;
        exit(0);
    }
    if(options.count("input-fastq")==0 || options.count("output-fastq")==0)
    {
        clog<<"Error: Fastq file must be specified\n"<<endl;
        clog<<options.help()<<endl;
        exit(1);
    }

    if((options.count("paired-input-fastq")==0) ^ (options.count("paired-output-fastq")==0))
    {
        clog<<"Error: Input and output files for paired fastq file must be specified\n"<<endl;
        clog<<options.help()<<endl;
        exit(1);
    }

    // formal parameters
    const int kVerboseFlag = options["verbose"].as<bool>() ? 1 : 0;
    const size_t kBatchSize = options["batch-size"].as<size_t>();
    const string kInputFastq1 = options["input-fastq"].as<string>();
    const string kOutputFastq1 = options["output-fastq"].as<string>();
    const string kInputFastq2 = options["paired-input-fastq"].as<string>();
    const string kOutputFastq2 = options["paired-output-fastq"].as<string>();

    // set CLASSPATH
    if(getenv("CLASSPATH")==nullptr)
    {
        FILE* fp = popen("hadoop classpath --glob", "r");
        if(fp==nullptr)
        {
            clog<<"Error: Cannot set CLASSPATH with `hadoop classpath --glob`\n";
            exit(1);
        }
        char buffer[1<<16];
        fgets(buffer, 1<<16, fp);
        setenv("CLASSPATH", buffer, 0);
        pclose(fp);
    }

    SplitFASTQ(kVerboseFlag, kBatchSize, kInputFastq1, kOutputFastq1, kInputFastq2, kOutputFastq2, options["hdfs-buffer-size"].as<int>(), options["hdfs-replication"].as<short>(), options["hdfs-block-size"].as<size_t>());

    return 0;
}
