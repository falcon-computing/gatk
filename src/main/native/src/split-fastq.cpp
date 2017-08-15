#include<cstring>

#include<iostream>
#include<string>
#include<thread>
#include<tuple>
#include<vector>

#include<fcntl.h>
#include<sys/types.h>
#include<sys/wait.h>
#include<unistd.h>

#include<zlib.h>
#include<hdfs.h>

using namespace std;

static void GetBatch(FILE* input, vector<tuple<char*, size_t> >* output, size_t* line_count_ptr)
{
    // output.size() is kBatchSize*4 
    size_t line_count = 0;
    for(auto& read: *output)
    {
        char* ptr = get<0>(read);
        size_t n = 0;
        ssize_t read_bytes = getline(&ptr, &n, input);
        if(read_bytes <= 0)
        {
            if(!feof(input))
            {
                clog<<"ERROR: Cannot read file: "<<strerror(errno)<<endl;
                exit(2);
            }
            break;
        }
        get<0>(read) = ptr;
        get<1>(read) = read_bytes;
        ++line_count;
    }
    if((line_count & 0x3) != 0)
    {
        clog<<"ERROR: Corrupted input fastq file"<<endl;
        exit(1);
    }
    *line_count_ptr = line_count;
}

static void WriteBatch(const vector<tuple<char*, size_t> >* output, size_t line_count, gzFile gzfile)
{
    size_t line = 0;
    int gzerr = 0;
    for(const auto& read: *output)
    {
        gzerr = gzwrite(gzfile, get<0>(read), get<1>(read));
        if(gzerr == 0)
        {
            clog<<"ERROR: Cannot write file: "<<gzerror(gzfile, &gzerr)<<endl;
            exit(2);
        }
        ++line;
        free(get<0>(read));
        if(line>=line_count)
        {
            break;
        }
    }
    gzclose(gzfile);
    delete output;
}

static void WriteHDFS(int read_fd, hdfsFS hdfs_fs, hdfsFile hdfs_file)
{
    const size_t kBufferSize = 1<<16;
    char* buffer = new char[kBufferSize];
    ssize_t read_bytes = 0;
    for(;;)
    {
        read_bytes = read(read_fd, buffer, kBufferSize);
        if(read_bytes==0)
        {
            break;
        }
        else if(read_bytes < 0)
        {
            clog<<"ERROR: Cannot read from pipe: "<<strerror(errno)<<endl;
            exit(2);
        }
        tSize write_bytes = hdfsWrite(hdfs_fs, hdfs_file, buffer, read_bytes);
        if(write_bytes < 0)
        {
            clog<<"ERROR: Cannot write to HDFS"<<endl;
            exit(2);
        }
    }
    close(read_fd);
    if(hdfsCloseFile(hdfs_fs, hdfs_file) < 0)
    {
        clog<<"ERROR: Cannot close HDFS file: "<<strerror(errno)<<endl;
        exit(2);
    }
    delete[] buffer;
}

void SplitFASTQ(const int kVerboseFlag, const size_t kBatchSize, const string& kInputFastq1, const string& kOutputFastq1, const string& kInputFastq2, const string& kOutputFastq2, const int kHDFSBufferSize = 0, const short kHDFSReplication = 0, const size_t kHDFSBlockSize = 0, const int8_t kCompressionLevel = 1)
{
    // local parameters
    vector<pid_t> children;
    vector<thread> threads;
    const bool kIsPaired = !kInputFastq2.empty();
    char gzip[] = "gzip";
    char decompress[] = "-cd";

    pid_t input_pid1 = -1;
    pid_t input_pid2 = -1;
    int input_pipe1[2] = {};
    int input_pipe2[2] = {};

    if(pipe(input_pipe1)!=0)
    {
        clog<<"ERROR: Cannot make pipe: "<<strerror(errno)<<endl;
        exit(2);
    }
    input_pid1 = fork();
    if(input_pid1==-1)
    {
        clog<<"ERROR: Cannot fork: "<<strerror(errno)<<endl;
        exit(2);
    }
    else if(input_pid1==0)
    {
        // child process
        dup2(input_pipe1[1], STDOUT_FILENO);
        close(input_pipe1[0]);
        close(input_pipe1[1]);
        char input_filename[1<<16] = {};
        strncpy(input_filename, kInputFastq1.c_str(), kInputFastq1.size());
        char* args[3] = {gzip, decompress, input_filename};
        execvp(gzip, args);
        clog<<"ERROR: Cannot exec gzip: "<<strerror(errno)<<endl;
        exit(2);
    }
    else
    {
        children.push_back(input_pid1);
        close(input_pipe1[1]);
    }

    FILE* input1 = fdopen(input_pipe1[0], "r");
    FILE* input2 = nullptr;
    vector<tuple<char*, size_t> >* output1;
    vector<tuple<char*, size_t> >* output2;

    if(kIsPaired)
    {
        if(pipe(input_pipe2)!=0)
        {
            clog<<"ERROR: Cannot make pipe: "<<strerror(errno)<<endl;
            exit(2);
        }
        input_pid2 = fork();
        if(input_pid2==-1)
        {
            clog<<"ERROR: Cannot fork: "<<strerror(errno)<<endl;
            exit(2);
        }
        else if(input_pid2==0)
        {
            // child process
            dup2(input_pipe2[1], STDOUT_FILENO);
            close(input_pipe2[0]);
            close(input_pipe2[1]);
            char input_filename[1<<16] = {};
            strncpy(input_filename, kInputFastq2.c_str(), kInputFastq2.size());
            char* args[3] = {gzip, decompress, input_filename};
            execvp(gzip, args);
            clog<<"ERROR: Cannot exec gzip: "<<strerror(errno)<<endl;
            exit(2);
        }
        else
        {
            children.push_back(input_pid2);
            close(input_pipe2[1]);
        }
        input2 = fdopen(input_pipe2[0], "r");
    }

    // initialize hdfs
    const string kHDFSProto = "hdfs://";
    const bool kOutput1IsHDFS = kOutputFastq1.substr(0, kHDFSProto.size()) == kHDFSProto;
    const bool kOutput2IsHDFS = kOutputFastq2.substr(0, kHDFSProto.size()) == kHDFSProto;
    hdfsFS hdfs_fs1 = nullptr;
    hdfsFS hdfs_fs2 = nullptr;
    string hdfs1_nodename;
    string hdfs2_nodename;
    string output_fastq1_hdfs_path;
    string output_fastq2_hdfs_path;

    if(kOutput1IsHDFS)
    {
        hdfs1_nodename = kOutputFastq1.substr(0, kOutputFastq1.find('/', kHDFSProto.size()));
        output_fastq1_hdfs_path = kOutputFastq1.substr(hdfs1_nodename.size(), string::npos);
        if(kVerboseFlag)
        {
            clog<<"INFO: Connect to "<<hdfs1_nodename<<endl;
        }
        hdfsBuilder* hdfs_builder = hdfsNewBuilder();   // freed by hdfsBuilderConnect
        hdfsBuilderSetNameNode(hdfs_builder, hdfs1_nodename.c_str());
        hdfs_fs1 = hdfsBuilderConnect(hdfs_builder);
        if(hdfs_fs1==nullptr)
        {
            clog<<"ERROR: Cannot connect to "<<hdfs1_nodename<<endl;
            exit(2);
        }
    }

    if(kOutput2IsHDFS)  // implies kIsPaired
    {
        hdfs2_nodename = kOutputFastq2.substr(0, kOutputFastq2.find('/', kHDFSProto.size()));
        output_fastq2_hdfs_path = kOutputFastq2.substr(hdfs2_nodename.size(), string::npos);
        if(hdfs1_nodename != hdfs2_nodename)
        {
            if(kVerboseFlag)
            {
                clog<<"INFO: Connect to "<<hdfs2_nodename<<endl;
            }
            hdfsBuilder* hdfs_builder = hdfsNewBuilder();   // freed by hdfsBuilderConnect
            hdfsBuilderSetNameNode(hdfs_builder, hdfs2_nodename.c_str());
            hdfs_fs2 = hdfsBuilderConnect(hdfs_builder);
            if(hdfs_fs2==nullptr)
            {
                clog<<"ERROR: Cannot connect to "<<hdfs2_nodename<<endl;
                exit(2);
            }
        }
        else
        {
            hdfs_fs2 = hdfs_fs1;
        }
    }

    int batch_id = 0;
    for(;;++batch_id)
    {
        size_t line_count1 = 0;
        size_t line_count2 = 0;
        output1 = new vector<tuple<char*, size_t> >(kBatchSize*4);
        if(kIsPaired)
        {
            output2 = new vector<tuple<char*, size_t> >(kBatchSize*4);
            thread thread1 = thread(GetBatch, input1, output1, &line_count1);
            thread thread2 = thread(GetBatch, input2, output2, &line_count2);
            thread1.join();
            thread2.join();
            if(line_count1!=line_count2)
            {
                clog<<"ERROR: Two input fastq files do not match"<<endl;
                exit(1);
            }
        }
        else
        {
            GetBatch(input1, output1, &line_count1);
        }
        if(line_count1==0)
        {
            break;
        }

        char gz_write_mode[3] = "w0";
        gz_write_mode[1] += (kCompressionLevel>=0 && kCompressionLevel<=9) ? kCompressionLevel : 1;
        // write the outputs
        if(kOutput1IsHDFS)
        {
            string output_filename_string = output_fastq1_hdfs_path+".part"+to_string(batch_id);
            int write_pipe[2];
            pipe(write_pipe);
            gzFile output_file1 = gzdopen(write_pipe[1], gz_write_mode);
            hdfsFile hdfs_file1 = hdfsOpenFile(hdfs_fs1, output_filename_string.c_str(), O_WRONLY,
                kHDFSBufferSize, kHDFSReplication, kHDFSBlockSize);
            if(hdfs_file1==nullptr)
            {
                clog<<"ERROR: Cannot open file on HDFS"<<endl;
                exit(2);
            }
            threads.push_back(thread(WriteBatch, output1, line_count1, output_file1));
            threads.push_back(thread(WriteHDFS, write_pipe[0], hdfs_fs1, hdfs_file1));
        }
        else
        {
            string output_filename_string = kOutputFastq1+".part"+to_string(batch_id);
            gzFile output_file1 = gzopen(output_filename_string.c_str(), gz_write_mode);
            threads.push_back(thread(WriteBatch, output1, line_count1, output_file1));
        }

        if(kIsPaired)
        {
            if(kOutput2IsHDFS)
            {
                string output_filename_string = output_fastq2_hdfs_path+".part"+to_string(batch_id);
                int write_pipe[2];
                pipe(write_pipe);
                gzFile output_file2 = gzdopen(write_pipe[1], gz_write_mode);
                hdfsFile hdfs_file2 = hdfsOpenFile(hdfs_fs2, output_filename_string.c_str(), O_WRONLY,
                    kHDFSBufferSize, kHDFSReplication, kHDFSBlockSize);
                if(hdfs_file2==nullptr)
                {
                    clog<<"ERROR: Cannot open file on HDFS"<<endl;
                    exit(2);
                }
                threads.push_back(thread(WriteBatch, output2, line_count2, output_file2));
                threads.push_back(thread(WriteHDFS, write_pipe[0], hdfs_fs2, hdfs_file2));
            }
            else
            {
                string output_filename_string = kOutputFastq2+".part"+to_string(batch_id);
                gzFile output_file2 = gzopen(output_filename_string.c_str(), gz_write_mode);
                threads.push_back(thread(WriteBatch, output2, line_count2, output_file2));
            }
        }

        if(kVerboseFlag)
        {
            clog<<"INFO: Get batch "<<batch_id<<" with "<<line_count1/4<<(kIsPaired?" paired":"")<<" reads"<<endl;
        }
    }

    // delete extra files
    for(;;++batch_id)
    {
        if(kOutput1IsHDFS)
        {
            string output_filename_string = output_fastq1_hdfs_path+".part"+to_string(batch_id);
            if(hdfsExists(hdfs_fs1, output_filename_string.c_str())==0)
            {
                if(hdfsDelete(hdfs_fs1, output_filename_string.c_str(), 0)!=0)
                {
                    clog<<"WARNING: Cannot remove extra file \""<<output_filename_string<<"\": "<<strerror(errno)<<endl;
                }
            }
            else
            {
                break;
            }
        }
        else
        {
            string output_filename_string = kOutputFastq1+".part"+to_string(batch_id);
            int err = remove(output_filename_string.c_str());
            if(err!=0)
            {
                if(errno!=ENOENT)
                {
                    clog<<"WARNING: Cannot remove extra file \""<<output_filename_string<<"\": "<<strerror(errno)<<endl;
                }
                break;
            }
        }
        if(kIsPaired)
        {
            if(kOutput2IsHDFS)
            {
                string output_filename_string = output_fastq2_hdfs_path+".part"+to_string(batch_id);
                if(hdfsExists(hdfs_fs2, output_filename_string.c_str())==0)
                {
                    if(hdfsDelete(hdfs_fs2, output_filename_string.c_str(), 0)!=0)
                    {
                        clog<<"WARNING: Cannot remove extra file \""<<output_filename_string<<"\": "<<strerror(errno)<<endl;
                    }
                }
                else
                {
                    break;
                }
            }
            else
            {
                string output_filename_string = kOutputFastq2+".part"+to_string(batch_id);
                int err = remove(output_filename_string.c_str());
                if(err!=0)
                {
                    if(errno!=ENOENT)
                    {
                        clog<<"WARNING: Cannot remove extra file \""<<output_filename_string<<"\": "<<strerror(errno)<<endl;
                    }
                    break;
                }
            }
        }
    }

    if(kVerboseFlag)
    {
        clog<<"INFO: Wait for children"<<endl;
    }

    for(size_t i = 0; i<children.size(); ++i)
    {
        wait(nullptr);
    }
    for(auto& t: threads)
    {
        t.join();
    }

    // clean up
    fclose(input1);
    if(kIsPaired)
    {
        fclose(input2);
    }

    if(kOutput1IsHDFS)
    {
        hdfsDisconnect(hdfs_fs1);
    }

    if(kOutput2IsHDFS && hdfs_fs1!=hdfs_fs2)
    {
        hdfsDisconnect(hdfs_fs2);
    }
}

