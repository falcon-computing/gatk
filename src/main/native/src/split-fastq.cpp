#include<cstring>

#include<iostream>
#include<map>
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
                clog<<"WARNING: Retry read from HDFS: "<<strerror(errno)<<endl;
                continue;
            }
            clog<<"ERROR: Cannot read from HDFS: "<<strerror(errno)<<endl;
            exit(2);
        }
        read_bytes = write(write_fd, buffer, read_bytes);
        if(read_bytes==0)
        {
            break;
        }
        else if(read_bytes < 0)
        {
            clog<<"ERROR: Cannot write to pipe: "<<strerror(errno)<<endl;
            exit(2);
        }
    }
    if(hdfsCloseFile(hdfs_fs, hdfs_file) < 0)
    {
        clog<<"ERROR: Cannot close HDFS file: "<<strerror(errno)<<endl;
        exit(2);
    }
    close(write_fd);
    delete[] buffer;
}

int SplitFASTQ(const int kVerboseFlag, const size_t kBatchSize, const string& kInputFastq1, const string& kOutputFastq1, const string& kInputFastq2, const string& kOutputFastq2, const int kHDFSBufferSize = 0, const short kHDFSReplication = 0, const size_t kHDFSBlockSize = 0, const int8_t kCompressionLevel = 1)
{
    // initialize hdfs
    const string kHDFSProto = "hdfs://";
    const bool kInput1IsHDFS = kInputFastq1.substr(0, kHDFSProto.size()) == kHDFSProto;
    const bool kInput2IsHDFS = kInputFastq2.substr(0, kHDFSProto.size()) == kHDFSProto;
    const bool kOutput1IsHDFS = kOutputFastq1.substr(0, kHDFSProto.size()) == kHDFSProto;
    const bool kOutput2IsHDFS = kOutputFastq2.substr(0, kHDFSProto.size()) == kHDFSProto;
    hdfsFS input_hdfs1 = nullptr;
    hdfsFS input_hdfs2 = nullptr;
    hdfsFS output_hdfs1 = nullptr;
    hdfsFS output_hdfs2 = nullptr;
    string input_hdfs1_nodename;
    string input_hdfs2_nodename;
    string output_hdfs1_nodename;
    string output_hdfs2_nodename;
    string output_fastq1_hdfs_path;
    string output_fastq2_hdfs_path;
    map<string, hdfsFS> hdfs_map;

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

    int gz_read_pipe1[2] = {};
    int gz_read_pipe2[2] = {};

    if(kInput1IsHDFS)
    {
        // connect to HDFS if not already connected
        input_hdfs1_nodename = kInputFastq1.substr(0, kInputFastq1.find('/', kHDFSProto.size()));
        if(hdfs_map.count(input_hdfs1_nodename)==0)
        {
            if(kVerboseFlag)
            {
                clog<<"INFO: Connect to "<<input_hdfs1_nodename<<endl;
            }
            hdfsBuilder* hdfs_builder = hdfsNewBuilder();   // freed by hdfsBuilderConnect
            hdfsBuilderSetNameNode(hdfs_builder, input_hdfs1_nodename.c_str());
            input_hdfs1 = hdfsBuilderConnect(hdfs_builder);
            if(input_hdfs1==nullptr)
            {
                clog<<"ERROR: Cannot connect to "<<input_hdfs1_nodename<<endl;
                exit(2);
            }
            hdfs_map[input_hdfs1_nodename] = input_hdfs1;
        }
        else
        {
            input_hdfs1 = hdfs_map[input_hdfs1_nodename];
        }

        // open HDFS file
        if(pipe(gz_read_pipe1)!=0)
        {
            clog<<"ERROR: Cannot make pipe: "<<strerror(errno)<<endl;
            exit(2);
        }
        hdfsFile hdfs_file1 = hdfsOpenFile(input_hdfs1, kInputFastq1.c_str(), O_RDONLY,
            kHDFSBufferSize, kHDFSReplication, kHDFSBlockSize);
        if(hdfs_file1==nullptr)
        {
            clog<<"ERROR: Cannot open file on HDFS"<<endl;
            exit(2);
        }
        threads.push_back(thread(ReadHDFS, input_hdfs1, hdfs_file1, gz_read_pipe1[1]));
    }
    else
    {
        gz_read_pipe1[0] = open(kInputFastq1.c_str(), O_RDONLY);
    }

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
        dup2(gz_read_pipe1[0], STDIN_FILENO);
        close(gz_read_pipe1[0]);
        if(kInput1IsHDFS)
        {
            close(gz_read_pipe1[1]);
        }
        char* args[3] = {gzip, decompress, nullptr};
        execvp(gzip, args);
        clog<<"ERROR: Cannot exec gzip: "<<strerror(errno)<<endl;
        exit(2);
    }
    else
    {
        children.push_back(input_pid1);
        close(input_pipe1[1]);
        if(!kInput1IsHDFS)
        {
            close(gz_read_pipe1[0]);
        }
    }

    FILE* input1 = fdopen(input_pipe1[0], "r");
    FILE* input2 = nullptr;
    vector<tuple<char*, size_t> >* output1;
    vector<tuple<char*, size_t> >* output2;

    if(kIsPaired)
    {
        if(kInput2IsHDFS)
        {
            // connect to HDFS if not already connected
            input_hdfs2_nodename = kInputFastq2.substr(0, kInputFastq2.find('/', kHDFSProto.size()));
            if(hdfs_map.count(input_hdfs2_nodename)==0)
            {
                if(kVerboseFlag)
                {
                    clog<<"INFO: Connect to "<<input_hdfs2_nodename<<endl;
                }
                hdfsBuilder* hdfs_builder = hdfsNewBuilder();   // freed by hdfsBuilderConnect
                hdfsBuilderSetNameNode(hdfs_builder, input_hdfs2_nodename.c_str());
                input_hdfs2 = hdfsBuilderConnect(hdfs_builder);
                if(input_hdfs2==nullptr)
                {
                    clog<<"ERROR: Cannot connect to "<<input_hdfs2_nodename<<endl;
                    exit(2);
                }
                hdfs_map[input_hdfs2_nodename] = input_hdfs2;
            }
            else
            {
                input_hdfs2 = hdfs_map[input_hdfs2_nodename];
            }

            // open HDFS file
            if(pipe(gz_read_pipe2)!=0)
            {
                clog<<"ERROR: Cannot make pipe: "<<strerror(errno)<<endl;
                exit(2);
            }
            hdfsFile hdfs_file2 = hdfsOpenFile(input_hdfs2, kInputFastq2.c_str(), O_RDONLY,
                kHDFSBufferSize, kHDFSReplication, kHDFSBlockSize);
            if(hdfs_file2==nullptr)
            {
                clog<<"ERROR: Cannot open file on HDFS"<<endl;
                exit(2);
            }
            threads.push_back(thread(ReadHDFS, input_hdfs2, hdfs_file2, gz_read_pipe2[1]));
        }
        else
        {
            gz_read_pipe2[0] = open(kInputFastq2.c_str(), O_RDONLY);
        }
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
            dup2(gz_read_pipe2[0], STDIN_FILENO);
            close(gz_read_pipe2[0]);
            if(kInput2IsHDFS)
            {
                close(gz_read_pipe2[1]);
            }
            char* args[3] = {gzip, decompress, nullptr};
            execvp(gzip, args);
            clog<<"ERROR: Cannot exec gzip: "<<strerror(errno)<<endl;
            exit(2);
        }
        else
        {
            children.push_back(input_pid2);
            close(input_pipe2[1]);
            if(!kInput2IsHDFS)
            {
                close(gz_read_pipe2[0]);
            }
        }
        input2 = fdopen(input_pipe2[0], "r");
    }

    if(kOutput1IsHDFS)
    {
        // connect to HDFS if not already connected
        output_hdfs1_nodename = kOutputFastq1.substr(0, kOutputFastq1.find('/', kHDFSProto.size()));
        output_fastq1_hdfs_path = kOutputFastq1.substr(output_hdfs1_nodename.size(), string::npos);
        if(hdfs_map.count(output_hdfs1_nodename)==0)
        {
            if(kVerboseFlag)
            {
                clog<<"INFO: Connect to "<<output_hdfs1_nodename<<endl;
            }
            hdfsBuilder* hdfs_builder = hdfsNewBuilder();   // freed by hdfsBuilderConnect
            hdfsBuilderSetNameNode(hdfs_builder, output_hdfs1_nodename.c_str());
            output_hdfs1 = hdfsBuilderConnect(hdfs_builder);
            if(output_hdfs1==nullptr)
            {
                clog<<"ERROR: Cannot connect to "<<output_hdfs1_nodename<<endl;
                exit(2);
            }
            hdfs_map[output_hdfs1_nodename] = output_hdfs1;
        }
        else
        {
            output_hdfs1 = hdfs_map[output_hdfs1_nodename];
        }
    }

    if(kOutput2IsHDFS)  // implies kIsPaired
    {
        output_hdfs2_nodename = kOutputFastq2.substr(0, kOutputFastq2.find('/', kHDFSProto.size()));
        output_fastq2_hdfs_path = kOutputFastq2.substr(output_hdfs2_nodename.size(), string::npos);
        if(hdfs_map.count(output_hdfs2_nodename)==0)
        {
            if(kVerboseFlag)
            {
                clog<<"INFO: Connect to "<<output_hdfs2_nodename<<endl;
            }
            hdfsBuilder* hdfs_builder = hdfsNewBuilder();   // freed by hdfsBuilderConnect
            hdfsBuilderSetNameNode(hdfs_builder, output_hdfs2_nodename.c_str());
            output_hdfs2 = hdfsBuilderConnect(hdfs_builder);
            if(output_hdfs2==nullptr)
            {
                clog<<"ERROR: Cannot connect to "<<output_hdfs2_nodename<<endl;
                exit(2);
            }
        }
        else
        {
            output_hdfs2 = hdfs_map[output_hdfs2_nodename];
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
            hdfsFile hdfs_file1 = hdfsOpenFile(output_hdfs1, output_filename_string.c_str(), O_WRONLY,
                kHDFSBufferSize, kHDFSReplication, kHDFSBlockSize);
            if(hdfs_file1==nullptr)
            {
                clog<<"ERROR: Cannot open file on HDFS"<<endl;
                exit(2);
            }
            threads.push_back(thread(WriteBatch, output1, line_count1, output_file1));
            threads.push_back(thread(WriteHDFS, write_pipe[0], output_hdfs1, hdfs_file1));
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
                hdfsFile hdfs_file2 = hdfsOpenFile(output_hdfs2, output_filename_string.c_str(), O_WRONLY,
                    kHDFSBufferSize, kHDFSReplication, kHDFSBlockSize);
                if(hdfs_file2==nullptr)
                {
                    clog<<"ERROR: Cannot open file on HDFS"<<endl;
                    exit(2);
                }
                threads.push_back(thread(WriteBatch, output2, line_count2, output_file2));
                threads.push_back(thread(WriteHDFS, write_pipe[0], output_hdfs2, hdfs_file2));
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

    int num_splits = batch_id;

    // delete extra files
    for(;;++batch_id)
    {
        if(kOutput1IsHDFS)
        {
            string output_filename_string = output_fastq1_hdfs_path+".part"+to_string(batch_id);
            if(hdfsExists(output_hdfs1, output_filename_string.c_str())==0)
            {
                if(hdfsDelete(output_hdfs1, output_filename_string.c_str(), 0)!=0)
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
                if(hdfsExists(output_hdfs2, output_filename_string.c_str())==0)
                {
                    if(hdfsDelete(output_hdfs2, output_filename_string.c_str(), 0)!=0)
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
    if(!kInput1IsHDFS)
    {
        close(gz_read_pipe1[0]);
    }
    if(!kInput2IsHDFS)
    {
        close(gz_read_pipe2[0]);
    }

    for(auto i: hdfs_map)
    {
        hdfsDisconnect(i.second);
    }

    return int(num_splits);
}

