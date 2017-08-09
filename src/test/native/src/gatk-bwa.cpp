#include<iostream>
#include<string>

#include<bwa_wrapper.h>

using namespace std;

void AlignSeqs(bseq1_t* seqs, int read_count, uint64_t start_idx, const bool kIsPaired);

extern ktp_aux_t* aux;

int main(int argc, char* argv[])
{
    aux = new ktp_aux_t;
    memset(aux, 0, sizeof(ktp_aux_t));
    pre_process(argc, argv, aux, false);
    aux->out = sam_open("/dev/stdout", "w");
    const bool kIsPaired = aux->opt->flag&MEM_F_PE;
    int read_count = 0;
    uint64_t start_idx = 0;
    bool all_right = true;
    for(;;)
    {
        bseq1_t* seqs = bseq_read(aux->actual_chunk_size, &read_count, aux->ks, aux->ks2);
        if(seqs==nullptr)
        {
            break;
        }
        clog<<"gatkbwa:INFO Get "<<read_count<<" sequences\n";

        AlignSeqs(seqs, read_count, start_idx, kIsPaired);
        
        clog<<"gatkbwa:INFO Aligned "<<read_count<<" sequences\n";

        start_idx += read_count;

        for(int i = 0; i < read_count; ++i)
        {
            for(int j = 0; j < seqs[i].bams->l && all_right; ++j)
            {
                if(sam_write1(aux->out, aux->h, seqs[i].bams->bams[j])==-1)
                {
                    clog<<"gatk-bwa:ERROR Cannot write bam file\n";
                }
            }
            bams_destroy(seqs[i].bams);
        }
        free(seqs);
    }

    sam_close(aux->out);

    return 0;
}

