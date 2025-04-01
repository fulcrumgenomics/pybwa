#include "bntseq.h"
#include "bwt.h"
#include "bwtaln.h"
#include "htslib/kstring.h"
#include "bwase.h"
#include "htslib/sam.h"
#include "libbwaaln_utils.h"

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

void bwa_cal_pac_pos_with_bwt(const bntseq_t *bns, bwt_t *bwt, bwa_seq_t *p, int max_mm, float fnr)
{
    int j, strand, n_multi;
    bwa_cal_pac_pos_core(bns, bwt, p, max_mm, fnr);
    for (j = n_multi = 0; j < p->n_multi; ++j) {
        bwt_multi1_t *q = p->multi + j;
        q->pos = bwa_sa2pos(bns, bwt, q->pos, p->len + q->ref_shift, &strand);
        q->strand = strand;
        if (q->pos != p->pos && q->pos != (bwtint_t)-1)
            p->multi[n_multi++] = *q;
    }
    p->n_multi = n_multi;
}

// Copy as-is from bwtaln.c
#ifdef HAVE_PTHREAD
typedef struct {
    int tid;
    const bntseq_t *bns;
    bwt_t *bwt;
    uint8_t *pac;
    sam_hdr_t *h;
    int n_seqs;
    bwa_seq_t *seqs;
    const gap_opt_t *opt;
    int max_hits;
    int with_md;
    bam1_t **bams;
} thread_aux_t;

static int
 _bwa_aln_core(const bntseq_t *bns, bwt_t *bwt, uint8_t *pac, sam_hdr_t *h, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, int max_hits, int with_md, bam1_t **bams, int tid)
{
    int i, j, i_increment;
    kstring_t *kstr = (kstring_t*)calloc(1, sizeof(kstring_t));

    // NB: this _complements_ the read
    bwa_cal_sa_reg_gap(tid, bwt, n_seqs, seqs, opt, 0);

    i = (tid == -1) ? 0 : tid;
    i_increment = (tid == -1) ? 1 : opt->n_threads;
    for (; i < n_seqs; i += i_increment) {
        // undo the complement done by bwa_cal_sa_reg_gap
        for (j = 0; j < seqs[i].full_len; j++) {
            seqs[i].seq[j] = (seqs[i].seq[j] > 3) ? seqs[i].seq[j] : (3 - seqs[i].seq[j]);
        }

        // Find the hits
        bwa_aln2seq_core(seqs[i].n_aln, seqs[i].aln, &seqs[i], 1, max_hits);

        // calculate the genomic position given the suffix array offsite
        bwa_cal_pac_pos_with_bwt(bns, bwt, &seqs[i], opt->max_diff, opt->fnr);

        // refine gapped alignment
        bwa_refine_gapped(bns, 1, &seqs[i], pac, with_md);

        // create the htslib record
        bams[i] = bwa_print_sam1(bns, &seqs[i], NULL, opt->mode, opt->max_top2, kstr, h);
    }

    free(kstr->s);
    free(kstr);
    return 0;
}

static void *worker(void *data)
{
    thread_aux_t *d = (thread_aux_t*)data;
    _bwa_aln_core(d->bns, d->bwt, d->pac, d->h, d->n_seqs, d->seqs, d->opt, d->max_hits, d->with_md, d->bams, d->tid);
    return 0;
}
#endif

bam1_t **bwa_aln_and_samse(const bntseq_t *bns, bwt_t *const bwt, uint8_t *pac, sam_hdr_t *h, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, int max_hits, int with_md)
{
    bam1_t **bams = (bam1_t**)calloc(n_seqs, sizeof(bam1_t*));
#ifdef HAVE_PTHREAD
    if (opt->n_threads <= 1) { // no multi-threading at all
        _bwa_aln_core(bns, bwt, pac, h, n_seqs, seqs, opt, max_hits, with_md, bams, -1);
    } else {
        pthread_t *tid;
        pthread_attr_t attr;
        thread_aux_t *data;
        int j;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
        tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
        for (j = 0; j < opt->n_threads; ++j) {
            data[j].tid = j; data[j].bns = bns; data[j].bwt = bwt;
            data[j].pac = pac; data[j].h = h;
            data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = opt;
            data[j].max_hits = max_hits; data[j].with_md = with_md;
            data[j].bams = bams;
            pthread_create(&tid[j], &attr, worker, data + j);
        }
        for (j = 0; j < opt->n_threads; ++j)  {
            pthread_join(tid[j], 0);
        }
        free(data); free(tid);
    }
#else
    _bwa_aln_core(bns, bwt, pac, h, n_seqs, seqs, opt, max_hits, with_md, bams, -1);
#endif
    return bams;
}