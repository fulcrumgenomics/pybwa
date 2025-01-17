#include "bntseq.h"
#include "bwt.h"
#include "bwtaln.h"
#include "kstring.h"
#include "bwase.h"
#include "libbwaaln_utils.h"

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

void bwa_cal_pac_pos_with_bwt(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, int max_mm, float fnr, bwt_t *bwt)
{
    int i, j, strand, n_multi;
    char str[1024];
    for (i = 0; i != n_seqs; ++i) {
        bwa_seq_t *p = &seqs[i];
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
}

// Copy as-is from bwtaln.c
#ifdef HAVE_PTHREAD
typedef struct {
    int tid;
    bwt_t *bwt;
    int n_seqs;
    bwa_seq_t *seqs;
    const gap_opt_t *opt;
} thread_aux_t;

static void *worker(void *data)
{
    thread_aux_t *d = (thread_aux_t*)data;
    bwa_cal_sa_reg_gap(d->tid, d->bwt, d->n_seqs, d->seqs, d->opt);
    return 0;
}
#endif

void bwa_cal_sa_reg_gap_threaded(int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt)
{
#ifdef HAVE_PTHREAD
    if (opt->n_threads <= 1) { // no multi-threading at all
        bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
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
            data[j].tid = j; data[j].bwt = bwt;
            data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = opt;
            pthread_create(&tid[j], &attr, worker, data + j);
        }
        for (j = 0; j < opt->n_threads; ++j)  {
            pthread_join(tid[j], 0);
        }
        free(data); free(tid);
    }
#else
    bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
#endif
}