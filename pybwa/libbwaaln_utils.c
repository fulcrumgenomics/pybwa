#include "bntseq.h"
#include "bwt.h"
#include "bwtaln.h"
#include "kstring.h"
#include "bwase.h"
#include "libbwaaln_utils.h"


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