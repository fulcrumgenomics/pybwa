#include <limits.h>
#include "bntseq.h"
#include "bwt.h"
#include "kstring.h"
#include "kvec.h"
#include "bwamem.h"
#include "libbwamem_utils.h"
#include "utils.h"

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

// Port of worker_t that adds the alns field.
typedef struct {
    const mem_opt_t *opt;
    const bwt_t *bwt;
    const bntseq_t *bns;
    const uint8_t *pac;
    const mem_pestat_t *pes;
    smem_aux_t **aux;
    bseq1_t *seqs;
    mem_alnreg_v *regs;
    int64_t n_processed;
    mem_alns_t* alns;
} worker_alt_t;

// Port of mem_reg2sam from "bwamem.c", where we do not print the SAM output, but instead return the vector of alignments
mem_alns_t mem_reg2sam_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m)
{
    extern char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, mem_alnreg_v *a, int l_query, const char *query);
    kstring_t str;
    mem_alns_t aa;
    int k, l;
    char **XA = 0;

    if (!(opt->flag & MEM_F_ALL))
        XA = mem_gen_alt(opt, bns, pac, a, s->l_seq, s->seq);
    kv_init(aa);
    str.l = str.m = 0; str.s = 0;
    for (k = l = 0; k < a->n; ++k) {
        mem_alnreg_t *p = &a->a[k];
        mem_aln_t *q;
        if (p->score < opt->T) continue;
        if (p->secondary >= 0 && (p->is_alt || !(opt->flag&MEM_F_ALL))) continue;
        if (p->secondary >= 0 && p->secondary < INT_MAX && p->score < a->a[p->secondary].score * opt->drop_ratio) continue;
        q = kv_pushp(mem_aln_t, aa);
        *q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, p);
        assert(q->rid >= 0); // this should not happen with the new code
        q->XA = XA? XA[k] : 0;
        q->flag |= extra_flag; // flag secondary
        if (p->secondary >= 0) q->sub = -1; // don't output sub-optimal score
        if (l && p->secondary < 0) // if supplementary
            q->flag |= (opt->flag&MEM_F_NO_MULTI)? 0x10000 : 0x800;
        if (!(opt->flag & MEM_F_KEEP_SUPP_MAPQ) && l && !p->is_alt && q->mapq > aa.a[0].mapq)
            q->mapq = aa.a[0].mapq; // lower mapq for supplementary mappings, unless -5 or -q is applied
        ++l;
    }
    if (aa.n == 0) { // no alignments good enough; then create an unaligned record
        mem_aln_t *q;
        q = kv_pushp(mem_aln_t, aa);
        *q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, 0);
        q->flag |= extra_flag;
    }
    return aa;
}

// Port of worker2 from "bwamem.c" so that we can return (and not print) the alignment
static void worker2_alt(void *data, int i, int tid)
{
    extern int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2]);
    extern void mem_reg2ovlp(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a);
    worker_alt_t *w = (worker_alt_t*)data;
    // NB: paired end not supported
    if (bwa_verbose >= 4) printf("=====> Finalizing read '%s' <=====\n", w->seqs[i].name);
    mem_mark_primary_se(w->opt, w->regs[i].n, w->regs[i].a, w->n_processed + i);
    if (w->opt->flag & MEM_F_PRIMARY5) mem_reorder_primary5(w->opt->T, &w->regs[i]);
    w->alns[i] = mem_reg2sam_alt(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
    free(w->regs[i].a);
}

// Port of "mem_process_seqs" from "bwamem.c" so that we can align multi-threaded and return the vector of alignments
// (one vector of alignments per input sequence
mem_alns_t* mem_process_seqs_alt(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac,
                      int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0)
{
    extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
    worker_t w;
    worker_alt_t w_alt;
    mem_pestat_t pes[4];
    double ctime, rtime;
    int i;

    ctime = cputime(); rtime = realtime();
    w.regs = malloc(n * sizeof(mem_alnreg_v));
    w.opt = opt; w.bwt = bwt; w.bns = bns; w.pac = pac;
    w.seqs = seqs; w.n_processed = n_processed;
    w.pes = &pes[0];
    w.aux = malloc(opt->n_threads * sizeof(smem_aux_t));
    for (i = 0; i < opt->n_threads; ++i)
        w.aux[i] = smem_aux_init();
    kt_for(opt->n_threads, worker1, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // find mapping positions
    for (i = 0; i < opt->n_threads; ++i)
        smem_aux_destroy(w.aux[i]);
    free(w.aux);
    // Paired-end not supported
//    if (opt->flag&MEM_F_PE) { // infer insert sizes if not provided
//        if (pes0) memcpy(pes, pes0, 4 * sizeof(mem_pestat_t)); // if pes0 != NULL, set the insert-size distribution as pes0
//        else mem_pestat(opt, bns->l_pac, n, w.regs, pes); // otherwise, infer the insert size distribution from data
//    }
    w_alt.alns = malloc(n * sizeof(mem_alns_t));
    // copy worker_t to worker_alt_t
    w_alt.opt = w.opt; w_alt.bwt = bwt; w_alt.bns = w.bns; w_alt.pac = w.pac;
    w_alt.pes = w.pes; w_alt.aux = w.aux; w_alt.seqs = w.seqs; w_alt.regs = w.regs; w_alt.n_processed = w.n_processed;
    kt_for(opt->n_threads, worker2_alt, &w_alt, (opt->flag&MEM_F_PE)? n>>1 : n); // generate alignment
    free(w.regs);
    if (bwa_verbose >= 3)
        fprintf(stderr, "[M::%s] Processed %d reads in %.3f CPU sec, %.3f real sec\n", __func__, n, cputime() - ctime, realtime() - rtime);

    return w_alt.alns;
}
