#include "bntseq.h"
#include "bwt.h"
#include "bwtaln.h"
#include "kstring.h"


// Version of "bwa_cal_pac_pos" bwase.c that can use an already loaded forward suffix array (BWT).
void bwa_cal_pac_pos_with_bwt(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, int max_mm, float fnr, bwt_t *bwt);

// Port of running "bwa_cal_sa_reg_gap" in the "bwa_aln_core" method in bwtaln.c so we can support multi-threading.
void bwa_cal_sa_reg_gap_threaded(int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt);
