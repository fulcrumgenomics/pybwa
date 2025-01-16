#include "bntseq.h"
#include "bwt.h"
#include "bwtaln.h"
#include "kstring.h"


// Version of "bwa_cal_pac_pos" bwase.c that can use an already loaded forward suffix array (BWT)
void bwa_cal_pac_pos_with_bwt(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, int max_mm, float fnr, bwt_t *bwt);
