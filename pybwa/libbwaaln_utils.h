#include "bntseq.h"
#include "bwt.h"
#include "bwtaln.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"

bam1_t **bwa_aln_and_samse(const bntseq_t *bns, bwt_t *const bwt, uint8_t *pac, sam_hdr_t *h, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, int max_hits, int with_md);