#include "bntseq.h"
#include "bwt.h"
#include "bwamem.h"
#include "kstring.h"
#include "kvec.h"

typedef kvec_t(mem_aln_t) mem_alns_t;

mem_alns_t* mem_process_seqs_alt(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0);
