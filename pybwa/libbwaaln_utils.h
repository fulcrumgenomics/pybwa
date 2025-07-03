#include "bntseq.h"
#include "bwt.h"
#include "bwtaln.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"

typedef struct {
    char *s;
    uint32_t s_len;
    char *refname;
    uint32_t refname_len;
    uint32_t pos;
    char negative;
    uint32_t *cigar;
    uint32_t n_cigar;
    uint32_t edits;
    char *md;
    uint32_t md_len;
    char *rest;
    uint32_t rest_len;
} bwa_hit_t;

typedef struct {
    uint32_t n;
    bwa_hit_t *hits;
} bwa_hits_t;

bam1_t **bwa_aln_and_samse(const bntseq_t *bns, bwt_t *const bwt, uint8_t *pac, sam_hdr_t *h, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, int max_hits, int with_md);
bwa_hit_t* parse_xa_entry(char *value, int32_t len, bwa_hit_t *hit);
bwa_hits_t* parse_xa(char *value, uint32_t max_hits);