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

#include <errno.h>
#include <unistd.h>

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

/** Finds the index of the given delimiter in the given string with a given length.  If
 *  found, the delimiter is replaced with the NULL terminator.  Returns the index of the
 *  delimiter, otherwise sets the errno to EINVAL and returns -1 (i.e. cannot be found).
 */
static inline
int32_t _next_index(char *value, int32_t start, int32_t len, char delim) {
    errno = 0;
    int32_t end = start;
    while (end < len && value[end] != '\0' && value[end] != delim) {
        end++;
    }
    if (end < len) {
        value[end] = '\0';
        return end;
    }
    errno = EINVAL;
    return -1;
}

/** Little helper function to copy the source string to the destination, after allocation. */
static inline
char* _malloc_and_strncopy(char *dest, char *src, size_t len) {
    dest = (char*)malloc(sizeof(char)*(len+1));
    if (dest == NULL) {
        errno = ENOMEM;
        return NULL;
    }
    dest = strncpy(dest, src, len);
    dest[len] = '\0';
    return dest;
}

/** Little helper function to parse an integer from a string, with some error checking */
static inline
long int _strtol_helper(const char *value) {
    errno = 0;
    char *endptr = NULL;
    long int retval = strtol(value, &endptr, 0);
    if (endptr == value) { // no digits found
        errno = EINVAL;
    }
    return retval;
}

/** Parse the XA tag e.g. XA:Z:chr4,-97592047,24M,3 */
bwa_hit_t* parse_xa_entry(char *value, int32_t len, bwa_hit_t *hit) {
    int32_t start = 0, end = 0;
    errno = 0;

    hit->s = value;
    hit->s_len = len;
    hit->refname = NULL;
    hit->cigar = NULL;
    hit->md = NULL;
    hit->md_len = 0;
    hit->rest = NULL;
    hit->rest_len = 0;

    if (len == 0) {
        errno = EINVAL;
        return NULL;
    }

    value[len] = '\0';

    // refname
    end = _next_index(value, start, len, ',');
    if (errno != 0) goto fail_at_refname;
    hit->refname_len = end - start;
    hit->refname = _malloc_and_strncopy(hit->refname, value + start, hit->refname_len);
    if (hit->refname == NULL) goto fail_at_refname;
    value[end] = ',';

    // strand and pos
    start = end + 1;
    end = _next_index(value, start, len, ',');
    if (errno != 0) goto fail_at_strand;
    hit->negative = value[start];
    hit->pos = _strtol_helper(value + start + 1);
    if (errno != 0) goto fail_at_pos;
    value[end] = ',';

    // strand and cigar
    start = end + 1;
    end = _next_index(value, start, len, ',');
    if (errno != 0) goto fail_at_strand;
    size_t cigar_mem = 0;
    ssize_t n_cigar = 0;
    if ((n_cigar = sam_parse_cigar(value + start , NULL, &hit->cigar, &cigar_mem)) < 0) {
        errno = EINVAL;
        goto fail_at_cigar;
    }
    hit->n_cigar = n_cigar;
    value[end] = ',';

    // edits
    start = end + 1;
    end = _next_index(value, start, len, ',');
    if (errno == 0) value[end] = ',';
    else { // read to the end of the string, which is OK
        end = len;
        errno = 0;
    }
    hit->edits = _strtol_helper(value + start);
    if (errno != 0) goto fail_at_edits;

    // md and the rest
    if (end != len) {
        start = end + 1;
        end = _next_index(value, start, len, ',');
        if (errno == 0) value[end] = ',';
        else { // read to the end of the string, which is OK
            end = len;
            errno = 0;
        }
        hit->md_len = end - start;
        hit->md = _malloc_and_strncopy(hit->md, value + start, hit->md_len);
        if (hit->md == NULL) goto fail_at_md;
        if (end < len) {
            start = end + 1;
            hit->rest_len = len - start;
            hit->rest = _malloc_and_strncopy(hit->rest, value + start, hit->rest_len);
            if (hit->rest == NULL) goto fail_at_rest;
        }
    }

    value[len] = ';';

    return hit;

fail_at_rest:
    free(hit->md);
    hit->md = NULL;
    hit->md_len = 0;
fail_at_md:
fail_at_edits:
    free(hit->cigar);
    hit->cigar = NULL;
    hit->n_cigar = 0;
fail_at_cigar:
fail_at_pos:
fail_at_strand:
    free(hit->refname);
    hit->refname = NULL;
    hit->refname_len = 0;
fail_at_refname:
    return NULL;
}

bwa_hits_t* parse_xa(char *value) {
    errno = 0;
    if (value == NULL) {
        errno = EINVAL;
        return NULL;
    }
    int32_t end = 0;
    uint32_t num_hits = 0;

    while (value[end] != '\0') {
        if (end > 0 && value[end] == ';') num_hits++;
        end++;
    }
    if (num_hits == 0) {
        errno = EINVAL;
        return NULL;
    }

    bwa_hits_t *hits = (bwa_hits_t*)calloc(sizeof(bwa_hits_t), 1);
    if (hits == NULL) {
        errno = ENOMEM;
        return NULL;
    }
    hits->hits = (bwa_hit_t*)calloc(sizeof(bwa_hit_t), num_hits);
    if (hits->hits == NULL) {
        errno = ENOMEM;
        return NULL;
    }

    end = 0;
    hits->n = 0;
    while (value[end] != '\0') {
        int32_t start = end;
        while (value[end] != '\0' && value[end] != ';') {
            end++;
        }
        if (NULL == parse_xa_entry(value + start, end - start, &hits->hits[hits->n])) {
            int32_t i = 0;
            while (i < hits->n) {
                free(hits->hits[i].refname);
                free(hits->hits[i].cigar);
                free(hits->hits[i].md);
                free(hits->hits[i].rest);
                i++;
            }
            free(hits->hits);
            free(hits);
            errno = EINVAL;
            return NULL;
        }
        end++; // move past the delimiter
        hits->n++;
    }
    return hits;
}
