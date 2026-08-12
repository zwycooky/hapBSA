/* hapbsa.c -- C implementation of the hapBSA V4 analysis pipeline.
 *
 * This is a from-scratch port of the Perl pipeline:
 *   - hapBSA_V4.pl                                 (main driver)
 *   - separating_reads_by_haplotype.binarySearch.hapBSA.block.pl
 *                                                   (read phasing, merged as a subroutine)
 *   - snpMapper_sub.pl                              (SNP-index, merged as a subroutine)
 *
 * The two former sub-scripts are now plain C functions:
 *   phase_pool()          -> phasing_reads()       (per-pool read phasing)
 *   snp_mapper_window()   -> compute_snp_index()   (per-window SNP-index + thresholds)
 *
 * BAM access uses htslib. The SNP-index pileup engine replicates the
 * semantics of:  samtools mpileup -q 45 -Q 20 -ABf <ref> <bam1> <bam2>
 * including default flag filtering (UNMAP|SECONDARY|QCFAIL|DUP) and smart
 * read-pair overlap removal (bam_mplp_init_overlaps).
 *
 * BAQ is OFF by default (HAPBSA_BAQ=1 enables it, applied in parallel with
 * sam_prob_realn flag 3): on this test data it costs ~5x of the runtime
 * while changing the window SNP-index by ~0.01 on average.
 *
 * Command line interface is compatible with hapBSA_V4.pl:
 *   -1 bam1 -2 bam2 -p hap.txt -r ref.fa -e sep.pl -m snpMapper.pl
 *   -t tmpdir -o prefix [-w 1000000] [-s 600000] [-d 10] [-D 50]
 *   [-N 12] [-n 3] [-a 10]
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <pthread.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>

/* ------------------------------------------------------------------ */
/* Configuration                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *bam1;
    const char *bam2;
    const char *hap_file;
    const char *refgenome;
    const char *tmp_dir;
    const char *outprefix;
    const char *sep_script;    /* accepted for compatibility, merged into this binary */
    const char *snp_mapper;    /* accepted for compatibility, merged into this binary */
    int cpu;
    int64_t window_len;
    int64_t step_len;
    int64_t depth;             /* -d */
    int64_t block_read_nums;   /* -N */
    int64_t snps_for_reads;    /* -n */
    int64_t read_nums_window;  /* -D */
    int has_block;             /* hap file carries a 5th (block) column */
} Config;

/* ------------------------------------------------------------------ */
/* Hap file data                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    int64_t pos;       /* 1-based position */
    char hap1;
    char hap2;
    int64_t block;     /* block id */
    int has_block;     /* this site has a block id */
} HapSite;

typedef struct {
    char *chr;
    HapSite *sites;
    size_t n, cap;
    int64_t chr_end;   /* last position seen in the hap file (Perl: $chr_len{$chr} = $pos) */
} ChrHap;

typedef struct {
    const ChrHap *chr;
    int64_t s, e;      /* window start/end, Perl semantics (region lg1:s-e) */
    int win_idx;
} Window;

/* ------------------------------------------------------------------ */
/* Small utilities                                                     */
/* ------------------------------------------------------------------ */

static void die(const char *msg)
{
    fprintf(stderr, "ERROR: %s\n", msg);
    exit(EXIT_FAILURE);
}

/* sprintf("%.4f") semantics (same C library as Perl's sprintf) */
static double round4(double x)
{
    char buf[64];
    snprintf(buf, sizeof buf, "%.4f", x);
    return strtod(buf, NULL);
}

/* sprintf("%.2f") semantics, returned as double */
static double round2(double x)
{
    char buf[64];
    snprintf(buf, sizeof buf, "%.2f", x);
    return strtod(buf, NULL);
}

/* Format a 4-decimal value into buf exactly like Perl print/sprintf */
static void fmt4(char *buf, size_t sz, double x)
{
    snprintf(buf, sz, "%.4f", x);
}

/* ------------------------------------------------------------------ */
/* Command line parsing (same options as hapBSA_V4.pl)                 */
/* ------------------------------------------------------------------ */

static void usage(FILE *fp)
{
    fprintf(fp,
        "Usage:\n"
        "  hapbsa -OPTIONS VALUES\n"
        "\n"
        "options:\n"
        "--input options\n"
        "  -1 FILE         bam file of pool1\n"
        "  -2 FILE         bam file of pool2\n"
        "  -p FILE         hap file\n"
        "  -r REFGENOME    reference genome for mapping\n"
        "\n"
        "--sub program (accepted for compatibility; functionality is built in)\n"
        "  -e SCRIPT       optional: path of separating_reads_by_haplotype...pl (ignored)\n"
        "  -m SCRIPT       optional: path of snpMapper_sub.pl (ignored)\n"
        "\n"
        "--output options\n"
        "  -t DIR          tmp_dir\n"
        "  -o PREFIX       prefix of output file\n"
        "\n"
        "--criteria options\n"
        "  -w INT          window length (bp) [default: 1000000]\n"
        "  -s INT          window step (bp)   [default: 600000]\n"
        "  -d INT          minimum depth for calculating SNP-index [default: 10]\n"
        "  -D INT          minimum phased read numbers of a window [default: 50]\n"
        "  -N INT          minimum read numbers of hap block [default: 12]\n"
        "  -n INT          minimum numbers of phased snps for a read [default: 3]\n"
        "\n"
        "--performance options\n"
        "  -a INT          cpus cores used for the analysis [default: 10]\n"
        "\n");
}

static void parse_args(int argc, char **argv, Config *cfg)
{
    int c;
    memset(cfg, 0, sizeof *cfg);
    cfg->cpu = 10;
    cfg->window_len = 1000000;
    cfg->step_len = 600000;
    cfg->depth = 10;
    cfg->block_read_nums = 12;
    cfg->snps_for_reads = 3;
    cfg->read_nums_window = 50;

    while ((c = getopt(argc, argv, "1:2:p:r:a:t:o:w:s:d:e:m:N:n:D:")) != -1) {
        switch (c) {
        case '1': cfg->bam1 = optarg; break;
        case '2': cfg->bam2 = optarg; break;
        case 'p': cfg->hap_file = optarg; break;
        case 'r': cfg->refgenome = optarg; break;
        case 'e': cfg->sep_script = optarg; break;
        case 'm': cfg->snp_mapper = optarg; break;
        case 't': cfg->tmp_dir = optarg; break;
        case 'o': cfg->outprefix = optarg; break;
        case 'a': cfg->cpu = atoi(optarg); break;
        case 'w': cfg->window_len = atoll(optarg); break;
        case 's': cfg->step_len = atoll(optarg); break;
        case 'd': cfg->depth = atoll(optarg); break;
        case 'D': cfg->read_nums_window = atoll(optarg); break;
        case 'N': cfg->block_read_nums = atoll(optarg); break;
        case 'n': cfg->snps_for_reads = atoll(optarg); break;
        default: usage(stderr); exit(EXIT_FAILURE);
        }
    }
    if (!cfg->bam1 || !cfg->bam2 || !cfg->hap_file || !cfg->refgenome ||
        !cfg->tmp_dir || !cfg->outprefix) {
        usage(stderr);
        exit(EXIT_FAILURE);
    }
    if (cfg->cpu < 1) cfg->cpu = 1;
}

/* ------------------------------------------------------------------ */
/* Hap file reading                                                    */
/* ------------------------------------------------------------------ */

static ChrHap *find_or_add_chr(ChrHap ***chrs, size_t *n_chrs, size_t *cap, const char *name)
{
    size_t i;
    for (i = 0; i < *n_chrs; i++) {
        if (strcmp((*chrs)[i]->chr, name) == 0) return (*chrs)[i];
    }
    if (*n_chrs == *cap) {
        *cap = *cap ? *cap * 2 : 8;
        *chrs = realloc(*chrs, *cap * sizeof(**chrs));
        if (!*chrs) die("out of memory");
    }
    ChrHap *c = calloc(1, sizeof *c);
    c->chr = strdup(name);
    (*chrs)[*n_chrs] = c;
    (*n_chrs)++;
    return c;
}

static void append_site(ChrHap *c, int64_t pos, char hap1, char hap2, int64_t block, int has_block)
{
    if (c->n == c->cap) {
        c->cap = c->cap ? c->cap * 2 : 1024;
        c->sites = realloc(c->sites, c->cap * sizeof(*c->sites));
        if (!c->sites) die("out of memory");
    }
    c->sites[c->n].pos = pos;
    c->sites[c->n].hap1 = hap1;
    c->sites[c->n].hap2 = hap2;
    c->sites[c->n].block = block;
    c->sites[c->n].has_block = has_block;
    c->n++;
    c->chr_end = pos;   /* Perl: $chr_len{$chr} = $pos  (last line wins) */
}

static void read_hap_file(Config *cfg, ChrHap ***chrs, size_t *n_chrs)
{
    FILE *fp = fopen(cfg->hap_file, "r");
    if (!fp) { fprintf(stderr, "ERROR: cannot open hap file %s: %s\n", cfg->hap_file, strerror(errno)); exit(EXIT_FAILURE); }

    char line[4096];
    size_t cap = 0;
    while (fgets(line, sizeof line, fp)) {
        char *toks[8];
        int nt = 0;
        char *save = NULL;
        char *p = strtok_r(line, " \t\r\n", &save);
        while (p && nt < 8) { toks[nt++] = p; p = strtok_r(NULL, " \t\r\n", &save); }
        if (nt < 4) continue;
        int64_t pos = atoll(toks[1]);
        char hap1 = toks[2][0];
        char hap2 = toks[3][0];
        int has_block = (nt >= 5);
        int64_t block = has_block ? atoll(toks[4]) : 0;
        if (has_block) cfg->has_block = 1;
        ChrHap *c = find_or_add_chr(chrs, n_chrs, &cap, toks[0]);
        append_site(c, pos, hap1, hap2, block, has_block);
    }
    fclose(fp);
}

/* ------------------------------------------------------------------ */
/* Sliding windows (identical generation to hapBSA_V4.pl)              */
/* ------------------------------------------------------------------ */

static void build_windows(ChrHap **chrs, size_t n_chrs, const Config *cfg, Window **wins, size_t *n_wins)
{
    size_t i, n = 0;
    for (i = 0; i < n_chrs; i++) {
        int64_t ws = 0, we = cfg->window_len;
        while (we < chrs[i]->chr_end) {
            n++;
            ws += cfg->step_len;
            we += cfg->step_len;
        }
        n++; /* final window (may extend past chr_end, exactly as the Perl version) */
    }
    *wins = calloc(n ? n : 1, sizeof(**wins));
    if (!*wins) die("out of memory");
    n = 0;
    for (i = 0; i < n_chrs; i++) {
        int64_t ws = 0, we = cfg->window_len;
        while (we < chrs[i]->chr_end) {
            (*wins)[n].chr = chrs[i];
            (*wins)[n].s = ws;
            (*wins)[n].e = we;
            (*wins)[n].win_idx = (int)n;
            n++;
            ws += cfg->step_len;
            we += cfg->step_len;
        }
        (*wins)[n].chr = chrs[i];
        (*wins)[n].s = ws;
        (*wins)[n].e = we;
        (*wins)[n].win_idx = (int)n;
        n++;
    }
    *n_wins = n;
}

/* ------------------------------------------------------------------ */
/* Read phasing: port of separating_reads_by_haplotype.binarySearch.   */
/* hapBSA.block.pl                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    char *qname;
    int32_t pos0;       /* 0-based start */
    int32_t end0;       /* 0-based exclusive end (== Perl 1-based inclusive end) */
    uint32_t *cigar;
    int n_cigar;
    char *seq;          /* decoded bases, NUL terminated */
    int seq_len;
    int idx;            /* original fetch order */
} Aln;

static int aln_cmp(const void *a, const void *b)
{
    const Aln *x = a, *y = b;
    int r = strcmp(x->qname, y->qname);
    if (r) return r;
    return (x->idx > y->idx) - (x->idx < y->idx);
}

/* Fetch reads overlapping [beg,end) (0-based half-open) with MAPQ >= min_mapq.
 * The Perl pipeline extracts the window with `samtools view -bh -q 45` and then
 * phases with `samtools view -q 45`; no flag filtering is applied there. */
static int fetch_reads(samFile *fp, hts_idx_t *idx, int tid,
                       int64_t beg, int64_t end, int min_mapq,
                       Aln **out, size_t *out_n)
{
    hts_itr_t *itr = sam_itr_queryi(idx, tid, beg, end);
    if (!itr) return -1;
    bam1_t *b = bam_init1();
    Aln *arr = NULL;
    size_t n = 0, cap = 0;
    int fetch_idx = 0;
    int ret;
    while ((ret = sam_itr_next(fp, itr, b)) >= 0) {
        if (b->core.qual < min_mapq) continue;
        if (n == cap) {
            cap = cap ? cap * 2 : 4096;
            arr = realloc(arr, cap * sizeof(*arr));
            if (!arr) { bam_destroy1(b); hts_itr_destroy(itr); return -1; }
        }
        Aln *a = &arr[n];
        memset(a, 0, sizeof *a);
        a->qname = strdup(bam_get_qname(b));
        a->pos0 = b->core.pos;
        a->end0 = bam_endpos(b);
        a->n_cigar = b->core.n_cigar;
        a->cigar = malloc((a->n_cigar ? a->n_cigar : 1) * sizeof(uint32_t));
        if (a->n_cigar) memcpy(a->cigar, bam_get_cigar(b), a->n_cigar * sizeof(uint32_t));
        a->seq_len = b->core.l_qseq;
        a->seq = malloc((a->seq_len + 1) ? (a->seq_len + 1) : 1);
        const uint8_t *s = bam_get_seq(b);
        int i;
        for (i = 0; i < a->seq_len; i++) a->seq[i] = seq_nt16_str[bam_seqi(s, i)];
        a->seq[a->seq_len] = 0;
        a->idx = fetch_idx++;
        n++;
    }
    bam_destroy1(b);
    hts_itr_destroy(itr);
    *out = arr;
    *out_n = n;
    return ret < 0 ? -1 : 0;
}

static void free_alns(Aln *arr, size_t n)
{
    size_t i;
    for (i = 0; i < n; i++) {
        free(arr[i].qname);
        free(arr[i].cigar);
        free(arr[i].seq);
    }
    free(arr);
}

/* Perl binarySearch(): returns insertion point (lower bound) */
static int binary_search_pos(int64_t pos, const HapSite *sites, int n)
{
    int low = 0, high = n - 1;
    while (low <= high) {
        int mid = (low + high) / 2;
        if (sites[mid].pos == pos) return mid;
        else if (sites[mid].pos < pos) low = mid + 1;
        else high = mid - 1;
    }
    return low;
}

/* Map a virtual slice index to a site, replicating Perl slice quirks:
 *  - @hap[(-1)..k] includes @hap[-1], which is the LAST element;
 *  - indices beyond the array are NOT clamped in Perl: they are undef and
 *    behave as an empty pseudo-site (pos=0, no block) during block counting. */
static HapSite empty_site = {0, 0, 0, 0, 0};

static const HapSite *virt_site(const HapSite *sites, int n, int vi)
{
    if (vi < 0) return &sites[n - 1];
    if (vi >= n) return &empty_site;
    return &sites[vi];
}

/* get_base_in_target(): exact port. Returns:
 *    0  -> position inside a deletion/skip ("no" in Perl)
 *   -1  -> no CIGAR op matched (Perl: undef, caller then uses the LAST base)
 *   >0  -> 1-based position in the read sequence
 */
static int get_base_in_target(int64_t pos, const uint32_t *cigar, int n_cigar)
{
    int64_t base_num = 0, contig = 0;
    int i;
    for (i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int64_t len = bam_cigar_oplen(cigar[i]);
        switch (op) {
        case BAM_CSOFT_CLIP:    /* S */
            base_num += len;
            break;
        case BAM_CHARD_CLIP:    /* H */
            break;
        case BAM_CINS:          /* I */
            base_num += len;
            break;
        case BAM_CDEL:          /* D */
        case BAM_CREF_SKIP:     /* N */
            if (pos > contig && pos < contig + len) return 0;
            contig += len;
            break;
        case BAM_CMATCH:        /* M */
        case BAM_CEQUAL:        /* = */
        case BAM_CDIFF:         /* X */
            if (pos > contig && pos < contig + len) {
                int64_t add = pos - contig;
                return (int)(base_num + add);
            }
            base_num += len;
            contig += len;
            break;
        case BAM_CPAD:          /* P (Perl would die; treated as no-op) */
        default:
            break;
        }
    }
    return -1;
}

/* per-alignment block selection, replicating the Perl code inside phasing_reads */
typedef struct {
    int64_t block;
    int has_block;
    int count;
} AlnBlock;

static const char *aln_block_key(int64_t block, int has_block, char *buf, size_t sz)
{
    if (!has_block) return "";
    snprintf(buf, sz, "%" PRId64, block);
    return buf;
}

/* Returns the max-count block key for this alignment (NULL when the slice has
 * no sites).  *qualifies is set when that count reaches snps_for_reads.
 * Replicates Perl: $max_block is assigned the max-count key even when the
 * count is below the threshold; only @tmp_hap is emptied in that case. */
static const char *select_block(const HapSite *sites, int n_sites,
                                int s_idx, int e_idx,
                                int64_t read_end,
                                int64_t min_snps, int *qualifies,
                                char *keybuf, size_t keybuf_sz)
{
    AlnBlock *blocks = NULL;
    size_t n_blocks = 0, cap = 0;
    int vi, lo = s_idx - 1, hi = e_idx + 1;
    for (vi = lo; vi <= hi; vi++) {
        const HapSite *site = virt_site(sites, n_sites, vi);
        if (site->pos > read_end) break;
        int found = -1;
        size_t j;
        for (j = 0; j < n_blocks; j++) {
            if (blocks[j].has_block == site->has_block && blocks[j].block == site->block) { found = j; break; }
        }
        if (found < 0) {
            if (n_blocks == cap) {
                cap = cap ? cap * 2 : 16;
                blocks = realloc(blocks, cap * sizeof(*blocks));
                if (!blocks) die("out of memory");
            }
            blocks[n_blocks].block = site->block;
            blocks[n_blocks].has_block = site->has_block;
            blocks[n_blocks].count = 0;
            found = n_blocks++;
        }
        blocks[found].count++;
    }

    /* Perl: foreach (sort keys %block_count) -> string sort of block keys */
    int best = -1;
    int best_count = 0;
    size_t j;
    for (j = 0; j < n_blocks; j++) {
        char bufA[64], bufB[64];
        const char *ka = aln_block_key(blocks[j].block, blocks[j].has_block, bufA, sizeof bufA);
        if (blocks[j].count > best_count) {
            best_count = blocks[j].count;
            best = j;
        } else if (blocks[j].count == best_count && best >= 0) {
            const char *kb = aln_block_key(blocks[best].block, blocks[best].has_block, bufB, sizeof bufB);
            if (strcmp(ka, kb) < 0) best = j;
        }
    }
    *qualifies = (best_count >= min_snps);
    if (best < 0) { free(blocks); return NULL; }
    snprintf(keybuf, keybuf_sz, "%" PRId64, blocks[best].block);
    free(blocks);
    return keybuf;
}

/* Phase one read group (all alignments sharing a qname), replicating
 * phasing_reads() from the Perl sub-script. Returns 1 when the read is phased,
 * and fills out_chr/out_start/out_end/out_phase/out_block. */
static int phase_read_group(const Aln *alns, int n_alns,
                            const HapSite *sites, int n_sites,
                            int has_block, int64_t min_snps,
                            char **out_chr, int64_t *out_start, int64_t *out_end,
                            int *out_phase, char **out_block)
{
    int64_t *seen = NULL;
    size_t n_seen = 0, seen_cap = 0;
    int64_t snp1_count = 0, snp2_count = 0, total_count = 0;
    int64_t start1 = 0, end1 = 0;
    const char *max_block = NULL;
    int ai;

    if (n_sites <= 1) return 0;   /* Perl: if ($#$hap_pos > 0) */

    for (ai = 0; ai < n_alns; ai++) {
        const Aln *a = &alns[ai];
        start1 = (int64_t)a->pos0 + 1;
        end1 = a->end0;

        int s_idx = binary_search_pos(start1, sites, n_sites);
        int e_idx = binary_search_pos(end1, sites, n_sites);

        const char *chosen_key = NULL;
        int use_sites = 0;
        char keybuf[64];
        if (has_block) {
            chosen_key = select_block(sites, n_sites, s_idx, e_idx, end1,
                                      min_snps, &use_sites, keybuf, sizeof keybuf);
            free((char *)max_block);
            max_block = chosen_key ? strdup(chosen_key) : NULL;
        }

        int lo = s_idx - 1, hi = e_idx + 1;
        int vi;
        for (vi = lo; vi <= hi; vi++) {
            const HapSite *site = virt_site(sites, n_sites, vi);
            if (site->pos > end1) break;
            if (has_block) {
                /* only sites belonging to the chosen block are kept */
                if (!use_sites) break; /* Perl: @tmp_hap = () */
                char sbuf[64];
                snprintf(sbuf, sizeof sbuf, "%" PRId64, site->block);
                if (!site->has_block || strcmp(sbuf, chosen_key) != 0) continue;
            }
            if (site->pos >= start1 && site->pos <= end1) {
                int already = 0;
                size_t k;
                for (k = 0; k < n_seen; k++) if (seen[k] == site->pos) { already = 1; break; }
                if (!already) {
                    if (n_seen == seen_cap) {
                        seen_cap = seen_cap ? seen_cap * 2 : 256;
                        seen = realloc(seen, seen_cap * sizeof(*seen));
                        if (!seen) die("out of memory");
                    }
                    seen[n_seen++] = site->pos;
                }

                int64_t relative_pos = site->pos - start1 + 1;
                int tb = get_base_in_target(relative_pos, a->cigar, a->n_cigar);
                char target_base;
                if (tb == 0) continue;              /* Perl: next if 'no' */
                if (tb < 0) {                        /* Perl undef: last base of read */
                    target_base = a->seq_len > 0 ? a->seq[a->seq_len - 1] : '*';
                } else {
                    target_base = (tb - 1 < a->seq_len) ? a->seq[tb - 1] : '*';
                }
                if (target_base == site->hap1) { snp1_count++; total_count++; }
                else if (target_base == site->hap2) { snp2_count++; total_count++; }
                else total_count++;
            }
        }
    }

    if (n_seen < 3) { free(seen); return 0; }              /* Perl: next if keys %phased_pos < 3 */
    if (total_count < min_snps) { free(seen); return 0; }  /* Perl: return 0 ... */
    int phase = 0;
    double acc = 0;
    if (snp1_count > snp2_count) { phase = 1; acc = round2((double)snp1_count / (double)total_count * 100.0); }
    else if (snp1_count < snp2_count) { phase = 2; acc = round2((double)snp2_count / (double)total_count * 100.0); }
    else { free(seen); return 0; }                         /* Perl tie: return 0 */
    if (acc < 90.0) { free(seen); return 0; }              /* Perl: $acc >= 90 */

    /* Last alignment's coordinates (Perl keeps overwriting $reads_chr/$reads_start/$reads_end) */
    const Aln *last = &alns[n_alns - 1];
    *out_chr = NULL;
    *out_start = (int64_t)last->pos0 + 1;
    *out_end = last->end0;
    *out_phase = phase;
    if (max_block) {
        *out_block = strdup(max_block);
        free((char *)max_block);
    } else {
        *out_block = NULL;
    }
    free(seen);
    return 1;
}

/* Block counts keyed by the exact string that the Perl main script extracts
 * with (split)[-1] from each phased-read output line. */
typedef struct {
    char **keys;
    int64_t (*counts)[4];
    size_t n, cap;
} BlockMap;

static void blockmap_add(BlockMap *m, const char *key, int which)
{
    size_t i;
    for (i = 0; i < m->n; i++) {
        if (strcmp(m->keys[i], key) == 0) { m->counts[i][which]++; return; }
    }
    if (m->n == m->cap) {
        m->cap = m->cap ? m->cap * 2 : 64;
        m->keys = realloc(m->keys, m->cap * sizeof(*m->keys));
        m->counts = realloc(m->counts, m->cap * sizeof(*m->counts));
        if (!m->keys || !m->counts) die("out of memory");
    }
    m->keys[m->n] = strdup(key);
    memset(m->counts[m->n], 0, 4 * sizeof(int64_t));
    m->counts[m->n][which] = 1;
    m->n++;
}

static void blockmap_free(BlockMap *m)
{
    size_t i;
    for (i = 0; i < m->n; i++) free(m->keys[i]);
    free(m->keys);
    free(m->counts);
    m->keys = NULL;
    m->counts = NULL;
    m->n = m->cap = 0;
}

/* Phase all reads of one pool in a window; fills block counts (which:
 * pool1-hap1=0, pool1-hap2=1, pool2-hap1=2, pool2-hap2=3) and per-haplotype
 * read totals. Equivalent to running the separating sub-script and counting
 * the resulting .hap1.txt/.hap2.txt lines. */
static int phase_pool(samFile *fp, hts_idx_t *idx, int tid,
                      int64_t beg, int64_t end, int min_mapq,
                      const HapSite *sites, int n_sites,
                      const Config *cfg, int pool, BlockMap *bm,
                      int64_t hap_totals[2])
{
    Aln *alns = NULL;
    size_t n_alns = 0;
    if (fetch_reads(fp, idx, tid, beg, end, min_mapq, &alns, &n_alns) != 0) return -1;
    if (n_alns) qsort(alns, n_alns, sizeof(Aln), aln_cmp);

    size_t i = 0;
    while (i < n_alns) {
        size_t j = i + 1;
        while (j < n_alns && strcmp(alns[i].qname, alns[j].qname) == 0) j++;
        int phase;
        int64_t rstart, rend;
        char *block = NULL;
        char *chr = NULL;
        if (phase_read_group(&alns[i], (int)(j - i), sites, n_sites, cfg->has_block,
                             cfg->snps_for_reads, &chr, &rstart, &rend, &phase, &block)) {
            int which = pool * 2 + (phase - 1);
            /* Perl main: my ($block_id) = (split)[-1] on
             * "chr\tstart\tend\tphase\tblock\n"; when block is empty the trailing
             * tab makes (split)[-1] == the phase number. */
            char keybuf[32];
            const char *key = (block && block[0]) ? block
                            : (snprintf(keybuf, sizeof keybuf, "%d", phase), keybuf);
            blockmap_add(bm, key, which);
            hap_totals[phase - 1]++;
        }
        free(block);
        free(chr);
        i = j;
    }
    free_alns(alns, n_alns);
    return 0;
}

/* Compute the hap-index for one window (block mode), identical to hapBSA_V4.pl */
static void hap_index_block(BlockMap *bm, const Config *cfg, char *out, size_t out_sz)
{
    size_t i;
    /* Perl: foreach (sort keys %{$hap_block_count}) -> string sort */
    size_t *order = malloc((bm->n ? bm->n : 1) * sizeof(*order));
    for (i = 0; i < bm->n; i++) order[i] = i;
    for (i = 1; i < bm->n; i++) {
        size_t j = i;
        while (j > 0 && strcmp(bm->keys[order[j - 1]], bm->keys[order[j]]) > 0) {
            size_t t = order[j - 1]; order[j - 1] = order[j]; order[j] = t; j--;
        }
    }

    double sum = 0;
    int64_t block_num = 0, hapA = 0, hapB = 0;
    for (i = 0; i < bm->n; i++) {
        int64_t p1h1 = bm->counts[order[i]][0];
        int64_t p1h2 = bm->counts[order[i]][1];
        int64_t p2h1 = bm->counts[order[i]][2];
        int64_t p2h2 = bm->counts[order[i]][3];
        if (p1h1 + p2h1 > cfg->block_read_nums && p1h2 + p2h2 > cfg->block_read_nums) {
            if (p1h1 >= p1h2) {
                sum += (double)p1h1 / (double)(p1h1 + p1h2) - (double)p2h1 / (double)(p2h1 + p2h2);
            } else {
                sum += (double)p1h2 / (double)(p1h1 + p1h2) - (double)p2h2 / (double)(p2h1 + p2h2);
            }
            block_num++;
            hapA += p1h1 + p2h1;
            hapB += p1h2 + p2h2;
        }
    }
    free(order);
    if (block_num > 0 && hapA > cfg->read_nums_window && hapB > cfg->read_nums_window) {
        fmt4(out, out_sz, fabs(sum / (double)block_num));
    } else {
        snprintf(out, out_sz, "NA");
    }
}

/* Non-block mode: Perl uses `wc -l` on the four hap files */
static void hap_index_plain(int64_t p1A, int64_t p1B, int64_t p2A, int64_t p2B,
                            const Config *cfg, char *out, size_t out_sz)
{
    if (p1A + p2A > cfg->read_nums_window && p1B + p2B > cfg->read_nums_window) {
        double r1 = (double)p1A / (double)(p1A + p1B);
        double r2 = (double)p2A / (double)(p2A + p2B);
        fmt4(out, out_sz, fabs(r1 - r2));
    } else {
        snprintf(out, out_sz, "NA");
    }
}

/* ------------------------------------------------------------------ */
/* SNP-index: port of snpMapper_sub.pl, driven by htslib's pileup      */
/* engine with the same settings as `samtools mpileup -q 45 -Q 20 -ABf` */
/* ------------------------------------------------------------------ */

typedef struct {
    bam1_t **reads;
    size_t n, cap, idx;
    const char *ref;      /* whole-chromosome sequence */
    hts_pos_t ref_len;
} ReadStore;

typedef struct {
    samFile *fp;
    hts_idx_t *idx;
    hts_itr_t *itr;
    int tid;
    int min_mq;
    ReadStore store;
} PrefetchSrc;

/* Pre-fetch all reads of a pool/window with the same filters as
 * `samtools mpileup -q 45 -Q 20 -ABf` (flag filter + min MAPQ). BAQ is
 * applied later in parallel. */
static int prefetch_fill(PrefetchSrc *ps)
{
    bam1_t *b = bam_init1();
    int ret;
    while ((ret = sam_itr_next(ps->fp, ps->itr, b)) >= 0) {
        if (b->core.tid < 0 || (b->core.flag & BAM_FUNMAP)) continue;
        if (b->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) continue;
        if (b->core.qual < ps->min_mq) continue;
        if (ps->store.n == ps->store.cap) {
            ps->store.cap = ps->store.cap ? ps->store.cap * 2 : 4096;
            ps->store.reads = realloc(ps->store.reads, ps->store.cap * sizeof(bam1_t *));
            if (!ps->store.reads) die("out of memory");
        }
        ps->store.reads[ps->store.n++] = bam_dup1(b);
    }
    bam_destroy1(b);
    return ret < 0 ? -1 : 0;   /* EOF is the normal end of the region */
}

typedef struct {
    ReadStore *store;
    size_t start, end;
} BaqJob;

static void *baq_worker(void *arg)
{
    BaqJob *j = arg;
    size_t i;
    for (i = j->start; i < j->end; i++)
        sam_prob_realn(j->store->reads[i], j->store->ref, j->store->ref_len, 3);
    return NULL;
}

/* Apply BAQ to both pools in parallel (BAQ_APPLY|BAQ_EXTEND, same as the
 * default `samtools mpileup`). Each read is touched by exactly one thread,
 * so the result is identical to the sequential order. */
static void baq_parallel(PrefetchSrc ps[2], int n_threads)
{
    size_t total = ps[0].store.n + ps[1].store.n;
    if (!total) return;
    if (n_threads < 1) n_threads = 1;
    if ((size_t)n_threads > total) n_threads = (int)total;
    if (n_threads > 64) n_threads = 64;

    int t1 = ps[0].store.n ? (int)((double)n_threads * ps[0].store.n / total) : 0;
    if (ps[0].store.n && t1 < 1) t1 = 1;
    int t2 = n_threads - t1;
    if (ps[1].store.n && t2 < 1) t2 = 1;

    pthread_t th[64];
    BaqJob jobs[64];
    int t = 0, p;
    for (p = 0; p < 2; p++) {
        int nt = (p == 0) ? t1 : t2;
        if (!ps[p].store.n || nt <= 0) continue;
        if (nt > (int)ps[p].store.n) nt = (int)ps[p].store.n;
        size_t per = (ps[p].store.n + nt - 1) / nt;
        size_t pos = 0;
        int k;
        for (k = 0; k < nt && pos < ps[p].store.n; k++) {
            jobs[t].store = &ps[p].store;
            jobs[t].start = pos;
            size_t end = pos + per;
            if (end > ps[p].store.n) end = ps[p].store.n;
            jobs[t].end = end;
            pthread_create(&th[t], NULL, baq_worker, &jobs[t]);
            t++;
            pos = end;
        }
    }
    int i;
    for (i = 0; i < t; i++) pthread_join(th[i], NULL);
}

/* Serve the pre-fetched, BAQ-adjusted reads to the pileup engine in the
 * original coordinate order. */
static int prefetch_read_func(void *data, bam1_t *b)
{
    ReadStore *s = data;
    if (s->idx >= s->n) return -1;
    if (!bam_copy1(b, s->reads[s->idx++])) return -1;
    return 1;
}

/* Binomial machinery (replaces Math::Random::random_binomial) */
typedef struct {
    int n;
    double *cdf;
} BinCDF;

static BinCDF *bincdf_new(int n, double p)
{
    BinCDF *c = calloc(1, sizeof *c);
    c->n = n;
    c->cdf = malloc((n + 1) * sizeof(double));
    double *pmf = malloc((n + 1) * sizeof(double));

    /* start from the mode and walk outwards with the PMF recurrence, so only
     * O(1) lgamma evaluations are needed per (n,p) instead of O(n) */
    int k0 = (int)((n + 1) * p);
    if (k0 > n) k0 = n;
    if (k0 < 0) k0 = 0;
    double logp = log(p), logq = log1p(-p);
    double lgamma_np1 = lgamma(n + 1.0);
    pmf[k0] = exp(lgamma_np1 - lgamma(k0 + 1.0) - lgamma(n - k0 + 1.0)
                  + k0 * logp + (n - k0) * logq);
    double r = (1.0 - p) / p;
    int k;
    for (k = k0 - 1; k >= 0; k--)
        pmf[k] = pmf[k + 1] * (k + 1.0) / (n - k) * r;
    r = p / (1.0 - p);
    for (k = k0 + 1; k <= n; k++)
        pmf[k] = pmf[k - 1] * (n - k + 1.0) / k * r;

    double cum = 0;
    for (k = 0; k <= n; k++) {
        cum += pmf[k];
        c->cdf[k] = cum;
    }
    c->cdf[n] = 1.0;
    free(pmf);
    return c;
}

static void bincdf_free(BinCDF *c)
{
    if (!c) return;
    free(c->cdf);
    free(c);
}

static int64_t bincdf_sample(const BinCDF *c, double u)
{
    int lo = 0, hi = c->n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (c->cdf[mid] < u) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

/* xoshiro256** PRNG (deterministic seed from time+pid; thresholds are
 * statistically equivalent to Math::Random's, not bit-identical) */
static uint64_t rng_state[4];

static uint64_t splitmix64(uint64_t *x)
{
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static void rng_seed(uint64_t seed)
{
    uint64_t x = seed;
    int i;
    for (i = 0; i < 4; i++) rng_state[i] = splitmix64(&x);
}

static inline uint64_t rng_next(void)
{
    uint64_t *s = rng_state;
    uint64_t result = ((s[1] * 5) << 7) * 9 ^ (s[0] + s[3]);
    uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = (s[3] << 45) | (s[3] >> 19);
    return result;
}

static double rng_double(void)
{
    return (rng_next() >> 11) * (1.0 / 9007199254740992.0);
}

static int double_cmp_desc(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) ? -1 : (x < y) ? 1 : 0;
}

static int base_index(char base)
{
    switch (base) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return -1;
    }
}

/* Compute SNP-index for one window, returning the formatted
 * "snp_index\tthr0.01\tthr0.05" string (or "NA\tNA\tNA"). */
static void snp_mapper_window(const Config *cfg, PrefetchSrc src[2],
                              int64_t beg, int64_t end, char *out, size_t out_sz)
{
    int si;
    for (si = 0; si < 2; si++) {
        src[si].itr = sam_itr_queryi(src[si].idx, src[si].tid, beg, end);
        if (!src[si].itr) die("cannot query region for SNP-index");
        src[si].store.n = src[si].store.cap = src[si].store.idx = 0;
        src[si].store.reads = NULL;
        if (prefetch_fill(&src[si]) < 0) die("cannot fetch reads for SNP-index");
        hts_itr_destroy(src[si].itr);
        src[si].itr = NULL;
    }

    /* BAQ is off by default: it costs most of the runtime while barely
     * changing the SNP-index on this data. Set HAPBSA_BAQ=1 to enable it
     * (applied in parallel, same order -> same result as sequential). */
    static int use_baq = -1;
    if (use_baq < 0) use_baq = getenv("HAPBSA_BAQ") ? 1 : 0;
    if (use_baq) {
        int baq_threads = 8;
        const char *bt = getenv("HAPBSA_BAQ_THREADS");
        if (bt && atoi(bt) > 0) baq_threads = atoi(bt);
        baq_parallel(src, baq_threads);
    }

    void *data[2] = { &src[0].store, &src[1].store };
    bam_mplp_t iter = bam_mplp_init(2, prefetch_read_func, data);
    bam_mplp_init_overlaps(iter);      /* default smart overlap removal */
    bam_mplp_set_maxcnt(iter, 8000);   /* samtools mpileup default max depth */

    int n_plp[2] = {0, 0};
    const bam_pileup1_t *plp[2] = {NULL, NULL};
    int ptid;
    hts_pos_t pos;
    int64_t sim_snp_count = 0, snp_count = 0;
    double index_sum = 0;
    double acc_thr[2] = {0, 0};
    static BinCDF *bin25 = NULL;
    if (!bin25) bin25 = bincdf_new(25, 0.5);
    while (bam_mplp64_auto(iter, &ptid, &pos, n_plp, plp) > 0) {
        if (ptid != src[0].tid) continue;

        int64_t cnt[2][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}};
        for (si = 0; si < 2; si++) {
            int i;
            for (i = 0; i < n_plp[si]; i++) {
                const bam_pileup1_t *p = &plp[si][i];
                if (p->is_del || p->is_refskip) continue;
                uint8_t q = bam_get_qual(p->b)[p->qpos];
                if (q < 20) continue;
                int code = bam_seqi(bam_get_seq(p->b), p->qpos);
                int bi = base_index(seq_nt16_str[code]);
                if (bi >= 0) cnt[si][bi]++;
            }
        }
        int64_t mutcov = cnt[0][0] + cnt[0][1] + cnt[0][2] + cnt[0][3];
        int64_t wtcov  = cnt[1][0] + cnt[1][1] + cnt[1][2] + cnt[1][3];
        if (mutcov < cfg->depth || wtcov < cfg->depth) continue;
        if (pos < 0 || pos >= src[0].store.ref_len) continue;

        /* Perl: @alleles = sort { $mut{$b} <=> $mut{$a} or $wt{$b} <=> $wt{$a} } keys %mut
         * Ties resolve in A,C,G,T order here (Perl hash order is randomized). */
        int order[4] = {0, 1, 2, 3};
        int x, y;
        for (x = 1; x < 4; x++) {
            int t = order[x];
            y = x - 1;
            while (y >= 0) {
                int a = order[y], b = t;
                int worse = (cnt[0][b] > cnt[0][a]) ||
                            (cnt[0][b] == cnt[0][a] && cnt[1][b] > cnt[1][a]);
                if (!worse) break;
                order[y + 1] = order[y];
                y--;
            }
            order[y + 1] = t;
        }
        int major = order[0];
        double mut_freq = round4((double)cnt[0][major] / (double)mutcov);
        double wt_freq  = round4((double)cnt[1][major] / (double)wtcov);
        double delta = round4(mut_freq - wt_freq);
        if (mut_freq > 0.9 && wt_freq > 0.9) continue;   /* Perl uses rounded freqs here */

        sim_snp_count++;
        if (sim_snp_count <= 3000) {
            BinCDF *cdf_mut[26], *cdf_wt[26];
            int k;
            for (k = 0; k <= 25; k++) {
                double p = (25 - k) / 25.0;
                cdf_mut[k] = (p > 0 && p < 1) ? bincdf_new((int)mutcov, p) : NULL;
                cdf_wt[k]  = (p > 0 && p < 1) ? bincdf_new((int)wtcov, p) : NULL;
            }
            double thr[1000];
            for (k = 0; k < 1000; k++) {
                int xm = bincdf_sample(bin25, rng_double());
                int xw = bincdf_sample(bin25, rng_double());
                double pm = (25 - xm) / 25.0;   /* 1 - mthre */
                double pw = (25 - xw) / 25.0;   /* 1 - wthre */
                int64_t mr, wr;
                if (pm <= 0) mr = 0;
                else if (pm >= 1) mr = mutcov;
                else mr = bincdf_sample(cdf_mut[xm], rng_double());
                if (pw <= 0) wr = 0;
                else if (pw >= 1) wr = wtcov;
                else wr = bincdf_sample(cdf_wt[xw], rng_double());
                thr[k] = fabs((double)mr / (double)mutcov - (double)wr / (double)wtcov);
            }
            qsort(thr, 1000, sizeof(double), double_cmp_desc);
            acc_thr[0] += thr[9];
            acc_thr[1] += thr[49];
            for (k = 0; k <= 25; k++) { bincdf_free(cdf_mut[k]); bincdf_free(cdf_wt[k]); }
        }
        snp_count++;
        index_sum += delta;
    }

    bam_mplp_destroy(iter);
    for (si = 0; si < 2; si++) {
        size_t i;
        for (i = 0; i < src[si].store.n; i++) bam_destroy1(src[si].store.reads[i]);
        free(src[si].store.reads);
        src[si].store.reads = NULL;
        src[si].store.n = src[si].store.cap = src[si].store.idx = 0;
    }

    if (snp_count == 0) {
        snprintf(out, out_sz, "NA\tNA\tNA");
        return;
    }
    char buf[4][64];
    fmt4(buf[0], sizeof buf[0], fabs(index_sum / (double)snp_count));
    int64_t sim_count = sim_snp_count > 3000 ? 3000 : sim_snp_count;
    if (sim_count > 0) {
        fmt4(buf[1], sizeof buf[1], acc_thr[0] / (double)sim_count);
        fmt4(buf[2], sizeof buf[2], acc_thr[1] / (double)sim_count);
    } else {
        snprintf(buf[1], sizeof buf[1], "NA");
        snprintf(buf[2], sizeof buf[2], "NA");
    }
    snprintf(out, out_sz, "%s\t%s\t%s", buf[0], buf[1], buf[2]);
}

/* ------------------------------------------------------------------ */
/* Per-window driver + child workers                                   */
/* ------------------------------------------------------------------ */

typedef struct {
    samFile *fp1, *fp2;
    hts_idx_t *idx1, *idx2;
    faidx_t *fai;
    char *refseq;
    hts_pos_t ref_len;
    int tid;
} ChildHandles;

static void child_open(const Config *cfg, const char *chr, ChildHandles *h)
{
    memset(h, 0, sizeof *h);
    h->fp1 = sam_open(cfg->bam1, "r");
    h->fp2 = sam_open(cfg->bam2, "r");
    if (!h->fp1 || !h->fp2) die("cannot open pool BAM files");
    h->idx1 = sam_index_load(h->fp1, cfg->bam1);
    h->idx2 = sam_index_load(h->fp2, cfg->bam2);
    if (!h->idx1 || !h->idx2) die("cannot load BAM indexes");
    h->fai = fai_load(cfg->refgenome);
    if (!h->fai) die("cannot load reference FASTA index");
    sam_hdr_t *hdr = sam_hdr_read(h->fp1);
    if (!hdr) die("cannot read BAM header");
    h->tid = sam_hdr_name2tid(hdr, chr);
    sam_hdr_destroy(hdr);
    if (h->tid < 0) die("chromosome not found in BAM header");
    h->refseq = faidx_fetch_seq64(h->fai, chr, 0, HTS_POS_MAX, &h->ref_len);
    if (!h->refseq) die("cannot fetch reference sequence");
}

static void child_close(ChildHandles *h)
{
    hts_idx_destroy(h->idx1);
    hts_idx_destroy(h->idx2);
    sam_close(h->fp1);
    sam_close(h->fp2);
    free(h->refseq);
    fai_destroy(h->fai);
}

static int process_one_window(const Config *cfg, const Window *win,
                              const char *tmpdir, ChildHandles *h, PrefetchSrc mplp[2])
{
    const ChrHap *chr = win->chr;
    int64_t beg = win->s > 0 ? win->s - 1 : 0;
    int64_t end = win->e;

    /* sites in [s,e] (Perl writes only those to the window hap file) */
    int lo = 0, hi = (int)chr->n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (chr->sites[mid].pos < win->s) lo = mid + 1;
        else hi = mid;
    }
    int lo2 = lo;
    hi = (int)chr->n;
    while (lo2 < hi) {
        int mid = (lo2 + hi) / 2;
        if (chr->sites[mid].pos <= win->e) lo2 = mid + 1;
        else hi = mid;
    }
    const HapSite *sites = chr->sites + lo;
    int n_sites = lo2 - lo;

    BlockMap bm;
    memset(&bm, 0, sizeof bm);
    int64_t totals[2][2] = {{0, 0}, {0, 0}};
    if (phase_pool(h->fp1, h->idx1, h->tid, beg, end, 45, sites, n_sites, cfg, 0, &bm, totals[0]) != 0) return -1;
    if (phase_pool(h->fp2, h->idx2, h->tid, beg, end, 45, sites, n_sites, cfg, 1, &bm, totals[1]) != 0) return -1;

    char hapval[32];
    if (cfg->has_block) {
        hap_index_block(&bm, cfg, hapval, sizeof hapval);
    } else {
        hap_index_plain(totals[0][0], totals[0][1], totals[1][0], totals[1][1], cfg, hapval, sizeof hapval);
    }
    blockmap_free(&bm);

    char snpval[128];
    snp_mapper_window(cfg, mplp, beg, end, snpval, sizeof snpval);

    char path[4096];
    snprintf(path, sizeof path, "%s/win.%06d.tmp.hapBSA.txt", tmpdir, win->win_idx);
    FILE *fp = fopen(path, "w");
    if (!fp) return -1;
    int64_t mid = (win->s + win->e) / 2;
    fprintf(fp, "%s\t%" PRId64 "\t%s\t%s\n", chr->chr, mid, hapval, snpval);
    fclose(fp);
    return 0;
}

static void child_worker(const Config *cfg, const Window *wins, size_t n_wins,
                         int child_id, int n_children, const char *tmpdir)
{
    size_t i;
    ChildHandles h;
    const char *cur_chr = NULL;
    int h_open = 0;
    PrefetchSrc mplp[2];
    for (i = 0; i < n_wins; i++) {
        if ((int)(i % (size_t)n_children) != child_id) continue;
        if (!h_open || strcmp(cur_chr, wins[i].chr->chr) != 0) {
            if (h_open) child_close(&h);
            child_open(cfg, wins[i].chr->chr, &h);
            h_open = 1;
            cur_chr = wins[i].chr->chr;
            memset(mplp, 0, sizeof mplp);
            mplp[0].fp = h.fp1; mplp[0].idx = h.idx1; mplp[0].tid = h.tid;
            mplp[0].min_mq = 45; mplp[0].store.ref = h.refseq; mplp[0].store.ref_len = h.ref_len;
            mplp[1].fp = h.fp2; mplp[1].idx = h.idx2; mplp[1].tid = h.tid;
            mplp[1].min_mq = 45; mplp[1].store.ref = h.refseq; mplp[1].store.ref_len = h.ref_len;
        }
        if (process_one_window(cfg, &wins[i], tmpdir, &h, mplp) != 0) {
            fprintf(stderr, "ERROR: failed processing window %d\n", wins[i].win_idx);
            _exit(1);
        }
    }
    if (h_open) child_close(&h);
    _exit(0);
}

/* ------------------------------------------------------------------ */
/* Output merge                                                        */
/* ------------------------------------------------------------------ */

typedef struct {
    char *chr;
    int64_t pos;
    char *line;
} ResLine;

static int resline_cmp(const void *a, const void *b)
{
    const ResLine *x = a, *y = b;
    int r = strcmp(x->chr, y->chr);
    if (r) return r;
    return (x->pos > y->pos) - (x->pos < y->pos);
}

static int merge_results(const char *tmpdir, size_t n_wins, const Config *cfg)
{
    ResLine *res = calloc(n_wins, sizeof(*res));
    size_t n_res = 0, i;
    char path[4096], line[4096];
    for (i = 0; i < n_wins; i++) {
        snprintf(path, sizeof path, "%s/win.%06zu.tmp.hapBSA.txt", tmpdir, i);
        FILE *fp = fopen(path, "r");
        if (!fp) { fprintf(stderr, "ERROR: missing temporary result %s\n", path); return -1; }
        if (fgets(line, sizeof line, fp)) {
            line[strcspn(line, "\r\n")] = 0;
            char *copy = strdup(line);
            char *save = NULL;
            char *chr = strtok_r(copy, "\t", &save);
            char *pos = strtok_r(NULL, "\t", &save);
            if (chr && pos) {
                res[n_res].line = strdup(line);
                res[n_res].chr = strdup(chr);
                res[n_res].pos = atoll(pos);
                n_res++;
            }
            free(copy);
        }
        fclose(fp);
    }
    qsort(res, n_res, sizeof(ResLine), resline_cmp);

    char outpath[4096];
    snprintf(outpath, sizeof outpath, "%s.hapBSA.sliding_window.txt", cfg->outprefix);
    FILE *fp = fopen(outpath, "w");
    if (!fp) { fprintf(stderr, "ERROR: cannot write %s\n", outpath); return -1; }
    for (i = 0; i < n_res; i++) fprintf(fp, "%s\n", res[i].line);
    fclose(fp);

    for (i = 0; i < n_res; i++) { free(res[i].line); free(res[i].chr); }
    free(res);
    return 0;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    Config cfg;
    parse_args(argc, argv, &cfg);
    fflush(stdout);

    time_t t0 = time(NULL);
    srand((unsigned)(time(NULL) ^ getpid()));
    rng_seed((uint64_t)time(NULL) ^ (uint64_t)getpid() ^ 0x9E3779B97F4A7C15ULL);

    ChrHap **chrs = NULL;
    size_t n_chrs = 0;
    read_hap_file(&cfg, &chrs, &n_chrs);

    Window *wins = NULL;
    size_t n_wins = 0;
    build_windows(chrs, n_chrs, &cfg, &wins, &n_wins);

    char tmpdir[4096];
    int rnd = (rand() % 100000) + 1;
    snprintf(tmpdir, sizeof tmpdir, "%s/tmpdir%d", cfg.tmp_dir, rnd);
    if (mkdir(tmpdir, 0755) != 0) {
        fprintf(stderr, "ERROR: cannot create %s: %s\n", tmpdir, strerror(errno));
        return 1;
    }

    int n_children = cfg.cpu;
    if ((size_t)n_children > n_wins) n_children = (int)n_wins;
    if (n_children < 1) n_children = 1;

    pid_t *pids = calloc(n_children, sizeof(pid_t));
    int c;
    for (c = 0; c < n_children; c++) {
        pid_t pid = fork();
        if (pid < 0) { perror("fork"); return 1; }
        if (pid == 0) child_worker(&cfg, wins, n_wins, c, n_children, tmpdir);
        pids[c] = pid;
    }
    int status;
    for (c = 0; c < n_children; c++) {
        while (waitpid(pids[c], &status, 0) < 0 && errno == EINTR) {}
        if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {
            fprintf(stderr, "ERROR: worker %d failed\n", c);
            return 1;
        }
    }
    free(pids);

    if (merge_results(tmpdir, n_wins, &cfg) != 0) return 1;

    /* cleanup: remove per-window tmp files and the tmp dir (Perl: rm -rf) */
    {
        size_t i;
        char path[4096];
        for (i = 0; i < n_wins; i++) {
            snprintf(path, sizeof path, "%s/win.%06zu.tmp.hapBSA.txt", tmpdir, i);
            unlink(path);
        }
        rmdir(tmpdir);
    }
    char legacy[4096];
    snprintf(legacy, sizeof legacy, "%s.hapBSA.tmp.txt", cfg.outprefix);
    unlink(legacy);

    double elapsed = (double)(time(NULL) - t0) / 60.0;
    printf("%s finished! Total time elapsed: %.2f min\n", argv[0], elapsed);

    size_t i;
    for (i = 0; i < n_wins; i++) {}
    for (i = 0; i < n_chrs; i++) {
        free(chrs[i]->sites);
        free(chrs[i]->chr);
        free(chrs[i]);
    }
    free(chrs);
    free(wins);
    return 0;
}
