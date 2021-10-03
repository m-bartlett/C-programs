/* https://www.chiark.greenend.org.uk/~sgtatham/filigram/filigram.c */

/*
 * This single standalone C program is not the original source form of
 * this program. For convenience of distribution, it's been bolted
 * together into this form by a Perl script, from a collection of more
 * normally organised modules.
 *
 * If you'd rather see the original version, you can find it at
 *   https://git.tartarus.org/simon/pix.git
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdbool.h>
#include <assert.h>
#include <limits.h>

/* -------------------- originally from lz77.h -------------------- */

/*
 * Header file for an LZ77 implementation.
 */

typedef struct LZ77 LZ77;

/*
 * Set up an LZ77 context. The user supplies two function pointer
 * parameters: `literal', called when the compressor outputs a
 * literal byte, and `match', called when the compressor outputs a
 * match. Each of these functions, in addition to parameters
 * describing the data being output, also receives a void * context
 * parameter which is passed in to lz77_new().
 */
LZ77 *lz77_new(void (*literal)(void *ctx, unsigned char c),
               void (*match)(void *ctx, int distance, int len),
               void *ctx);

/*
 * Supply data to be compressed. This produces output by calling
 * the literal() and match() functions passed to lz77_new().
 * 
 * This function buffers data internally, so in any given call it
 * will not necessarily output literals and matches which cover all
 * the data passed in to it. To force a buffer flush, use
 * lz77_flush().
 */
void lz77_compress(LZ77 *lz, const void *data, int len);

/*
 * Force the LZ77 compressor to output literals and matches which
 * cover all the data it has received as input until now.
 *
 * After calling this function, the LZ77 stream is still active:
 * future input data will still be matched against data which was
 * provided before the flush. Therefore, this function can be used
 * to produce guaranteed record boundaries in a single data stream,
 * but can _not_ be used to reuse the same LZ77 context for two
 * entirely separate data streams (otherwise the second one might
 * make backward references beyond the start of its data).
 * 
 * Calling this in the middle of a data stream can decrease
 * compression quality. If maximum compression is the only concern,
 * do not call this function until the very end of the stream.
 */
void lz77_flush(LZ77 *lz);

/*
 * Free an LZ77 context when it's finished with. This does not
 * automatically flush buffered data; call lz77_flush explicitly to
 * do that.
 */
void lz77_free(LZ77 *lz);

/* -------------------- originally from deflate.h -------------------- */

/*
 * Header file for my independent implementation of Deflate
 * (RFC1951) compression.
 */

/*
 * Types of Deflate data stream.
 * 
 * DEFLATE_TYPE_BARE represents the basic Deflate data format, as
 * defined in RFC 1951. It has no checksum to detect errors and no
 * magic-number header for ease of recognition, but it does have
 * internal EOF indication.
 * 
 * DEFLATE_TYPE_ZLIB represents the zlib container format, as
 * defined in RFC 1950. It has a two-byte header, and a four-byte
 * Adler32 checksum at the end to verify correct decoding, but
 * apart from those six bytes it's exactly equivalent to
 * DEFLATE_TYPE_BARE.
 * 
 * DEFLATE_TYPE_GZIP represents the gzip compressed file format, as
 * defined in RFC 1952. This is a more full-featured format, with a
 * magic number, a CRC checksum of the compressed data, and various
 * header features including storing the original filename. This
 * implementation accepts but ignores all of those features on
 * input except the checksum, and outputs them in the most trivial
 * fashion. Also, this implementation will not decode multiple
 * concatenated gzip members (permitted by the RFC).
 */
enum {
    DEFLATE_TYPE_BARE,
    DEFLATE_TYPE_ZLIB,
    DEFLATE_TYPE_GZIP
};

/* ----------------------------------------------------------------------
 * Compression functions. Create a compression context with
 * deflate_compress_new(); feed it data with repeated calls to
 * deflate_compress_data(); destroy it with
 * deflate_compress_free().
 */

typedef struct deflate_compress_ctx deflate_compress_ctx;

/*
 * Create a new compression context. `type' indicates whether it's
 * bare Deflate (as used in, say, zip files) or Zlib (as used in,
 * say, PDF).
 */
deflate_compress_ctx *deflate_compress_new(int type);

/*
 * Free a compression context previously allocated by
 * deflate_compress_new().
 */
void deflate_compress_free(deflate_compress_ctx *ctx);

/*
 * Give the compression context some data to compress. The input
 * data is passed in `inblock', and has length `inlen'. This
 * function may or may not produce some output data; if so, it is
 * written to a dynamically allocated chunk of memory, a pointer to
 * that memory is stored in `outblock', and the length of output
 * data is stored in `outlen'. It is common for no data to be
 * output, if the input data has merely been stored in internal
 * buffers.
 * 
 * `flushtype' indicates whether you want to force buffered data to
 * be output. It can be one of the following values:
 * 
 *  - DEFLATE_NO_FLUSH: nothing is output if the compressor would
 *    rather not. Use this when the best compression is desired
 *    (i.e. most of the time).
 *
 *  - DEFLATE_SYNC_FLUSH: all the buffered data is output, but the
 *    compressed data stream remains open and ready to continue
 *    compressing data. Use this in interactive protocols when a
 *    single compressed data stream is split across several network
 *    packets.
 * 
 *  - DEFLATE_END_OF_DATA: all the buffered data is output and the
 *    compressed data stream is cleaned up. Any checksums required
 *    at the end of the stream are also output.
 */
void deflate_compress_data(deflate_compress_ctx *ctx,
                           const void *inblock, int inlen, int flushtype,
                           void **outblock, int *outlen);

enum {
    DEFLATE_NO_FLUSH,
    DEFLATE_SYNC_FLUSH,
    DEFLATE_END_OF_DATA
};

/* ----------------------------------------------------------------------
 * Decompression functions. Create a decompression context with
 * deflate_decompress_new(); feed it data with repeated calls to
 * deflate_decompress_data(); destroy it with
 * deflate_decompress_free().
 */

typedef struct deflate_decompress_ctx deflate_decompress_ctx;

/*
 * Create a new decompression context. `type' means the same as it
 * does in deflate_compress_new().
 */
deflate_decompress_ctx *deflate_decompress_new(int type);

/*
 * Free a decompression context previously allocated by
 * deflate_decompress_new().
 */
void deflate_decompress_free(deflate_decompress_ctx *ctx);

/*
 * Give the decompression context some data to decompress. The
 * input data is passed in `inblock', and has length `inlen'. This
 * function may or may not produce some output data; if so, it is
 * written to a dynamically allocated chunk of memory, a pointer to
 * that memory is stored in `outblock', and the length of output
 * data is stored in `outlen'.
 *
 * Returns 0 on success, or a non-zero error code if there was a
 * decoding error. In case of an error return, the data decoded
 * before the error is still returned as well. The possible errors
 * are listed below.
 * 
 * If you want to check that the compressed data stream was
 * correctly terminated, you can call this function with inlen==0
 * to signal input EOF and see if an error comes back. If you don't
 * care, don't bother.
 */
int deflate_decompress_data(deflate_decompress_ctx *ctx,
                            const void *inblock, int inlen,
                            void **outblock, int *outlen);

/*
 * Enumeration of error codes. The strange macro is so that I can
 * define description arrays in the accompanying source.
 */
#define DEFLATE_ERRORLIST(A) \
    A(DEFLATE_NO_ERR, "success"), \
    A(DEFLATE_ERR_ZLIB_HEADER, "invalid zlib header"), \
    A(DEFLATE_ERR_ZLIB_WRONGCOMP, "zlib header specifies non-deflate compression"), \
    A(DEFLATE_ERR_GZIP_HEADER, "invalid gzip header"), \
    A(DEFLATE_ERR_GZIP_WRONGCOMP, "gzip header specifies non-deflate compression"), \
    A(DEFLATE_ERR_GZIP_FHCRC, "gzip header specifies disputed FHCRC flag"), \
    A(DEFLATE_ERR_SMALL_HUFTABLE, "under-committed Huffman code space"), \
    A(DEFLATE_ERR_LARGE_HUFTABLE, "over-committed Huffman code space"), \
    A(DEFLATE_ERR_UNCOMP_HDR, "wrongly formatted header in uncompressed block"), \
    A(DEFLATE_ERR_NODISTTABLE, "backward copy encoded in block without distances table"), \
    A(DEFLATE_ERR_BADDISTCODE, "invalid distance code 30 or 31 found in block"), \
    A(DEFLATE_ERR_CHECKSUM, "incorrect data checksum"), \
    A(DEFLATE_ERR_INLEN, "incorrect data length"), \
    A(DEFLATE_ERR_UNEXPECTED_EOF, "unexpected end of data")
#define DEFLATE_ENUM_DEF(x,y) x
enum { DEFLATE_ERRORLIST(DEFLATE_ENUM_DEF), DEFLATE_NUM_ERRORS };
#undef DEFLATE_ENUM_DEF

/*
 * Arrays mapping the above error codes to, respectively, a text
 * error string and a textual representation of the symbolic error
 * code.
 */
extern const char *const deflate_error_msg[DEFLATE_NUM_ERRORS];
extern const char *const deflate_error_sym[DEFLATE_NUM_ERRORS];

/* -------------------- originally from crc32.h -------------------- */

unsigned long crc32_update(unsigned long crcword, const void *vdata, int len);
/* -------------------- originally from pngout.h -------------------- */


void *pngwrite(int w, int h, int channels, int depth,
               const unsigned short *bitmap, size_t *size);

/* -------------------- originally from cmdline.h -------------------- */

/* ----------------------------------------------------------------------
 * Declarations for generic command-line parsing functions.
 */

bool parsestr(char *string, void *ret);
bool parseint(char *string, void *ret);
bool parsesignedint(char *string, void *ret);
bool parsesize(char *string, void *ret);
bool parseflt(char *string, void *ret);
bool parsebool(char *string, void *ret);
bool parsecol(char *string, void *ret);
bool incrementint(char *string, void *vret);

struct Cmdline {
    int nlongopts;
    char *longopt;                     /* `nlongopts' NUL-separated strings */
    char shortopt;
    char *arghelp;
    char *deschelp;
    char *valname;                     /* for use in `cannot parse' messages */
    bool (*parse)(char *string, void *ret);
    int parse_ret_off;                 /* offset into options structure */
    int gotflag_off;                   /* and another one */
};

void parse_cmdline(char const *programname, int argc, char **argv,
                   const struct Cmdline *options, int noptions, void *optdata);

void usage_message(char const *usageline,
                   const struct Cmdline *options, int noptions,
                   char **extratext, int nextra);

/* -------------------- originally from bmpwrite.h -------------------- */

/* ----------------------------------------------------------------------
 * Declarations for functions to write out .BMP files.
 */

typedef enum { BMP, PPM, PNG } imagetype;

imagetype infer_type(const char *progname, const char *type,
                     const char *filename);

struct Bitmap;
struct Bitmap *bmpinit(char const *filename, int width, int height,
                       imagetype type);
void bmppixel(struct Bitmap *bm,
              unsigned char r, unsigned char g, unsigned char b);
void bmpendrow(struct Bitmap *bm);
void bmpclose(struct Bitmap *bm);

/* -------------------- originally from misc.h -------------------- */

/* ----------------------------------------------------------------------
 * Generally useful declarations.
 */

#define lenof(x) (sizeof ((x)) / sizeof ( *(x) ))

struct Size {
    int w, h;
};

struct RGB {
    double r, g, b;
};

/* -------------------- originally from lz77.c -------------------- */

/*
 * LZ77 compressor. In the probably vain hope that this will be the
 * _last time_ I have to write one.
 */

/*
 * Modifiable parameters. These are currently tuned for Deflate
 * (RFC1951) compression, and deflate.c expects this when linking
 * against it.
 */
#define MAXMATCHDIST 32768             /* maximum backward distance */
#define MAXMATCHLEN 258                /* maximum length of a match */
#define HASHMAX 2039                   /* one more than max hash value */
#define HASHCHARS 3                    /* how many chars make a hash */
#define MAXLAZY 3                      /* limit of lazy matching */

/*
 * Derived parameters.
 */
#define MAXREWIND (HASHCHARS * 2)
#define WINSIZE (MAXMATCHDIST + HASHCHARS + MAXREWIND)
#define NMATCHES (MAXLAZY + MAXLAZY * MAXLAZY)

struct LZ77 {
    /*
     * Administrative data passed in from lz77_new().
     */
    void *ctx;
    void (*literal)(void *ctx, unsigned char c);
    void (*match)(void *ctx, int distance, int len);

    /*
     * The actual data in the sliding window, stored as a circular
     * buffer. `winpos' marks the head of the buffer. So
     * data[winpos] is the most recent character; data[winpos+1] is
     * the most recent but one; and data[winpos-1] (all mod
     * WINSIZE, of course) is the least recent character still in
     * the buffer.
     * 
     * We store slightly more than one windowful of data: the first
     * HASHCHARS bytes of the window are duplicated at the far end.
     * This permits us to retrieve HASHCHARS bytes starting at any
     * window position without having to do expensive mod
     * operations.
     */
    unsigned char data[WINSIZE+HASHCHARS];
    int winpos;

    /*
     * This variable indicates the number of input characters we
     * have accepted into the window which haven't been output yet.
     * It can exceed WINSIZE if we're tracking a particularly long
     * match.
     */
    int k;

    /*
     * This variable indicates the number of valid bytes in the
     * sliding window. It starts at zero at the beginning of the
     * data stream, then rapidly rises to WINSIZE; it can also
     * decrease by small amounts during the compression when we
     * need to rewind the state by a few bytes.
     */
    int nvalid;

    /*
     * This variable indicates that the compressor state has been
     * rewound by a few bytes, which causes the main loop to reuse
     * characters already in the sliding window instead of taking
     * more from the input data.
     */
    int rewound;

    /*
     * This array links the window positions into disjoint lists
     * chained from the hash table. hashnext[pos] gives the
     * position of the next most recent sequence of three
     * characters with the same hash as this one. Indices in this
     * table have the same meaning as indices in `data', and mark
     * the _most recent_ of the three characters involved.
     * 
     * The list need only be singly linked: we detect the end of a
     * hash chain by observing that pos and hashnext[pos] are on
     * opposite sides of the maximum match distance.
     * 
     * If hashnext[pos] < 0, that also indicates the end of a hash
     * chain. This only comes up at the start of the algorithm
     * before the hash chain has been used.
     */
    int hashnext[WINSIZE];

    /*
     * The hash table. For each hash value, this stores the
     * location within the window of the most recent sequence of
     * three characters with that hash value, or -1 if there hasn't
     * been one yet.
     */
    int hashhead[HASHMAX];

    /*
     * This list links together the set of backward distances which
     * currently constitute valid matches. Unlike the other arrays
     * of size WINSIZE, this one is indexed by absolute backward
     * distance: its indices are not dependent on winpos.
     *
     * This list can contain multiple disjoint linked lists when
     * we're doing lazy matching. (Two candidate matches starting
     * at positions which differ by less than HASHCHARS cannot both
     * have the same backward distance, because if they did then
     * they'd be part of the same longer match and obviously the
     * one starting earlier would be preferable; so there can be no
     * collision within this array.)
     * 
     * matchnext[pos] is
     *  - less than 0 if pos is the last element in such a list
     *  - exactly 0 if pos is not part of such a list at all (this
     *    does not clash with a real backward distance because real
     *    backward distances must be at least 1!)
     *  - otherwise gives the next shortest backward distance of a
     *    so-far valid match.
     */
    int matchnext[WINSIZE];

    /*
     * These variables track the heads of linked lists in
     * nextmatch. Each value matchhead[i] points to the head of a
     * list of possible matches, or -1 if that list is currently
     * empty. Indices up to and including HASHCHARS-1 always
     * indicate a match starting i bytes after the last output
     * byte; indices after that are allocated dynamically by the
     * `matchlater' array.
     */
    int matchhead[NMATCHES];

    /*
     * These variables track the current length, and best backward
     * distance, of the matches listed in matchhead. Entries in
     * these two arrays can persist after matchhead[i] has become
     * -1, because only after we know how long all the matches are
     * going to be can we decide which one to output.
     */
    int matchlen[NMATCHES], matchdist[NMATCHES];

    /*
     * This array stores dynamically allocated indices in the
     * `matchhead', `matchlen' and `matchdist' arrays, starting
     * from HASHCHARS. These indices are used to track matches
     * which start after other matches have finished, in order to
     * make the best available choice when deciding on a lazy
     * matching strategy.
     * 
     * matchlater[p*HASHCHARS+q] stores the index of a match which
     * begins q bytes after the end of the one in matchhead[p], or
     * is zero if no such match is known.
     * 
     * Also, `nextindex' stores the next unused index in the
     * matchhead array and friends.
     */
    int matchlater[HASHCHARS*HASHCHARS];
    int nextindex;

    /*
     * These are the literals saved during lazy matching: if we
     * output a match starting at k+i, we need i literals from
     * position k.
     */
    unsigned char literals[HASHCHARS * (HASHCHARS+1)];
};

static int lz77_hash(const unsigned char *data) {
    return (257*data[0] + 263*data[1] + 269*data[2]) % HASHMAX;
}

LZ77 *lz77_new(void (*literal)(void *ctx, unsigned char c),
               void (*match)(void *ctx, int distance, int len),
               void *ctx)
{
    LZ77 *lz;
    int i;

    lz = (LZ77 *)malloc(sizeof(LZ77));
    if (!lz)
        return NULL;

    lz->ctx = ctx;
    lz->literal = literal;
    lz->match = match;

    memset(lz->data, 0, sizeof(lz->data));

    for (i = 0; i < WINSIZE; i++) {
        lz->hashnext[i] = -1;
        lz->matchnext[i] = 0;
    }

    for (i = 0; i < HASHMAX; i++) {
        lz->hashhead[i] = -1;
    }

    for (i = 0; i < NMATCHES; i++) {
        lz->matchhead[i] = -1;
        lz->matchlen[i] = 0;
        lz->matchdist[i] = 0;
    }

    for (i = 0; i < HASHCHARS * (HASHCHARS+1); i++) {
        lz->literals[i] = 0;
    }

    for (i = 0; i < HASHCHARS * HASHCHARS; i++) {
        lz->matchlater[i] = 0;
    }

    lz->nextindex = HASHCHARS;
    lz->winpos = 0;
    lz->k = 0;
    lz->nvalid = 0;
    lz->rewound = 0;

    return lz;
}

void lz77_free(LZ77 *lz)
{
    free(lz);
}

static void lz77_hashsearch(LZ77 *lz, int hash, unsigned char *currchars,
                            int index)
{
    int pos, nextpos, matchlimit;
    int posval, nextposval;

    assert(lz->matchhead[index] < 0);
    assert(lz->matchdist[index] == 0);
    assert(lz->matchlen[index] == 0);

    /*
     * Find the most distant place in the window where we could
     * viably start a match. This is limited by the maximum match
     * distance, and it's also limited by the current amount of
     * valid data in the window.
     */
    matchlimit = (lz->nvalid < MAXMATCHDIST ? lz->nvalid : MAXMATCHDIST);
    matchlimit = (matchlimit + lz->winpos) % WINSIZE;

    pos = lz->hashhead[hash];
    posval = (matchlimit + WINSIZE - pos) % WINSIZE;
    while (pos != -1) {
        int bdist;
        int i;

        /*
         * Compute the backward distance.
         */
        bdist = (pos + WINSIZE - lz->winpos) % WINSIZE;

        /*
         * If there's already a match being tracked at this
         * distance, don't bother trying this one. Also, it's
         * possible this match may be in the future (because if we
         * rewound the compressor, some entries at the head of the
         * hash chain will not be valid again _yet_), so check that
         * too.
         */
        if (bdist > 0 && bdist <= MAXMATCHDIST && !lz->matchnext[bdist]) {

            /*
             * Make sure the characters at pos really do match the
             * ones in currchars.
             */
            for (i = 0; i < HASHCHARS; i++)
                if (lz->data[pos + i] != currchars[i])
                    break;
            if (i == HASHCHARS) {
                /*
                 * They did, so this is a valid match position. Add
                 * it to the list of possible match locations.
                 */
                lz->matchnext[bdist] = lz->matchhead[index];
                lz->matchhead[index] = bdist;
                /*
                 * We're working through the window from most to
                 * least recent, and we want to prefer matching
                 * against more recent data; so here we only set
                 * matchdist on the first (i.e. most recent) match
                 * we find.
                 */
                if (!lz->matchlen[index]) {
                    lz->matchlen[index] = HASHCHARS;
                    lz->matchdist[index] = bdist;
                }
            }
        }

        /*
         * Step along the hash chain to try the next candidate
         * location. If we've gone past matchlimit, stop.
         * 
         * Special case: we could also link to ourself here. This
         * occurs when this window position is repeating a previous
         * hash and nothing else has had the same hash in between.
         */
        nextpos = lz->hashnext[pos];
        nextposval = (matchlimit + WINSIZE - nextpos) % WINSIZE;
        if (nextposval >= posval) {
            break;
        } else {
            pos = nextpos;
            posval = nextposval;
        }
    }
}

static void lz77_outputmatch(LZ77 *lz)
{
    int besti = -1, bestval = -1;
    int i;
    int matchadjust[HASHCHARS*HASHCHARS], matchstart[NMATCHES];

    /*
     * Before we do anything else, we need to fiddle with the
     * `matchlater' array. It's possible that a match we need to be
     * aware of in a particular position won't be listed there. An
     * example is if we track three possible matches at k, k+1 and
     * k+2, the middle one ends first and a new match begins
     * immediately, then the match at k ends, and the match after
     * the one at k+1 carries on. In this situation the
     * second-order match would work just as well after k as it
     * would after k+1, and this may be vital to know. So here we
     * go through and copy entries of `matchlater' into each other.
     */
    memset(matchstart, 0, sizeof(matchstart));
    for (i = 0; i < MAXLAZY; i++) {
        /* First tag each match with its starting point relative to k. */
        int j;

        for (j = 0; j < HASHCHARS; j++) {
            int m = lz->matchlater[i*HASHCHARS+j];
            if (m)
                matchstart[m] = i + lz->matchlen[i] + j;
        }
    }
    /* Now go through the matches and replace each one with an earlier-
     * starting one if it can do better that way. */
    memset(matchadjust, 0, sizeof(matchadjust));
    for (i = 0; i < HASHCHARS; i++) {
        int k;
        for (k = 0; k < HASHCHARS; k++) {
            int ii = lz->matchlater[i*HASHCHARS+k];
            int j;

            if (!ii)
                continue;

            for (j = HASHCHARS; j < ii; j++) {
                /*
                 * To be worth replacing the current match in this
                 * slot, the replacement match must start at least
                 * as early, finish strictly later, and be at least
                 * the minimum length.
                 */
                if ((matchstart[j] <= matchstart[ii]) &&
                    (matchstart[j] + lz->matchlen[j] >
                     matchstart[ii] + lz->matchlen[ii]) &&
                    (matchstart[j] + lz->matchlen[j] >
                     matchstart[ii] + HASHCHARS)) {
                    /*
                     * Replace it, and remember that its length isn't
                     * quite what we'd have expected it to be.
                     */
                    lz->matchlater[i*HASHCHARS+k] = j;
                    matchadjust[i*HASHCHARS+k] =
                        matchstart[ii] - matchstart[j];
                    /*
                     * Now don't do any more adjustment of this i.
                     * (Once we know that a match of length n can
                     * be found straight after another match, we
                     * don't also need to know that a match of n-1
                     * at the same distance can be found one byte
                     * later.)
                     */
                    goto nexti;
                }
            }
        }
        nexti:;
    }

    /*
     * Decide which match to output. The basic rule is simply that
     * we pick the longest match we can find, but break ties in
     * favour of earlier matches (since a literal we know we have
     * to output before the match is worse than a literal after the
     * match that we _might_ get away without).
     *
     * However, there's another consideration: if an early-starting
     * match has a subsequent match beginning shortly after its end
     * which lasts until at least MAXLAZY-1 bytes after the end of
     * the _longest_ available match, then we must consider the
     * number of literals output in each case and potentially
     * disqualify later-starting matches on that basis.
     */
    for (i = 0; i < MAXLAZY; i++) {
        int j, k;

        /*
         * See if there's a second-order match which disqualifies
         * match[i] from consideration.
         */
        for (j = 0; j < i; j++) {    /* j indexes earlier 1st-order matches */
            for (k = 0; k < HASHCHARS; k++) {   /* k indexes subsequent posn */
                int mi = j*HASHCHARS+k;
                if (lz->matchlater[mi]) {
                    int m = lz->matchlater[mi];
                    int ma = matchadjust[mi];

                    /*
                     * The match must last until at least MAXLAZY-1
                     * bytes after the end of this one, _and_ it
                     * must end within HASHCHARS bytes of the
                     * current window position.
                     */
                    int totallen = j+lz->matchlen[j] + k + lz->matchlen[m]-ma;
                    if (lz->k - totallen <= HASHCHARS &&
                        totallen >= i + lz->matchlen[i] + MAXLAZY-1)
                        goto disqualified;   /* high-order `continue;' */
                }
            }
        }

        /*
         * If we reach here, there wasn't, so compare matches in
         * the usual way.
         */
        if (lz->matchlen[i] > bestval) {
            bestval = lz->matchlen[i];
            besti = i;
        }

        disqualified:;
    }

    /*
     * Output the match plus its pre-literals (if any).
     */
    for (i = 0; i < besti; i++)
        lz->literal(lz->ctx, lz->literals[i]);
    lz->match(lz->ctx, lz->matchdist[besti], lz->matchlen[besti]);

    /*
     * Decrease k by the amount of data we've just output.
     */
    lz->k -= besti + lz->matchlen[besti];

    /*
     * Make sure the first-order component of the match arrays has
     * been cleaned out. (In normal use this will happen naturally,
     * because this function won't be called until all matches have
     * terminated; but if we're called from lz77_flush() then we
     * might still have live candidate matches.)
     */
    for (i = 0; i < HASHCHARS; i++) {
        lz->matchdist[i] = lz->matchlen[i] = 0;
        while (lz->matchhead[i] >= 0) {
            int tmp = lz->matchhead[i];

            lz->matchhead[i] = lz->matchnext[tmp];
            lz->matchnext[tmp] = 0;
        }
    }

    /*
     * If k is HASHCHARS or more, see if there were second-order
     * match possibilities, and if so move them into first-order
     * position.
     *
     * If there weren't, we must rewind the compressor state by a
     * few bytes so that we can reprocess those bytes and look for
     * new matches starting there. Because the window size is a few
     * bytes more than MAXMATCHDIST, this does not impact our
     * ability to find matches at the maximum distance even after
     * this rewinding.
     */
    if (lz->k >= HASHCHARS) {
        bool got_one = false;

        for (i = 0; i < HASHCHARS; i++) {
            int mi = besti*HASHCHARS+i;
            int m = lz->matchlater[mi];
            lz->literals[i] = lz->literals[HASHCHARS+besti*HASHCHARS+i];
            if (m > 0) {
                lz->matchhead[i] = lz->matchhead[m];
                lz->matchdist[i] = lz->matchdist[m];
                lz->matchlen[i] = lz->matchlen[m] - matchadjust[mi];
                lz->matchhead[m] = -1;
                lz->matchdist[m] = 0;
                lz->matchlen[m] = 0;
                if (lz->matchlen[i] > 0)
                    got_one = true;
            }
        }

        if (!got_one) {
            int rw = lz->k - (HASHCHARS-1);

            assert(lz->rewound + rw < MAXREWIND);

            lz->winpos = (lz->winpos + rw) % WINSIZE;
            lz->rewound += rw;
            lz->k -= rw;
            lz->nvalid -= rw;
        }
    }

    /*
     * Clean out the second-order part of the matchdist and
     * matchlen arrays, wipe any remaining lists out of
     * matchhead/matchnext, and reset the second-order match index
     * allocation.
     */
    for (i = HASHCHARS; i < NMATCHES; i++) {
        lz->matchdist[i] = lz->matchlen[i] = 0;
        while (lz->matchhead[i] >= 0) {
            int tmp = lz->matchhead[i];

            lz->matchhead[i] = lz->matchnext[tmp];
            lz->matchnext[tmp] = 0;
        }
    }
    for (i = 0; i < HASHCHARS*HASHCHARS; i++)
        lz->matchlater[i] = 0;
    lz->nextindex = HASHCHARS;
}

void lz77_flush(LZ77 *lz)
{
    while (lz->k >= HASHCHARS) {
        /*
         * Output a match.
         */
        lz77_outputmatch(lz);

        /*
         * In case we've rewound, call the main compress function
         * with length zero to reprocess the rewound data. Then we
         * might go back round this loop if that spotted a final
         * match after the one we output above.
         */
        lz77_compress(lz, NULL, 0);
    }

    assert(lz->k < HASHCHARS);
    /*
     * Now just output literals until we reduce k to zero.
     */
    while (lz->k > 0) {
        lz->k--;
        lz->literal(lz->ctx, lz->data[(lz->winpos + lz->k) % WINSIZE]);
    }
}

void lz77_compress(LZ77 *lz, const void *vdata, int len)
{
    const unsigned char *data = (const unsigned char *)vdata;

    while (lz->rewound > 0 || len > 0) {
        unsigned char currchars[HASHCHARS];
        int hash;
        bool rewound_this_time;
        int thissearchindex = -1;

        /*
         * First thing we do with every character we see: add it to
         * the sliding window, and shift the window round.
         */
        lz->winpos = (lz->winpos + WINSIZE-1) % WINSIZE;
        if (lz->rewound == 0) {
            len--;
            lz->data[lz->winpos] = *data++;
            if (lz->winpos < HASHCHARS)
                lz->data[lz->winpos + WINSIZE] = lz->data[lz->winpos];
            rewound_this_time = false;
        } else {
            lz->rewound--;
            rewound_this_time = true;
        }
        lz->k++;
        if (lz->nvalid < WINSIZE)
            lz->nvalid++;

        /*
         * If we're still within one hash length of the start of
         * the input data, there really is nothing more we can do
         * at all just yet.
         */
        if (lz->nvalid < HASHCHARS)
            continue;

        /*
         * Hash the most recent three characters.
         */
        {
            int i;
            for (i = 0; i < HASHCHARS; i++)
                currchars[i] = lz->data[lz->winpos + i];
            hash = lz77_hash(currchars);
        }

        /*
         * If k is within the range [HASHCHARS, HASHCHARS+MAXLAZY),
         * make a list of matches starting HASHCHARS bytes before
         * the end of the window, and store the head of the list in
         * matchhead[k-HASHCHARS].
         */
        if (lz->k >= HASHCHARS && lz->k < HASHCHARS+MAXLAZY) {
            thissearchindex = lz->k - HASHCHARS;
            lz77_hashsearch(lz, hash, currchars, thissearchindex);
        }

        if (lz->k >= HASHCHARS && lz->k < HASHCHARS+MAXLAZY-1) {
            /*
             * Save the literal at this position, in case we need it.
             */
            lz->literals[lz->k - HASHCHARS] =
                lz->data[lz->winpos + HASHCHARS-1];
        }

        /*
         * Alternatively, if we're currently close to the end of a
         * match we've been tracking, see if there's a match
         * starting here so we can take that into account when
         * making our lazy matching decision.
         */
        {
            int i, k;
            bool do_search = false;

            for (i = 0; i < MAXLAZY; i++) {
                if (lz->matchhead[i] < 0 && lz->matchlen[i] > 0) {
                    int distafter = lz->k - HASHCHARS - i - lz->matchlen[i];
                    if (distafter >= 0 && distafter < MAXLAZY-1-i) {
                        assert(lz->nextindex < NMATCHES);
                        assert(lz->matchlater[i*HASHCHARS+distafter] == 0);
                        lz->matchlater[i*HASHCHARS+distafter] = lz->nextindex;
                        do_search = true;
                    }
                    /*
                     * Also this is a good time to save any
                     * literals that come before this match.
                     */
                    if (distafter >= 0 && distafter < HASHCHARS) {
                        for (k = 0; k <= distafter; k++)
                            lz->literals[HASHCHARS + i*HASHCHARS+distafter-k] =
                            lz->data[(lz->winpos + HASHCHARS-1 + k) % WINSIZE];
                    }
                }
            }

            if (do_search) {
                thissearchindex = lz->nextindex++;
                lz77_hashsearch(lz, hash, currchars, thissearchindex);
            }
        }

        /*
         * We never let k get bigger than HASHCHARS unless there's
         * a match starting at it. So if k==HASHCHARS and there
         * isn't such a match, output a literal and decrement k.
         */
        if (lz->k == HASHCHARS && lz->matchhead[0] < 0) {
            lz->literal(lz->ctx, currchars[HASHCHARS-1]);
            lz->k--;
        }

        /*
         * For each list of matches we're currently tracking which
         * we _haven't_ only just started, go through the list and
         * winnow it to only those which have continued to match.
         */
        if (lz->k > HASHCHARS) {
            int i;
            for (i = 0; i < lz->nextindex; i++) {
                int outpos = -1, inpos = lz->matchhead[i], nextinpos;
                int bestdist = 0; 

                if (i == thissearchindex)
                    continue;          /* this is the one we've just started */

                while (inpos >= 0) {
                    nextinpos = lz->matchnext[inpos];

                    /*
                     * See if this match is still going.
                     */
                    if (
#ifdef MAXMATCHLEN
                        lz->matchlen[i] < MAXMATCHLEN &&
#endif

                        lz->data[(lz->winpos + inpos) % WINSIZE] ==
                        lz->data[lz->winpos]) {
                        /*
                         * It is; put it back on the winnowed list.
                         */
                        if (outpos < 0)
                            lz->matchhead[i] = inpos;
                        else
                            lz->matchnext[outpos] = inpos;
                        outpos = inpos;
                        /*
                         * Because we built up the match list in
                         * reverse order compared to the hash
                         * chain, we must prefer _later_ entries in
                         * the match list in order to prefer
                         * matching against the most recent data.
                         */
                        bestdist = inpos;
                    } else {
                        /*
                         * It isn't; mark it as unused in the
                         * matchnext array.
                         */
                        lz->matchnext[inpos] = 0;
                    }

                    inpos = nextinpos;
                }

                /*
                 * Terminate the new list.
                 */
                if (outpos < 0)
                    lz->matchhead[i] = -1;
                else
                    lz->matchnext[outpos] = -1;

                /*
                 * And update the distance/length tracker for this
                 * match.
                 */
                if (bestdist) {
                    lz->matchdist[i] = bestdist;
                    lz->matchlen[i]++;
                }
            }
        }

        /*
         * Add the current window position to the head of the
         * appropriate hash chain.
         * 
         * (If the compressor is still recovering from a rewind,
         * this position will _already_ be on its hash chain, so we
         * don't do this.)
         */
        if (!rewound_this_time) {
            lz->hashnext[lz->winpos] = lz->hashhead[hash];
            lz->hashhead[hash] = lz->winpos;
        }

        /*
         * If k >= HASHCHARS+MAXLAZY-1 (meaning that we've had a
         * chance to try a match at every viable starting point)
         * and all of the first-order ongoing matches end
         * HASHCHARS-1 bytes ago or more, we must output something.
         */
        if (lz->k >= HASHCHARS+MAXLAZY-1) {
            int i;
            for (i = 0; i < HASHCHARS; i++)
                if (lz->matchhead[i] >= 0 ||
                    lz->k-i-lz->matchlen[i] < HASHCHARS-1)
                    break;
            if (i == HASHCHARS)
                lz77_outputmatch(lz);
        }
    }
}

#ifdef LZ77_TESTMODE

/*
 * These tests are mostly constructed to check specific
 * match-finding situations. Strings of capital letters are
 * arranged so that they never repeat any sequence of three
 * characters, and hence they are incompressible by this algorithm.
 */

const char *const tests[] = {
    /*
     * Basics: an incompressible string to make sure we don't
     * compress it by mistake.
     */
    "AAABAACAADAAEAAFAAGAAHAAIAAJAAKAALAAMAANAAOAAPAAQAARAASAATAAUAAVAAWAAXAA",

    /*
     * Simple repeated string. Repeating the string three times
     * rather than two also checks that we prefer to match against
     * more recent data when there's no length difference.
     */
    "ZAabcBBABCABabcDABEABFabcABGABHA",
    "ZAabcdBBABCABabcdDABEABFabcdABGABHA",
    "ZAabcdeBBABCABabcdeDABEABFabcdeABGABHA",

    /*
     * Self-overlapping match.
     */
    "RABSABTAspoonspoonspoonspoonspoonspoonspoonspoonBUABVAB",

    /*
     * Lazy matching, one step on. `flapping' should be rendered as
     * a literal `f' plus a match against `lapping', rather than as
     * a match against `flap' followed by one against `ping'.
     */
    "ACCACDACflapEACFACGlappingACHACIAflappingCJACK",

    /*
     * Lazy matching, two steps on. (Transmitting `fl' as literals
     * is still superior to transmitting two matches.) Then
     * gradually reduce the length of the later match until it's no
     * longer profitable to use it.
     */
    "ACCACDACflapEACFACGapplianceACHACIAflapplianceCJACK",
    "ACCACDACflapEACFACGappliancACHACIAflapplianceCJACK",
    "ACCACDACflapEACFACGapplianACHACIAflapplianceCJACK",
    "ACCACDACflapEACFACGappliaACHACIAflapplianceCJACK",
    "ACCACDACflapEACFACGappliACHACIAflapplianceCJACK",
    "ACCACDACflapEACFACGapplACHACIAflapplianceCJACK",

    /*
     * Non-lazy matching, three steps on. (Transmitting the `fla'
     * as literals is _not_ superior to transmitting it as a match;
     * and in fact it's best to transmit the longest initial match
     * we can - `flap' - and then cover `athetic' with the next
     * match.)
     */
    "ACCACDACflapEACFACGpatheticACHACIAflapatheticCJACK",

    /*
     * Test that various kinds of match correctly find an
     * immediately following match in all circumstances.
     */
    "WAabcdeCXfghijACYACabcdefghijZAD",
    "WAabcdeCXbcdefACYACfghijZADabcdefghijBADCA",
    "WAabcdeCXbcdefgACYACfghijZADabcdefghijBADCA",
    "0WAabcdeCXbcdefACcdefghYACfghijklmZADabcdefghijklmBADCA",
    "2WAabcdeCXbcdefACcdefghiYACfghijklmZADabcdefghijklmBADCA",
    "1WAabcdeCXbcdefgACcdefghYACfghijklmZADabcdefghijklmBADCA",
    "0WAabcdefCXbcdefACcdefgYACfghijklmZADabcdefghijklmBADCA",
    "0WAabcdefCXbcdefgACcdefgYACfghijklmZADabcdefghijklmBADCA",
    "0WAabcdefCXbcdefACcdefghYACfghijklmZADabcdefghijklmBADCA",
    "1WAabcdeCXbcdefgACcdefghYACfghijklmZADabcdefghijklmBADCA",
    "0WAabcdefCXbcdefgACcdefghYACfghijklmZADabcdefghijklmBADCA",
    "WAabcdeCXbcdefgACcdefghiYACfghijZADabcdefghijBADCA",

    /*
     * Lazy matching: nasty cases in which it can be marginally
     * better _not_ to lazily match. In some of these cases,
     * choosing the superficially longer lazy match eliminates an
     * opportunity to render the entire final lower-case section
     * using one more match and at least three fewer literals.
     */
    "0WAabcdeCXbcdefghijklACYACfghijklmnoZADabcdefghijklmnoBADCA",
    "[01]WAabcdeCXbcdefghijklACYACghijklmnoZADabcdefghijklmnoBADCA",
    "0WAabcdeCXbcdefghijklACYACfghijklmnZADabcdefghijklmnBADCA",
    "1WAabcdeCXbcdefghijklACYACghijklmnZADabcdefghijklmnBADCA",
    "1WAabcdeCXbcdefghijklACYACfghijklmZADabcdefghijklmBADCA",
    "1WAabcdeCXbcdefghijklACYACghijklmZADabcdefghijklmBADCA",
    "1WAabcdCXbcdefACcdefghijkYACghijklmnZADabcdefghijklmnBADCA",
    "1WAabcdCXbcdefACcdefghijkYACfghijklmnZADabcdefghijklmnBADCA",
    "0WAabcdCXbcdefACcdefghijkYACefghijklmnZADabcdefghijklmnBADCA",
    "1WAabcdCXbcdefACcdefghijkYACghijklmZADabcdefghijklmBADCA",
    "[01]WAabcdCXbcdefACcdefghijkYACfghijklmZADabcdefghijklmBADCA",
    "0WAabcdCXbcdefACcdefghijkYACefghijklmZADabcdefghijklmBADCA",
    "2WAabcdCXbcdefACcdefghijkYACghijklmZADabcdefghijklBADCA",
    "2WAabcdCXbcdefACcdefghijkYACfghijklmZADabcdefghijklBADCA",
    "0WAabcdCXbcdefACcdefghijkYACefghijklmZADabcdefghijklBADCA",

    /*
     * Regression tests against specific things I've seen go wrong
     * in the past. All I really ask of these cases is that they
     * don't fail assertions; optimal compression is not critical.
     */
    "AabcBcdefgCdefDefghiEhijklFabcdefghijklG",
    "AabcBcdeCefgDfghijkEFabcdefghijklG",
    "AabcBbcdefgCcdefghiDhijklEjklmnopFabcdefghijklmnopqrstG",
    "AabcBbcdefghCcdefghijklmnopqrstuvDijklmnopqrstuvwxyzEdefghijklmnoF"
        "abcdefghijklmnopqrstuvwxyzG",
    "AabcdefBbcdeCcdefghijklmDfghijklmnoEFGHIJabcdefghijklmnoK",
    "AabcdBbcdCcdefDefgEabcdefgF",
    "AabcdeBbcdCcdefgDefgEabcdefghiF",

    /*
     * Fun final test.
     */
    "Pease porridge hot, pease porridge cold, pease porridge"
        " in the pot, nine days old.",
};

struct testctx {
    const char *data;
    int len, ptr;
};

void match(void *vctx, int distance, int len) {
    struct testctx *ctx = (struct testctx *)vctx;

    assert(distance > 0);
    assert(distance <= ctx->ptr);
    assert(len >= HASHCHARS);
    assert(len <= ctx->len - ctx->ptr);
    assert(len <= MAXMATCHLEN);
    assert(!memcmp(ctx->data + ctx->ptr, ctx->data + ctx->ptr - distance,
                   len));

    printf("<%d,%d>", distance, len);
    fflush(stdout);

    ctx->ptr += len;
}

void literal(void *vctx, unsigned char c) {
    struct testctx *ctx = (struct testctx *)vctx;

    assert(ctx->ptr < ctx->len);
    assert(c == (unsigned char)(ctx->data[ctx->ptr]));

    fputc(c, stdout);
    fflush(stdout);

    ctx->ptr++;
}

void dotest(const void *data, int len, int step)
{
    struct testctx t;
    LZ77 *lz;
    int j;

    t.data = data;
    t.len = len;
    t.ptr = 0;
    lz = lz77_new(literal, match, &t);
    for (j = 0; j < t.len; j += step)
        lz77_compress(lz, t.data + j, (t.len - j < step ? t.len - j : step));
    lz77_flush(lz);
    lz77_free(lz);
    assert(t.len == t.ptr);
    printf("\n");
}

int main(int argc, char **argv) {
    int i, len, truncate = 0;
    int step;
    char *filename = NULL;

    step = 48000;                      /* big step by default */

    while (--argc) {
        char *p = *++argv;
        if (!strcmp(p, "-t")) {
            truncate = 1;
        } else if (p[0] == '-' && p[1] == 'b') {
            step = atoi(p+2);          /* -bN sets block size to N */
        } else if (p[0] != '-') {
            filename = p;
        }
    }

    if (filename) {
        char *data = NULL;
        int datalen = 0, datasize = 0;
        int c;
        FILE *fp = fopen(filename, "rb");

        while ( (c = fgetc(fp)) != EOF) {
            if (datalen >= datasize) {
                datasize = (datalen * 3 / 2) + 512;
                data = realloc(data, datasize);
            }
            data[datalen++] = c;
        }

        fclose(fp);
        dotest(data, datalen, step);
        
    } else {
        for (i = 0; i < lenof(tests); i++) {
            for (len = (truncate ? 0 : strlen(tests[i]));
                 len <= strlen(tests[i]); len++) {
                dotest(tests[i], len, step);
            }
        }
    }

    return 0;
}

#endif
/* -------------------- originally from deflate.c -------------------- */

/*
 * Reimplementation of Deflate (RFC1951) compression. Adapted from
 * the version in PuTTY, and extended to write dynamic Huffman
 * trees and choose block boundaries usefully.
 * 
 * This source file is not complete: you also need to link it
 * against my lz77.c.
 */

/*
 * TODO:
 * 
 *  - Feature: could do with forms of flush other than SYNC_FLUSH.
 *    I'm not sure exactly how those work when you don't know in
 *    advance that your next block will be static (as we did in
 *    PuTTY). And remember the 9-bit limitation of zlib.
 *     + also, zlib has FULL_FLUSH which clears the LZ77 state as
 *       well, for random access.
 *
 *  - Compression quality: chooseblock() appears to be computing
 *    wildly inaccurate block size estimates. Possible resolutions:
 *     + find and fix some trivial bug I haven't spotted yet
 *     + abandon the entropic approximation and go with trial
 *       Huffman runs
 *
 *  - Compression quality: see if increasing SYMLIMIT causes
 *    dynamic blocks to start being consistently smaller than it.
 *     + actually we seem to be there already, but check on a
 *       larger corpus.
 *
 *  - Compression quality: we ought to be able to fall right back
 *    to actual uncompressed blocks if really necessary, though
 *    it's not clear what the criterion for doing so would be.
 */

/*
 * This software is copyright 2000-2006 Simon Tatham.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
 * IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#define snew(type) ( (type *) malloc(sizeof(type)) )
#define snewn(n, type) ( (type *) malloc((n) * sizeof(type)) )
#define sresize(x, n, type) ( (type *) realloc((x), (n) * sizeof(type)) )
#define sfree(x) ( free((x)) )

/* ----------------------------------------------------------------------
 * This file can be compiled in a number of modes.
 * 
 * With -DSTANDALONE, it builds a self-contained deflate tool which
 * can compress, decompress, and also analyse a deflated file to
 * print out the sequence of literals and copy commands it
 * contains.
 * 
 * With -DTESTMODE, it builds a test application which is given a
 * file on standard input, both compresses and decompresses it, and
 * outputs the re-decompressed result so it can be conveniently
 * diffed against the original. Define -DTESTDBG as well for lots
 * of diagnostics.
 */

#if defined TESTDBG
/* gcc-specific diagnostic macro */
#define debug_int(x...) ( fprintf(stderr, x) )
#define debug(x) ( debug_int x )
#else
#define debug(x) ((void)0)
#endif

#ifdef STANDALONE
#define ANALYSIS
#endif

#ifdef ANALYSIS
int analyse_level = 0;
#endif

/* ----------------------------------------------------------------------
 * Deflate functionality common to both compression and decompression.
 */

static const unsigned char lenlenmap[] = {
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
};

#define MAXCODELEN 16

/*
 * Given a sequence of Huffman code lengths, compute the actual
 * codes, in the final form suitable for feeding to outbits (i.e.
 * already bit-mirrored).
 *
 * Returns the maximum code length found. Can also return -1 to
 * indicate the table was overcommitted (too many or too short
 * codes to exactly cover the possible space), or -2 to indicate it
 * was undercommitted (too few or too long codes).
 */
static int hufcodes(const unsigned char *lengths, int *codes, int nsyms)
{
    int count[MAXCODELEN], startcode[MAXCODELEN];
    int code, maxlen;
    int i, j;

    /* Count the codes of each length. */
    maxlen = 0;
    for (i = 1; i < MAXCODELEN; i++)
        count[i] = 0;
    for (i = 0; i < nsyms; i++) {
        count[lengths[i]]++;
        if (maxlen < lengths[i])
            maxlen = lengths[i];
    }
    /* Determine the starting code for each length block. */
    code = 0;
    for (i = 1; i < MAXCODELEN; i++) {
        startcode[i] = code;
        code += count[i];
        if (code > (1 << i))
            maxlen = -1;               /* overcommitted */
        code <<= 1;
    }
    if (code < (1 << MAXCODELEN))
        maxlen = -2;                   /* undercommitted */
    /* Determine the code for each symbol. Mirrored, of course. */
    for (i = 0; i < nsyms; i++) {
        code = startcode[lengths[i]]++;
        codes[i] = 0;
        for (j = 0; j < lengths[i]; j++) {
            codes[i] = (codes[i] << 1) | (code & 1);
            code >>= 1;
        }
    }

    return maxlen;
}

/*
 * Adler32 checksum function.
 */
static unsigned long adler32_update(unsigned long s,
                                    const unsigned char *data, int len)
{
    unsigned s1 = s & 0xFFFF, s2 = (s >> 16) & 0xFFFF;
    int i;

    for (i = 0; i < len; i++) {
        s1 += data[i];
        s2 += s1;
        if (!(i & 0xFFF)) {
            s1 %= 65521;
            s2 %= 65521;
        }
    }

    return ((s2 % 65521) << 16) | (s1 % 65521);
}

typedef struct {
    short code, extrabits;
    int min, max;
} coderecord;

static const coderecord lencodes[] = {
    {257, 0, 3, 3},
    {258, 0, 4, 4},
    {259, 0, 5, 5},
    {260, 0, 6, 6},
    {261, 0, 7, 7},
    {262, 0, 8, 8},
    {263, 0, 9, 9},
    {264, 0, 10, 10},
    {265, 1, 11, 12},
    {266, 1, 13, 14},
    {267, 1, 15, 16},
    {268, 1, 17, 18},
    {269, 2, 19, 22},
    {270, 2, 23, 26},
    {271, 2, 27, 30},
    {272, 2, 31, 34},
    {273, 3, 35, 42},
    {274, 3, 43, 50},
    {275, 3, 51, 58},
    {276, 3, 59, 66},
    {277, 4, 67, 82},
    {278, 4, 83, 98},
    {279, 4, 99, 114},
    {280, 4, 115, 130},
    {281, 5, 131, 162},
    {282, 5, 163, 194},
    {283, 5, 195, 226},
    {284, 5, 227, 257},
    {285, 0, 258, 258},
};

static const coderecord distcodes[] = {
    {0, 0, 1, 1},
    {1, 0, 2, 2},
    {2, 0, 3, 3},
    {3, 0, 4, 4},
    {4, 1, 5, 6},
    {5, 1, 7, 8},
    {6, 2, 9, 12},
    {7, 2, 13, 16},
    {8, 3, 17, 24},
    {9, 3, 25, 32},
    {10, 4, 33, 48},
    {11, 4, 49, 64},
    {12, 5, 65, 96},
    {13, 5, 97, 128},
    {14, 6, 129, 192},
    {15, 6, 193, 256},
    {16, 7, 257, 384},
    {17, 7, 385, 512},
    {18, 8, 513, 768},
    {19, 8, 769, 1024},
    {20, 9, 1025, 1536},
    {21, 9, 1537, 2048},
    {22, 10, 2049, 3072},
    {23, 10, 3073, 4096},
    {24, 11, 4097, 6144},
    {25, 11, 6145, 8192},
    {26, 12, 8193, 12288},
    {27, 12, 12289, 16384},
    {28, 13, 16385, 24576},
    {29, 13, 24577, 32768},
};

/* ----------------------------------------------------------------------
 * Deflate compression.
 */

#define SYMLIMIT 65536
#define SYMPFX_LITLEN    0x00000000U
#define SYMPFX_DIST      0x40000000U
#define SYMPFX_EXTRABITS 0x80000000U
#define SYMPFX_CODELEN   0xC0000000U
#define SYMPFX_MASK      0xC0000000U

#define SYM_EXTRABITS_MASK 0x3C000000U
#define SYM_EXTRABITS_SHIFT 26

struct huftrees {
    unsigned char *len_litlen;
    int *code_litlen;
    unsigned char *len_dist;
    int *code_dist;
    unsigned char *len_codelen;
    int *code_codelen;
};

struct deflate_compress_ctx {
    LZ77 *lz;
    unsigned char *outbuf;
    int outlen, outsize;
    unsigned long outbits;
    int noutbits;
    bool firstblock;
    unsigned long *syms;
    int symstart, nsyms;
    int type;
    unsigned long checksum;
    unsigned long datasize;
    bool lastblock;
    bool finished;
    unsigned char static_len1[288], static_len2[30];
    int static_code1[288], static_code2[30];
    struct huftrees sht;
#ifdef STATISTICS
    unsigned long bitcount;
#endif
};

static void outbits(deflate_compress_ctx *out,
                    unsigned long bits, int nbits)
{
    assert(out->noutbits + nbits <= 32);
    out->outbits |= bits << out->noutbits;
    out->noutbits += nbits;
    while (out->noutbits >= 8) {
        if (out->outlen >= out->outsize) {
            out->outsize = out->outlen + 64;
            out->outbuf = sresize(out->outbuf, out->outsize, unsigned char);
        }
        out->outbuf[out->outlen++] = (unsigned char) (out->outbits & 0xFF);
        out->outbits >>= 8;
        out->noutbits -= 8;
    }
#ifdef STATISTICS
    out->bitcount += nbits;
#endif
}

/*
 * Binary heap functions used by buildhuf(). Each one assumes the
 * heap to be stored in an array of ints, with two ints per node
 * (user data and key). They take in the old heap length, and
 * return the new one.
 */
#define HEAPPARENT(x) (((x)-2)/4*2)
#define HEAPLEFT(x) ((x)*2+2)
#define HEAPRIGHT(x) ((x)*2+4)
static int addheap(int *heap, int len, int userdata, int key)
{
    int me, dad, tmp;

    me = len;
    heap[len++] = userdata;
    heap[len++] = key;

    while (me > 0) {
        dad = HEAPPARENT(me);
        if (heap[me+1] < heap[dad+1]) {
            tmp = heap[me]; heap[me] = heap[dad]; heap[dad] = tmp;
            tmp = heap[me+1]; heap[me+1] = heap[dad+1]; heap[dad+1] = tmp;
            me = dad;
        } else
            break;
    }

    return len;
}
static int rmheap(int *heap, int len, int *userdata, int *key)
{
    int me, lc, rc, c, tmp;

    len -= 2;
    *userdata = heap[0];
    *key = heap[1];
    heap[0] = heap[len];
    heap[1] = heap[len+1];

    me = 0;

    while (1) {
        lc = HEAPLEFT(me);
        rc = HEAPRIGHT(me);
        if (lc >= len)
            break;
        else if (rc >= len || heap[lc+1] < heap[rc+1])
            c = lc;
        else
            c = rc;
        if (heap[me+1] > heap[c+1]) {
            tmp = heap[me]; heap[me] = heap[c]; heap[c] = tmp;
            tmp = heap[me+1]; heap[me+1] = heap[c+1]; heap[c+1] = tmp;
        } else
            break;
        me = c;
    }

    return len;
}

/*
 * The core of the Huffman algorithm: takes an input array of
 * symbol frequencies, and produces an output array of code
 * lengths.
 *
 * This is basically a generic Huffman implementation, but it has
 * one zlib-related quirk which is that it caps the output code
 * lengths to fit in an unsigned char (which is safe since Deflate
 * will reject anything longer than 15 anyway). Anyone wanting to
 * rip it out and use it in another context should find that easy
 * to remove.
 */
#define HUFMAX 286
static void buildhuf(const int *freqs, unsigned char *lengths, int nsyms)
{
    int parent[2*HUFMAX-1];
    int length[2*HUFMAX-1];
    int heap[2*HUFMAX];
    int heapsize;
    int i, j, n;
    int si, sj;

    assert(nsyms <= HUFMAX);

    memset(parent, 0, sizeof(parent));

    /*
     * Begin by building the heap.
     */
    heapsize = 0;
    for (i = 0; i < nsyms; i++)
        if (freqs[i] > 0)              /* leave unused symbols out totally */
            heapsize = addheap(heap, heapsize, i, freqs[i]);

    /*
     * Now repeatedly take two elements off the heap and merge
     * them.
     */
    n = HUFMAX;
    while (heapsize > 2) {
        heapsize = rmheap(heap, heapsize, &i, &si);
        heapsize = rmheap(heap, heapsize, &j, &sj);
        parent[i] = n;
        parent[j] = n;
        heapsize = addheap(heap, heapsize, n, si + sj);
        n++;
    }

    /*
     * Now we have our tree, in the form of a link from each node
     * to the index of its parent. Count back down the tree to
     * determine the code lengths.
     */
    memset(length, 0, sizeof(length));
    /* The tree root has length 0 after that, which is correct. */
    for (i = n-1; i-- ;)
        if (parent[i] > 0)
            length[i] = 1 + length[parent[i]];

    /*
     * And that's it. (Simple, wasn't it?) Copy the lengths into
     * the output array and leave.
     * 
     * Here we cap lengths to fit in unsigned char.
     */
    for (i = 0; i < nsyms; i++)
        lengths[i] = (length[i] > 255 ? 255 : length[i]);
}

/*
 * Wrapper around buildhuf() which enforces the Deflate restriction
 * that no code length may exceed 15 bits, or 7 for the auxiliary
 * code length alphabet. This function has the same calling
 * semantics as buildhuf(), except that it might modify the freqs
 * array.
 */
static void deflate_buildhuf(int *freqs, unsigned char *lengths,
                             int nsyms, int limit)
{
    int smallestfreq, totalfreq, nactivesyms;
    int num, denom, adjust;
    int i;
    int maxprob;

    /*
     * Nasty special case: if the frequency table has fewer than
     * two non-zero elements, we must invent some, because we can't
     * have fewer than one bit encoding a symbol.
     */
    assert(nsyms >= 2);
    {
        int count = 0;
        for (i = 0; i < nsyms; i++)
            if (freqs[i] > 0)
                count++;
        if (count < 2) {
            for (i = 0; i < nsyms && count > 0; i++)
                if (freqs[i] == 0) {
                    freqs[i] = 1;
                    count--;
                }
        }
    }

    /*
     * First, try building the Huffman table the normal way. If
     * this works, it's optimal, so we don't want to mess with it.
     */
    buildhuf(freqs, lengths, nsyms);

    for (i = 0; i < nsyms; i++)
        if (lengths[i] > limit)
            break;

    if (i == nsyms)
        return;                        /* OK */

    /*
     * The Huffman algorithm can only ever generate a code length
     * of N bits or more if there is a symbol whose probability is
     * less than the reciprocal of the (N+2)th Fibonacci number
     * (counting from F_0=0 and F_1=1), i.e. 1/2584 for N=16, or
     * 1/55 for N=8. (This is a necessary though not sufficient
     * condition.)
     *
     * Why is this? Well, consider the input symbol with the
     * smallest probability. Let that probability be x. In order
     * for this symbol to have a code length of at least 1, the
     * Huffman algorithm will have to merge it with some other
     * node; and since x is the smallest probability, the node it
     * gets merged with must be at least x. Thus, the probability
     * of the resulting combined node will be at least 2x. Now in
     * order for our node to reach depth 2, this 2x-node must be
     * merged again. But what with? We can't assume the node it
     * merges with is at least 2x, because this one might only be
     * the _second_ smallest remaining node. But we do know the
     * node it merges with must be at least x, so our order-2
     * internal node is at least 3x.
     *
     * How small a node can merge with _that_ to get an order-3
     * internal node? Well, it must be at least 2x, because if it
     * was smaller than that then it would have been one of the two
     * smallest nodes in the previous step and been merged at that
     * point. So at least 3x, plus at least 2x, comes to at least
     * 5x for an order-3 node.
     *
     * And so it goes on: at every stage we must merge our current
     * node with a node at least as big as the bigger of this one's
     * two parents, and from this starting point that gives rise to
     * the Fibonacci sequence. So we find that in order to have a
     * node n levels deep (i.e. a maximum code length of n), the
     * overall probability of the root of the entire tree must be
     * at least F_{n+2} times the probability of the rarest symbol.
     * In other words, since the overall probability is 1, it is a
     * necessary condition for a code length of 16 or more that
     * there must be at least one symbol with probability <=
     * 1/F_18.
     *
     * (To demonstrate that a probability this big really can give
     * rise to a code length of 16, consider the set of input
     * frequencies { 1-epsilon, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55,
     * 89, 144, 233, 377, 610, 987 }, for arbitrarily small
     * epsilon.)
     *
     * So here buildhuf() has returned us an overlong code. So to
     * ensure it doesn't do it again, we add a constant to all the
     * (non-zero) symbol frequencies, causing them to become more
     * balanced and removing the danger. We can then feed the
     * results back to the standard buildhuf() and be
     * assert()-level confident that the resulting code lengths
     * contain nothing outside the permitted range.
     */
    assert(limit == 15 || limit == 7);
    maxprob = (limit == 15 ? 2584 : 55);   /* no point in computing full F_n */
    totalfreq = nactivesyms = 0;
    smallestfreq = -1;
    for (i = 0; i < nsyms; i++) {
        if (freqs[i] == 0)
            continue;
        if (smallestfreq < 0 || smallestfreq > freqs[i])
            smallestfreq = freqs[i];
        totalfreq += freqs[i];
        nactivesyms++;
    }
    assert(smallestfreq <= totalfreq / maxprob);

    /*
     * We want to find the smallest integer `adjust' such that
     * (totalfreq + nactivesyms * adjust) / (smallestfreq +
     * adjust) is less than maxprob. A bit of algebra tells us
     * that the threshold value is equal to
     *
     *   totalfreq - maxprob * smallestfreq
     *   ----------------------------------
     *          maxprob - nactivesyms
     *
     * rounded up, of course. And we'll only even be trying
     * this if
     */
    num = totalfreq - smallestfreq * maxprob;
    denom = maxprob - nactivesyms;
    adjust = (num + denom - 1) / denom;

    /*
     * Now add `adjust' to all the input symbol frequencies.
     */
    for (i = 0; i < nsyms; i++)
        if (freqs[i] != 0)
            freqs[i] += adjust;

    /*
     * Rebuild the Huffman tree...
     */
    buildhuf(freqs, lengths, nsyms);

    /*
     * ... and this time it ought to be OK.
     */
    for (i = 0; i < nsyms; i++)
        assert(lengths[i] <= limit);
}

/*
 * Compute the bit length of a symbol, given the three Huffman
 * trees.
 */
static int symsize(unsigned sym, const struct huftrees *trees)
{
    unsigned basesym = sym &~ SYMPFX_MASK;

    switch (sym & SYMPFX_MASK) {
      case SYMPFX_LITLEN:
        return trees->len_litlen[basesym];
      case SYMPFX_DIST:
        return trees->len_dist[basesym];
      case SYMPFX_CODELEN:
        return trees->len_codelen[basesym];
      default /*case SYMPFX_EXTRABITS*/:
        return basesym >> SYM_EXTRABITS_SHIFT;
    }
}

/*
 * Write out a single symbol, given the three Huffman trees.
 */
static void writesym(deflate_compress_ctx *out,
                     unsigned sym, const struct huftrees *trees)
{
    unsigned basesym = sym &~ SYMPFX_MASK;
    int i;

    switch (sym & SYMPFX_MASK) {
      case SYMPFX_LITLEN:
        debug(("send: litlen %d\n", basesym));
        outbits(out, trees->code_litlen[basesym], trees->len_litlen[basesym]);
        break;
      case SYMPFX_DIST:
        debug(("send: dist %d\n", basesym));
        outbits(out, trees->code_dist[basesym], trees->len_dist[basesym]);
        break;
      case SYMPFX_CODELEN:
        debug(("send: codelen %d\n", basesym));
        outbits(out, trees->code_codelen[basesym],trees->len_codelen[basesym]);
        break;
      case SYMPFX_EXTRABITS:
        i = basesym >> SYM_EXTRABITS_SHIFT;
        basesym &= ~SYM_EXTRABITS_MASK;
        debug(("send: extrabits %d/%d\n", basesym, i));
        outbits(out, basesym, i);
        break;
    }
}

/*
 * outblock() must output _either_ a dynamic block of length
 * `dynamic_len', _or_ a static block of length `static_len', but
 * it gets to choose which.
 */
static void outblock(deflate_compress_ctx *out,
                     int dynamic_len, int static_len)
{
    int freqs1[286], freqs2[30], freqs3[19];
    unsigned char len1[286], len2[30], len3[19];
    int code1[286], code2[30], code3[19];
    int hlit, hdist, hclen, bfinal, btype;
    int treesrc[286 + 30];
    int treesyms[286 + 30];
    int codelen[19];
    int i, ntreesrc, ntreesyms;
    bool dynamic;
    int blklen;
    struct huftrees dht;
    const struct huftrees *ht;
#ifdef STATISTICS
    unsigned long bitcount_before;
#endif

    dht.len_litlen = len1;
    dht.len_dist = len2;
    dht.len_codelen = len3;
    dht.code_litlen = code1;
    dht.code_dist = code2;
    dht.code_codelen = code3;

    /*
     * We make our choice of block to output by doing all the
     * detailed work to determine the exact length of each possible
     * block. Then we choose the one which has fewest output bits
     * per symbol.
     */

    /*
     * First build the two main Huffman trees for the dynamic
     * block.
     */

    /*
     * Count up the frequency tables.
     */
    memset(freqs1, 0, sizeof(freqs1));
    memset(freqs2, 0, sizeof(freqs2));
    freqs1[256] = 1;           /* we're bound to need one EOB */
    for (i = 0; i < dynamic_len; i++) {
        unsigned sym = out->syms[(out->symstart + i) % SYMLIMIT];

        /*
         * Increment the occurrence counter for this symbol, if
         * it's in one of the Huffman alphabets and isn't extra
         * bits.
         */
        if ((sym & SYMPFX_MASK) == SYMPFX_LITLEN) {
            sym &= ~SYMPFX_MASK;
            assert(sym < lenof(freqs1));
            freqs1[sym]++;
        } else if ((sym & SYMPFX_MASK) == SYMPFX_DIST) {
            sym &= ~SYMPFX_MASK;
            assert(sym < lenof(freqs2));
            freqs2[sym]++;
        }
    }
    deflate_buildhuf(freqs1, len1, lenof(freqs1), 15);
    deflate_buildhuf(freqs2, len2, lenof(freqs2), 15);
    hufcodes(len1, code1, lenof(freqs1));
    hufcodes(len2, code2, lenof(freqs2));

    /*
     * Determine HLIT and HDIST.
     */
    for (hlit = 286; hlit > 257 && len1[hlit-1] == 0; hlit--);
    for (hdist = 30; hdist > 1 && len2[hdist-1] == 0; hdist--);

    /*
     * Write out the list of symbols used to transmit the
     * trees.
     */
    ntreesrc = 0;
    for (i = 0; i < hlit; i++)
        treesrc[ntreesrc++] = len1[i];
    for (i = 0; i < hdist; i++)
        treesrc[ntreesrc++] = len2[i];
    ntreesyms = 0;
    for (i = 0; i < ntreesrc ;) {
        int j = 1;
        int k;

        /* Find length of run of the same length code. */
        while (i+j < ntreesrc && treesrc[i+j] == treesrc[i])
            j++;

        /* Encode that run as economically as we can. */
        k = j;
        if (treesrc[i] == 0) {
            /*
             * Zero code length: we can output run codes for
             * 3-138 zeroes. So if we have fewer than 3 zeroes,
             * we just output literals. Otherwise, we output
             * nothing but run codes, and tweak their lengths
             * to make sure we aren't left with under 3 at the
             * end.
             */
            if (k < 3) {
                while (k--)
                    treesyms[ntreesyms++] = 0 | SYMPFX_CODELEN;
            } else {
                while (k > 0) {
                    int rpt = (k < 138 ? k : 138);
                    if (rpt > k-3 && rpt < k)
                        rpt = k-3;
                    assert(rpt >= 3 && rpt <= 138);
                    if (rpt < 11) {
                        treesyms[ntreesyms++] = 17 | SYMPFX_CODELEN;
                        treesyms[ntreesyms++] =
                            (SYMPFX_EXTRABITS | (rpt - 3) |
                             (3 << SYM_EXTRABITS_SHIFT));
                    } else {
                        treesyms[ntreesyms++] = 18 | SYMPFX_CODELEN;
                        treesyms[ntreesyms++] =
                            (SYMPFX_EXTRABITS | (rpt - 11) |
                             (7 << SYM_EXTRABITS_SHIFT));
                    }
                    k -= rpt;
                }
            }
        } else {
            /*
             * Non-zero code length: we must output the first
             * one explicitly, then we can output a copy code
             * for 3-6 repeats. So if we have fewer than 4
             * repeats, we _just_ output literals. Otherwise,
             * we output one literal plus at least one copy
             * code, and tweak the copy codes to make sure we
             * aren't left with under 3 at the end.
             */
            assert(treesrc[i] < 16);
            treesyms[ntreesyms++] = treesrc[i] | SYMPFX_CODELEN;
            k--;
            if (k < 3) {
                while (k--)
                    treesyms[ntreesyms++] = treesrc[i] | SYMPFX_CODELEN;
            } else {
                while (k > 0) {
                    int rpt = (k < 6 ? k : 6);
                    if (rpt > k-3 && rpt < k)
                        rpt = k-3;
                    assert(rpt >= 3 && rpt <= 6);
                    treesyms[ntreesyms++] = 16 | SYMPFX_CODELEN;
                    treesyms[ntreesyms++] = (SYMPFX_EXTRABITS | (rpt - 3) |
                                             (2 << SYM_EXTRABITS_SHIFT));
                    k -= rpt;
                }
            }
        }

        i += j;
    }
    assert((unsigned)ntreesyms < lenof(treesyms));

    /*
     * Count up the frequency table for the tree-transmission
     * symbols, and build the auxiliary Huffman tree for that.
     */
    memset(freqs3, 0, sizeof(freqs3));
    for (i = 0; i < ntreesyms; i++) {
        unsigned sym = treesyms[i];

        /*
         * Increment the occurrence counter for this symbol, if
         * it's the Huffman alphabet and isn't extra bits.
         */
        if ((sym & SYMPFX_MASK) == SYMPFX_CODELEN) {
            sym &= ~SYMPFX_MASK;
            assert(sym < lenof(freqs3));
            freqs3[sym]++;
        }
    }
    deflate_buildhuf(freqs3, len3, lenof(freqs3), 7);
    hufcodes(len3, code3, lenof(freqs3));

    /*
     * Reorder the code length codes into transmission order, and
     * determine HCLEN.
     */
    for (i = 0; i < 19; i++)
        codelen[i] = len3[lenlenmap[i]];
    for (hclen = 19; hclen > 4 && codelen[hclen-1] == 0; hclen--)
        /* empty loop body */;

    /*
     * Now work out the exact size of both the dynamic and the
     * static block, in bits.
     */
    {
        int ssize, dsize;

        /*
         * First the dynamic block.
         */
        dsize = 3 + 5 + 5 + 4;         /* 3-bit header, HLIT, HDIST, HCLEN */
        dsize += 3 * hclen;            /* code-length-alphabet code lengths */
        /* Code lengths */
        for (i = 0; i < ntreesyms; i++)
            dsize += symsize(treesyms[i], &dht);
        /* The actual block data */
        for (i = 0; i < dynamic_len; i++) {
            unsigned sym = out->syms[(out->symstart + i) % SYMLIMIT];
            dsize += symsize(sym, &dht);
        }
        /* And the end-of-data symbol. */
        dsize += symsize(SYMPFX_LITLEN | 256, &dht);

        /*
         * Now the static block.
         */
        ssize = 3;                     /* 3-bit block header */
        /* The actual block data */
        for (i = 0; i < static_len; i++) {
            unsigned sym = out->syms[(out->symstart + i) % SYMLIMIT];
            ssize += symsize(sym, &out->sht);
        }
        /* And the end-of-data symbol. */
        ssize += symsize(SYMPFX_LITLEN | 256, &out->sht);

        /*
         * Compare the two and decide which to output. We break
         * exact ties in favour of the static block, because of the
         * special case in which that block has zero length.
         */
        dynamic = ((double)ssize * dynamic_len > (double)dsize * static_len);
        ht = dynamic ? &dht : &out->sht;
        blklen = dynamic ? dynamic_len : static_len;
    }

    /*
     * Actually transmit the block.
     */

    /* 3-bit block header */
    bfinal = (out->lastblock ? 1 : 0);
    btype = dynamic ? 2 : 1;
    debug(("send: bfinal=%d btype=%d\n", bfinal, btype));
    outbits(out, bfinal, 1);
    outbits(out, btype, 2);

#ifdef STATISTICS
    bitcount_before = out->bitcount;
#endif

    if (dynamic) {
        /* HLIT, HDIST and HCLEN */
        debug(("send: hlit=%d hdist=%d hclen=%d\n", hlit, hdist, hclen));
        outbits(out, hlit - 257, 5);
        outbits(out, hdist - 1, 5);
        outbits(out, hclen - 4, 4);

        /* Code lengths for the auxiliary tree */
        for (i = 0; i < hclen; i++) {
            debug(("send: lenlen %d\n", codelen[i]));
            outbits(out, codelen[i], 3);
        }

        /* Code lengths for the literal/length and distance trees */
        for (i = 0; i < ntreesyms; i++)
            writesym(out, treesyms[i], ht);
#ifdef STATISTICS
        fprintf(stderr, "total tree size %lu bits\n",
                out->bitcount - bitcount_before);
#endif
    }

    /* Output the actual symbols from the buffer */
    for (i = 0; i < blklen; i++) {
        unsigned sym = out->syms[(out->symstart + i) % SYMLIMIT];
        writesym(out, sym, ht);
    }

    /* Output the end-of-data symbol */
    writesym(out, SYMPFX_LITLEN | 256, ht);

    /*
     * Remove all the just-output symbols from the symbol buffer by
     * adjusting symstart and nsyms.
     */
    out->symstart = (out->symstart + blklen) % SYMLIMIT;
    out->nsyms -= blklen;
}

/*
 * Give the approximate log-base-2 of an input integer, measured in
 * 8ths of a bit. (I.e. this computes an integer approximation to
 * 8*logbase2(x).)
 */
static int approxlog2(unsigned x)
{
    int ret = 31*8;

    /*
     * Binary-search to get the top bit of x up to bit 31.
     */
    if (x < 0x00010000U) x <<= 16, ret -= 16*8;
    if (x < 0x01000000U) x <<=  8, ret -=  8*8;
    if (x < 0x10000000U) x <<=  4, ret -=  4*8;
    if (x < 0x40000000U) x <<=  2, ret -=  2*8;
    if (x < 0x80000000U) x <<=  1, ret -=  1*8;

    /*
     * Now we know the logarithm we want is in [ret,ret+1).
     * Determine the bottom three bits by checking against
     * threshold values.
     * 
     * (Each of these threshold values is 0x80000000 times an odd
     * power of 2^(1/16). Therefore, this function rounds to
     * nearest.)
     */
    if (x <= 0xAD583EEAU) {
        if (x <= 0x91C3D373U)
            ret += (x <= 0x85AAC367U ? 0 : 1);
        else
            ret += (x <= 0x9EF53260U ? 2 : 3);
    } else {
        if (x <= 0xCE248C15U)
            ret += (x <= 0xBD08A39FU ? 4 : 5);
        else
            ret += (x <= 0xE0CCDEECU ? 6 : x <= 0xF5257D15L ? 7 : 8);
    }

    return ret;
}

static void chooseblock(deflate_compress_ctx *out)
{
    int freqs1[286], freqs2[30];
    int i, len, bestlen, longestlen = 0;
    int total1, total2;
    int bestvfm;

    memset(freqs1, 0, sizeof(freqs1));
    memset(freqs2, 0, sizeof(freqs2));
    freqs1[256] = 1;                   /* we're bound to need one EOB */
    total1 = 1;
    total2 = 0;

    /*
     * Iterate over all possible block lengths, computing the
     * entropic coding approximation to the final length at every
     * stage. We divide the result by the number of symbols
     * encoded, to determine the `value for money' (overall
     * bits-per-symbol count) of a block of that length.
     */
    bestlen = -1;
    bestvfm = 0;

    len = 300 * 8;            /* very approximate size of the Huffman trees */

    for (i = 0; i < out->nsyms; i++) {
        unsigned sym = out->syms[(out->symstart + i) % SYMLIMIT];

        if (i > 0 && (sym & SYMPFX_MASK) == SYMPFX_LITLEN) {
            /*
             * This is a viable point at which to end the block.
             * Compute the value for money.
             */
            int vfm = i * 32768 / len;      /* symbols encoded per bit */

            if (bestlen < 0 || vfm > bestvfm) {
                bestlen = i;
                bestvfm = vfm;
            }

            longestlen = i;
        }

        /*
         * Increment the occurrence counter for this symbol, if
         * it's in one of the Huffman alphabets and isn't extra
         * bits.
         */
        if ((sym & SYMPFX_MASK) == SYMPFX_LITLEN) {
            sym &= ~SYMPFX_MASK;
            assert(sym < lenof(freqs1));
            len += freqs1[sym] * approxlog2(freqs1[sym]);
            len -= total1 * approxlog2(total1);
            freqs1[sym]++;
            total1++;
            len -= freqs1[sym] * approxlog2(freqs1[sym]);
            len += total1 * approxlog2(total1);
        } else if ((sym & SYMPFX_MASK) == SYMPFX_DIST) {
            sym &= ~SYMPFX_MASK;
            assert(sym < lenof(freqs2));
            len += freqs2[sym] * approxlog2(freqs2[sym]);
            len -= total2 * approxlog2(total2);
            freqs2[sym]++;
            total2++;
            len -= freqs2[sym] * approxlog2(freqs2[sym]);
            len += total2 * approxlog2(total2);
        } else if ((sym & SYMPFX_MASK) == SYMPFX_EXTRABITS) {
            len += 8 * ((sym &~ SYMPFX_MASK) >> SYM_EXTRABITS_SHIFT);
        }
    }

    assert(bestlen > 0);

    outblock(out, bestlen, longestlen);
}

/*
 * Force the current symbol buffer to be flushed out as a single
 * block.
 */
static void flushblock(deflate_compress_ctx *out)
{
    /*
     * No need to check that out->nsyms is a valid block length: we
     * know it has to be, because flushblock() is called in between
     * two matches/literals.
     */
    outblock(out, out->nsyms, out->nsyms);
    assert(out->nsyms == 0);
}

/*
 * Place a symbol into the symbols buffer.
 */
static void outsym(deflate_compress_ctx *out, unsigned long sym)
{
    assert(out->nsyms < SYMLIMIT);
    out->syms[(out->symstart + out->nsyms++) % SYMLIMIT] = sym;

    if (out->nsyms == SYMLIMIT)
        chooseblock(out);
}

static void literal(void *vctx, unsigned char c)
{
    deflate_compress_ctx *out = (deflate_compress_ctx *)vctx;

    outsym(out, SYMPFX_LITLEN | c);
}

static void match(void *vctx, int distance, int len)
{
    const coderecord *d, *l;
    int i, j, k;
    deflate_compress_ctx *out = (deflate_compress_ctx *)vctx;

    assert(len >= 3 && len <= 258);

    /*
     * Binary-search to find which length code we're
     * transmitting.
     */
    i = -1;
    j = sizeof(lencodes) / sizeof(*lencodes);
    while (1) {
        assert(j - i >= 2);
        k = (j + i) / 2;
        if (len < lencodes[k].min)
            j = k;
        else if (len > lencodes[k].max)
            i = k;
        else {
            l = &lencodes[k];
            break;                     /* found it! */
        }
    }

    /*
     * Transmit the length code.
     */
    outsym(out, SYMPFX_LITLEN | l->code);

    /*
     * Transmit the extra bits.
     */
    if (l->extrabits) {
        outsym(out, (SYMPFX_EXTRABITS | (len - l->min) |
                     (l->extrabits << SYM_EXTRABITS_SHIFT)));
    }

    /*
     * Binary-search to find which distance code we're
     * transmitting.
     */
    i = -1;
    j = sizeof(distcodes) / sizeof(*distcodes);
    while (1) {
        assert(j - i >= 2);
        k = (j + i) / 2;
        if (distance < distcodes[k].min)
            j = k;
        else if (distance > distcodes[k].max)
            i = k;
        else {
            d = &distcodes[k];
            break;                     /* found it! */
        }
    }

    /*
     * Write the distance code.
     */
    outsym(out, SYMPFX_DIST | d->code);

    /*
     * Transmit the extra bits.
     */
    if (d->extrabits) {
        outsym(out, (SYMPFX_EXTRABITS | (distance - d->min) |
                     (d->extrabits << SYM_EXTRABITS_SHIFT)));
    }
}

deflate_compress_ctx *deflate_compress_new(int type)
{
    deflate_compress_ctx *out;

    out = snew(deflate_compress_ctx);
    out->type = type;
    out->outbits = out->noutbits = 0;
    out->firstblock = true;
#ifdef STATISTICS
    out->bitcount = 0;
#endif

    out->syms = snewn(SYMLIMIT, unsigned long);
    out->symstart = out->nsyms = 0;

    out->checksum = (type == DEFLATE_TYPE_ZLIB ? 1 : 0);
    out->datasize = 0;
    out->lastblock = false;
    out->finished = false;

    out->lz = lz77_new(literal, match, out);

    /*
     * Build the static Huffman tables now, so we'll have them
     * available every time outblock() is called.
     */
    {
        int i;

        for (i = 0; i < (int)lenof(out->static_len1); i++)
            out->static_len1[i] = (i < 144 ? 8 :
                                   i < 256 ? 9 :
                                   i < 280 ? 7 : 8);
        for (i = 0; i < (int)lenof(out->static_len2); i++)
            out->static_len2[i] = 5;
    }
    hufcodes(out->static_len1, out->static_code1, lenof(out->static_code1));
    hufcodes(out->static_len2, out->static_code2, lenof(out->static_code2));
    out->sht.len_litlen = out->static_len1;
    out->sht.len_dist = out->static_len2;
    out->sht.len_codelen = NULL;
    out->sht.code_litlen = out->static_code1;
    out->sht.code_dist = out->static_code2;
    out->sht.code_codelen = NULL;

    return out;
}

void deflate_compress_free(deflate_compress_ctx *out)
{
    sfree(out->syms);
    lz77_free(out->lz);
    sfree(out);
}

void deflate_compress_data(deflate_compress_ctx *out,
                           const void *vblock, int len, int flushtype,
                           void **outblock, int *outlen)
{
    const unsigned char *block = (const unsigned char *)vblock;

    assert(!out->finished);

    out->outbuf = NULL;
    out->outlen = out->outsize = 0;

    /*
     * If this is the first block, output the header.
     */
    if (out->firstblock) {
        switch (out->type) {
          case DEFLATE_TYPE_BARE:
            break;                     /* no header */
          case DEFLATE_TYPE_ZLIB:
            /*
             * zlib (RFC1950) header bytes: 78 9C. (Deflate
             * compression, 32K window size, default algorithm.)
             */
            outbits(out, 0x9C78, 16);
            break;
          case DEFLATE_TYPE_GZIP:
            /*
             * Minimal gzip (RFC1952) header:
             * 
             *  - basic header of 1F 8B
             *  - compression method byte (8 = deflate)
             *  - flags byte (zero: we use no optional features)
             *  - modification time (zero: no time stamp available)
             *  - extra flags byte (2: we use maximum compression
             *    always)
             *  - operating system byte (255: we do not specify)
             */
            outbits(out, 0x00088B1F, 32);   /* header, CM, flags */
            outbits(out, 0, 32);       /* mtime */
            outbits(out, 0xFF02, 16);  /* xflags, OS */
            break;
        }
        out->firstblock = false;
    }

    /*
     * Feed our data to the LZ77 compression phase.
     */
    lz77_compress(out->lz, block, len);

    /*
     * Update checksums and counters.
     */
    switch (out->type) {
      case DEFLATE_TYPE_ZLIB:
        out->checksum = adler32_update(out->checksum, block, len);
        break;
      case DEFLATE_TYPE_GZIP:
        out->checksum = crc32_update(out->checksum, block, len);
        break;
    }
    out->datasize += len;

    switch (flushtype) {
        /*
         * FIXME: what other flush types are available and useful?
         * In PuTTY, it was clear that we generally wanted to be in
         * a static block so it was safe to open one. Here, we
         * probably prefer to be _outside_ a block if we can. Think
         * about this.
         */
      case DEFLATE_NO_FLUSH:
        break;                         /* don't flush any data at all (duh) */
      case DEFLATE_SYNC_FLUSH:
        /*
         * Flush the LZ77 compressor.
         */
        lz77_flush(out->lz);

        /*
         * Close the current block.
         */
        flushblock(out);

        /*
         * Then output an empty _uncompressed_ block: send 000,
         * then sync to byte boundary, then send bytes 00 00 FF
         * FF.
         */
        outbits(out, 0, 3);
        if (out->noutbits)
            outbits(out, 0, 8 - out->noutbits);
        outbits(out, 0, 16);
        outbits(out, 0xFFFF, 16);
        break;
      case DEFLATE_END_OF_DATA:
        /*
         * Flush the LZ77 compressor.
         */
        lz77_flush(out->lz);

        /*
         * Output a block with BFINAL set.
         */
        out->lastblock = true;
        flushblock(out);

        /*
         * Sync to byte boundary, flushing out the final byte.
         */
        if (out->noutbits)
            outbits(out, 0, 8 - out->noutbits);

        /*
         * Format-specific trailer data.
         */
        switch (out->type) {
          case DEFLATE_TYPE_ZLIB:
            /*
             * Just write out the Adler32 checksum.
             */
            outbits(out, (out->checksum >> 24) & 0xFF, 8);
            outbits(out, (out->checksum >> 16) & 0xFF, 8);
            outbits(out, (out->checksum >>  8) & 0xFF, 8);
            outbits(out, (out->checksum >>  0) & 0xFF, 8);
            break;
          case DEFLATE_TYPE_GZIP:
            /*
             * Write out the CRC32 checksum and the data length.
             */
            outbits(out, out->checksum, 32);
            outbits(out, out->datasize, 32);
            break;
        }

        out->finished = true;
        break;
    }

    /*
     * Return any data that we've generated.
     */
    *outblock = (void *)out->outbuf;
    *outlen = out->outlen;
}

/* ----------------------------------------------------------------------
 * Deflate decompression.
 */

/*
 * The way we work the Huffman decode is to have a table lookup on
 * the first N bits of the input stream (in the order they arrive,
 * of course, i.e. the first bit of the Huffman code is in bit 0).
 * Each table entry lists the number of bits to consume, plus
 * either an output code or a pointer to a secondary table.
 */
struct table;
struct tableentry;

struct tableentry {
    unsigned char nbits;
    short code;
    struct table *nexttable;
};

struct table {
    int mask;                          /* mask applied to input bit stream */
    struct tableentry *table;
};

#define MAXSYMS 288

#define DWINSIZE 32768

/*
 * Build a single-level decode table for elements
 * [minlength,maxlength) of the provided code/length tables, and
 * recurse to build subtables.
 */
static struct table *mkonetab(int *codes, unsigned char *lengths, int nsyms,
                              int pfx, int pfxbits, int bits)
{
    struct table *tab = snew(struct table);
    int pfxmask = (1 << pfxbits) - 1;
    int nbits, i, j, code;

    tab->table = snewn(1 << bits, struct tableentry);
    tab->mask = (1 << bits) - 1;

    for (code = 0; code <= tab->mask; code++) {
        tab->table[code].code = -1;
        tab->table[code].nbits = 0;
        tab->table[code].nexttable = NULL;
    }

    for (i = 0; i < nsyms; i++) {
        if (lengths[i] <= pfxbits || (codes[i] & pfxmask) != pfx)
            continue;
        code = (codes[i] >> pfxbits) & tab->mask;
        for (j = code; j <= tab->mask; j += 1 << (lengths[i] - pfxbits)) {
            tab->table[j].code = i;
            nbits = lengths[i] - pfxbits;
            if (tab->table[j].nbits < nbits)
                tab->table[j].nbits = nbits;
        }
    }
    for (code = 0; code <= tab->mask; code++) {
        if (tab->table[code].nbits <= bits)
            continue;
        /* Generate a subtable. */
        tab->table[code].code = -1;
        nbits = tab->table[code].nbits - bits;
        if (nbits > 7)
            nbits = 7;
        tab->table[code].nbits = bits;
        tab->table[code].nexttable = mkonetab(codes, lengths, nsyms,
                                              pfx | (code << pfxbits),
                                              pfxbits + bits, nbits);
    }

    return tab;
}

/*
 * Build a decode table, given a set of Huffman tree lengths.
 */
static struct table *mktable(unsigned char *lengths, int nlengths,
#ifdef ANALYSIS
                             const char *alphabet,
#endif
                             int *error)
{
    int codes[MAXSYMS];
    int maxlen;

#ifdef ANALYSIS
    if (alphabet && analyse_level > 1) {
        int i, col = 0;
        printf("code lengths for %s alphabet:\n", alphabet);
        for (i = 0; i < nlengths; i++) {
            col += printf("%3d", lengths[i]);
            if (col > 72) {
                putchar('\n');
                col = 0;
            }
        }
        if (col > 0)
            putchar('\n');
    }
#endif

    maxlen = hufcodes(lengths, codes, nlengths);

    if (maxlen < 0) {
        *error = (maxlen == -1 ? DEFLATE_ERR_LARGE_HUFTABLE :
                  DEFLATE_ERR_SMALL_HUFTABLE);
        return NULL;
    }

    /*
     * Now we have the complete list of Huffman codes. Build a
     * table.
     */
    return mkonetab(codes, lengths, nlengths, 0, 0, maxlen < 9 ? maxlen : 9);
}

static int freetable(struct table **ztab)
{
    struct table *tab;
    int code;

    if (ztab == NULL)
        return -1;

    if (*ztab == NULL)
        return 0;

    tab = *ztab;

    for (code = 0; code <= tab->mask; code++)
        if (tab->table[code].nexttable != NULL)
            freetable(&tab->table[code].nexttable);

    sfree(tab->table);
    tab->table = NULL;

    sfree(tab);
    *ztab = NULL;

    return (0);
}

struct deflate_decompress_ctx {
    struct table *staticlentable, *staticdisttable;
    struct table *currlentable, *currdisttable, *lenlentable;
    enum {
        ZLIBSTART,
        GZIPSTART, GZIPMETHFLAGS, GZIPIGNORE1, GZIPIGNORE2, GZIPIGNORE3,
        GZIPEXTRA, GZIPFNAME, GZIPCOMMENT,
        OUTSIDEBLK, TREES_HDR, TREES_LENLEN, TREES_LEN, TREES_LENREP,
        INBLK, GOTLENSYM, GOTLEN, GOTDISTSYM,
        UNCOMP_LEN, UNCOMP_NLEN, UNCOMP_DATA,
        END,
        ADLER1, ADLER2,
        CRC1, CRC2, ILEN1, ILEN2,
        FINALSPIN
    } state;
    int sym, hlit, hdist, hclen, lenptr, lenextrabits, lenaddon, len,
        lenrep;
    bool lastblock;
    int uncomplen;
    unsigned char lenlen[19];
    unsigned char lengths[286 + 32];
    unsigned long bits;
    int nbits;
    unsigned char window[DWINSIZE];
    int winpos;
    unsigned char *outblk;
    int outlen, outsize;
    int type;
    unsigned long checksum;
    unsigned long bytesout;
    int gzflags, gzextralen;
#ifdef ANALYSIS
    int bytesread;
    int bitcount_before;
#define BITCOUNT(dctx) ( (dctx)->bytesread * 8 - (dctx)->nbits )
#endif
};

deflate_decompress_ctx *deflate_decompress_new(int type)
{
    deflate_decompress_ctx *dctx = snew(deflate_decompress_ctx);
    unsigned char lengths[288];

    memset(lengths, 8, 144);
    memset(lengths + 144, 9, 256 - 144);
    memset(lengths + 256, 7, 280 - 256);
    memset(lengths + 280, 8, 288 - 280);
    dctx->staticlentable = mktable(lengths, 288,
#ifdef ANALYSIS
                                   NULL,
#endif
                                   NULL);
    assert(dctx->staticlentable);
    memset(lengths, 5, 32);
    dctx->staticdisttable = mktable(lengths, 32,
#ifdef ANALYSIS
                                    NULL,
#endif
                                    NULL);
    assert(dctx->staticdisttable);
    dctx->state = (type == DEFLATE_TYPE_ZLIB ? ZLIBSTART :
                   type == DEFLATE_TYPE_GZIP ? GZIPSTART :
                   OUTSIDEBLK);
    dctx->currlentable = dctx->currdisttable = dctx->lenlentable = NULL;
    dctx->bits = 0;
    dctx->nbits = 0;
    dctx->winpos = 0;
    dctx->type = type;
    dctx->lastblock = false;
    dctx->checksum = (type == DEFLATE_TYPE_ZLIB ? 1 : 0);
    dctx->bytesout = 0;
    dctx->gzflags = dctx->gzextralen = 0;
#ifdef ANALYSIS
    dctx->bytesread = dctx->bitcount_before = 0;
#endif

    return dctx;
}

void deflate_decompress_free(deflate_decompress_ctx *dctx)
{
    if (dctx->currlentable && dctx->currlentable != dctx->staticlentable)
        freetable(&dctx->currlentable);
    if (dctx->currdisttable && dctx->currdisttable != dctx->staticdisttable)
        freetable(&dctx->currdisttable);
    if (dctx->lenlentable)
        freetable(&dctx->lenlentable);
    freetable(&dctx->staticlentable);
    freetable(&dctx->staticdisttable);
    sfree(dctx);
}

static int huflookup(unsigned long *bitsp, int *nbitsp, struct table *tab)
{
    unsigned long bits = *bitsp;
    int nbits = *nbitsp;
    while (1) {
        struct tableentry *ent;
        ent = &tab->table[bits & tab->mask];
        if (ent->nbits > nbits)
            return -1;                 /* not enough data */
        bits >>= ent->nbits;
        nbits -= ent->nbits;
        if (ent->code == -1)
            tab = ent->nexttable;
        else {
            *bitsp = bits;
            *nbitsp = nbits;
            return ent->code;
        }

        /*
         * If we reach here with `tab' null, it can only be because
         * there was a missing entry in the Huffman table. This
         * should never occur even with invalid input data, because
         * we enforce at mktable time that the Huffman codes should
         * precisely cover the code space; so we can enforce this
         * by assertion.
         */
        assert(tab);
    }
}

static void emit_char(deflate_decompress_ctx *dctx, int c)
{
    dctx->window[dctx->winpos] = c;
    dctx->winpos = (dctx->winpos + 1) & (DWINSIZE - 1);
    if (dctx->outlen >= dctx->outsize) {
        dctx->outsize = dctx->outlen * 3 / 2 + 512;
        dctx->outblk = sresize(dctx->outblk, dctx->outsize, unsigned char);
    }
    if (dctx->type == DEFLATE_TYPE_ZLIB) {
        unsigned char uc = c;
        dctx->checksum = adler32_update(dctx->checksum, &uc, 1);
    } else if (dctx->type == DEFLATE_TYPE_GZIP) {
        unsigned char uc = c;
        dctx->checksum = crc32_update(dctx->checksum, &uc, 1);
    }
    dctx->outblk[dctx->outlen++] = c;
    dctx->bytesout++;
}

#define EATBITS(n) ( dctx->nbits -= (n), dctx->bits >>= (n) )

int deflate_decompress_data(deflate_decompress_ctx *dctx,
                            const void *vblock, int len,
                            void **outblock, int *outlen)
{
    const coderecord *rec;
    const unsigned char *block = (const unsigned char *)vblock;
    int code, bfinal, btype, rep, dist, nlen, header;
    unsigned long cksum;
    int error = 0;

    if (len == 0) {
        *outblock = NULL;
        *outlen = 0;
        if (dctx->state != FINALSPIN)
            return DEFLATE_ERR_UNEXPECTED_EOF;
        else
            return 0;
    }

    dctx->outblk = NULL;
    dctx->outsize = 0;
    dctx->outlen = 0;

    while (len > 0 || dctx->nbits > 0) {
        while (dctx->nbits < 24 && len > 0) {
            dctx->bits |= (*block++) << dctx->nbits;
            dctx->nbits += 8;
            len--;
#ifdef ANALYSIS
            dctx->bytesread++;
#endif
        }
        switch (dctx->state) {
          case ZLIBSTART:
            /* Expect 16-bit zlib header. */
            if (dctx->nbits < 16)
                goto finished;         /* done all we can */

            /*
             * The header is stored as a big-endian 16-bit integer,
             * in contrast to the general little-endian policy in
             * the rest of the format :-(
             */
            header = (((dctx->bits & 0xFF00) >> 8) |
                      ((dctx->bits & 0x00FF) << 8));
            EATBITS(16);

            /*
             * Check the header:
             *
             *  - bits 8-11 should be 1000 (Deflate/RFC1951)
             *  - bits 12-15 should be at most 0111 (window size)
             *  - bit 5 should be zero (no dictionary present)
             *  - we don't care about bits 6-7 (compression rate)
             *  - bits 0-4 should be set up to make the whole thing
             *    a multiple of 31 (checksum).
             */
            if ((header & 0xF000) >  0x7000 ||
                (header & 0x0020) != 0x0000 ||
                (header % 31) != 0) {
                error = DEFLATE_ERR_ZLIB_HEADER;
                goto finished;
            }
            if ((header & 0x0F00) != 0x0800) {
                error = DEFLATE_ERR_ZLIB_WRONGCOMP;
                goto finished;
            }
            dctx->state = OUTSIDEBLK;
            break;
          case GZIPSTART:
            /* Expect 16-bit gzip header. */
            if (dctx->nbits < 16)
                goto finished;
            header = dctx->bits & 0xFFFF;
            EATBITS(16);
            if (header != 0x8B1F) {
                error = DEFLATE_ERR_GZIP_HEADER;
                goto finished;
            }
            dctx->state = GZIPMETHFLAGS;
            break;
          case GZIPMETHFLAGS:
            /* Expect gzip compression method and flags bytes. */
            if (dctx->nbits < 16)
                goto finished;
            header = dctx->bits & 0xFF;
            EATBITS(8);
            if (header != 8) {
                error = DEFLATE_ERR_GZIP_WRONGCOMP;
                goto finished;
            }
            dctx->gzflags = dctx->bits & 0xFF;
            if (dctx->gzflags & 2) {
                /*
                 * The FHCRC flag is slightly confusing. RFC1952
                 * documents it as indicating the presence of a
                 * two-byte CRC16 of the gzip header, occurring
                 * just before the beginning of the Deflate stream.
                 * However, gzip itself (as of 1.3.5) appears to
                 * believe it indicates that the current gzip
                 * `member' is not the final one, i.e. that the
                 * stream is composed of multiple gzip members
                 * concatenated together, and furthermore gzip will
                 * refuse to decode any file that has it set.
                 * 
                 * For this reason, I label it as `disputed' and
                 * also refuse to decode anything that has it set.
                 * I don't expect this to be a problem in practice.
                 */
                error = DEFLATE_ERR_GZIP_FHCRC;
                goto finished;
            }
            EATBITS(8);
            dctx->state = GZIPIGNORE1;
            break;
          case GZIPIGNORE1:
          case GZIPIGNORE2:
          case GZIPIGNORE3:
            /* Expect two bytes of gzip timestamp/XFL/OS, which we ignore. */
            if (dctx->nbits < 16)
                goto finished;
            EATBITS(16);
            if (dctx->state == GZIPIGNORE3) {
                dctx->state = GZIPEXTRA;
            } else
                dctx->state++;         /* maps IGNORE1 -> IGNORE2 -> IGNORE3 */
            break;
          case GZIPEXTRA:
            if (dctx->gzflags & 4) {
                /* Expect two bytes of extra-length count, then that many
                 * extra bytes of header data, all of which we ignore. */
                if (!dctx->gzextralen) {
                    if (dctx->nbits < 16)
                        goto finished;
                    dctx->gzextralen = dctx->bits & 0xFFFF;
                    EATBITS(16);
                    break;
                } else if (dctx->gzextralen > 0) {
                    if (dctx->nbits < 8)
                        goto finished;
                    EATBITS(8);
                    if (--dctx->gzextralen > 0)
                        break;
                }
            }
            dctx->state = GZIPFNAME;
            break;
          case GZIPFNAME:
            if (dctx->gzflags & 8) {
                /*
                 * Expect a NUL-terminated filename.
                 */
                if (dctx->nbits < 8)
                    goto finished;
                code = dctx->bits & 0xFF;
                EATBITS(8);
            } else
                code = 0;
            if (code == 0)
                dctx->state = GZIPCOMMENT;
            break;
          case GZIPCOMMENT:
            if (dctx->gzflags & 16) {
                /*
                 * Expect a NUL-terminated filename.
                 */
                if (dctx->nbits < 8)
                    goto finished;
                code = dctx->bits & 0xFF;
                EATBITS(8);
            } else
                code = 0;
            if (code == 0)
                dctx->state = OUTSIDEBLK;
            break;
          case OUTSIDEBLK:
            /* Expect 3-bit block header. */
            if (dctx->nbits < 3)
                goto finished;         /* done all we can */
            bfinal = dctx->bits & 1;
            if (bfinal)
                dctx->lastblock = true;
            EATBITS(1);
            btype = dctx->bits & 3;
            EATBITS(2);
            if (btype == 0) {
                int to_eat = dctx->nbits & 7;
                dctx->state = UNCOMP_LEN;
                EATBITS(to_eat);       /* align to byte boundary */
            } else if (btype == 1) {
                dctx->currlentable = dctx->staticlentable;
                dctx->currdisttable = dctx->staticdisttable;
                dctx->state = INBLK;
            } else if (btype == 2) {
                dctx->state = TREES_HDR;
            }
            debug(("recv: bfinal=%d btype=%d\n", bfinal, btype));
#ifdef ANALYSIS
            if (analyse_level > 1) {
                static const char *const btypes[] = {
                    "uncompressed", "static", "dynamic", "type 3 (unknown)"
                };
                printf("new block, %sfinal, %s\n",
                       bfinal ? "" : "not ",
                       btypes[btype]);
            }
#endif
            break;
          case TREES_HDR:
            /*
             * Dynamic block header. Five bits of HLIT, five of
             * HDIST, four of HCLEN.
             */
            if (dctx->nbits < 5 + 5 + 4)
                goto finished;         /* done all we can */
            dctx->hlit = 257 + (dctx->bits & 31);
            EATBITS(5);
            dctx->hdist = 1 + (dctx->bits & 31);
            EATBITS(5);
            dctx->hclen = 4 + (dctx->bits & 15);
            EATBITS(4);
            debug(("recv: hlit=%d hdist=%d hclen=%d\n", dctx->hlit,
                   dctx->hdist, dctx->hclen));
#ifdef ANALYSIS
            if (analyse_level > 1)
                printf("hlit=%d, hdist=%d, hclen=%d\n",
                        dctx->hlit, dctx->hdist, dctx->hclen);
#endif
            dctx->lenptr = 0;
            dctx->state = TREES_LENLEN;
            memset(dctx->lenlen, 0, sizeof(dctx->lenlen));
            break;
          case TREES_LENLEN:
            if (dctx->nbits < 3)
                goto finished;
            while (dctx->lenptr < dctx->hclen && dctx->nbits >= 3) {
                dctx->lenlen[lenlenmap[dctx->lenptr++]] =
                    (unsigned char) (dctx->bits & 7);
                debug(("recv: lenlen %d\n", (unsigned char) (dctx->bits & 7)));
                EATBITS(3);
            }
            if (dctx->lenptr == dctx->hclen) {
                dctx->lenlentable = mktable(dctx->lenlen, 19,
#ifdef ANALYSIS
                                            "code length",
#endif
                                            &error);
                if (!dctx->lenlentable)
                    goto finished;     /* error code set up by mktable */
                dctx->state = TREES_LEN;
                dctx->lenptr = 0;
            }
            break;
          case TREES_LEN:
            if (dctx->lenptr >= dctx->hlit + dctx->hdist) {
                dctx->currlentable = mktable(dctx->lengths, dctx->hlit,
#ifdef ANALYSIS
                                             "literal/length",
#endif
                                             &error);
                if (!dctx->currlentable)
                    goto finished;     /* error code set up by mktable */
                if (dctx->hdist == 1 && dctx->lengths[dctx->hlit] == 0) {
                    /*
                     * Special case: if the code length list for the
                     * backward-distance table contains a single zero
                     * entry, it means this block will never encode a
                     * backward distance at all (i.e. it's all
                     * literals).
                     */
                    dctx->currdisttable = NULL;
                } else {
                    dctx->currdisttable = mktable(dctx->lengths + dctx->hlit,
                                                  dctx->hdist,
#ifdef ANALYSIS
                                                  "distance",
#endif
                                                  &error);
                    if (!dctx->currdisttable)
                        goto finished;     /* error code set up by mktable */
                }
                freetable(&dctx->lenlentable);
                dctx->lenlentable = NULL;
                dctx->state = INBLK;
                break;
            }
            code = huflookup(&dctx->bits, &dctx->nbits, dctx->lenlentable);
            debug(("recv: codelen %d\n", code));
            if (code == -1)
                goto finished;
            if (code < 16) {
#ifdef ANALYSIS
                if (analyse_level > 1)
                    printf("code-length %d\n", code);
#endif
                dctx->lengths[dctx->lenptr++] = code;
            } else {
                dctx->lenextrabits = (code == 16 ? 2 : code == 17 ? 3 : 7);
                dctx->lenaddon = (code == 18 ? 11 : 3);
                dctx->lenrep = (code == 16 && dctx->lenptr > 0 ?
                                dctx->lengths[dctx->lenptr - 1] : 0);
                dctx->state = TREES_LENREP;
            }
            break;
          case TREES_LENREP:
            if (dctx->nbits < dctx->lenextrabits)
                goto finished;
            rep =
                dctx->lenaddon +
                (dctx->bits & ((1 << dctx->lenextrabits) - 1));
            EATBITS(dctx->lenextrabits);
            if (dctx->lenextrabits)
                debug(("recv: codelen-extrabits %d/%d\n", rep - dctx->lenaddon,
                       dctx->lenextrabits));
#ifdef ANALYSIS
            if (analyse_level > 1)
                printf("code-length-repeat: %d copies of %d\n", rep,
                       dctx->lenrep);
#endif
            while (rep > 0 && dctx->lenptr < dctx->hlit + dctx->hdist) {
                dctx->lengths[dctx->lenptr] = dctx->lenrep;
                dctx->lenptr++;
                rep--;
            }
            dctx->state = TREES_LEN;
            break;
          case INBLK:
#ifdef ANALYSIS
            dctx->bitcount_before = BITCOUNT(dctx);
#endif
            code = huflookup(&dctx->bits, &dctx->nbits, dctx->currlentable);
            debug(("recv: litlen %d\n", code));
            if (code == -1)
                goto finished;
            if (code < 256) {
#ifdef ANALYSIS
                if (analyse_level > 0)
                    printf("%lu: literal %d [%d]\n", dctx->bytesout, code,
                           BITCOUNT(dctx) - dctx->bitcount_before);
#endif
                emit_char(dctx, code);
            } else if (code == 256) {
                if (dctx->lastblock)
                    dctx->state = END;
                else
                    dctx->state = OUTSIDEBLK;
                if (dctx->currlentable != dctx->staticlentable) {
                    freetable(&dctx->currlentable);
                    dctx->currlentable = NULL;
                }
                if (dctx->currdisttable &&
                    dctx->currdisttable != dctx->staticdisttable) {
                    freetable(&dctx->currdisttable);
                    dctx->currdisttable = NULL;
                }
            } else if (code < 286) {   /* static tree can give >285; ignore */
                dctx->state = GOTLENSYM;
                dctx->sym = code;
            }
            break;
          case GOTLENSYM:
            rec = &lencodes[dctx->sym - 257];
            if (dctx->nbits < rec->extrabits)
                goto finished;
            dctx->len =
                rec->min + (dctx->bits & ((1 << rec->extrabits) - 1));
            if (rec->extrabits)
                debug(("recv: litlen-extrabits %d/%d\n",
                       dctx->len - rec->min, rec->extrabits));
            EATBITS(rec->extrabits);
            dctx->state = GOTLEN;
            break;
          case GOTLEN:
            if (!dctx->currdisttable) {
                error = DEFLATE_ERR_NODISTTABLE;
                goto finished;
            }
            code = huflookup(&dctx->bits, &dctx->nbits, dctx->currdisttable);
            debug(("recv: dist %d\n", code));
            if (code == -1)
                goto finished;
            if (code >= 30) {
                error = DEFLATE_ERR_BADDISTCODE;
                goto finished;
            }
            dctx->state = GOTDISTSYM;
            dctx->sym = code;
            break;
          case GOTDISTSYM:
            rec = &distcodes[dctx->sym];
            if (dctx->nbits < rec->extrabits)
                goto finished;
            dist = rec->min + (dctx->bits & ((1 << rec->extrabits) - 1));
            if (rec->extrabits)
                debug(("recv: dist-extrabits %d/%d\n",
                       dist - rec->min, rec->extrabits));
            EATBITS(rec->extrabits);
            dctx->state = INBLK;
#ifdef ANALYSIS
            if (analyse_level > 0)
                printf("%lu: copy len=%d dist=%d [%d]\n", dctx->bytesout,
                       dctx->len, dist,
                       BITCOUNT(dctx) - dctx->bitcount_before);
#endif
            while (dctx->len--)
                emit_char(dctx, dctx->window[(dctx->winpos - dist) &
                                             (DWINSIZE - 1)]);
            break;
          case UNCOMP_LEN:
            /*
             * Uncompressed block. We expect to see a 16-bit LEN.
             */
            if (dctx->nbits < 16)
                goto finished;
            dctx->uncomplen = dctx->bits & 0xFFFF;
            EATBITS(16);
            dctx->state = UNCOMP_NLEN;
            break;
          case UNCOMP_NLEN:
            /*
             * Uncompressed block. We expect to see a 16-bit NLEN,
             * which should be the one's complement of the previous
             * LEN.
             */
            if (dctx->nbits < 16)
                goto finished;
            nlen = dctx->bits & 0xFFFF;
            EATBITS(16);
            if (dctx->uncomplen != (nlen ^ 0xFFFF)) {
                error = DEFLATE_ERR_UNCOMP_HDR;
                goto finished;
            }
            if (dctx->uncomplen == 0) {/* block is empty */
                if (dctx->lastblock)
                    dctx->state = END;
                else
                    dctx->state = OUTSIDEBLK;
            } else
                dctx->state = UNCOMP_DATA;
            break;
          case UNCOMP_DATA:
            if (dctx->nbits < 8)
                goto finished;
#ifdef ANALYSIS
            if (analyse_level > 0)
                printf("%lu: uncompressed %d [8]\n", dctx->bytesout,
                       (int)(dctx->bits & 0xFF));
#endif
            emit_char(dctx, dctx->bits & 0xFF);
            EATBITS(8);
            if (--dctx->uncomplen == 0) {       /* end of uncompressed block */
                if (dctx->lastblock)
                    dctx->state = END;
                else
                    dctx->state = OUTSIDEBLK;
            }
            break;
          case END:
            /*
             * End of compressed data. We align to a byte boundary,
             * and then look for format-specific trailer data.
             */
            {
                int to_eat = dctx->nbits & 7;
                EATBITS(to_eat);
            }
            if (dctx->type == DEFLATE_TYPE_ZLIB)
                dctx->state = ADLER1;
            else if (dctx->type == DEFLATE_TYPE_GZIP)
                dctx->state = CRC1;
            else
                dctx->state = FINALSPIN;
            break;
          case ADLER1:
            if (dctx->nbits < 16)
                goto finished;
            cksum = (dctx->bits & 0xFF) << 8;
            EATBITS(8);
            cksum |= (dctx->bits & 0xFF);
            EATBITS(8);
            if (cksum != ((dctx->checksum >> 16) & 0xFFFF)) {
                error = DEFLATE_ERR_CHECKSUM;
                goto finished;
            }
            dctx->state = ADLER2;
            break;
          case ADLER2:
            if (dctx->nbits < 16)
                goto finished;
            cksum = (dctx->bits & 0xFF) << 8;
            EATBITS(8);
            cksum |= (dctx->bits & 0xFF);
            EATBITS(8);
            if (cksum != (dctx->checksum & 0xFFFF)) {
                error = DEFLATE_ERR_CHECKSUM;
                goto finished;
            }
            dctx->state = FINALSPIN;
            break;
          case CRC1:
            if (dctx->nbits < 16)
                goto finished;
            cksum = dctx->bits & 0xFFFF;
            EATBITS(16);
            if (cksum != (dctx->checksum & 0xFFFF)) {
                error = DEFLATE_ERR_CHECKSUM;
                goto finished;
            }
            dctx->state = CRC2;
            break;
          case CRC2:
            if (dctx->nbits < 16)
                goto finished;
            cksum = dctx->bits & 0xFFFF;
            EATBITS(16);
            if (cksum != ((dctx->checksum >> 16) & 0xFFFF)) {
                error = DEFLATE_ERR_CHECKSUM;
                goto finished;
            }
            dctx->state = ILEN1;
            break;
          case ILEN1:
            if (dctx->nbits < 16)
                goto finished;
            cksum = dctx->bits & 0xFFFF;
            EATBITS(16);
            if (cksum != (dctx->bytesout & 0xFFFF)) {
                error = DEFLATE_ERR_INLEN;
                goto finished;
            }
            dctx->state = ILEN2;
            break;
          case ILEN2:
            if (dctx->nbits < 16)
                goto finished;
            cksum = dctx->bits & 0xFFFF;
            EATBITS(16);
            if (cksum != ((dctx->bytesout >> 16) & 0xFFFF)) {
                error = DEFLATE_ERR_INLEN;
                goto finished;
            }
            dctx->state = FINALSPIN;
            break;
          case FINALSPIN:
            /* Just ignore any trailing garbage on the data stream. */
            /* (We could alternatively throw an error here, if we wanted
             * to detect and complain about trailing garbage.) */
            EATBITS(dctx->nbits);
            break;
        }
    }

    finished:
    *outblock = dctx->outblk;
    *outlen = dctx->outlen;
    return error;
}

#define A(code,str) str
const char *const deflate_error_msg[DEFLATE_NUM_ERRORS] = {
    DEFLATE_ERRORLIST(A)
};
#undef A

#define A(code,str) #code
const char *const deflate_error_sym[DEFLATE_NUM_ERRORS] = {
    DEFLATE_ERRORLIST(A)
};
#undef A

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64)
#define WINDOWS_IO
#endif

#if defined(WINDOWS_IO) && (defined(STANDALONE) || defined(TESTMODE))
#endif

#ifdef STANDALONE

enum { DEFLATE_TYPE_AUTO = 3 };

int main(int argc, char **argv)
{
    unsigned char buf[65536];
    void *outbuf;
    int ret, err, outlen;
    deflate_decompress_ctx *dhandle;
    deflate_compress_ctx *chandle;
    int type = DEFLATE_TYPE_AUTO;
    bool opts = true, compress = false, decompress = false, got_arg = false;
    char *filename = NULL;
    FILE *fp;

    while (--argc) {
        char *p = *++argv;

        got_arg = TRUE;

        if (p[0] == '-' && opts) {
            if (!strcmp(p, "-b"))
                type = DEFLATE_TYPE_BARE;
            else if (!strcmp(p, "-g"))
                type = DEFLATE_TYPE_GZIP;
            else if (!strcmp(p, "-z"))
                type = DEFLATE_TYPE_ZLIB;
            else if (!strcmp(p, "-c"))
                compress = TRUE;
            else if (!strcmp(p, "-d"))
                decompress = TRUE;
            else if (!strcmp(p, "-a"))
                analyse_level++, decompress = true;
            else if (!strcmp(p, "--"))
                opts = false;          /* next thing is filename */
            else {
                fprintf(stderr, "unknown command line option '%s'\n", p);
                return 1;
            }
        } else if (!filename) {
            filename = p;
        } else {
            fprintf(stderr, "can only handle one filename\n");
            return 1;
        }
    }

    if (!compress && !decompress) {
        fprintf(stderr, "usage: deflate [ -c | -d | -a ] [ -b | -g | -z ]"
                " [filename]\n");
        return (got_arg ? 1 : 0);
    }

    if (compress && decompress) {
        fprintf(stderr, "please do not specify both compression and"
                " decompression\n");
        return (got_arg ? 1 : 0);
    }

    if (compress && type == DEFLATE_TYPE_AUTO) {
        fprintf(stderr, "please specify a file format for compression "
                "(-b, -g, -z)\n");
        return (got_arg ? 1 : 0);
    }

    if (compress) {
        chandle = deflate_compress_new(type);
        dhandle = NULL;
    } else {
        dhandle = NULL;
        chandle = NULL;
    }

    if (filename)
        fp = fopen(filename, "rb");
    else
        fp = stdin;

    if (!fp) {
        assert(filename);
        fprintf(stderr, "unable to open '%s'\n", filename);
        return 1;
    }
    
#ifdef WINDOWS_IO   
    if(_setmode(_fileno(stdout), _O_BINARY ) == -1)
    {
        fprintf(stderr, "Can't set stdout to binary mode\n");
        return 1;
    }
#endif

    do {
        ret = fread(buf, 1, sizeof(buf), fp);
        if (decompress && !dhandle) {
            if (type == DEFLATE_TYPE_AUTO) {
                /*
                 * Attempt to autodetect the input file type.
                 */
                if (ret >= 2 && buf[0] == 0x1F && buf[1] == 0x8B)
                    type = DEFLATE_TYPE_GZIP;
                else if (ret >= 2 && buf[0] == 0x78 &&
                         (buf[1] & 0x20) == 0 &&
                         (buf[0]*256+buf[1]) % 31 == 0)
                    type = DEFLATE_TYPE_ZLIB;
                else
                    type = DEFLATE_TYPE_BARE;
            }
            dhandle = deflate_decompress_new(type);
        }
        outbuf = NULL;
        if (dhandle) {
            if (ret > 0)
                err = deflate_decompress_data(dhandle, buf, ret,
                                              (void **)&outbuf, &outlen);
            else
                err = deflate_decompress_data(dhandle, NULL, 0,
                                              (void **)&outbuf, &outlen);
        } else {
            if (ret > 0)
                deflate_compress_data(chandle, buf, ret, DEFLATE_NO_FLUSH,
                                      (void **)&outbuf, &outlen);
            else
                deflate_compress_data(chandle, buf, ret, DEFLATE_END_OF_DATA,
                                      (void **)&outbuf, &outlen);
            err = 0;
        }
        if (outbuf) {
            if (!analyse_level && outlen)
                fwrite(outbuf, 1, outlen, stdout);
            sfree(outbuf);
        }
        if (err > 0) {
            fprintf(stderr, "decoding error: %s\n", deflate_error_msg[err]);
            return 1;
        }
    } while (ret > 0);

    if (dhandle)
        deflate_decompress_free(dhandle);
    if (chandle)
        deflate_compress_free(chandle);

    if (filename)
        fclose(fp);

    return 0;
}

#endif

#ifdef TESTMODE

int main(int argc, char **argv)
{
    char *filename = NULL;
    FILE *fp;
    deflate_compress_ctx *chandle;
    deflate_decompress_ctx *dhandle;
    unsigned char buf[65536], *outbuf, *outbuf2;
    int ret, err, outlen, outlen2;
    int dlen = 0, clen = 0;
    bool opts = true;

    while (--argc) {
        char *p = *++argv;

        if (p[0] == '-' && opts) {
            if (!strcmp(p, "--"))
                opts = false;          /* next thing is filename */
            else {
                fprintf(stderr, "unknown command line option '%s'\n", p);
                return 1;
            }
        } else if (!filename) {
            filename = p;
        } else {
            fprintf(stderr, "can only handle one filename\n");
            return 1;
        }
    }

    if (filename)
        fp = fopen(filename, "rb");
    else
        fp = stdin;

    if (!fp) {
        assert(filename);
        fprintf(stderr, "unable to open '%s'\n", filename);
        return 1;
    }

    chandle = deflate_compress_new(DEFLATE_TYPE_ZLIB);
    dhandle = deflate_decompress_new(DEFLATE_TYPE_ZLIB);
    
#ifdef WINDOWS_IO   
    if(_setmode(_fileno(stdout), _O_BINARY ) == -1)
    {
        fprintf(stderr, "Can't set stdout to binary mode\n");
        return 1;
    }
#endif

    do {
        ret = fread(buf, 1, sizeof(buf), fp);
        if (ret <= 0) {
            deflate_compress_data(chandle, NULL, 0, DEFLATE_END_OF_DATA,
                                  (void **)&outbuf, &outlen);
        } else {
            dlen += ret;
            deflate_compress_data(chandle, buf, ret, DEFLATE_NO_FLUSH,
                                  (void **)&outbuf, &outlen);
        }
        if (outbuf) {
            clen += outlen;
            err = deflate_decompress_data(dhandle, outbuf, outlen,
                                          (void **)&outbuf2, &outlen2);
            sfree(outbuf);
            if (outbuf2) {
                if (outlen2)
                    fwrite(outbuf2, 1, outlen2, stdout);
                sfree(outbuf2);
            }
            if (!err && ret <= 0) {
                /*
                 * signal EOF
                 */
                err = deflate_decompress_data(dhandle, NULL, 0,
                                              (void **)&outbuf2, &outlen2);
                assert(outbuf2 == NULL);
            }
            if (err) {
                fprintf(stderr, "decoding error: %s\n",
                        deflate_error_msg[err]);
                return 1;
            }
        }
    } while (ret > 0);

    fprintf(stderr, "%d plaintext -> %d compressed\n", dlen, clen);

    return 0;
}

#endif
/* -------------------- originally from crc32.c -------------------- */

/*
 * CRC32 checksum function.
 */

unsigned long crc32_update(unsigned long crcword, const void *vdata, int len)
{
    static const unsigned long crc32_table[256] = {
        0x00000000L, 0x77073096L, 0xEE0E612CL, 0x990951BAL,
        0x076DC419L, 0x706AF48FL, 0xE963A535L, 0x9E6495A3L,
        0x0EDB8832L, 0x79DCB8A4L, 0xE0D5E91EL, 0x97D2D988L,
        0x09B64C2BL, 0x7EB17CBDL, 0xE7B82D07L, 0x90BF1D91L,
        0x1DB71064L, 0x6AB020F2L, 0xF3B97148L, 0x84BE41DEL,
        0x1ADAD47DL, 0x6DDDE4EBL, 0xF4D4B551L, 0x83D385C7L,
        0x136C9856L, 0x646BA8C0L, 0xFD62F97AL, 0x8A65C9ECL,
        0x14015C4FL, 0x63066CD9L, 0xFA0F3D63L, 0x8D080DF5L,
        0x3B6E20C8L, 0x4C69105EL, 0xD56041E4L, 0xA2677172L,
        0x3C03E4D1L, 0x4B04D447L, 0xD20D85FDL, 0xA50AB56BL,
        0x35B5A8FAL, 0x42B2986CL, 0xDBBBC9D6L, 0xACBCF940L,
        0x32D86CE3L, 0x45DF5C75L, 0xDCD60DCFL, 0xABD13D59L,
        0x26D930ACL, 0x51DE003AL, 0xC8D75180L, 0xBFD06116L,
        0x21B4F4B5L, 0x56B3C423L, 0xCFBA9599L, 0xB8BDA50FL,
        0x2802B89EL, 0x5F058808L, 0xC60CD9B2L, 0xB10BE924L,
        0x2F6F7C87L, 0x58684C11L, 0xC1611DABL, 0xB6662D3DL,
        0x76DC4190L, 0x01DB7106L, 0x98D220BCL, 0xEFD5102AL,
        0x71B18589L, 0x06B6B51FL, 0x9FBFE4A5L, 0xE8B8D433L,
        0x7807C9A2L, 0x0F00F934L, 0x9609A88EL, 0xE10E9818L,
        0x7F6A0DBBL, 0x086D3D2DL, 0x91646C97L, 0xE6635C01L,
        0x6B6B51F4L, 0x1C6C6162L, 0x856530D8L, 0xF262004EL,
        0x6C0695EDL, 0x1B01A57BL, 0x8208F4C1L, 0xF50FC457L,
        0x65B0D9C6L, 0x12B7E950L, 0x8BBEB8EAL, 0xFCB9887CL,
        0x62DD1DDFL, 0x15DA2D49L, 0x8CD37CF3L, 0xFBD44C65L,
        0x4DB26158L, 0x3AB551CEL, 0xA3BC0074L, 0xD4BB30E2L,
        0x4ADFA541L, 0x3DD895D7L, 0xA4D1C46DL, 0xD3D6F4FBL,
        0x4369E96AL, 0x346ED9FCL, 0xAD678846L, 0xDA60B8D0L,
        0x44042D73L, 0x33031DE5L, 0xAA0A4C5FL, 0xDD0D7CC9L,
        0x5005713CL, 0x270241AAL, 0xBE0B1010L, 0xC90C2086L,
        0x5768B525L, 0x206F85B3L, 0xB966D409L, 0xCE61E49FL,
        0x5EDEF90EL, 0x29D9C998L, 0xB0D09822L, 0xC7D7A8B4L,
        0x59B33D17L, 0x2EB40D81L, 0xB7BD5C3BL, 0xC0BA6CADL,
        0xEDB88320L, 0x9ABFB3B6L, 0x03B6E20CL, 0x74B1D29AL,
        0xEAD54739L, 0x9DD277AFL, 0x04DB2615L, 0x73DC1683L,
        0xE3630B12L, 0x94643B84L, 0x0D6D6A3EL, 0x7A6A5AA8L,
        0xE40ECF0BL, 0x9309FF9DL, 0x0A00AE27L, 0x7D079EB1L,
        0xF00F9344L, 0x8708A3D2L, 0x1E01F268L, 0x6906C2FEL,
        0xF762575DL, 0x806567CBL, 0x196C3671L, 0x6E6B06E7L,
        0xFED41B76L, 0x89D32BE0L, 0x10DA7A5AL, 0x67DD4ACCL,
        0xF9B9DF6FL, 0x8EBEEFF9L, 0x17B7BE43L, 0x60B08ED5L,
        0xD6D6A3E8L, 0xA1D1937EL, 0x38D8C2C4L, 0x4FDFF252L,
        0xD1BB67F1L, 0xA6BC5767L, 0x3FB506DDL, 0x48B2364BL,
        0xD80D2BDAL, 0xAF0A1B4CL, 0x36034AF6L, 0x41047A60L,
        0xDF60EFC3L, 0xA867DF55L, 0x316E8EEFL, 0x4669BE79L,
        0xCB61B38CL, 0xBC66831AL, 0x256FD2A0L, 0x5268E236L,
        0xCC0C7795L, 0xBB0B4703L, 0x220216B9L, 0x5505262FL,
        0xC5BA3BBEL, 0xB2BD0B28L, 0x2BB45A92L, 0x5CB36A04L,
        0xC2D7FFA7L, 0xB5D0CF31L, 0x2CD99E8BL, 0x5BDEAE1DL,
        0x9B64C2B0L, 0xEC63F226L, 0x756AA39CL, 0x026D930AL,
        0x9C0906A9L, 0xEB0E363FL, 0x72076785L, 0x05005713L,
        0x95BF4A82L, 0xE2B87A14L, 0x7BB12BAEL, 0x0CB61B38L,
        0x92D28E9BL, 0xE5D5BE0DL, 0x7CDCEFB7L, 0x0BDBDF21L,
        0x86D3D2D4L, 0xF1D4E242L, 0x68DDB3F8L, 0x1FDA836EL,
        0x81BE16CDL, 0xF6B9265BL, 0x6FB077E1L, 0x18B74777L,
        0x88085AE6L, 0xFF0F6A70L, 0x66063BCAL, 0x11010B5CL,
        0x8F659EFFL, 0xF862AE69L, 0x616BFFD3L, 0x166CCF45L,
        0xA00AE278L, 0xD70DD2EEL, 0x4E048354L, 0x3903B3C2L,
        0xA7672661L, 0xD06016F7L, 0x4969474DL, 0x3E6E77DBL,
        0xAED16A4AL, 0xD9D65ADCL, 0x40DF0B66L, 0x37D83BF0L,
        0xA9BCAE53L, 0xDEBB9EC5L, 0x47B2CF7FL, 0x30B5FFE9L,
        0xBDBDF21CL, 0xCABAC28AL, 0x53B39330L, 0x24B4A3A6L,
        0xBAD03605L, 0xCDD70693L, 0x54DE5729L, 0x23D967BFL,
        0xB3667A2EL, 0xC4614AB8L, 0x5D681B02L, 0x2A6F2B94L,
        0xB40BBE37L, 0xC30C8EA1L, 0x5A05DF1BL, 0x2D02EF8DL
    };
    const unsigned char *data = (const unsigned char *)vdata;
    crcword ^= 0xFFFFFFFFL;
    while (len--) {
        unsigned long newbyte = *data++;
        newbyte ^= crcword & 0xFFL;
        crcword = (crcword >> 8) ^ crc32_table[newbyte];
    }
    return crcword ^ 0xFFFFFFFFL;
}
/* -------------------- originally from pngout.c -------------------- */

/*
 * PNG encoder.
 */

/*
 * TODO:
 *
 *  - Test thoroughly.
 *
 *  - tRNS can be shortened if the last pixels are non-transparent,
 *    so it might be nice to sort the palette to take account of
 *    that?
 */

#define PUT_32BIT_MSB_FIRST(cp, value) ( \
  (cp)[0] = (unsigned char)((value) >> 24), \
  (cp)[1] = (unsigned char)((value) >> 16), \
  (cp)[2] = (unsigned char)((value) >> 8), \
  (cp)[3] = (unsigned char)(value) )

struct pixel {
    unsigned short r, g, b, a;
};

struct interlace_pass {
    int xstart, xstep, ystart, ystep;
};

/*
 * With no interlacing, there is one pass over the image which
 * takes in all its pixels.
 */
#if 0 /* we don't actually have an option to use this at the moment */
static const struct interlace_pass null_interlace_passes[] = {
    {0, 1, 0, 1},
};
#endif

/*
 * Adam7 interlacing consists of seven passes which between them
 * specify every pixel exactly once.
 */
static const struct interlace_pass adam7_interlace_passes[] = {
    {0, 8, 0, 8},
    {4, 8, 0, 8},
    {0, 4, 4, 8},
    {2, 4, 0, 4},
    {0, 2, 2, 4},
    {1, 2, 0, 2},
    {0, 1, 1, 2},
};

static struct pixel read_pixel(const unsigned short **bitmap, int channels,
                               unsigned short mul)
{
    struct pixel ret;
    ret.r = *(*bitmap)++ * mul;
    if (channels < 3) {
        ret.g = ret.b = ret.r;
    } else {
        ret.g = *(*bitmap)++ * mul;
        ret.b = *(*bitmap)++ * mul;
    }
    if (channels & 1) {
        ret.a = 0xFFFF;
    } else {
        ret.a = *(*bitmap)++ * mul;
    }
    return ret;
}

static unsigned uabs(unsigned u)
{
    /* Return the absolute value of u, if u were interpreted as signed. */
    return u > (unsigned)INT_MAX ? -u : u;
}

/*
 * inputs:
 *  - w,h == size of image
 *  - channels = 1, 2, 3 or 4 (greyscale, grey+alpha, rgb,
 *    rgb+alpha)
 *  - depth = number of significant bits in each sample; must be 1,
 *    2, 4, 8 or 16
 *  - bitmap = array of h rows of w columns of 'channels' samples
 *    each ranging from 0 to (1<<depth)-1
 *
 * returns:
 *  - return value points to a dynamically allocated PNG in memory
 *  - *size gives the number of bytes in that PNG data
 *  - NULL return means failure to allocate memory
 */
void *pngwrite(int w, int h, int channels, int depth,
               const unsigned short *bitmap, size_t *size)
{
    struct pixel palette[256];
    int npalette;
    unsigned short mul = 0xFFFFU / ((1 << depth)-1);
    unsigned short d8, d4, d2, d1;
    const struct interlace_pass *passes;
    int npasses;
    int ochannels;
    unsigned filter_offset;
    unsigned pass_width[7], pass_height[7];
    unsigned max_width;
    size_t scandata_size, filesize;
    unsigned char *scandata, *p, *q, *filtered[5], *output;
    void *compressed;
    int complen;
    int i, j, colour, alpha, coltype, cdepth, bitdepth;
    const unsigned short *bp;
    deflate_compress_ctx *def;
    unsigned long crc;

    /*
     * Count the distinct pixel values, and build up a palette. At
     * the end of this loop we've either gone through the whole
     * image and found at most 256 pixel values which are contained
     * in 'palette', or else we've found more than 256 and have to
     * do a true-colour image.
     *
     * We also determine in this loop whether the image is
     * greyscale, whether it has a non-trivial alpha channel, and
     * (for use if we need to do true colour) what bit depth is
     * required.
     */
    colour = 0;
    alpha = 0;
    npalette = 0;
    d8 = d4 = d2 = d1 = 0;
    for (i = 0, bp = bitmap; i < w*h; i++) {
        /*
         * Read a pixel.
         */
        struct pixel pix = read_pixel(&bp, channels, mul);
        struct pixel p1, p2;

        /*
         * Check the static conditions: detect non-greyscale pixels
         * (their r, g, b components are not all identical), pixels
         * with interesting alpha values (i.e. not 0xFFFF), and
         * count the bit depth required to represent any given pixel
         * accurately.
         *
         * The last of these jobs - bit depth determination - is
         * done using some modular bit-twiddling trickery to save
         * time. The obvious approach is to observe that any 8-bit
         * sample converted to 16 bits will be an exact multiple of
         * 0x0101, so therefore we can determine if a 16-bit RGBA
         * tuple can be represented losslessly as an 8-bit one by
         * dividing each sample by 0x0101 and making sure the
         * remainder is zero. However, that involves four divisions,
         * and another four to check each of the 4-, 2- and 1-bit
         * depths, so that's all a bit unpleasant.
         *
         * A nicer approach is to use modular arithmetic. An 8-bit
         * sample converted to 16 bits will be the result of
         * multiplying a number in the range [0,255] by 0x0101
         * modulo 2^16. Multiplication by 0x0101 mod 2^16 is an
         * invertible operation (because 0x0101 is odd), and the
         * inverse operation is to multiply by 0xFF01 mod 2^16
         * (because 0xFF01 * 0x0101 = 0x1000001 which is congruent
         * to 1 mod 2^16). So we can take our sample, multiply by
         * 0xFF01, and check that no bits above the bottom 8 are
         * set. Likewise for the other depths: the multipliers for
         * 8-, 4-, 2- and 1- bit samples are respectively 0x0101,
         * 0x1111, 0x5555 and 0xFFFF, and their modular inverses are
         * 0xFF01, 0xFFF1, 0xFFFD and 0xFFFF respectively - also
         * known as -255, -15, -3 and -1. (Of course, for the last
         * of those, an actual multiplication would be overkill.)
         *
         * In fact, we streamline the process even further, by ORing
         * together the results of those multiplications across all
         * samples of all pixels, and then only checking for
         * unwanted set bits once at the very end of the loop. So
         * our variables d8, d4, d2 and d1 accumulate the bitwise
         * ORs of all those products.
         */
        if (!colour && (pix.r != pix.g || pix.r != pix.b))
            colour = 1;
        if (!alpha && pix.a != 0xFFFFU)
            alpha = 1;
        d8 |= ((pix.r * -255) | (pix.g * -255) |
               (pix.b * -255) | (pix.a * -255));
        d4 |= ((pix.r * -15) | (pix.g * -15) |
               (pix.b * -15) | (pix.a * -15));
        d2 |= ((pix.r * -3) | (pix.g * -3) |
               (pix.b * -3) | (pix.a * -3));
        d1 |= ((-pix.r) | (-pix.g) | (-pix.b) | (-pix.a));

        /*
         * Find this pixel value in the palette. We also shuffle the
         * palette upwards as we search, so that it's constantly
         * being move-to-front encoded. This hopefully arranges that
         * in typical images we take significantly less than
         * O(number of pixels * palette size) time to do this scan.
         */
        if (npalette <= 256) {
            p1 = pix;
            for (j = 0; j < npalette; j++) {
                p2 = palette[j];
                palette[j] = p1;
                if (p2.r == pix.r && p2.g == pix.g &&
                    p2.b == pix.b && p2.a == pix.a) {
                    break;
                }
                p1 = p2;
            }
            if (j == npalette) {
                if (npalette < 256)
                    palette[npalette] = p1;
                npalette++;
            }
        }
    }
    colour &= 1;                       /* normalise 'maybe' to 'no' */
    alpha &= 1;
    cdepth = (d8 & 0xFF00 ? 16 :
              d4 & 0xFFF0 ? 8 :
              d2 & 0xFFFC ? 4 :
              d1 & 0xFFFE ? 2 : 1);

    /*
     * Choose the PNG-level bit depth and colour type we'll use for
     * the real output.
     *
     * (Note that cdepth <= 8 is a second condition for using a
     * palette, because while PNG in general supports 16
     * bits/channel, its palette chunks only go up to 8 bits.)
     */
    if (npalette <= 256 && cdepth <= 8) {
        coltype = 3;                   /* paletted */
        bitdepth = (npalette <= 2 ? 1 :
                    npalette <= 4 ? 2 :
                    npalette <= 16 ? 4 : 8);
        /*
         * Special case: if we can do a pure greyscale image without
         * increasing the bit depth, do that instead, and save
         * storing a palette.
         */
        if (bitdepth == cdepth && !colour && !alpha)
            coltype = 0;
    } else {
        /* Colour types 0,2,4,6 mean grey,rgb,grey+alpha,rgba */
        coltype = 2*colour + 4*alpha;
        /* Only pure greyscale can go below 8 bit depth */
        if (coltype != 0 && cdepth < 8)
            bitdepth = 8;
        else
            bitdepth = cdepth;
    }
    ochannels = (coltype == 2 ? 3 :
                 coltype == 4 ? 2 :
                 coltype == 6 ? 4 : 1);

    /*
     * Decide on the division into scanlines (i.e. are we
     * interlacing?) and encode the unfiltered scanline data.
     *
     * Currently I fix on always interlacing, though it would be
     * easy enough to change that decision configurably.
     */
    passes = adam7_interlace_passes;
    npasses = lenof(adam7_interlace_passes);
    assert(npasses <= lenof(pass_width));
    /*
     * Compute the size of the unfiltered data. Each scanline in
     * each pass consists of the appropriate number of pixels,
     * multiplied by the number of bits per pixel, then rounded up
     * to the nearest byte, plus one prefixed byte for the filter
     * type.
     */
    scandata_size = 0;
    max_width = 0;
    for (i = 0; i < npasses; i++) {
        pass_width[i] =
            (w - passes[i].xstart + passes[i].xstep - 1) / passes[i].xstep;
        pass_width[i] *= ochannels * bitdepth;
        pass_width[i] += 7;
        pass_width[i] >>= 3;

        pass_height[i] =
            (h - passes[i].ystart + passes[i].ystep - 1) / passes[i].ystep;

        /* Empty passes don't even have the filter type byte */
        if (pass_width[i] && pass_height[i]) {
            pass_width[i]++;
            scandata_size += pass_width[i] * pass_height[i];
            if (max_width < pass_width[i])
                max_width = pass_width[i];
        }
    }
    scandata = (unsigned char *)malloc(scandata_size + 4 * max_width);
    if (!scandata)
        return NULL;                   /* couldn't allocate memory */
    /*
     * Now write out the unfiltered scanline data.
     */
    p = scandata;
    for (i = 0; i < npasses; i++) {
        int shift = 16 - bitdepth;
#ifndef NDEBUG
        unsigned char *q = p;          /* only used for assertion */
#endif
        int x, y;

        if (!pass_width[i] || !pass_height[i])
            continue;                  /* skip this pass completely */

        for (y = passes[i].ystart; y < h; y += passes[i].ystep) {
            unsigned int bits = 0, nbits = 0;
#ifndef NDEBUG
            unsigned char *r = p;      /* only used for assertion */
#endif
            *p++ = 0;                /* filter type 0: currently unfiltered */
            for (x = passes[i].xstart; x < w; x += passes[i].xstep) {
                struct pixel pix;
                unsigned short samples[4];

                bp = bitmap + ((size_t)y * w + x) * channels;
                pix = read_pixel(&bp, channels, mul);

                /*
                 * Convert our pixel into a list of samples.
                 */
                switch (coltype) {
                  case 0:
                    samples[0] = pix.r >> shift;
                    break;
                  case 2:
                    samples[0] = pix.r >> shift;
                    samples[1] = pix.g >> shift;
                    samples[2] = pix.b >> shift;
                    break;
                  case 4:
                    samples[0] = pix.r >> shift;
                    samples[1] = pix.a >> shift;
                    break;
                  case 6:
                    samples[0] = pix.r >> shift;
                    samples[1] = pix.g >> shift;
                    samples[2] = pix.b >> shift;
                    samples[3] = pix.a >> shift;
                    break;
                  case 3:
                    for (j = 0; j < npalette; j++) {
                        struct pixel p1 = palette[j];
                        if (p1.r == pix.r && p1.g == pix.g &&
                            p1.b == pix.b && p1.a == pix.a) {
                            break;
                        }
                    }
                    assert(j < npalette);
                    samples[0] = j;
                    break;
                }
                /*
                 * Write the samples out into the bit stream.
                 */
                for (j = 0; j < ochannels; j++) {
                    bits = (bits << bitdepth) | samples[j];
                    nbits += bitdepth;
                    while (nbits >= 8) {
                        *p++ = bits >> (nbits-8);
                        nbits -= 8;
                    }
                }
            }

            /*
             * Output a partial byte if we have one.
             */
            if (nbits > 0) {
                *p++ = bits << (8-nbits);
            }

            assert(p - r == pass_width[i]);
        }
        assert(p - q == pass_width[i] * pass_height[i]);
    }
    assert(p - scandata == scandata_size);

    /*
     * Filter the scanline data.
     *
     * We work backwards along the data we've just written (so that
     * when we filter each scanline its predecessor is still
     * available in its untransformed state). For each scanline, we
     * compute all five filtered forms, and choose the one which
     * gives the smallest mean absolute value of all image bytes.
     */
    filtered[0] = p;
    for (i = 1; i < 4; i++)
        filtered[i] = filtered[i-1] + max_width;
    filter_offset = (ochannels * bitdepth) >> 3;
    if (filter_offset == 0)
        filter_offset = 1;
    for (i = npasses; i-- > 0 ;) {
        unsigned x, y;

        if (!pass_width[i] || !pass_height[i])
            continue;                  /* skip this pass completely */

        for (y = passes[i].ystart; y < h; y += passes[i].ystep) {
            unsigned char *scan = (p -= pass_width[i]);
            unsigned char *prev = scan - pass_width[i];
            unsigned sumabs[5];
            int filt;

            /*
             * For the topmost scanline of the pass, there is no
             * previous scanline to use in the filtering..
             */
            if (y + passes[i].ystep >= h)
                prev = NULL;

            for (j = 0; j < 4; j++)
                filtered[j][0] = j+1;
            for (j = 0; j < 5; j++)
                sumabs[j] = 0;

            for (x = 1; x < pass_width[i]; x++) {
                unsigned orig = scan[x];
                unsigned a = (x >= filter_offset ?
                              scan[x-filter_offset] : 0);
                unsigned b = (prev ? prev[x] : 0);
                unsigned c = (prev && x >= filter_offset ?
                              prev[x-filter_offset] : 0);
                unsigned p, pa, pb, pc, pr;

                filtered[0][x] = orig - a;
                filtered[1][x] = orig - b;
                filtered[2][x] = orig - ((a + b) >> 1);
                p = a + b - c;
                pa = uabs(p - a);
                pb = uabs(p - b);
                pc = uabs(p - c);
                if (pa <= pb && pa <= pc)
                    pr = a;
                else if (pb <= pc)
                    pr = b;
                else
                    pr = c;
                filtered[3][x] = orig - pr;

                sumabs[0] += abs((signed char)scan[x]);
                sumabs[1] += abs((signed char)filtered[0][x]);
                sumabs[2] += abs((signed char)filtered[1][x]);
                sumabs[3] += abs((signed char)filtered[2][x]);
                sumabs[4] += abs((signed char)filtered[3][x]);
            }

            filt = 0;
            for (j = 0; j < 4; j++)
                if (sumabs[j+1] < sumabs[filt])
                    filt = j+1;
            if (filt > 0)
                memcpy(scan, filtered[filt-1], pass_width[i]);
        }
    }
    assert(p == scandata);

    /*
     * Compress the bitmap data.
     */
    def = deflate_compress_new(DEFLATE_TYPE_ZLIB);
    deflate_compress_data(def, scandata, scandata_size, DEFLATE_END_OF_DATA,
                          &compressed, &complen);
    deflate_compress_free(def);
    free(scandata);

    /*
     * Construct the full output file.
     */
    filesize = 8;                      /* PNG header */
    filesize += 12+13;                 /* IHDR chunk */
    if (coltype == 3) {
        filesize += 12 + 3 * npalette; /* PLTE chunk */
        if (alpha) {
            filesize += 12 + npalette; /* tRNS chunk */
        }
    }
    filesize += 12+complen;            /* IDAT chunk */
    filesize += 12;                    /* IEND chunk */
    output = (unsigned char *)malloc(filesize);
    if (!output) {
        free(compressed);
        return NULL;
    }

    memcpy(output, "\x89PNG\x0D\x0A\x1A\x0A", 8);   /* PNG header */
    p = output + 8;

    p += 8; q = p;
    memcpy(q-4, "IHDR", 4);
    PUT_32BIT_MSB_FIRST(p, w);
    PUT_32BIT_MSB_FIRST(p+4, h);
    p[8] = bitdepth;
    p[9] = coltype;
    p[10] = 0;                         /* fixed compression method: zlib */
    p[11] = 0;                         /* fixed filter method */
    p[12] = 1;                         /* fixed interlace method: adam7 */
    p += 13;
    crc = crc32_update(0, q-4, p-q+4);
    PUT_32BIT_MSB_FIRST(q-8, p-q);
    PUT_32BIT_MSB_FIRST(p, crc);
    p += 4;

    if (coltype == 3) {
        p += 8; q = p;
        memcpy(p-4, "PLTE", 4);
        for (i = 0; i < npalette; i++) {
            *p++ = palette[i].r >> 8;
            *p++ = palette[i].g >> 8;
            *p++ = palette[i].b >> 8;
        }
        crc = crc32_update(0, q-4, p-q+4);
        PUT_32BIT_MSB_FIRST(q-8, p-q);
        PUT_32BIT_MSB_FIRST(p, crc);
        p += 4;

        if (alpha) {
            p += 8; q = p;
            memcpy(p-4, "tRNS", 4);
            for (i = 0; i < npalette; i++) {
                *p++ = palette[i].a >> 8;
            }
            crc = crc32_update(0, q-4, p-q+4);
            PUT_32BIT_MSB_FIRST(q-8, p-q);
            PUT_32BIT_MSB_FIRST(p, crc);
            p += 4;
        }
    }

    p += 8; q = p;
    PUT_32BIT_MSB_FIRST(p, 13);
    memcpy(p-4, "IDAT", 4);
    memcpy(p, compressed, complen);
    free(compressed);
    p += complen;
    crc = crc32_update(0, q-4, p-q+4);
    PUT_32BIT_MSB_FIRST(q-8, p-q);
    PUT_32BIT_MSB_FIRST(p, crc);
    p += 4;

    p += 8; q = p;
    PUT_32BIT_MSB_FIRST(p, 13);
    memcpy(p-4, "IEND", 4);
    /* this chunk is empty */
    crc = crc32_update(0, q-4, p-q+4);
    PUT_32BIT_MSB_FIRST(q-8, p-q);
    PUT_32BIT_MSB_FIRST(p, crc);
    p += 4;

    assert(p - output == filesize);
    *size = filesize;
    return output;
}

#ifdef TEST_PNGOUT

unsigned short test[256*256*4];

void write(void *data, size_t size, char *name)
{
    FILE *fp = fopen(name, "w");
    fwrite(data, size, 1, fp);
    fclose(fp);
}

int main(void)
{
    int x, y;
    unsigned short *p;
    size_t size;
    void *data;

    /* Paletted, greyscale, no alpha. */
    p = test;
    for (y = 0; y < 256; y++)
        for (x = 0; x < 256; x++)
            *p++ = ((x >> 5) + (y >> 5)) & 1;
    data = pngwrite(256, 256, 1, 1, test, &size);
    write(data, size, "test1.png");

    /* Paletted, colour, no alpha. */
    p = test;
    for (y = 0; y < 256; y++)
        for (x = 0; x < 256; x++) {
            *p++ = ((x >> 5) + (y >> 5)) & 1;
            *p++ = ~((x >> 5) + (y >> 5)) & 1;
            *p++ = 0;
        }
    data = pngwrite(256, 256, 3, 1, test, &size);
    write(data, size, "test2.png");

    /* Paletted, greyscale, alpha. */
    p = test;
    for (y = 0; y < 256; y++)
        for (x = 0; x < 256; x++) {
            *p++ = (((x >> 5) + (y >> 5)) & 1) << 1;
            *p++ = 1;
        }
    data = pngwrite(256, 256, 2, 2, test, &size);
    write(data, size, "test3.png");

    /* Paletted, colour, alpha. */
    p = test;
    for (y = 0; y < 256; y++)
        for (x = 0; x < 256; x++) {
            *p++ = (((x >> 5) + (y >> 5)) & 1) << 1;
            *p++ = (~((x >> 5) + (y >> 5)) & 1) << 1;
            *p++ = 0;
            *p++ = 1;
        }
    data = pngwrite(256, 256, 4, 2, test, &size);
    write(data, size, "test4.png");

    /* True-colour, greyscale, no alpha. */
    p = test;
    for (y = 0; y < 256; y++)
        for (x = 0; x < 256; x++)
            *p++ = x*256+y;
    data = pngwrite(256, 256, 1, 16, test, &size);
    write(data, size, "test5.png");

    /* True-colour, colour, no alpha. */
    p = test;
    for (y = 0; y < 256; y++)
        for (x = 0; x < 256; x++) {
            *p++ = x;
            *p++ = y;
            *p++ = 0;
        }
    data = pngwrite(256, 256, 3, 8, test, &size);
    write(data, size, "test6.png");

    /* True-colour, greyscale, alpha. */
    p = test;
    for (y = 0; y < 256; y++)
        for (x = 0; x < 256; x++) {
            *p++ = x*256+y;
            *p++ = 0x8000;
        }
    data = pngwrite(256, 256, 2, 16, test, &size);
    write(data, size, "test7.png");

    /* True-colour, colour, alpha. */
    p = test;
    for (y = 0; y < 256; y++)
        for (x = 0; x < 256; x++) {
            *p++ = x;
            *p++ = y;
            *p++ = 0;
            *p++ = 0x80;
        }
    data = pngwrite(256, 256, 4, 8, test, &size);
    write(data, size, "test8.png");

    return 0;
}

#endif
/* -------------------- originally from cmdline.c -------------------- */

/* ----------------------------------------------------------------------
 * Generic command-line parsing functions.
 */

bool parsestr(char *string, void *vret) {
    char **ret = (char **)vret;
    *ret = string;
    return true;
}

bool parseint(char *string, void *vret) {
    int *ret = (int *)vret;
    int i;
    i = strspn(string, "0123456789");
    if (i > 0) {
        *ret = atoi(string);
        string += i;
    } else
        return false;
    if (*string)
        return false;
    return true;
}

bool parsesignedint(char *string, void *vret) {
    int *ret = (int *)vret;
    int i, j;
    j = 0;
    if (string[j] == '+' || string[j] == '-')
        j++;
    i = strspn(string + j, "0123456789");
    if (i > 0) {
        *ret = atoi(string);
        string += i+j;
    } else
        return false;
    if (*string)
        return false;
    return true;
}

bool parsesize(char *string, void *vret) {
    struct Size *ret = (struct Size *)vret;
    int i;
    i = strspn(string, "0123456789");
    if (i > 0) {
        ret->w = atoi(string);
        string += i;
    } else
        return false;
    if (*string++ != 'x')
        return false;
    i = strspn(string, "0123456789");
    if (i > 0) {
        ret->h = atoi(string);
        string += i;
    } else
        return false;
    if (*string)
        return false;
    return true;
}

bool parseflt(char *string, void *vret) {
    double *d = (double *)vret;
    char *endp;
    *d = strtod(string, &endp);
    if (endp && *endp) {
        return false;
    } else
        return true;
}

bool parsebool(char *string, void *vret) {
    bool *d = (bool *)vret;
    if (!strcmp(string, "yes") ||
        !strcmp(string, "true") ||
        !strcmp(string, "verily")) {
        *d = true;
        return true;
    } else if (!strcmp(string, "no") ||
               !strcmp(string, "false") ||
               !strcmp(string, "nowise")) {
        *d = false;
        return true;
    }
    return false;
}

/*
 * Read a colour into an RGB structure.
 */
bool parsecol(char *string, void *vret) {
    struct RGB *ret = (struct RGB *)vret;
    char *q;

    ret->r = strtod(string, &q);
    string = q;
    if (!*string) {
        return false;
    } else
        string++;
    ret->g = strtod(string, &q);
    string = q;
    if (!*string) {
        return false;
    } else
        string++;
    ret->b = strtod(string, &q);

    return true;
}

/*
 * 'Parsing' function used with flag options, to increment an integer.
 * Used for things like -v which increase verbosity more if you
 * specify them more times.
 */
bool incrementint(char *string, void *vret)
{
    int *pi = (int *)vret;
    assert(!string);
    (*pi)++;
    return true;
}

static void process_option(char const *programname,
                           const struct Cmdline *option,
                           char *arg, void *optdata)
{
    assert((arg != NULL) == (option->arghelp != NULL));

    if (!option->parse(arg, (char *)optdata + option->parse_ret_off)) {
        fprintf(stderr, "%s: unable to parse %s `%s'\n", programname,
                option->valname, arg);
        exit(EXIT_FAILURE);
    }
    if (option->gotflag_off >= 0)
        *(bool *)((char *)optdata + option->gotflag_off) = true;
}

void parse_cmdline(char const *programname, int argc, char **argv,
                   const struct Cmdline *options, int noptions, void *optdata)
{
    bool doing_options = true;
    int i;
    const struct Cmdline *argopt = NULL;

    for (i = 0; i < noptions; i++) {
        if (options[i].gotflag_off >= 0)
            *(bool *)((char *)optdata + options[i].gotflag_off) = false;
    }

    while (--argc > 0) {
        char *arg = *++argv;

        if (!strcmp(arg, "--")) {
            doing_options = false;
        } else if (!doing_options || arg[0] != '-') {
            if (!argopt) {
                for (i = 0; i < noptions; i++) {
                    if (!options[i].nlongopts && !options[i].shortopt) {
                        argopt = &options[i];
                        break;
                    }
                }

                if (!argopt) {
                    fprintf(stderr, "%s: no argument words expected\n",
                            programname);
                    exit(EXIT_FAILURE);
                }
            }

            process_option(programname, argopt, arg, optdata);
        } else {
            char c = arg[1];
            char *val = arg+2;
            int i, j, len;
            bool done = false;
            for (i = 0; i < noptions; i++) {
                /*
                 * Try a long option.
                 */
                if (c == '-') {
                    char *opt = options[i].longopt;
                    for (j = 0; j < options[i].nlongopts;
                         j++, opt += 1 + strlen(opt)) {
                         len = strlen(opt);
                        if (!strncmp(arg, opt, len) &&
                            len == strcspn(arg, "=")) {
                            if (options[i].arghelp) {
                                if (arg[len] == '=') {
                                    process_option(programname, &options[i],
                                                   arg + len + 1, optdata);
                                } else if (--argc > 0) {
                                    process_option(programname, &options[i],
                                                   *++argv, optdata);
                                } else {
                                    fprintf(stderr, "%s: option `%s' requires"
                                            " an argument\n", programname,
                                            arg);
                                    exit(EXIT_FAILURE);
                                }
                            } else {
                                process_option(programname, &options[i],
                                               NULL, optdata);
                            }
                            done = true;
                        }
                    }
                    if (!done) {
                        if (!strcmp(arg, "--version")) {
#ifdef VERSION
                            printf("%s version %s\n", programname, VERSION);
#else
                            printf("%s, unknown version\n", programname);
#endif
                            exit(EXIT_SUCCESS);
                        }
                    }
                } else if (c == options[i].shortopt) {
                    if (options[i].arghelp) {
                        if (*val) {
                            process_option(programname, &options[i], val,
                                           optdata);
                        } else if (--argc > 0) {
                            process_option(programname, &options[i], *++argv,
                                           optdata);
                        } else {
                            fprintf(stderr, "%s: option `%s' requires an"
                                    " argument\n", programname, arg);
                            exit(EXIT_FAILURE);
                        }
                    } else {
                        process_option(programname, &options[i], NULL,
                                       optdata);
                    }
                    done = true;
                }
            }
            if (!done && c == '-') {
                fprintf(stderr, "%s: unknown option `%.*s'\n",
                        programname, (int)strcspn(arg, "="), arg);
                exit(EXIT_FAILURE);
            }
        }
    }
}

void usage_message(char const *usageline,
                   const struct Cmdline *options, int noptions,
                   char **extratext, int nextra)
{
    int i, maxoptlen = 0;
    char *prefix;

    printf("usage: %s\n", usageline);

    /*
     * Work out the alignment for the help display.
     */
    for (i = 0; i < noptions; i++) {
        int optlen = 0;

        if (options[i].shortopt)
            optlen += 2;               /* "-X" */

        if (options[i].nlongopts) {
            if (optlen > 0)
                optlen += 2;           /* ", " */
            optlen += strlen(options[i].longopt);
        }

        if (options[i].arghelp) {
            if (optlen > 0)
                optlen++;              /* " " */
            optlen += strlen(options[i].arghelp);
        }

        if (maxoptlen < optlen)
            maxoptlen = optlen;
    }

    /*
     * Display the option help.
     */
    prefix = "where: ";
    for (i = 0; i < noptions; i++) {
        int optlen = 0;

        printf("%s", prefix);

        if (options[i].shortopt)
            optlen += printf("-%c", options[i].shortopt);

        if (options[i].nlongopts) {
            if (optlen > 0)
                optlen += printf(", ");
            optlen += printf("%s", options[i].longopt);
        }

        if (options[i].arghelp) {
            if (optlen > 0)
                optlen += printf(" ");
            optlen += printf("%s", options[i].arghelp);
        }

        printf("%*s  %s\n", maxoptlen-optlen, "", options[i].deschelp);

        prefix = "       ";
    }

    for (i = 0; i < nextra; i++) {
        printf("%s\n", extratext[i]);
    }

    exit(EXIT_FAILURE);
}
/* -------------------- originally from bmpwrite.c -------------------- */

/* ----------------------------------------------------------------------
 * Functions to write out .BMP files. Also supports PPM, because it
 * turns out to come in useful in some situations.
 */

static void fput32(unsigned long val, FILE *fp);
static void fput16(unsigned val, FILE *fp);

struct Bitmap {
    FILE *fp;
    bool close;
    int type;
    int width, height;
    unsigned long padding;
    unsigned short *data, *dataptr;
};

imagetype infer_type(const char *progname, const char *type,
                     const char *filename)
{
    if (!type) {
        /* Automatically infer type from file extension. */
        const char *dot = strrchr(filename, '.');
        if (dot && !strcmp(dot, ".png"))
            return PNG;
        if (dot && !strcmp(dot, ".ppm"))
            return PPM;
        if (dot && !strcmp(dot, ".bmp"))
            return BMP;
        fprintf(stderr, "%s: unable to infer image type from filename '%s';"
                " use --type={png,ppm,bmp}\n", progname, filename);
        exit(1);
    }
    if (!strcmp(type, "png"))
        return PNG;
    if (!strcmp(type, "ppm"))
        return PPM;
    if (!strcmp(type, "bmp"))
        return BMP;
    fprintf(stderr, "%s: unrecognised image type '%s';"
                " use --type={png,ppm,bmp}\n", progname, type);
    exit(1);
}

static void fput32(unsigned long val, FILE *fp) {
    fputc((val >>  0) & 0xFF, fp);
    fputc((val >>  8) & 0xFF, fp);
    fputc((val >> 16) & 0xFF, fp);
    fputc((val >> 24) & 0xFF, fp);
}
static void fput16(unsigned val, FILE *fp) {
    fputc((val >>  0) & 0xFF, fp);
    fputc((val >>  8) & 0xFF, fp);
}

struct Bitmap *bmpinit(char const *filename, int width, int height,
                       imagetype type)
{
    struct Bitmap *bm;

    bm = (struct Bitmap *)malloc(sizeof(struct Bitmap));
    bm->type = type;
    bm->width = width;
    bm->height = height;

    if (filename[0] == '-' && !filename[1]) {
        bm->fp = stdout;
        bm->close = false;
    } else {
        bm->fp = fopen(filename, "wb");
        bm->close = true;
    }

    if (type == BMP) {
        /*
         * BMP File format is:
         *
         * 2char "BM"
         * 32bit total file size
         * 16bit zero (reserved)
         * 16bit zero (reserved)
         * 32bit 0x36 (offset from start of file to image data)
         * 32bit 0x28 (size of following BITMAPINFOHEADER)
         * 32bit width
         * 32bit height
         * 16bit 0x01 (planes=1)
         * 16bit 0x18 (bitcount=24)
         * 32bit zero (no compression)
         * 32bit size of image data (total file size minus 0x36)
         * 32bit 0xB6D (XPelsPerMeter)
         * 32bit 0xB6D (YPelsPerMeter)
         * 32bit zero (palette colours used)
         * 32bit zero (palette colours important)
         *
         * then bitmap data, BGRBGRBGR... with padding zeros to bring
         * scan line to a multiple of 4 bytes. Padding zeros DO happen
         * after final scan line. Scan lines work from bottom upwards.
         */
        unsigned long scanlen, bitsize;

        scanlen = 3 * width;
        bm->padding = ((scanlen+3)&~3) - scanlen;
        bitsize = (scanlen + bm->padding) * height;

        fputs("BM", bm->fp);
        fput32(0x36 + bitsize, bm->fp);
        fput16(0, bm->fp); fput16(0, bm->fp);
        fput32(0x36, bm->fp); fput32(0x28, bm->fp);
        fput32(width, bm->fp); fput32(height, bm->fp);
        fput16(1, bm->fp); fput16(24, bm->fp);
        fput32(0, bm->fp); fput32(bitsize, bm->fp);
        fput32(0xB6D, bm->fp); fput32(0xB6D, bm->fp);
        fput32(0, bm->fp); fput32(0, bm->fp);
    } else if (type == PPM) {
        /*
         * PPM file format is:
         *
         * 2char "P6"
         * arbitrary whitespace
         * ASCII decimal width
         * arbitrary whitespace
         * ASCII decimal height
         * arbitrary whitespace
         * ASCII decimal maximum RGB value (here 255, for one-byte-per-sample)
         * one whitespace character
         *
         * then simply bitmap data, RGBRGBRGB...
         */
        fprintf(bm->fp, "P6 %d %d 255\n", width, height);

        bm->padding = 0;
    } else if (type == PNG) {
        bm->data = (unsigned short *)malloc(3 * width * height *
                                            sizeof(unsigned short));
        bm->dataptr = bm->data;
        bm->padding = 0;
    } else {
        assert(0 && "bad type");
    }

    return bm;
}

void bmppixel(struct Bitmap *bm,
              unsigned char r, unsigned char g, unsigned char b) {
    if (bm->type == BMP) {
        putc(b, bm->fp);
        putc(g, bm->fp);
        putc(r, bm->fp);
    } else if (bm->type == PPM) {
        putc(r, bm->fp);
        putc(g, bm->fp);
        putc(b, bm->fp);
    } else if (bm->type == PNG) {
        *bm->dataptr++ = r;
        *bm->dataptr++ = g;
        *bm->dataptr++ = b;
    } else {
        assert(0 && "bad type");
    }
}

void bmpendrow(struct Bitmap *bm) {
    int j;
    for (j = 0; j < bm->padding; j++)
        putc(0, bm->fp);
}

void bmpclose(struct Bitmap *bm) {
    if (bm->type == PNG) {
        size_t size;
        void *pngdata = pngwrite(bm->width, bm->height, 3, 8, bm->data, &size);
        fwrite(pngdata, 1, size, bm->fp);
        free(pngdata);
        free(bm->data);
    }
    if (bm->close)
        fclose(bm->fp);
    free(bm);
}
/* -------------------- originally from filigram.c -------------------- */

/* filigram.c - draw `filigram' pictures, based on the idea of plotting
 * the fractional part of a polynomial in x and y over a pixel grid
 *
 * This program is copyright 2000 Simon Tatham.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT.  IN NO EVENT SHALL SIMON TATHAM BE LIABLE FOR
 * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#define VERSION "20210403.80fbac2"

/* ----------------------------------------------------------------------
 * Function prototypes and structure type predeclarations.
 */

struct Poly;
static struct Poly *polyread(char const *string);
static struct Poly *polypdiff(struct Poly *input, bool wrtx);
static double polyeval(struct Poly *p, double x, double y);
static void polyfree(struct Poly *p);

struct Colours;
static struct Colours *colread(char const *string);
static void colfree(struct Colours *cols);
static struct RGB colfind(struct Colours *cols, int xg, int yg);

/* ----------------------------------------------------------------------
 * Routines for handling polynomials.
 */

struct PolyTerm {
    int constant;
    int xpower;
    int ypower;
    struct PolyTerm *next;
};

/*
 * A polynomial is stored as a 2-D array of PolyTerms, going up to
 * `maxdegree' in both cases. The nonzero terms are linked as a
 * list once the polynomial is finalised.
 */
struct Poly {
    int deg;
    struct PolyTerm *head;
    struct PolyTerm *terms;
};

/*
 * Build the linked list in a Poly structure.
 */
static void polymklist(struct Poly *p) {
    struct PolyTerm **tail;
    int i;
    tail = &p->head;
    for (i = 0; i < p->deg * p->deg; i++)
        if (p->terms[i].constant != 0) {
            *tail = &p->terms[i];
            tail = &p->terms[i].next;
        }
    *tail = NULL;
}

/*
 * Parse a string into a polynomial.
 *
 * String syntax must be: a number of terms, separated by + or -.
 * (Whitespace is fine, and is skipped.) Each term consists of an
 * optional integer, followed by an optional x power (x followed by
 * an integer), followed by an optional y power. For example:
 * 
 *   x2+2xy+y2
 * 
 * TODO if I ever get the energy: improve this to be a general
 * expression evaluator.
 */
static struct Poly *polyread(char const *string) {
    struct Poly *output;
    int sign, factor, xpower, ypower;
    int i, j, k;

    output = (struct Poly *)malloc(sizeof(struct Poly));
    output->deg = 0;
    output->terms = NULL;

    while (*string) {
        sign = +1;
        factor = 1;
        xpower = 0;
        ypower = 0;
        if (*string == '-') {
            sign = -1;
            string++;
        }
        i = strspn(string, "0123456789");
        if (i) {
            factor = atoi(string);
            string += i;
        }
        if (*string == 'x') {
            string++;
            xpower = 1;
            i = strspn(string, "0123456789");
            if (i) {
                xpower = atoi(string);
                string += i;
            }
        }
        if (*string == 'y') {
            string++;
            ypower = 1;
            i = strspn(string, "0123456789");
            if (i) {
                ypower = atoi(string);
                string += i;
            }
        }

        i = xpower;
        if (i < ypower) i = ypower;
        i++;
        if (i > output->deg) {
            int old;
            struct PolyTerm *newterms;

            old = output->deg;
            newterms = (struct PolyTerm *)malloc(i*i*sizeof(struct PolyTerm));
            for (j = 0; j < i; j++)
                for (k = 0; k < i; k++) {
                    newterms[k*i+j].xpower = j;
                    newterms[k*i+j].ypower = k;
                    if (j<old && k<old)
                        newterms[k*i+j].constant =
                        output->terms[k*old+j].constant;
                    else
                        newterms[k*i+j].constant = 0;
                }
            if (output->terms)
                free(output->terms);
            output->terms = newterms;
            output->deg = i;
        }
        output->terms[ypower*output->deg+xpower].constant += sign*factor;

        if (*string == '-') {
            /* do nothing */
        } else if (*string == '+') {
            string++;                  /* skip the plus */
        } else if (*string) {
            polyfree(output);
            return NULL;               /* error! */
        }
    }

    polymklist(output);
    return output;
}

/*
 * Partially differentiate a polynomial with respect to x (if
 * wrtx is true) or to y (if wrtx is false).
 */
static struct Poly *polypdiff(struct Poly *input, bool wrtx) {
    struct Poly *output;
    int deg, yp, xp, ydp, xdp, factor, constant;

    output = (struct Poly *)malloc(sizeof(struct Poly));
    deg = output->deg = input->deg;
    output->terms = (struct PolyTerm *)malloc(output->deg * output->deg *
                                              sizeof(struct PolyTerm));

    for (yp = 0; yp < deg; yp++)
        for (xp = 0; xp < deg; xp++) {
            xdp = xp; ydp = yp;
            factor = wrtx ? ++xdp : ++ydp;
            if (xdp >= deg || ydp >= deg)
                constant = 0;
            else
                constant = input->terms[ydp*deg+xdp].constant;
            output->terms[yp*deg+xp].xpower = xp;
            output->terms[yp*deg+xp].ypower = yp;
            output->terms[yp*deg+xp].constant = constant * factor;
        }

    polymklist(output);
    return output;
}

/*
 * Evaluate a polynomial.
 * 
 * TODO: make this more efficient by pre-computing the powers once
 * instead of doing them all, linearly, per term.
 */
static double polyeval(struct Poly *p, double x, double y) {
    struct PolyTerm *pt;
    double result, term;
    int i;

    result = 0.0;
    for (pt = p->head; pt; pt = pt->next) {
        term = pt->constant;
        for (i = 0; i < pt->xpower; i++) term *= x;
        for (i = 0; i < pt->ypower; i++) term *= y;
        result += term;
    }
    return result;
}

/*
 * Free a polynomial when finished.
 */
static void polyfree(struct Poly *p) {
    free(p->terms);
    free(p);
}

/* ----------------------------------------------------------------------
 * Routines for handling lists of colours.
 */

struct Colours {
    int xmod;
    int ncolours;
    struct RGB *colours;
};

/*
 * Read a colour list into a Colours structure.
 */
static struct Colours *colread(char const *string) {
    struct Colours *ret;
    int i;
    struct RGB c, *clist;
    int nc, csize;
    char *q;

    ret = (struct Colours *)malloc(sizeof(struct Colours));
    ret->xmod = ret->ncolours = 0;
    ret->colours = NULL;
    nc = csize = 0;
    clist = NULL;

    /* If there's an `integer+' prefix, set xmod. */
    i = strspn(string, "0123456789");
    if (string[i] == '+') {
        ret->xmod = atoi(string);
        string += i + 1;
    }

    while (*string) {
        /* Now we expect triplets of doubles separated by any punctuation. */
        c.r = strtod(string, &q); string = q; if (!*string) break; string++;
        c.g = strtod(string, &q); string = q; if (!*string) break; string++;
        c.b = strtod(string, &q); string = q;
        if (nc >= csize) {
            csize = nc + 16;
            clist = realloc(clist, csize * sizeof(struct RGB));
        }
        clist[nc++] = c;
        if (!*string)
            break;
        string++;
    }
    ret->ncolours = nc;
    ret->colours = clist;
    return ret;
}

/*
 * Free a Colours structure when done.
 */
static void colfree(struct Colours *cols) {
    free(cols->colours);
    free(cols);
}

/*
 * Look up a colour in a Colours structure.
 */
static struct RGB colfind(struct Colours *cols, int xg, int yg) {
    int n;
    if (cols->xmod) {
        xg %= cols->xmod;
        if (xg < 0) xg += cols->xmod;
        yg *= cols->xmod;
    }
    n = (xg+yg) % cols->ncolours;
    if (n < 0) n += cols->ncolours;
    return cols->colours[n];
}

/* ----------------------------------------------------------------------
 * The code to do the actual plotting.
 */

struct Params {
    int width, height;
    double x0, x1, y0, y1;
    double xscale, yscale, oscale;
    bool fading;
    char const *filename;
    int outtype;
    struct Poly *poly;
    struct Colours *colours;
};

static int toint(double d) {
    int i = (int)d;
    if (i <= 0) {
        i = (int)(d + -i + 1) - (-i + 1);
    }
    return i;
}

static bool plot(struct Params params) {
    struct Poly *dfdx, *dfdy;
    struct Bitmap *bm;

    double xstep, ystep;
    double x, xfrac, y, yfrac;
    double dzdx, dzdy, dxscale, dyscale;
    double z, xfade, yfade, fade;
    struct RGB c;
    int ii, i, j, xg, yg;

    dfdx = polypdiff(params.poly, true);
    dfdy = polypdiff(params.poly, false);

    bm = bmpinit(params.filename, params.width, params.height, params.outtype);

    xstep = (params.x1 - params.x0) / params.width;
    ystep = (params.y1 - params.y0) / params.height;
    dxscale = params.xscale * params.oscale;
    dyscale = params.yscale * params.oscale;
    for (ii = 0; ii < params.height; ii++) {
        if (params.outtype == BMP)
            i = ii;
        else
            i = params.height-1 - ii;
        y = params.y0 + ystep * i;
        yfrac = y / params.yscale; yfrac -= toint(yfrac);

        for (j = 0; j < params.width; j++) {
            x = params.x0 + xstep * j;
            xfrac = x / params.xscale; xfrac -= toint(xfrac);

            dzdx = polyeval(dfdx, x, y) * dxscale;
            dzdy = polyeval(dfdy, x, y) * dyscale;
            xg = toint(dzdx + 0.5); dzdx -= xg;
            yg = toint(dzdy + 0.5); dzdy -= yg;

            z = polyeval(params.poly, x, y) * params.oscale;
            z -= xg * xfrac;
            z -= yg * yfrac;
            z -= toint(z);

            xfade = dzdx; if (xfade < 0) xfade = -xfade;
            yfade = dzdy; if (yfade < 0) yfade = -yfade;
            fade = 1.0 - (xfade < yfade ? yfade : xfade) * 2;

            c = colfind(params.colours, xg, yg);

            if (params.fading) z *= fade;
            z *= 256.0;

            bmppixel(bm, toint(c.r*z), toint(c.g*z), toint(c.b*z));
        }
        bmpendrow(bm);
    }

    bmpclose(bm);
    polyfree(dfdy);
    polyfree(dfdx);
    return true;
}

/* ----------------------------------------------------------------------
 * Main program: parse the command line and call plot() when satisfied.
 */

bool parsepoly(char *string, void *vret) {
    struct Poly **ret = (struct Poly **)vret;
    *ret = polyread(string);
    if (!*ret)
        return false;
    return true;
}

bool parsecols(char *string, void *vret) {
    struct Colours **ret = (struct Colours **)vret;
    *ret = colread(string);
    if (!*ret)
        return false;
    return true;
}

int main(int argc, char **argv) {
    struct Params par;
    int i;
    double aspect;

    struct options {
        char *outfile;
        char *outtype;
        struct Size imagesize, basesize;
        double xcentre, ycentre, xrange, yrange, iscale, oscale;
        bool gotxcentre, gotycentre, gotxrange, gotyrange;
        bool gotiscale, gotoscale;
        bool fade, isbase, verbose;
        struct Poly *poly;
        struct Colours *colours;
    } optdata = {
        NULL, NULL, {0,0}, {0,0},
        0, 0, 0, 0, 0, 0,
        false, false, false, false, false, false,
        false, false, false, NULL, NULL
    };

    static const struct Cmdline options[] = {
        {1, "--output", 'o', "file.png", "output bitmap name",
                "filename", parsestr, offsetof(struct options, outfile), -1},
        {1, "--type", 0, "png|ppm|bmp", "set output file type",
                "image type", parsestr, offsetof(struct options, outtype), -1},
        {1, "--size", 's', "NNNxNNN", "output bitmap size",
                "output bitmap size", parsesize, offsetof(struct options, imagesize), -1},
        {1, "--xrange", 'x', "NNN", "mathematical x range",
                "x range", parseflt, offsetof(struct options, xrange), offsetof(struct options, gotxrange)},
        {1, "--yrange", 'y', "NNN", "mathematical y range",
                "y range", parseflt, offsetof(struct options, yrange), offsetof(struct options, gotyrange)},
        {2, "--xcentre\0--xcenter", 'X', "NNN", "mathematical x centre point",
                "x centre", parseflt, offsetof(struct options, xcentre), offsetof(struct options, gotxcentre)},
        {2, "--ycentre\0--ycenter", 'Y', "NNN", "mathematical y centre point",
                "y centre", parseflt, offsetof(struct options, ycentre), offsetof(struct options, gotycentre)},
        {1, "--basesize", 'b', "NNNxNNN", "base image size",
                "base image size", parsesize, offsetof(struct options, basesize), -1},
        {1, "--base", 'B', NULL, "this *is* the base image (default)",
                NULL, NULL, -1, offsetof(struct options, isbase)},
        {1, "--iscale", 'I', "NNN", "input scale (ie base pixel spacing)",
                "input scale", parseflt, offsetof(struct options, iscale), offsetof(struct options, gotiscale)},
        {1, "--oscale", 'O', "NNN", "output scale (ie modulus)",
                "output scale", parseflt, offsetof(struct options, oscale), offsetof(struct options, gotoscale)},
        {1, "--fade", 'f', NULL, "turn on fading at edges of patches",
                NULL, NULL, -1, offsetof(struct options, fade)},
        {1, "--poly", 'p', "x2+2xy+y2", "polynomial to plot",
                "polynomial", parsepoly, offsetof(struct options, poly), -1},
        {2, "--colours\0--colors", 'c', "1,0,0:0,1,0", "colours for patches",
                "colour specification", parsecols, offsetof(struct options, colours), -1},
        {1, "--verbose", 'v', NULL, "report details of what is done",
                NULL, NULL, -1, offsetof(struct options, verbose)},
    };

    char *usageextra[] = {
        " - at most one of -b, -B and -I should be used",
        " - can give 0 as one dimension in -s or -b and that dimension will",
        "   be computed so as to preserve the aspect ratio",
        " - if only one of -x and -y specified, will compute the other so",
        "   as to preserve the aspect ratio",
        " - colours can be prefixed with `N+' to get two-dimensional colour",
        "   selection (this needs to be better explained :-)"
    };

    parse_cmdline("filigram", argc, argv, options, lenof(options), &optdata);

    if (argc < 2)
        usage_message("filigram [options]",
                      options, lenof(options),
                      usageextra, lenof(usageextra));

    /*
     * Having read the arguments, now process them.
     */

    /* If no output scale, default to 1. */
    if (!optdata.gotoscale)
        par.oscale = 1.0;
    else
        par.oscale = optdata.oscale;

    /* If no output file, complain. */
    if (!optdata.outfile) {
        fprintf(stderr, "filigram: no output file specified: "
                "use something like `-o file.bmp'\n");
        return EXIT_FAILURE;
    } else
        par.filename = optdata.outfile;
    par.outtype = infer_type("bubbles", optdata.outtype, optdata.outfile);

    /* If no polynomial, complain. */
    if (!optdata.poly) {
        fprintf(stderr, "filigram: no polynomial specified: "
                "use something like `-p x2+2xy+y2'\n");
        return EXIT_FAILURE;
    } else
        par.poly = optdata.poly;

    /* If no colours, default to 1,1,1. */
    if (!optdata.colours)
        par.colours = colread("1,1,1");   /* assume this will succeed */
    else
        par.colours = optdata.colours;

    par.fading = optdata.fade;

    /*
     * If precisely one explicit aspect ratio specified, use it
     * to fill in blanks in other sizes.
     */
    aspect = 0;
    if (optdata.imagesize.w && optdata.imagesize.h) {
        aspect = (double)optdata.imagesize.w / optdata.imagesize.h;
    }
    if (optdata.basesize.w && optdata.basesize.h) {
        double newaspect = (double)optdata.basesize.w / optdata.basesize.h;
        if (newaspect != aspect)
            aspect = -1;
        else
            aspect = newaspect;
    }
    if (optdata.gotxrange && optdata.gotyrange) {
        double newaspect = optdata.xrange / optdata.yrange;
        if (newaspect != aspect)
            aspect = -1;
        else
            aspect = newaspect;
    }

    if (aspect > 0) {
        if (optdata.imagesize.w && !optdata.imagesize.h) optdata.imagesize.h = optdata.imagesize.w / aspect;
        if (!optdata.imagesize.w && optdata.imagesize.h) optdata.imagesize.w = optdata.imagesize.h * aspect;
        if (optdata.basesize.w && !optdata.basesize.h) optdata.basesize.h = optdata.basesize.w / aspect;
        if (!optdata.basesize.w && optdata.basesize.h) optdata.basesize.w = optdata.basesize.h * aspect;
        if (optdata.gotxrange && !optdata.gotyrange)
            optdata.yrange = optdata.xrange / aspect, optdata.gotyrange = true;
        if (!optdata.gotxrange && optdata.gotyrange)
            optdata.xrange = optdata.yrange * aspect, optdata.gotxrange = true;
    }

    /*
     * Now complain if no output image size was specified.
     */
    if (!optdata.imagesize.w || !optdata.imagesize.h) {
        fprintf(stderr, "filigram: no output size specified: "
                "use something like `-s 640x480'\n");
        return EXIT_FAILURE;
    } else {
        par.width = optdata.imagesize.w;
        par.height = optdata.imagesize.h;
    }

    /*
     * Also complain if no input mathematical extent specified.
     */
    if (!optdata.gotxrange || !optdata.gotyrange) {
        fprintf(stderr, "filigram: no image extent specified: "
                "use something like `-x 20 -y 15'\n");
        return EXIT_FAILURE;
    } else {
        if (!optdata.gotxcentre) optdata.xcentre = 0.0;
        if (!optdata.gotycentre) optdata.ycentre = 0.0;
        par.x0 = optdata.xcentre - optdata.xrange;
        par.x1 = optdata.xcentre + optdata.xrange;
        par.y0 = optdata.ycentre - optdata.yrange;
        par.y1 = optdata.ycentre + optdata.yrange;
    }

    /*
     * All that's left to set up is xscale and yscale. At this
     * stage we expect to see at most one of
     *   `optdata.isbase' true
     *   `optdata.basesize.w' and `optdata.basesize.h' both nonzero
     *   `optdata.gotiscale' true
     */
    i = (optdata.isbase) + (optdata.basesize.w && optdata.basesize.h) + (!!optdata.gotiscale);
    if (i > 1) {
        fprintf(stderr, "filigram: "
                "expected at most one of `-b', `-B', `-I'\n");
        return EXIT_FAILURE;
    } else if (i == 0) {
        /* Default: optdata.isbase is true. */
        optdata.isbase = true;
    }
    if (optdata.isbase) {
        optdata.basesize = optdata.imagesize;
    }
    if (optdata.gotiscale)
        par.xscale = par.yscale = optdata.iscale;
    else {
        par.xscale = (par.x1 - par.x0) / optdata.basesize.w;
        par.yscale = (par.y1 - par.y0) / optdata.basesize.h;
    }

    /*
     * If we're in verbose mode, regurgitate the final
     * parameters.
     */
    if (optdata.verbose) {
        printf("Output file `%s', %d x %d\n",
               par.filename, par.width, par.height);
        printf("Mathematical extent [%g,%g] x [%g,%g]\n",
               par.x0, par.x1, par.y0, par.y1);
        printf("Base pixel spacing %g (horiz), %g (vert)\n",
               par.xscale, par.yscale);
        printf("Output modulus %g\n",
               par.oscale);
    }

    i = plot(par) ? EXIT_SUCCESS : EXIT_FAILURE;

    colfree(par.colours);
    polyfree(par.poly);
    return i;
}
