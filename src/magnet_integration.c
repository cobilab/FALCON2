
#include "magnet_integration.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>   // access, unlink, mkstemps, close
#include <time.h>
#include <math.h>     // isnan

#include "common.h"

extern int32_t magnet_entry(int argc, char *argv[]);

static inline void push2(char **argv, int *argc, const char *flag, char *val) {
    argv[(*argc)++] = (char *)flag;
    argv[(*argc)++] = val;
}

static int run_magnet_like_cli(const char *reads,
                               const char *filter,
                               const char *out_filtered,   // -o
                               /* FALCON → MAGNET flags: */
                               int  magnetVerbose,         // -v
                               int  magnetInvert,          // -i
                               double magnetThreshold,     // -t
                               unsigned magnetLevel,       // -l
                               unsigned magnetPortion,     // -p
                               unsigned nThreads)          // -n
{
    char *argv[24];
    int argc = 0;

    argv[argc++] = (char *)"MAGNET";  // program name placeholder

    // Force overwrite like your previous external flow (safe to include)
    argv[argc++] = (char *)"-F";
    if (magnetVerbose) argv[argc++] = (char *)"-v";
    if (magnetInvert)  argv[argc++] = (char *)"-i";

    char tbuf[32], lbuf[16], pbuf[16], nbuf[16];

    // Only pass flags Falcon actually exposes:
    if (!isnan(magnetThreshold)) {
        snprintf(tbuf, sizeof tbuf, "%.6f", magnetThreshold);
        push2(argv, &argc, "-t", tbuf);
    }
    if (magnetLevel) {
        snprintf(lbuf, sizeof lbuf, "%u", magnetLevel);
        push2(argv, &argc, "-l", lbuf);
    }
    if (magnetPortion) {
        snprintf(pbuf, sizeof pbuf, "%u", magnetPortion);
        push2(argv, &argc, "-p", pbuf);
    }
    if (nThreads) {
        snprintf(nbuf, sizeof nbuf, "%u", nThreads);
        push2(argv, &argc, "-n", nbuf);
    }

    if (out_filtered) push2(argv, &argc, "-o", (char *)out_filtered);

    // Positional args: <reference> <reads> (this is MAGNET's expected order)
    argv[argc++] = (char *)filter;
    argv[argc++] = (char *)reads;

    return magnet_entry(argc, argv);
}

// Writes to a temp file, returns a FILE* to read it; auto-deletes on fclose()
FILE *RunMagnetPipe(const char *inputReads, const char *filterReference,
                    double threshold, U32 level, U8 invert, U8 verbose,
                    U32 portion, U32 nThreads)
{
    (void)portion; // passed through above; kept for signature parity

    char outTpl[] = "/tmp/falcon_magnet_pipe_XXXXXX.fq";
    int fd = mkstemps(outTpl, 3);  // keep ".fq"
    if (fd < 0) {
        perror("mkstemps");
        return NULL;
    }
    close(fd);

    int rc = run_magnet_like_cli(inputReads, filterReference, outTpl,
                                 verbose, invert, threshold, level, portion, nThreads);
    if (rc != 0) {
        unlink(outTpl);
        return NULL;
    }

    FILE *fp = fopen(outTpl, "r");
    if (!fp) {
        unlink(outTpl);
        return NULL;
    }
    // Remove on close: file persists while this FILE* is open
    unlink(outTpl);
    return fp;
}

int RunMagnet(const char *inputReads, const char *filterReference,
              double threshold, U32 level, U8 invert, U8 verbose,
              U32 portion, const char *outputFile, U32 nThreads)
{
    if (!outputFile) {
        fprintf(stderr, "Error: RunMagnet requires an output path.\n");
        return EXIT_FAILURE;
    }

    return run_magnet_like_cli(inputReads, filterReference, outputFile,
                               verbose, invert, threshold, level, portion, nThreads);
}