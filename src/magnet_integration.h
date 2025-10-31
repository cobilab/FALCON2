
#ifndef MAGNET_INTEGRATION_H
#define MAGNET_INTEGRATION_H

#include "defs.h"
#include <stdio.h>

/**
 * Run MAGNET and get its output as a FILE* stream
 * 
 * This function uses popen to create a pipe to MAGNET's output
 * so FALCON can process it directly without intermediate files.
 * 
 * @param inputReads The input FASTQ reads to filter
 * @param filterReference The reference FASTA to use for filtering 
 * @param threshold The similarity threshold [0.0-1.0]
 * @param level The sensitivity level [1-44]
 * @param invert Whether to invert the filtering (keep non-matching reads)
 * @param verbose Whether to display verbose information
 * @param portion The portion of acceptance parameter
 * @param nThreads Number of threads to use (0 for default)
 * @return FILE* handle to MAGNET's output stream or NULL on error
 */
FILE *RunMagnetPipe(const char *inputReads, const char *filterReference,
    double threshold, U32 level, U8 invert, U8 verbose,
    U32 portion, U32 nThreads);
    
/**
 * Run MAGNET to filter sequences and save output to a file
 * 
 * This function runs MAGNET as a system command and redirects output to a file.
 * 
 * @param inputReads The input file to filter
 * @param filterReference The reference FASTA to use for filtering
 * @param threshold The similarity threshold [0.0-1.0]
 * @param level The sensitivity level [1-44]
 * @param invert Whether to invert the filtering (keep non-matching reads)
 * @param verbose Whether to display verbose information
 * @param portion The portion of acceptance parameter
 * @param outputFile The path to save the filtered output to
 * @param nThreads Number of threads to use (0 for default)
 * @return 0 on success, non-zero on error
 */
int RunMagnet(const char *inputReads, const char *filterReference,
    double threshold, U32 level, U8 invert, U8 verbose,
    U32 portion, const char *outputFile, U32 nThreads);

#endif //MAGNET_INTEGRATION_H
