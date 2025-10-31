/**
 * @file /**
 * @file file_compression.c
 * @brief Implementation of compressed file I/O operations
 */

#include "file_compression.h"
#include "mem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

FILE* CFopen(const char* filename, const char* mode) {
    // Check if file has a compression extension (.gz)
    if (ends_with(filename, ".gz")) {
        // Create a command like "gzip -dc filename"
        char cmd[1024];
        sprintf(cmd, "gzip -dc %s", filename);
        
        // Open a pipe to the command
        FILE* f = popen(cmd, "r");
        if (!f) {
            fprintf(stderr, "Error opening compressed file: %s\n", filename);
            exit(1);
        }
        return f;
    }
    
    // For uncompressed files, use regular Fopen
    return Fopen(filename, mode);
}

int CFclose(const char *filename, FILE *f)
{
    if (!f) return 0;
    /* We only popen() for .gz reads in CFopen; close with pclose in that case. */
    if (ends_with(filename, ".gz")) return pclose(f);
    return fclose(f);
}

int ConcatWithCFopen(char *const *files, uint32_t nFiles, const char *outPath){
    FILE *OUT = Fopen(outPath, "wb");
    uint32_t i;
    size_t   k;
    uint8_t  buffer[ BUFFER_SIZE ];

    for(i = 0; i < nFiles; ++i){
        FILE *IN = CFopen(files[i], "r");
        fprintf(stderr, "      [+] Concatenating %u/%u: %s ... ", i + 1, nFiles, files[i]);
        while((k = fread(buffer, 1, sizeof buffer, IN)) > 0){
            if(fwrite(buffer, 1, k, OUT) != k){
                fprintf(stderr, "  [x] Error: fwrite failed while concatenating %s\n", files[i]);
                CFclose(files[i], IN);
                Fclose(OUT);
                return EXIT_FAILURE;
            }
        }
        if(ferror(IN)){
            fprintf(stderr, "  [x] Error: fread failed while reading %s\n", files[i]);
            CFclose(files[i], IN);
            Fclose(OUT);
            return EXIT_FAILURE;
        }
        CFclose(files[i], IN);
        fprintf(stderr, "Done!\n");
    }
    Fclose(OUT);
    return EXIT_SUCCESS;
}