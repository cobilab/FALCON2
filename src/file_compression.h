/**
 * @file file_compression.h
 * @brief Header file for compressed file I/O operations
 * 
 * This module provides functions to handle reading from both regular and 
 * compressed files (gzip) using a unified interface.
 */

#ifndef FILE_COMPRESSION_H
#define FILE_COMPRESSION_H
#include <stdint.h>
#include <stdio.h>

#endif //FILE_COMPRESSION_H

/**
 * @brief Opens a file, automatically detecting compression format
 * 
 * @param filename Path to the file
 * @param mode File opening mode ("r", "w", etc.)
 * @return File* Pointer to the opened file
 */
FILE *CFopen(const char *filename, const char *mode);

/**
 * @brief Closes a CFopen() stream correctly: pclose for .gz, fclose otherwise.
 *
 * @param filename Path to the file
 * @return File* Pointer to the opened file
 */
int CFclose(const char *filename, FILE *f);

/**
 * @brief Concatenate (possibly compressed) input files into a single output file.
 *        Each input is opened via CFopen (supports .gz), output via Fopen.
 *
 * @param files    Array of file paths
 * @param nFiles   Number of files
 * @param outPath  Output file path to write the concatenation into
 * @return 0 on success; non-zero error code on failure
 */
int ConcatWithCFopen(char *const *files, uint32_t nFiles, const char *outPath);