#ifndef READS_H_INCLUDED
#define READS_H_INCLUDED

#include "defs.h"

#define DEFAULT_HEADER_SIZE  100000
#define DEFAULT_READ_SIZE    100000
#define READ_LEFT_GUARD      32

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct{
  char      *header1[2];
  char      *bases;
  char      *header2;
  char      *scores;
  uint64_t  headerMaxSize;
  uint64_t  readMaxSize;
  uint64_t  header2Present;
  uint64_t  skipNs;
  char      lowestScore;
  }
Read;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Read    *CreateRead           (uint64_t, uint64_t);
Read    *GetRead              (FILE *, Read *);
void    PutRead               (Read *, FILE *);

#endif

