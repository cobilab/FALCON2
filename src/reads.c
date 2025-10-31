#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "defs.h"
#include "mem.h"
#include "reads.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Read *CreateRead(uint64_t headerMaxSize, uint64_t readMaxSize){

  Read *read           = (Read *)    Calloc(1, sizeof(Read));
  read->header1[0]     = (char *) Calloc(headerMaxSize + READ_LEFT_GUARD, 
                         sizeof(char));
  read->header1[1]     = (char *) Calloc(headerMaxSize + READ_LEFT_GUARD, 
                         sizeof(char));
  read->header2        = (char *) Calloc(headerMaxSize + READ_LEFT_GUARD, 
                         sizeof(char));
  read->headerMaxSize  = headerMaxSize;
  read->bases          = (char *) Calloc(readMaxSize + READ_LEFT_GUARD, 
                         sizeof(char));
  read->scores         = (char *) Calloc(readMaxSize + READ_LEFT_GUARD, 
                         sizeof(char));
  read->readMaxSize    = readMaxSize;
  read->header1[0]    += READ_LEFT_GUARD;
  read->header1[1]    += READ_LEFT_GUARD;
  read->header2       += READ_LEFT_GUARD;
  read->bases         += READ_LEFT_GUARD;
  read->scores        += READ_LEFT_GUARD;
  read->header2Present = 0;
  read->skipNs         = 0;
  read->lowestScore    = (uint8_t) 255;

  return read;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Read *GetRead(FILE *fp, Read *read){
  int n, c = fgetc(fp);

  if(c == EOF) return NULL;

  if(c != '@'){
    fprintf(stderr, "Error: failed to get the initial '@' character\n");
    exit(1);
    }

  if(!fgets((char *)read->header1[1], read->headerMaxSize, fp)){
    fprintf(stderr, "Error: unexpected end of file\n");
    exit(1);
    }

  if(!fgets((char *)read->bases, read->readMaxSize, fp)){
    fprintf(stderr, "Error: unexpected end of file\n");
    exit(1);
    }

  if(!fgets((char *)read->header2, read->headerMaxSize, fp)){
    fprintf(stderr, "Error: unexpected end of file\n");
    exit(1);
    }

  if(!fgets((char *)read->scores, read->readMaxSize, fp)){
    fprintf(stderr, "Error: unexpected end of file\n");
    exit(1);
    }

  return read;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PutRead(Read *read, FILE *fp){
  uint64_t n, nBases;

  fputc('@', fp);
  fputs((char *)read->header1[1], fp);
  fputs((char *)read->bases, fp);
  if(read->header2Present)
    fprintf(fp, "+%s", read->header1[1]);
  else
    fputs("+\n", fp);
  fputs((char *)read->scores, fp);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
