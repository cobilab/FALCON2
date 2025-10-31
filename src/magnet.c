#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <fcntl.h>
#include <unistd.h>

#include "mem.h"
#include "time.h"
#include "defs.h"
#include "param.h"
#include "msg.h"
#include "parser.h"
#include "reads.h"
#include "buffer.h"
#include "levels.h"
#include "common.h"
#include "file_compression.h"
#include "models.h"
#include "pmodels.h"
#include "stream.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - M O D E L S - - - - - - - - - - - - - - - -

CModel **ModelsMagnet;   // MEMORY SHARED BY THREADING
Parameters *PMagnet;

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - R E S E T   M O D E L S - - - - - - - - - - - -

void ResetModelsAndParamMagnet(CBUF *Buf, CModel **Shadow, CMWeight *CMW){
  uint32_t n;
  ResetCBuffer(Buf);
  for(n = 0 ; n < PMagnet->nModels ; ++n)
    ResetShadowModel(Shadow[n]);
  ResetWeightModel(CMW);
}

// - - - - - - - - - - - - - - C O M P R E S S I O N - - - - - - - - - - - - -

void CompressTargetMagnet(Threads T){
  FILE        *Reader = Fopen(PMagnet->base, "r");
  double      bits = 0;
  uint64_t    nBase = 0, r = 0, nSymbol, initNSymbol;
  uint32_t    n, k, idxPos, totModels, cModel;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     sym, *pos, conName[MAX_NAME];
  PModel      **pModel, *MX;
  CModel      **Shadow; // SHADOWS FOR SUPPORTING MODELS WITH THREADING
  FloatPModel *PT;
  CMWeight    *CMW;
  int         action;

  totModels = PMagnet->nModels; // EXTRA MODELS DERIVED FROM EDITS
  for(n = 0 ; n < PMagnet->nModels ; ++n)
    if(T.model[n].edits != 0)
      totModels += 1;

  Shadow      = (CModel **) Calloc(PMagnet->nModels, sizeof(CModel *));
  for(n = 0 ; n < PMagnet->nModels ; ++n)
    Shadow[n] = CreateShadowModel(ModelsMagnet[n]);

  pModel      = (PModel **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n] = CreatePModel(ALPHABET_SIZE);
  MX          = CreatePModel(ALPHABET_SIZE);
  PT          = CreateFloatPModel(ALPHABET_SIZE);
  CMW         = CreateWeightModel(totModels);

  FileType(PA, Reader);
  if(PA->type != 2){
    fprintf(stderr, "Error: input file must be in FASTQ format!\n");
    exit(1);
    }

  char name_o [10000];
  sprintf(name_o,  "%s.%u", PMagnet->output,  T.id);
  FILE *Writer  = Fopen(name_o,  "w");
  srand(T.id);
  Read *Read = CreateRead(50000, 100000);

  while((Read = GetRead(Reader, Read)) != NULL)
    {
    if(PA->nRead % PMagnet->nThreads == T.id)
      {
      nBase = strlen(Read->bases) - 1; // IT ALSO LOADS '\n' AT THE END
      bits  = 0;

      for(idxPos = 0 ; idxPos < nBase ; ++idxPos)
        {
        sym = Read->bases[idxPos];

	if((sym = DNASymToNum(sym)) == 4)
	  sym = rand() % 4;

        symBuf->buf[symBuf->idx] = sym;
        memset((void *)PT->freqs, 0, (ALPHABET_SIZE)* sizeof(double));
        n = 0;
        pos = &symBuf->buf[symBuf->idx-1];
        for(cModel = 0 ; cModel < PMagnet->nModels ; ++cModel){
          CModel *CM = Shadow[cModel];
          GetPModelIdx(pos, CM);
          ComputePModel(ModelsMagnet[cModel], pModel[n], CM->pModelIdx, CM->alphaDen);
          ComputeWeightedFreqs(CMW->weight[n], pModel[n], PT);
          if(CM->edits != 0){
            ++n;
            CM->SUBS.seq->buf[CM->SUBS.seq->idx] = sym;
            CM->SUBS.idx = GetPModelIdxCorr(CM->SUBS.seq->buf+CM->SUBS.seq->idx
            -1, CM, CM->SUBS.idx);
            ComputePModel(ModelsMagnet[cModel], pModel[n], CM->SUBS.idx, CM->SUBS.eDen);
            ComputeWeightedFreqs(CMW->weight[n], pModel[n], PT);
            }
          ++n;
          }

        ComputeMXProbs(PT, MX);
        bits += PModelSymbolLog(MX, sym);
        CalcDecayment(CMW, pModel, sym, PMagnet->gamma);
        RenormalizeWeights(CMW);
        CorrectXModels(Shadow, pModel, sym, PMagnet->nModels);
        UpdateCBuffer(symBuf);
        }

      if(BPBB(bits, nBase) < PMagnet->threshold){
        if(PMagnet->invert) fputc('0', Writer); // IGNORE READ
        else          fputc('1', Writer); // WRITE READ
        }
      else{
        if(PMagnet->invert) fputc('1', Writer); // WRITE READ
        else          fputc('0', Writer); // IGNORE READ
        }
      ResetModelsAndParamMagnet(symBuf, Shadow, CMW);
      }

    PA->nRead++;
    }

  DeleteWeightModel(CMW);
  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(pModel[n]);
  Free(pModel);
  RemovePModel(MX);
  RemoveFPModel(PT);
  for(n = 0 ; n < PMagnet->nModels ; ++n)
    FreeShadow(Shadow[n]);
  Free(Shadow);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Reader);
  fclose(Writer);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - F   T H R E A D I N G - - - - - - - - - - - - - - -

void *CompressThreadMagnet(void *Thr){
  Threads *T = (Threads *) Thr;
  CompressTargetMagnet(T[0]);
  pthread_exit(NULL);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E F E R E N C E - - - - - - - - - - - - -

void LoadReferenceMagnet(char *refName){
  FILE     *Reader = CFopen(refName, "r");
  uint32_t n;
  uint64_t idx = 0;
  uint64_t k, idxPos;
  PARSER   *PA = CreateParser();
  CBUF     *symBuf  = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t  *readBuf = Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t  sym, irSym = 0;
  FileType(PA, Reader);
  rewind(Reader);
  srand(0);

  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      if(ParseSym(PA, (sym = readBuf[idxPos])) == -1){
        idx = 0;
        continue;
        }

      if(sym == 'N') // WE CAN RAND HERE CAUSE IS ALWAYS IN ONE THREAD
        symBuf->buf[symBuf->idx] = sym = 0; //(rand() % 4);
      else
        symBuf->buf[symBuf->idx] = sym = DNASymToNum(sym);

      for(n = 0 ; n < PMagnet->nModels ; ++n){
        CModel *CM = ModelsMagnet[n];
        GetPModelIdx(symBuf->buf+symBuf->idx-1, CM);
        if(CM->ir == 1) // INVERTED REPEATS
          irSym = GetPModelIdxIR(symBuf->buf+symBuf->idx, CM);
        if(++idx > CM->ctx){
          UpdateCModelCounter(CM, sym, CM->pModelIdx);
          if(CM->ir == 1) // INVERTED REPEATS
            UpdateCModelCounter(CM, irSym, CM->pModelIdxIR);
          }
        }
      UpdateCBuffer(symBuf);
      }

  for(n = 0 ; n < PMagnet->nModels ; ++n)
    ResetCModelIdx(ModelsMagnet[n]);
  RemoveCBuffer(symBuf);
  Free(readBuf);
  RemoveParser(PA);
  fclose(Reader);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - C O M P R E S S O R   M A I N - - - - - - - - - - - -

void CompressActionMagnet(Threads *T, char *refName, char *baseName){
  pthread_t t[PMagnet->nThreads];
  uint32_t n;

  ModelsMagnet = (CModel **) Malloc(PMagnet->nModels * sizeof(CModel *));
  for(n = 0 ; n < PMagnet->nModels ; ++n)
    ModelsMagnet[n] = CreateCModel(T[0].model[n].ctx, T[0].model[n].den,
    T[0].model[n].ir, REFERENCE, PMagnet->col, T[0].model[n].edits,
    T[0].model[n].eDen);
  fprintf(stderr, "  [+] Loading reference file ....... ");
  LoadReferenceMagnet(refName);
  fprintf(stderr, "Done!\n");

  fprintf(stderr, "  [+] Filtering FASTQ reads ........ ");
  for(n = 0 ; n < PMagnet->nThreads ; ++n)
    pthread_create(&(t[n+1]), NULL, CompressThreadMagnet, (void *) &(T[n]));
  for(n = 0 ; n < PMagnet->nThreads ; ++n) // DO NOT JOIN FORS!
    pthread_join(t[n+1], NULL);
  fprintf(stderr, "Done!\n");

  fprintf(stderr, "  [+] Joinning streams ............. ");
  FILE *OUT   = Fopen(PMagnet->output,  "w");
  FILE *OUT2  = Fopen(PMagnet->output2, "w");
  FILE *IN    = CFopen(PMagnet->base,    "r");
  FILE **TMP  = (FILE **) Calloc(PMagnet->nThreads, sizeof(FILE *));
  for(n = 0 ; n < PMagnet->nThreads ; ++n){
    char name_o[MAX_NAME];
    sprintf(name_o, "%s.%u", PMagnet->output, n);
    TMP[n] = Fopen(name_o, "r");
    }

  Read *Read = CreateRead(10000, 40000);
  n = 0;
  while((Read = GetRead(IN, Read)) != NULL){
    if(fgetc(TMP[n++ % PMagnet->nThreads]) == '1')
      PutRead(Read, OUT);
    else
      PutRead(Read, OUT2);
    }

  for(n = 0 ; n < PMagnet->nThreads ; ++n){
    fclose(TMP[n]);
    char name_o[MAX_NAME];
    sprintf(name_o, "%s.%u", PMagnet->output, n);
    Fdelete(name_o);
    }
  fclose(IN);
  fclose(OUT);
  fclose(OUT2);
  fprintf(stderr, "Done!\n");
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


int32_t magnet_entry(int argc, char *argv[]){

  char     **p = *&argv, **xargv, *xpl = NULL;
  int32_t  xargc = 0;
  uint32_t n, k, col, ref;
  double   gamma;
  Threads  *T;

  PMagnet = (Parameters *) Malloc(1 * sizeof(Parameters));

  if(ArgsState(DEF_VERSION, p, argc, "-V", "--version")){
    PrintVersion();
    return EXIT_SUCCESS;
    }

  if(ArgsState(0, p, argc, "-s", "--show")){
    PrintLevels();
    return EXIT_SUCCESS;
    }

  PMagnet->verbose  = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v", "--verbose");
  PMagnet->force    = ArgsState  (DEFAULT_FORCE,   p, argc, "-F", "--force");
  PMagnet->invert   = ArgsState  (DEFAULT_INVERT,  p, argc, "-i", "--invert");
  PMagnet->sample   = ArgsNum    (DEFAULT_SAMPLE,  p, argc, "-p", MIN_SAP, MAX_SAP);
  PMagnet->level    = ArgsNum    (0,               p, argc, "-l", MIN_LEV, MAX_LEV);
  PMagnet->nThreads = ArgsNum    (DEFAULT_THREADS, p, argc, "-n", MIN_THREADS,
  MAX_THREADS);

  PMagnet->nModels = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-m") == 0)
      PMagnet->nModels += 1;

  if(PMagnet->nModels == 0 && PMagnet->level == 0)
    PMagnet->level = DEFAULT_LEVEL;

  if(PMagnet->level != 0){
    xpl = GetLevels(PMagnet->level);
    xargc = StrToArgv(xpl, &xargv);
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-m") == 0)
        PMagnet->nModels += 1;
    }

  gamma = DEFAULT_GAMMA;
  for(n = 1 ; n < xargc ; ++n)
    if(strcmp(xargv[n], "-g") == 0)
      gamma = atof(xargv[n+1]);

  col = MAX_COLLISIONS;
  for(n = 1 ; n < xargc ; ++n)
    if(strcmp(xargv[n], "-c") == 0)
      col = atoi(xargv[n+1]);

  PMagnet->col       = ArgsNum    (col,   p, argc, "-c", 1, 253);
  PMagnet->gamma     = ArgsDouble (gamma, p, argc, "-g");
  PMagnet->threshold = fabs(ArgsDouble (0.9,   p, argc, "-t"));
  PMagnet->gamma     = ((int)(PMagnet->gamma * 65536)) / 65536.0;
  PMagnet->output    = ArgsFileGen(p, argc, "-o", "filtered", ".fq");
  PMagnet->output2   = ArgsFileGen(p, argc, "-2", "not-filtered", ".fq");

  if(!PMagnet->force)
    FAccessWPerm(PMagnet->output);

  if(PMagnet->nModels == 0){
    fprintf(stderr, "Error: at least you need to use a context model!\n");
    return EXIT_FAILURE;
    }

  // READ MODEL PARAMETERS FROM XARGS & ARGS
  T = (Threads *) Calloc(PMagnet->nThreads, sizeof(Threads));
  for(ref = 0 ; ref < PMagnet->nThreads ; ++ref){
    T[ref].model = (ModelPar *) Calloc(PMagnet->nModels, sizeof(ModelPar));
    T[ref].id    = ref;
    k = 0;
    for(n = 1 ; n < argc ; ++n)
      if(strcmp(argv[n], "-m") == 0)
        T[ref].model[k++] = ArgsUniqModel(argv[n+1], 0);
    if(PMagnet->level != 0){
      for(n = 1 ; n < xargc ; ++n)
        if(strcmp(xargv[n], "-m") == 0)
          T[ref].model[k++] = ArgsUniqModel(xargv[n+1], 0);
      }
    }

  fprintf(stderr, "\n");
  if(PMagnet->verbose) PrintArgsMagnet(PMagnet, T[0], argv[argc-2], argv[argc-1]);

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  TIME *Time = CreateClock(clock());
  PMagnet->base = argv[argc-1];
  CompressActionMagnet(T, argv[argc-2], PMagnet->base);

  fprintf(stderr, "  [+] Freeing compression models ... ");
  for(n = 0 ; n < PMagnet->nModels ; ++n)
    FreeCModel(ModelsMagnet[n]);
  Free(ModelsMagnet);
  fprintf(stderr, "Done!\n");

  StopTimeNDRM(Time, clock());
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ RESULTS ]=======================\n");
  fprintf(stderr, "Filtered reads in   : %s\n", PMagnet->output);
  fprintf(stderr, "Unfiltered reads in : %s\n", PMagnet->output2);
  fprintf(stderr, "\n");

  if(PMagnet->verbose)
    {
    fprintf(stderr, "==[ STATISTICS ]====================\n");
    StopCalcAll(Time, clock());
    fprintf(stderr, "\n");
    }

  RemoveClock(Time);
  for(ref = 0 ; ref < PMagnet->nThreads ; ++ref)
    Free(T[ref].model);
  Free(T);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
