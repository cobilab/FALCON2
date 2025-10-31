#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <fcntl.h>
#include <regex.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/uio.h>
#include <sys/mman.h>

#include "mem.h"
#include "time.h"
#include "defs.h"
#include "keys.h"
#include "param.h"
#include "msg.h"
#include "top.h"
#include "parser.h"
#include "buffer.h"
#include "levels.h"
#include "common.h"
#include "file_compression.h"
#include "serialization.h"
#include "magnet_integration.h"
#include "filters.h"
#include "models.h"
#include "pmodels.h"
#include "kmodels.h"
#include "labels.h"
#include "paint.h"
#include "stream.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - M O D E L S   A N D   P A R A M E T E R S - - - - - - - - - -

CModel     **Models;   // MEMORY SHARED BY THREADING
KMODEL     **KModels;  // MEMORY SHARED BY THREADING
Parameters *P;
EYEPARAM   *PEYE;


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - R E S E T   M O D E L S - - - - - - - - - - - -

void ResetModelsAndParam(CBUF *Buf, CModel **Shadow, CMWeight *CMW){
  uint32_t n;
  ResetCBuffer(Buf);
  for(n = 0 ; n < P->nModels ; ++n)
    ResetShadowModel(Shadow[n]);
  ResetWeightModel(CMW);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - R E S E T   K M O D E L S - - - - - - - - - - - -

void ResetKModelsAndParam(CBUF *Buf, KMODEL **Shadow, CMWeight *CMW){
  uint32_t n;
  ResetCBuffer(Buf);
  for(n = 0 ; n < P->nModels ; ++n)
    ResetKShadowModel(Shadow[n]);
  ResetWeightModel(CMW);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - L O C A L   C O M P L E X I T Y - - - - - - - - - - - -

#ifdef LOCAL_SIMILARITY
void LocalComplexity(Threads T, TOP *Top, uint64_t topSize, FILE *OUT){
  FILE        *Reader = Fopen(P->base, "r");
  double      bits = 0, instant = 0;
  uint64_t    nBase = 0, entry;
  uint32_t    n, totModels, cModel;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t     *pos;
  PModel      **pModel, *MX;
  CModel      **Shadow; // SHADOWS FOR SUPPORTING MODELS WITH THREADING
  FloatPModel *PT;
  CMWeight    *CMW;
  int         sym;

  totModels = P->nModels; // EXTRA MODELS DERIVED FROM EDITS
  for(n = 0 ; n < P->nModels ; ++n) 
    if(T.model[n].edits != 0)
      totModels += 1;

  Shadow      = (CModel **) Calloc(P->nModels, sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    Shadow[n] = CreateShadowModel(Models[n]); 
  pModel      = (PModel **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n] = CreatePModel(ALPHABET_SIZE);
  MX          = CreatePModel(ALPHABET_SIZE);
  PT          = CreateFloatPModel(ALPHABET_SIZE);
  CMW         = CreateWeightModel(totModels);

  for(entry = 0 ; entry < topSize ; ++entry){
    if(Top->V[entry].size > 1){ 
      fprintf(stderr, "      [+] Running profile: %-5"PRIu64" ... ", entry + 1);

      // PRINT HEADER COMPLEXITY VALUE
      fprintf(OUT, "#\t%.5lf\t%"PRIu64"\t%s\n", (1.0-Top->V[entry].value)*100.0, 
      Top->V[entry].size, Top->V[entry].name);

      // MOVE POINTER FORWARD
      Fseeko(Reader, (off_t) Top->V[entry].iPos-1, SEEK_SET); 

      while((sym = fgetc(Reader)) != EOF){

        if(sym == '>'){ // FOUND HEADER & SKIP 
          while((sym = fgetc(Reader)) != '\n' && sym != EOF)
            ; // DO NOTHING

          if(sym == EOF) 
            break;     // END OF FILE: QUIT
          }

        if(nBase >= Top->V[entry].size) // IT PROCESSED ALL READ BASES: QUIT!
          break;

        if(sym == '\n') continue;  // SKIP '\n' IN FASTA

        if((sym = DNASymToNum(sym)) == 4){
          fprintf(OUT, "%c", PackByte(2.0, sym)); // PRINT COMPLEXITY & SYM IN1
          continue; // IT IGNORES EXTRA SYMBOLS
          }

        symBuf->buf[symBuf->idx] = sym;
        memset((void *)PT->freqs, 0, ALPHABET_SIZE * sizeof(double));
        n = 0;
        pos = &symBuf->buf[symBuf->idx-1];
        for(cModel = 0 ; cModel < P->nModels ; ++cModel){
          CModel *CM = Shadow[cModel];
          GetPModelIdx(pos, CM);
          ComputePModel(Models[cModel], pModel[n], CM->pModelIdx, CM->alphaDen);
          ComputeWeightedFreqs(CMW->weight[n], pModel[n], PT);
          if(CM->edits != 0){
            ++n;
            CM->SUBS.seq->buf[CM->SUBS.seq->idx] = sym;
            CM->SUBS.idx = GetPModelIdxCorr(CM->SUBS.seq->buf+CM->SUBS.seq->idx
            -1, CM, CM->SUBS.idx);
            ComputePModel(Models[cModel], pModel[n], CM->SUBS.idx, CM->SUBS.eDen);
            ComputeWeightedFreqs(CMW->weight[n], pModel[n], PT);
            }
          ++n;
          }

        ComputeMXProbs(PT, MX);
        instant = PModelSymbolLog(MX, sym);
        bits += instant;
        fprintf(OUT, "%c", PackByte(instant, sym)); // PRINT COMPLEX & SYM IN1
        ++nBase;
        CalcDecayment(CMW, pModel, sym, P->gamma);
        RenormalizeWeights(CMW);
        CorrectXModels(Shadow, pModel, sym, P->nModels);
        UpdateCBuffer(symBuf);
        }

      if(entry < topSize - 1){ // RESET MODELS & PROPERTIES
        ResetModelsAndParam(symBuf, Shadow, CMW);
        nBase = bits = 0;
        }

      fprintf(OUT, "\n");
      fprintf(stderr, "Done!\n");
      }
    } 

  DeleteWeightModel(CMW);
  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(pModel[n]);
  Free(pModel);
  RemovePModel(MX);
  RemoveFPModel(PT);
  for(n = 0 ; n < P->nModels ; ++n)
    FreeShadow(Shadow[n]);
  Free(Shadow);
  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Reader);
  }
#endif

#ifdef LOCAL_SIMILARITY
void LocalComplexityVariousDBs(Threads T, TOP *Top, uint64_t topSize, char **baseFiles, U8 fileCount, FILE *OUT){
  double      bits = 0, instant = 0;
  uint64_t    nBase = 0, entry;
  uint32_t    n, totModels, cModel;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t     *pos;
  PModel      **pModel, *MX;
  CModel      **Shadow; // SHADOWS FOR SUPPORTING MODELS WITH THREADING
  FloatPModel *PT;
  CMWeight    *CMW;
  int         sym;

  // Calculate total models including edits
  totModels = P->nModels; // EXTRA MODELS DERIVED FROM EDITS
  for(n = 0 ; n < P->nModels ; ++n)
    if(T.model[n].edits != 0)
      totModels += 1;

  // Initialize shadow models
  Shadow      = (CModel **) Calloc(P->nModels, sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    Shadow[n] = CreateShadowModel(Models[n]);

  // Initialize prediction models and weights
  pModel      = (PModel **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n] = CreatePModel(ALPHABET_SIZE);
  MX          = CreatePModel(ALPHABET_SIZE);
  PT          = CreateFloatPModel(ALPHABET_SIZE);
  CMW         = CreateWeightModel(totModels);

  // Loop over each input file
  for (uint32_t fIdx = 0; fIdx < fileCount; ++fIdx) {
    char *currentDb = baseFiles[fIdx];
    FILE        *Reader = Fopen(currentDb, "r");
    nBase = bits = 0;

    // Process each TOP entry
    for(entry = 0 ; entry < topSize ; ++entry){
      if(Top->V[entry].size > 1){
        fprintf(stderr, "      [+] Running profile: %-5"PRIu64" ... ", entry + 1);

        // PRINT HEADER COMPLEXITY VALUE
        fprintf(OUT, "#\t%.5lf\t%"PRIu64"\t%s\n", (1.0-Top->V[entry].value)*100.0,
        Top->V[entry].size, Top->V[entry].name);

        // MOVE POINTER FORWARD
        Fseeko(Reader, (off_t) Top->V[entry].iPos-1, SEEK_SET);

        while((sym = fgetc(Reader)) != EOF){

          if(sym == '>'){ // FOUND HEADER & SKIP
            while((sym = fgetc(Reader)) != '\n' && sym != EOF)
              ; // DO NOTHING

            if(sym == EOF)
              break;     // END OF FILE: QUIT
          }

          if(nBase >= Top->V[entry].size) // IT PROCESSED ALL READ BASES: QUIT!
            break;

          if(sym == '\n') continue;  // SKIP '\n' IN FASTA

          if((sym = DNASymToNum(sym)) == 4){
            fprintf(OUT, "%c", PackByte(2.0, sym)); // PRINT COMPLEXITY & SYM IN1
            continue; // IT IGNORES EXTRA SYMBOLS
          }

          symBuf->buf[symBuf->idx] = sym;
          memset((void *)PT->freqs, 0, ALPHABET_SIZE * sizeof(double));
          n = 0;
          pos = &symBuf->buf[symBuf->idx-1];
          for(cModel = 0 ; cModel < P->nModels ; ++cModel){
            CModel *CM = Shadow[cModel];
            GetPModelIdx(pos, CM);
            ComputePModel(Models[cModel], pModel[n], CM->pModelIdx, CM->alphaDen);
            ComputeWeightedFreqs(CMW->weight[n], pModel[n], PT);
            if(CM->edits != 0){
              ++n;
              CM->SUBS.seq->buf[CM->SUBS.seq->idx] = sym;
              CM->SUBS.idx = GetPModelIdxCorr(CM->SUBS.seq->buf+CM->SUBS.seq->idx
              -1, CM, CM->SUBS.idx);
              ComputePModel(Models[cModel], pModel[n], CM->SUBS.idx, CM->SUBS.eDen);
              ComputeWeightedFreqs(CMW->weight[n], pModel[n], PT);
            }
            ++n;
          }

          ComputeMXProbs(PT, MX);
          instant = PModelSymbolLog(MX, sym);
          bits += instant;
          fprintf(OUT, "%c", PackByte(instant, sym)); // PRINT COMPLEX & SYM IN1
          ++nBase;
          CalcDecayment(CMW, pModel, sym, P->gamma);
          RenormalizeWeights(CMW);
          CorrectXModels(Shadow, pModel, sym, P->nModels);
          UpdateCBuffer(symBuf);
        }

        if(entry < topSize - 1){ // RESET MODELS & PROPERTIES
          ResetModelsAndParam(symBuf, Shadow, CMW);
          nBase = bits = 0;
        }

        fprintf(OUT, "\n");
        fprintf(stderr, "Done!\n");
      }
    }
    fclose(Reader);
  }

  DeleteWeightModel(CMW);
  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(pModel[n]);
  Free(pModel);
  RemovePModel(MX);
  RemoveFPModel(PT);
  for(n = 0 ; n < P->nModels ; ++n)
    FreeShadow(Shadow[n]);
  Free(Shadow);
  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  }
#endif


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - P E R E G R I N E   C O M P R E S S I O N - - - - - - - -

void SamplingCompressTarget(Threads T){
  FILE        *Reader  = Fopen(P->base, "r");
  double      bits = 0;
  uint64_t    nBase = 0, r = 0, idx = 0, initNSymbol, nSymbol;
  uint32_t    k, idxPos;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t     sym, *pos, conName[MAX_NAME];
  PModel      *pModel;
  CModel      *Shadow; // SHADOW FOR SUPPORTING MODEL WITH THREADING
  int         action;

  Shadow = CreateShadowModel(Models[0]);
  pModel = CreatePModel(ALPHABET_SIZE);

  initNSymbol = nSymbol = 0;
  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      if((action = ParseMF(PA, (sym = readBuf[idxPos]))) < 0){
        switch(action){
          case -1: // IT IS THE BEGGINING OF THE HEADER
            if(PA->nRead>1 && ((PA->nRead-1) % P->nThreads == T.id) && nBase>1){
              #ifdef LOCAL_SIMILARITY
              if(P->local == 1){
                UpdateTopWP(BPBB(bits, nBase), conName, T.top, nBase,
                initNSymbol, nSymbol);
                }
              else
                UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
              #else
              UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
              #endif
              }
            #ifdef LOCAL_SIMILARITY
            initNSymbol = nSymbol;
            #endif
            // RESET MODELS 
            ResetCBuffer(symBuf);
            ResetShadowModel(Shadow);
            idx = r = nBase = bits = 0;
          break;
          case -2: conName[r] = '\0'; break; // IT IS THE '\n' HEADER END
          case -3: // IF IS A SYMBOL OF THE HEADER
            if(r >= MAX_NAME-1)
              conName[r] = '\0';
            else{
              if(sym == ' ' || sym < 32 || sym > 126){ // PROTECT INTERVAL
                if(r == 0) continue;
                else       sym = '_'; // PROTECT OUT SYM WITH UNDERL
                }
              conName[r++] = sym;
              }
          break;
          case -99: break; // IF IS A SIMPLE FORMAT BREAK
          default: exit(1);
          }
        continue; // GO TO NEXT SYMBOL
        }

      if(PA->nRead % P->nThreads == T.id){
        if((sym = DNASymToNum(sym)) == 4) continue; // IT IGNORES EXTRA SYMBOLS
        symBuf->buf[symBuf->idx] = sym;
        pos = &symBuf->buf[symBuf->idx-1];
        GetPModelIdx(pos, Shadow);

        if(idx++ % P->sample == 0 && idx > Shadow->ctx){
          ComputePModel(Models[0], pModel, Shadow->pModelIdx, Shadow->alphaDen);
          bits += PModelSymbolLog(pModel, sym);
          ++nBase;
          }

        UpdateCBuffer(symBuf);
        }
      }

  if(PA->nRead % P->nThreads == T.id){
    #ifdef LOCAL_SIMILARITY
    if(P->local == 1)
      UpdateTopWP(BPBB(bits, nBase), conName, T.top, nBase,
      initNSymbol, nSymbol);
    else
      UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
    #else
    UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
    #endif
    }

  RemovePModel(pModel);
  FreeShadow(Shadow);
  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Reader);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - F A L C O N   C O M P R E S S I O N - - - - - - - -

void FalconCompressTarget(Threads T){
  FILE        *Reader  = Fopen(P->base, "r");
  double      bits = 0;
  uint64_t    nBase = 0, r = 0, idx = 0, nSymbol, initNSymbol;
  uint32_t    k, idxPos;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t     sym, *pos, conName[MAX_NAME];
  PModel      *pModel;
  CModel      *Shadow; // SHADOW FOR SUPPORTING MODEL WITH THREADING
  int         action;

  Shadow = CreateShadowModel(Models[0]);
  pModel = CreatePModel(ALPHABET_SIZE);

  nSymbol = initNSymbol = 0;
  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      if((action = ParseMF(PA, (sym = readBuf[idxPos]))) < 0){
        switch(action){
          case -1: // IT IS THE BEGGINING OF THE HEADER
            if(PA->nRead>1 && ((PA->nRead-1) % P->nThreads == T.id && nBase>1)){
              #ifdef LOCAL_SIMILARITY
              if(P->local == 1){
                UpdateTopWP(BPBB(bits, nBase), conName, T.top, nBase,
                initNSymbol, nSymbol);
                }
              else
                UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
              #else
              UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
              #endif
              }
            #ifdef LOCAL_SIMILARITY
            initNSymbol = nSymbol;
            #endif
            // RESET MODELS 
            ResetCBuffer(symBuf);
            ResetShadowModel(Shadow);
            idx   = 0;
            r     = 0;
            nBase = 0; 
            bits  = 0;
          break;
          case -2: conName[r] = '\0'; break; // IT IS THE '\n' HEADER END
          case -3: // IF IS A SYMBOL OF THE HEADER
            if(r >= MAX_NAME-1)
              conName[r] = '\0';
            else{
              if(sym == ' ' || sym < 32 || sym > 126){ // PROTECT INTERVAL
                if(r == 0) continue;
                else       sym = '_'; // PROTECT OUT SYM WITH UNDERL
                }
              conName[r++] = sym;
              }
          break;
          case -99: break; // IF IS A SIMPLE FORMAT BREAK
          default: exit(1);
          }
        continue; // GO TO NEXT SYMBOL
        }

      if(PA->nRead % P->nThreads == T.id){
        if((sym = DNASymToNum(sym)) == 4) continue; // IT IGNORES EXTRA SYMBOLS
        symBuf->buf[symBuf->idx] = sym;
        pos = &symBuf->buf[symBuf->idx-1];
        GetPModelIdx(pos, Shadow);
        if(++idx >= Shadow->ctx){
          ComputePModel(Models[0], pModel, Shadow->pModelIdx, Shadow->alphaDen);
          bits += PModelSymbolLog(pModel, sym);
          ++nBase;
          }
        UpdateCBuffer(symBuf);
        }
      }

  if(PA->nRead % P->nThreads == T.id){
    #ifdef LOCAL_SIMILARITY
    if(P->local == 1)
      UpdateTopWP(BPBB(bits, nBase), conName, T.top, nBase,
      initNSymbol, nSymbol);
    else
      UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
    #else
    UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
    #endif
    }

  RemovePModel(pModel);
  FreeShadow(Shadow);
  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Reader);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - C O M P R E S S I O N   W I T H   K M O D E L S - - - - - - - -

void CompressTargetWKM(Threads T){
  FILE        *Reader = Fopen(P->base, "r");
  double      bits = 0;
  uint64_t    nBase = 0, r = 0, nSymbol, initNSymbol;
  uint32_t    n, k, idxPos, totModels, model;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t     sym, conName[MAX_NAME];
  PModel      **pModel, *MX;
  KMODEL      **Shadow; // SHADOWS FOR SUPPORTING MODELS WITH THREADING
  FloatPModel *PT;
  CMWeight    *CMW;
  int         action;

  totModels = P->nModels; // EXTRA MODELS DERIVED FROM EDITS
  for(n = 0 ; n < P->nModels ; ++n) 
    if(T.model[n].edits != 0)
      totModels += 1;

  Shadow      = (KMODEL **) Calloc(P->nModels, sizeof(KMODEL *));
  for(n = 0 ; n < P->nModels ; ++n)
    Shadow[n] = CreateKShadowModel(KModels[n]); 
  pModel      = (PModel **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n] = CreatePModel(ALPHABET_SIZE);
  MX          = CreatePModel(ALPHABET_SIZE);
  PT          = CreateFloatPModel(ALPHABET_SIZE);
  CMW         = CreateWeightModel(totModels);

  initNSymbol = nSymbol = 0;
  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      ++nSymbol;
      if((action = ParseMF(PA, (sym = readBuf[idxPos]))) < 0){
        switch(action){
          case -1: // IT IS THE BEGGINING OF THE HEADER
            if((PA->nRead-1) % P->nThreads == T.id && PA->nRead>1 && nBase>1){
              #ifdef LOCAL_SIMILARITY
              if(P->local == 1){
                UpdateTopWP(BPBB(bits, nBase), conName, T.top, nBase, 
                initNSymbol, nSymbol);
                }
              else
                UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
              #else
              UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
              #endif
              }
            #ifdef LOCAL_SIMILARITY
            initNSymbol = nSymbol; 
            #endif  
            ResetKModelsAndParam(symBuf, Shadow, CMW); // RESET MODELS
            r = nBase = bits = 0;
          break;
          case -2: conName[r] = '\0'; break; // IT IS THE '\n' HEADER END
          case -3: // IF IS A SYMBOL OF THE HEADER
            if(r >= MAX_NAME-1)
              conName[r] = '\0';
            else{
              if(sym == ' ' || sym < 32 || sym > 126){ // PROTECT INTERVAL
                if(r == 0) continue;
                else       sym = '_'; // PROTECT OUT SYM WITH UNDERL
                }
              conName[r++] = sym;
              }
          break; 
          case -99: break; // IF IS A SIMPLE FORMAT BREAK
          default: exit(1);
          }
        continue; // GO TO NEXT SYMBOL
        }

      if(PA->nRead % P->nThreads == T.id){

        if((sym = DNASymToNum(sym)) == 4){
          continue; // IT IGNORES EXTRA SYMBOLS
          }

        symBuf->buf[symBuf->idx] = sym;
        memset((void *)PT->freqs, 0, ALPHABET_SIZE * sizeof(double));
        n = 0;
        for(model = 0 ; model < P->nModels ; ++model){
          KMODEL *KM = Shadow[model];
          GetKIdx(symBuf->buf+symBuf->idx, KM);
          ComputeKPModel(KModels[model], pModel[n], KM->idx-sym, KM->alphaDen);
          ComputeWeightedFreqs(CMW->weight[n], pModel[n], PT);
          //if(KM->edits != 0){
          //  ++n;
          //  KM->SUBS.seq->buf[KM->SUBS.seq->idx] = sym;
          //  KM->SUBS.idx = GetPModelIdxCorr(KM->SUBS.seq->buf+KM->SUBS.seq->idx
          //  -1, KM, KM->SUBS.idx);
          //  ComputePModel(KModels[model], pModel[n], KM->SUBS.idx, KM->SUBS.eDen);
          //  ComputeWeightedFreqs(CMW->weight[n], pModel[n], PT);
          //  }
          ++n;
          }

        ComputeMXProbs(PT, MX);
        bits += PModelSymbolLog(MX, sym);
        ++nBase;
        CalcDecayment(CMW, pModel, sym, P->gamma);
        RenormalizeWeights(CMW);
        // CorrectXModels(Shadow, pModel, sym, P->nModels);
        UpdateCBuffer(symBuf);
        }
      }
        
  if(PA->nRead % P->nThreads == T.id){
    #ifdef LOCAL_SIMILARITY
    if(P->local == 1)
      UpdateTopWP(BPBB(bits, nBase), conName, T.top, nBase,
      initNSymbol, nSymbol);
    else
      UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
    #else
    UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
    #endif
    }

  DeleteWeightModel(CMW);
  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(pModel[n]);
  Free(pModel);
  RemovePModel(MX);
  RemoveFPModel(PT);
  for(n = 0 ; n < P->nModels ; ++n)
    FreeKShadow(Shadow[n]);
  Free(Shadow);
  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Reader);
  }
  

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S I O N - - - - - - - - - - - - - 

void CompressTarget(Threads T, char *dbFile){
  FILE        *Reader = CFopen(dbFile, "r");
  double      bits = 0;
  uint64_t    nBase = 0, r = 0, nSymbol, initNSymbol;
  uint32_t    n, k, idxPos, totModels, cModel;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t     sym, *pos, conName[MAX_NAME];
  PModel      **pModel, *MX;
  CModel      **Shadow; // SHADOWS FOR SUPPORTING MODELS WITH THREADING
  FloatPModel *PT;
  CMWeight    *CMW;
  int         action;

  totModels = P->nModels; // EXTRA MODELS DERIVED FROM EDITS
  for(n = 0 ; n < P->nModels ; ++n) 
    if(T.model[n].edits != 0)
      totModels += 1;

  Shadow      = (CModel **) Calloc(P->nModels, sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    Shadow[n] = CreateShadowModel(Models[n]); 
  pModel      = (PModel **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n] = CreatePModel(ALPHABET_SIZE);
  MX          = CreatePModel(ALPHABET_SIZE);
  PT          = CreateFloatPModel(ALPHABET_SIZE);
  CMW         = CreateWeightModel(totModels);

  initNSymbol = nSymbol = 0;
  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      ++nSymbol;
      if((action = ParseMF(PA, (sym = readBuf[idxPos]))) < 0){
        switch(action){
          case -1: // IT IS THE BEGGINING OF THE HEADER
            if((PA->nRead-1) % P->nThreads == T.id && PA->nRead>1 && nBase>1){
              #ifdef LOCAL_SIMILARITY
              if(P->local == 1){
                UpdateTopWPWithDb(BPBB(bits, nBase), conName, T.top, nBase,
                initNSymbol, nSymbol, P->currentDBIdx);
                }
              else
                UpdateTopWithDB(BPBB(bits, nBase), conName, T.top, nBase, P->currentDBIdx);
              #else
              UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
              #endif
              }
            #ifdef LOCAL_SIMILARITY
            initNSymbol = nSymbol; 
            #endif  
            ResetModelsAndParam(symBuf, Shadow, CMW); // RESET MODELS
            r = nBase = bits = 0;
          break;
          case -2: conName[r] = '\0'; break; // IT IS THE '\n' HEADER END
          case -3: // IF IS A SYMBOL OF THE HEADER
            if(r >= MAX_NAME-1)
              conName[r] = '\0';
            else{
              if(sym == ' ' || sym < 32 || sym > 126){ // PROTECT INTERVAL
                if(r == 0) continue;
                else       sym = '_'; // PROTECT OUT SYM WITH UNDERL
                }
              conName[r++] = sym;
              }
          break; 
          case -99: break; // IF IS A SIMPLE FORMAT BREAK
          default: exit(1);
          }
        continue; // GO TO NEXT SYMBOL
        }

      if(PA->nRead % P->nThreads == T.id){

        if((sym = DNASymToNum(sym)) == 4){
          continue; // IT IGNORES EXTRA SYMBOLS
          }

        symBuf->buf[symBuf->idx] = sym;
        memset((void *)PT->freqs, 0, ALPHABET_SIZE * sizeof(double));
        n = 0;
        pos = &symBuf->buf[symBuf->idx-1];
        for(cModel = 0 ; cModel < P->nModels ; ++cModel){
          CModel *CM = Shadow[cModel];
          GetPModelIdx(pos, CM);
          ComputePModel(Models[cModel], pModel[n], CM->pModelIdx, CM->alphaDen);
          ComputeWeightedFreqs(CMW->weight[n], pModel[n], PT);
          if(CM->edits != 0){
            ++n;
            CM->SUBS.seq->buf[CM->SUBS.seq->idx] = sym;
            CM->SUBS.idx = GetPModelIdxCorr(CM->SUBS.seq->buf+CM->SUBS.seq->idx
            -1, CM, CM->SUBS.idx);
            ComputePModel(Models[cModel], pModel[n], CM->SUBS.idx, CM->SUBS.eDen);
            ComputeWeightedFreqs(CMW->weight[n], pModel[n], PT);
            }
          ++n;
          }

        ComputeMXProbs(PT, MX);
        bits += PModelSymbolLog(MX, sym);
        ++nBase;
        CalcDecayment(CMW, pModel, sym, P->gamma);
        RenormalizeWeights(CMW);
        CorrectXModels(Shadow, pModel, sym, P->nModels);
        UpdateCBuffer(symBuf);
        }
      }
        
  if(PA->nRead % P->nThreads == T.id){
    #ifdef LOCAL_SIMILARITY
    if(P->local == 1)
      UpdateTopWPWithDb(BPBB(bits, nBase), conName, T.top, nBase,
      initNSymbol, nSymbol, P->currentDBIdx);
    else
      UpdateTopWithDB(BPBB(bits, nBase), conName, T.top, nBase, P->currentDBIdx);
    #else
    UpdateTop(BPBB(bits, nBase), conName, T.top, nBase);
    #endif
    }

  DeleteWeightModel(CMW);
  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(pModel[n]);
  Free(pModel);
  RemovePModel(MX);
  RemoveFPModel(PT);
  for(n = 0 ; n < P->nModels ; ++n)
    FreeShadow(Shadow[n]);
  Free(Shadow);
  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Reader);
  }

void CompressTargetInter(Threads T){
  FILE        *Reader  = Fopen(P->files[T.id], "r");
  double      *cModelWeight, cModelTotalWeight = 0, bits = 0, instance = 0;
  uint64_t    nBase = 0;
  uint32_t    n, k, idxPos, totModels, cModel;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t     sym, *pos;
  PModel      **pModel, *MX;
  CModel      **Shadow;
  FloatPModel *PT;

  totModels = P->nModels; // EXTRA MODELS DERIVED FROM EDITS
  for(n = 0 ; n < P->nModels ; ++n)
    if(T.model[n].edits != 0)
      totModels += 1;

  Shadow = (CModel **) Calloc(P->nModels, sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    Shadow[n] = CreateShadowModel(Models[n]);
  pModel        = (PModel  **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n]   = CreatePModel(ALPHABET_SIZE);
  MX            = CreatePModel(ALPHABET_SIZE);
  PT            = CreateFloatPModel(ALPHABET_SIZE);
  cModelWeight  = (double   *) Calloc(totModels, sizeof(double));

  for(n = 0 ; n < totModels ; ++n)
    cModelWeight[n] = 1.0 / totModels;

  FileType(PA, Reader);

  nBase = 0;
  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      if(ParseSym(PA, (sym = readBuf[idxPos])) == -1) continue;
      symBuf->buf[symBuf->idx] = sym = DNASymToNum(sym);

      memset((void *)PT->freqs, 0, ALPHABET_SIZE * sizeof(double));

      n = 0;
      pos = &symBuf->buf[symBuf->idx-1];
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        GetPModelIdx(pos, Shadow[cModel]);
        ComputePModel(Models[cModel], pModel[n], Shadow[cModel]->pModelIdx,
        Shadow[cModel]->alphaDen);
        ComputeWeightedFreqs(cModelWeight[n], pModel[n], PT);
        if(Shadow[cModel]->edits != 0){
          ++n;
          Shadow[cModel]->SUBS.seq->buf[Shadow[cModel]->SUBS.seq->idx] = sym;
          Shadow[cModel]->SUBS.idx = GetPModelIdxCorr(Shadow[cModel]->SUBS.
          seq->buf+Shadow[cModel]->SUBS.seq->idx-1, Shadow[cModel], Shadow
          [cModel]->SUBS.idx);
          ComputePModel(Models[cModel], pModel[n], Shadow[cModel]->SUBS.idx,
          Shadow[cModel]->SUBS.eDen);
          ComputeWeightedFreqs(cModelWeight[n], pModel[n], PT);
          }
        ++n;
        }

      MX->sum  = (MX->freqs[0] = 1 + (unsigned) (PT->freqs[0] * MX_PMODEL));
      MX->sum += (MX->freqs[1] = 1 + (unsigned) (PT->freqs[1] * MX_PMODEL));
      MX->sum += (MX->freqs[2] = 1 + (unsigned) (PT->freqs[2] * MX_PMODEL));
      MX->sum += (MX->freqs[3] = 1 + (unsigned) (PT->freqs[3] * MX_PMODEL));
      bits += (instance = PModelSymbolLog(MX, sym));
      nBase++;

      cModelTotalWeight = 0;
      for(n = 0 ; n < totModels ; ++n){
        cModelWeight[n] = Power(cModelWeight[n], P->gamma) * (double)
        pModel[n]->freqs[sym] / pModel[n]->sum;
        cModelTotalWeight += cModelWeight[n];
        }

      for(n = 0 ; n < totModels ; ++n)
        cModelWeight[n] /= cModelTotalWeight; // RENORMALIZE THE WEIGHTS

      n = 0;
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        if(Shadow[cModel]->edits != 0){
          CorrectCModelSUBS(Shadow[cModel], pModel[++n], sym);
          }
        ++n;
        }

      UpdateCBuffer(symBuf);
      }

  Free(cModelWeight);
  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(pModel[n]);
  Free(pModel);
  RemovePModel(MX);
  RemoveFPModel(PT);
  for(n = 0 ; n < P->nModels ; ++n)
    FreeShadow(Shadow[n]);
  Free(Shadow);
  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Reader);

  P->matrix[P->ref][T.id] = nBase == 0 ? 101 : bits / 2 / nBase; // 101 -> nan
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - F   T H R E A D I N G - - - - - - - - - - - - - - -

void *CompressThread(void *Thr){
  Threads *T = (Threads *) Thr;

  // Use the current database
  char *currentDb = P->dbFiles[P->currentDBIdx];
  
  //  if(P->nModels == 1 && T->model[0].edits == 0){
  //    if(P->sample > 1){
  //      SamplingCompressTarget(T[0]);
  //      }
  //    else{      
  //      FalconCompressTarget(T[0]);
  //      }
  //    }

  #ifdef KMODELSUSAGE
  CompressTargetWKM(T[0]);
  #else
  CompressTarget(T[0], currentDb);
  #endif

  pthread_exit(NULL);
  }

void *CompressThreadInter(void *Thr){
  Threads *T = (Threads *) Thr;
  CompressTargetInter(T[0]);
  pthread_exit(NULL);
}


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E F E R E N C E - - - - - - - - - - - - -

void LoadReference(char *refName){
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

  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){

      if(ParseSym(PA, (sym = readBuf[idxPos])) == -1){ 
        idx = 0; 
	continue; 
        }

      symBuf->buf[symBuf->idx] = sym = DNASymToNum(sym);

      for(n = 0 ; n < P->nModels ; ++n){
        CModel *CM = Models[n];
        GetPModelIdx(symBuf->buf+symBuf->idx-1, CM);
        if(CM->ir == 1) // INVERTED REPEATS
          irSym = GetPModelIdxIR(symBuf->buf+symBuf->idx, CM);
        if(idx >= CM->ctx){
          UpdateCModelCounter(CM, sym, CM->pModelIdx);
          if(CM->ir == 1) // INVERTED REPEATS
            UpdateCModelCounter(CM, irSym, CM->pModelIdxIR);
          }
        }

      ++idx;
      UpdateCBuffer(symBuf);
      }
 
  for(n = 0 ; n < P->nModels ; ++n)
    ResetCModelIdx(Models[n]);
  RemoveCBuffer(symBuf);
  Free(readBuf);
  RemoveParser(PA);
  fclose(Reader);
  }

void LoadReferenceInter(Threads T){
  FILE     *Reader = Fopen(P->files[T.id], "r");
  uint32_t n;
  uint64_t idx = 0;
  uint64_t k, idxPos;
  PARSER   *PA = CreateParser();
  CBUF     *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t  sym, irSym, *readBuf = Calloc(BUFFER_SIZE, sizeof(uint8_t));
  FileType(PA, Reader);
  rewind(Reader);

  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      if(ParseSym(PA, (sym = readBuf[idxPos])) == -1){ idx = 0; continue; }
      symBuf->buf[symBuf->idx] = sym = DNASymToNum(sym);
      for(n = 0 ; n < P->nModels ; ++n){
        CModel *CM = Models[n];
        GetPModelIdx(symBuf->buf+symBuf->idx-1, CM);
        if(++idx > CM->ctx){
          UpdateCModelCounter(CM, sym, CM->pModelIdx);
          if(CM->ir == 1){                         // INVERTED REPEATS
            irSym = GetPModelIdxIR(symBuf->buf+symBuf->idx, CM);
            UpdateCModelCounter(CM, irSym, CM->pModelIdxIR);
          }
        }
      }
      UpdateCBuffer(symBuf);
    }

  for(n = 0 ; n < P->nModels ; ++n)
    ResetCModelIdx(Models[n]);
  RemoveCBuffer(symBuf);
  Free(readBuf);
  RemoveParser(PA);
  fclose(Reader);
}

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - R E F E R E N C E   W I T H   K M O D E L S - - - - - - - -

void LoadReferenceWKM(char *refName){
  FILE     *Reader = Fopen(refName, "r");
  uint32_t n;
  uint64_t idx = 0, k, idxPos;
  PARSER   *PA = CreateParser();
  CBUF     *symBuf  = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t  *readBuf = Calloc(BUFFER_SIZE, sizeof(uint8_t)), sym;
  FileType(PA, Reader);
  rewind(Reader);

  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
 
      if(ParseSym(PA, (sym = readBuf[idxPos])) == -1){ 
        idx = 0; 
	continue; 
        }

      symBuf->buf[symBuf->idx] = sym = DNASymToNum(sym);

      for(n = 0 ; n < P->nModels ; ++n){
        KMODEL *KM = KModels[n];
        GetKIdx(symBuf->buf+symBuf->idx, KM);
        if(idx >= KM->ctx)
          UpdateKModelCounter(KM, sym, KM->idx);
        }

      ++idx;
      UpdateCBuffer(symBuf);
      }
 
  for(n = 0 ; n < P->nModels ; ++n)
    ResetKModelIdx(KModels[n]);
  RemoveCBuffer(symBuf);
  Free(readBuf);
  RemoveParser(PA);
  fclose(Reader);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - C O M P R E S S O R   M A I N - - - - - - - - - - - -

void CompressAction(Threads *T, char *refName, char *baseName) {
  pthread_t t[P->nThreads];
  uint32_t n, dbIdx;
  char     filteredFile[MAX_NAME]; // Enough space for filename
  char     concatFile[MAX_NAME]; // Enough space for filename
  int      useMagnetFilter = 0;

  if(P->loadModel) {
    // Load models from file instead of building them
    fprintf(stderr, "  [+] Loading models from file %s ... ", P->modelFile);
    int result = LoadModels(P->modelFile, &Models, &P->nModels, &P->col);
    if(result != 0) {
      fprintf(stderr, "Error loading models (code: %d)\n", result);
      exit(1);
    }
    fprintf(stderr, "Done!\n");
  } else {
    // Build models from reference files as usual
#ifdef KMODELSUSAGE
    KModels = (KMODEL **) Malloc(P->nModels * sizeof(KMODEL *));
    for(n = 0 ; n < P->nModels ; ++n)
      KModels[n] = CreateKModel(T[0].model[n].ctx, T[0].model[n].den,
      T[0].model[n].ir, REFERENCE, P->col, T[0].model[n].edits,
      T[0].model[n].eDen);
    fprintf(stderr, "  [+] Loading metagenomic file ..... ");
    LoadReferenceWKM(refName);
    fprintf(stderr, "Done!\n");
#else
    Models = (CModel **) Malloc(P->nModels * sizeof(CModel *));
    for(n = 0 ; n < P->nModels ; ++n)
      Models[n] = CreateCModel(T[0].model[n].ctx, T[0].model[n].den,
      T[0].model[n].ir, REFERENCE, P->col, T[0].model[n].edits,
      T[0].model[n].eDen);
    fprintf(stderr, "  [+] Loading %u metagenomic file(s):\n", P->nFiles);

    if(P->useMagnet) {
      useMagnetFilter = 1;
      strcpy(filteredFile, "falcon_magnet_filtered.fq");

      // Create a temporary filename for filtered reference different
      //snprintf(filteredFile, sizeof(filteredFile), "%s.magnet_filtered.fq", P->files[n]);

      fprintf(stderr, "      [+] Applying MAGNET filter to %u metagenomic file(s):\n", P->nFiles);

      if(P->nFiles > 1) {
        strcpy(concatFile, "falcon_concatenated.fq");
        ConcatWithCFopen(P->files, P->nFiles, concatFile);

        // Run MAGNET to filter the reference file
        int result = RunMagnet(
        concatFile,
        P->magnetFilter,
        P->magnetThreshold,
        P->magnetLevel,
        P->magnetInvert,
        P->magnetVerbose,
        P->magnetPortion,
        filteredFile,
        P->nThreads);

        if(result != 0) {
          fprintf(stderr, "      [+] MAGNET filtering failed with code %d.\n", result);
          for(n = 0 ; n < P->nFiles ; ++n){
            fprintf(stderr, "      [+] Loading original reference %s ... ", P->files[n]);
            fprintf(stderr, "      [+] Loading %u ... ", n+1);
            LoadReference(P->files[n]);
            fprintf(stderr, "Done! \n");
          }
        }
        else {
          fprintf(stderr, "      [+] Loading filtered reference %s ... ", filteredFile);
          // Load the filtered reference file
          LoadReference(filteredFile);
          fprintf(stderr, "Done!\n");
          useMagnetFilter = 1;
        }
      } else {
        // Run MAGNET to filter the reference file
        int result = RunMagnet(
        P->files[0],
        P->magnetFilter,
        P->magnetThreshold,
        P->magnetLevel,
        P->magnetInvert,
        P->magnetVerbose,
        P->magnetPortion,
        filteredFile,
        P->nThreads);

        if(result != 0) {
          fprintf(stderr, "      [+] MAGNET filtering failed with code %d.\n", result);
          fprintf(stderr, "      [+] Using original reference file %s instead.\n", P->files[0]);
          fprintf(stderr, "      [+] Loading original reference %s ... ", P->files[0]);
          // Load the original reference file
          LoadReference(P->files[0]);
          fprintf(stderr, "Done!\n");
        }
        else {
          fprintf(stderr, "      [+] Loading filtered reference %s ... ", filteredFile);
          // Load the filtered reference file
          LoadReference(filteredFile);
          fprintf(stderr, "Done!\n");
        }
      }
    } else {
      for(n = 0 ; n < P->nFiles ; ++n){
        fprintf(stderr, "      [+] Loading %u ... ", n+1);
        LoadReference(P->files[n]);
        fprintf(stderr, "Done! \n");
      }
    }
    fprintf(stderr, "  [+] Done! Learning phase complete!\n");

    // Save models if requested
    if(P->saveModel) {
      fprintf(stderr, "  [+] Saving models to file %s ... ", P->modelFile);
      int result = SaveModels(P->modelFile, Models, P->nModels, P->col);
      if(result != 0) {
        fprintf(stderr, "Error saving models (code: %d)\n", result);
        exit(1);
      }
      fprintf(stderr, "Done!\n");
    }
#endif

    fprintf(stderr, "  [+] Compressing database ......... %u file(s):\n", P->nDatabases);

    for(dbIdx = 0 ; dbIdx < P->nDatabases ; ++dbIdx){
      fprintf(stderr, "      [+] Loading %u ... ", dbIdx+1);

      // Set current database for threads
      P->currentDBIdx = dbIdx;

      for(n = 0 ; n < P->nThreads ; ++n)
        pthread_create(&(t[n+1]), NULL, CompressThread, (void *) &(T[n]));
      for(n = 0 ; n < P->nThreads ; ++n) // DO NOT JOIN FORS!
        pthread_join(t[n+1], NULL);
      fprintf(stderr, "Done!\n");

    }

    if(useMagnetFilter) {
      // Remove the filtered file if it was created
      if(remove(filteredFile) != 0) {
        fprintf(stderr, "Warning: Could not remove temporary file %s\n", filteredFile);
      }
    }
  }
}

void CompressActionTraining(Threads *T, char *refName){
  uint32_t n;
  char     filteredFile[MAX_NAME]; // Enough space for filename
  char     concatFile[MAX_NAME]; // Enough space for filename
  int      useMagnetFilter = 0;

#ifdef KMODELSUSAGE
    KModels = (KMODEL **) Malloc(P->nModels * sizeof(KMODEL *));
    for(n = 0 ; n < P->nModels ; ++n)
      KModels[n] = CreateKModel(T[0].model[n].ctx, T[0].model[n].den,
      T[0].model[n].ir, REFERENCE, P->col, T[0].model[n].edits,
      T[0].model[n].eDen);
    fprintf(stderr, "  [+] Loading metagenomic file ..... ");
    LoadReferenceWKM(refName);
    fprintf(stderr, "Done!\n");
#else
    Models = (CModel **) Malloc(P->nModels * sizeof(CModel *));
    for(n = 0 ; n < P->nModels ; ++n)
      Models[n] = CreateCModel(T[0].model[n].ctx, T[0].model[n].den,
      T[0].model[n].ir, REFERENCE, P->col, T[0].model[n].edits,
      T[0].model[n].eDen);
    fprintf(stderr, "  [+] Loading %u metagenomic file(s):\n", P->nFiles);

        if(P->useMagnet) {
      useMagnetFilter = 1;
      strcpy(filteredFile, "falcon_magnet_filtered.fq");

      // Create a temporary filename for filtered reference different
      //snprintf(filteredFile, sizeof(filteredFile), "%s.magnet_filtered.fq", P->files[n]);

      fprintf(stderr, "      [+] Applying MAGNET filter to %u metagenomic file(s):\n", P->nFiles);

      if(P->nFiles > 1) {
        strcpy(concatFile, "falcon_concatenated.fq");
        ConcatWithCFopen(P->files, P->nFiles, concatFile);

        // Run MAGNET to filter the reference file
        int result = RunMagnet(
        concatFile,
        P->magnetFilter,
        P->magnetThreshold,
        P->magnetLevel,
        P->magnetInvert,
        P->magnetVerbose,
        P->magnetPortion,
        filteredFile,
        P->nThreads);

        if(result != 0) {
          fprintf(stderr, "      [+] MAGNET filtering failed with code %d.\n", result);
          for(n = 0 ; n < P->nFiles ; ++n){
            fprintf(stderr, "      [+] Loading original reference %s ... ", P->files[n]);
            fprintf(stderr, "      [+] Loading %u ... ", n+1);
            LoadReference(P->files[n]);
            fprintf(stderr, "Done! \n");
          }
        }
        else {
          fprintf(stderr, "      [+] Loading filtered reference %s ... ", filteredFile);
          // Load the filtered reference file
          LoadReference(filteredFile);
          fprintf(stderr, "Done!\n");
          useMagnetFilter = 1;
        }
      } else {
        // Run MAGNET to filter the reference file
        int result = RunMagnet(
        P->files[0],
        P->magnetFilter,
        P->magnetThreshold,
        P->magnetLevel,
        P->magnetInvert,
        P->magnetVerbose,
        P->magnetPortion,
        filteredFile,
        P->nThreads);

        if(result != 0) {
          fprintf(stderr, "      [+] MAGNET filtering failed with code %d.\n", result);
          fprintf(stderr, "      [+] Using original reference file %s instead.\n", P->files[0]);
          fprintf(stderr, "      [+] Loading original reference %s ... ", P->files[0]);
          // Load the original reference file
          LoadReference(P->files[0]);
          fprintf(stderr, "Done!\n");
        }
        else {
          fprintf(stderr, "      [+] Loading filtered reference %s ... ", filteredFile);
          // Load the filtered reference file
          LoadReference(filteredFile);
          fprintf(stderr, "Done!\n");
        }
      }
    } else {
      for(n = 0 ; n < P->nFiles ; ++n){
        fprintf(stderr, "      [+] Loading %u ... ", n+1);
        LoadReference(P->files[n]);
        fprintf(stderr, "Done! \n");
      }
    }
    fprintf(stderr, "  [+] Done! Learning phase complete!\n");

    // Save models
    fprintf(stderr, "  [+] Saving models to file %s ... ", P->modelFile);
    int result = SaveModels(P->modelFile, Models, P->nModels, P->col);
    if(result != 0) {
      fprintf(stderr, "Error saving models (code: %d)\n", result);
      exit(1);
    }
    fprintf(stderr, "Done!\n");

#endif

  if(useMagnetFilter) {
    // Remove the filtered file if it was created
    if(remove(filteredFile) != 0) {
      fprintf(stderr, "Warning: Could not remove temporary file %s\n", filteredFile);
    }
  }
}

void CompressActionInter(Threads *T, uint32_t ref){
  uint32_t n, k;
  pthread_t t[P->nThreads];
  P->ref = ref;

  Models = (CModel **) Malloc(P->nModels * sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    Models[n] = CreateCModel(T[ref].model[n].ctx, T[ref].model[n].den,
    T[ref].model[n].ir, REFERENCE, P->col, T[ref].model[n].edits,
    T[ref].model[n].eDen);

  fprintf(stderr, "  [+] Loading reference %u ... ", ref+1);
  LoadReferenceInter(T[ref]);
  fprintf(stderr, "Done!\n");

  fprintf(stderr, "      [+] Compressing %u targets ... ", P->nFiles);
  ref = 0;
  do{
    for(n = 0 ; n < P->nThreads ; ++n)
      pthread_create(&(t[n+1]), NULL, CompressThreadInter, (void *) &(T[ref+n]));
    for(n = 0 ; n < P->nThreads ; ++n) // DO NOT JOIN FORS!
      pthread_join(t[n+1], NULL);
  }
  while((ref += P->nThreads) < P->nFiles && ref + P->nThreads <= P->nFiles);

  if(ref < P->nFiles){ // EXTRA - OUT OF THE MAIN LOOP
    for(n = ref, k = 0 ; n < P->nFiles ; ++n)
      pthread_create(&(t[++k]), NULL, CompressThreadInter, (void *) &(T[n]));
    for(n = ref, k = 0 ; n < P->nFiles ; ++n) // DO NOT JOIN FORS!
      pthread_join(t[++k], NULL);
  }
  fprintf(stderr, "Done!\n");

  for(n = 0 ; n < P->nModels ; ++n)
    FreeCModel(Models[n]);
  Free(Models);
}

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - R E A D   L A B E L S - - - - - - - - - - - - -

void ReadLabels(void){
  if(access(P->labels, F_OK ) == -1)
    return;
  FILE *F = Fopen(P->labels, "r");
  uint32_t n;
  for(n = 0 ; n < P->nFiles ; ++n){
    if(!fscanf(F, "%s", P->files[n])){
      fclose(F);
      return;
      }
    }
  fclose(F);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - M A T R I X   S I Z E - - - - - - - - - - - - -

void ReadMatrixSize(char *fn){
  FILE *F = Fopen(fn, "r");
  char c;
  P->nFiles = 0;
  while((c = fgetc(F)) != '\n'){
    if(c == EOF){
      fclose(F);
      break;
      }
    if(c == '\t')
      ++P->nFiles;
    }
  if(P->nFiles == 0){
    fprintf(stderr, "[x] Error: invalid file!\n");
    exit(1);
    }
  fprintf(stderr, "  [>] Matrix size: %u\n", P->nFiles);
  fclose(F);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - M A T R I X   S I Z E - - - - - - - - - - - - -

void ReadMatrix(char *fn){
  FILE *F = Fopen(fn, "r");
  uint32_t n, k;
  for(n = 0 ; n < P->nFiles ; ++n){
    for(k = 0 ; k < P->nFiles ; ++k){
      if(!fscanf(F, "%lf", &P->matrix[n][k])){
        fclose(F);
        return;
        }
      }
    }
  fclose(F);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - P A I N T - - - - - - - - - - - - - - - -

void PaintMatrix(double start, double rotations, double hue, double gamma,
double width, double space){
  FILE *Plot = Fopen(P->image, "w");
  Painter *Paint;
  COLORS  *CLR = (COLORS *) Calloc(1, sizeof(COLORS));
  uint32_t ref, tar;

  CLR->start     = start;
  CLR->rotations = rotations;
  CLR->hue       = hue;
  CLR->gamma     = gamma;

  Paint = CreateBasicPainter(DEFAULT_CX*2+((width+space)*P->nFiles)+5, width, space);

  PrintHead(Plot, (2 * DEFAULT_CX) + (((Paint->width + Paint->space) *
  P->nFiles) - Paint->space), Paint->size + EXTRA);

  Rect(Plot, (2 * DEFAULT_CX) + (((Paint->width + Paint->space) * P->nFiles) -
  Paint->space), Paint->size + EXTRA, 0, 0, "#ffffff");

  // PRINT HEATMAP SCALE
  uint32_t size = (Paint->width + Paint->space) * P->nFiles - Paint->space;
  for(ref = 0 ; ref < size ; ++ref){
    char color[12];
    Rect(Plot, Paint->width, 1, DEFAULT_CX - (Paint->width*2), Paint->cy + ref,
    HeatMapColor(((double) ref / size), color, CLR));
    }
  if(P->nFiles > 4)
    Text90d(Plot, -DEFAULT_CX - ((size/2)+46), DEFAULT_CX-(Paint->width*2 + 2),
    "SIMILARITY");
  Text   (Plot, DEFAULT_CX-(Paint->width*2 + 14), Paint->cy+13, "+");
  Text   (Plot, DEFAULT_CX-(Paint->width*2 + 14), Paint->cy+size, "-");

  for(ref = 0 ; ref < P->nFiles ; ++ref){
    //for(tar = P->nFiles ; tar-- ; ){ // INVERT LOOP TO INVERT HORIZONTAL ORDER
    for(tar = 0 ; tar < P->nFiles ; ++tar){
      char color[12];
      Rect(Plot, Paint->width, Paint->width, Paint->cx, Paint->cy,
      HeatMapColor(BoundDouble(0.0, P->matrix[ref][tar], 1.0), color, CLR));
      Paint->cx += Paint->width + Paint->space;
      }
    // TEXT HAS 16 PX -> CALCULATE AVERAGE POSITION
    Text   (Plot, Paint->cx + 4, (Paint->cy+Paint->width/2)+6, P->files[ref]);

    // USE THE FOLLOWING TO INVERT THE NAMES ORDER -> THEN CHANGE IN THE SQUARES
    //Text90d(Plot, 4-DEFAULT_CX, (Paint->cy+Paint->width/2)+10, P->files[P->nFiles-1-ref]);
    Text90d(Plot, 4-DEFAULT_CX, (Paint->cy+Paint->width/2)+10, P->files[ref]);

    Paint->cx =  DEFAULT_CX;
    Paint->cy += Paint->width + Paint->space;
    }

  PrintFinal(Plot);
  RemovePainter(Paint);
  }


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int32_t P_Falcon(char **argv, int argc){
  char     **p = *&argv, **xargv, *xpl = NULL;
  int32_t  xargc = 0;
  uint32_t n, k, col, ref, topSize;
  double   gamma;
  Threads  *T;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEFAULT_HELP, p, argc, "-h", "--help")) == 1 || argc < 2){
    PrintMenuFalcon();
    Free(P);
    return EXIT_SUCCESS;
  }

  if(ArgsState(DEF_VERSION, p, argc, "-V", "--version")){
    PrintVersion();
    Free(P);
    return EXIT_SUCCESS;
  }

  if(ArgsState(0, p, argc, "-s", "--show")){
    PrintLevels();
    Free(P);
    return EXIT_SUCCESS;
  }

  P->verbose  = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v", "--verbose");
  P->force    = ArgsState  (DEFAULT_FORCE,   p, argc, "-F", "--force");
  #ifdef LOCAL_SIMILARITY
  P->local    = ArgsState  (DEFAULT_LOCAL,   p, argc, "-Z", "--local");
  #endif
  P->sample   = ArgsNum    (DEFAULT_SAMPLE,  p, argc, "-p", MIN_SAP, MAX_SAP);
  P->level    = ArgsNum    (0,               p, argc, "-l", MIN_LEV, MAX_LEV);
  topSize     = ArgsNum    (DEF_TOP,         p, argc, "-t", MIN_TOP, MAX_TOP);
  P->nThreads = ArgsNum    (DEFAULT_THREADS, p, argc, "-n", MIN_THREADS,
  MAX_THREADS);
  
  // Magnet Integration Flags
  P->useMagnet       = ArgsState  (0, p, argc, "-mg", "--magnet");
  P->magnetVerbose   = ArgsState  (DEFAULT_VERBOSE, p, argc, "-mv", "--magnet-verbose");
  P->magnetFilter    = ArgsString (NULL, p, argc, "-mf", "--magnet-filter");
  P->magnetThreshold = ArgsDouble (0.9, p, argc, "-mt");
  P->magnetLevel     = ArgsNum    (36, p, argc, "-ml", MIN_LEV, 44);
  P->magnetInvert    = ArgsState  (0, p, argc, "-mi", "--magnet-invert");
  P->magnetPortion   = ArgsNum    (1, p, argc, "-mp", MIN_SAP, MAX_SAP);

  // Model Saving and Loading Flags
  P->saveModel  = ArgsState  (0, p, argc, "-S", "--save-model");
  P->loadModel  = ArgsState  (0, p, argc, "-L", "--load-model");
  P->modelInfo  = ArgsState  (0, p, argc, "-I", "--model-info");
  P->trainModel = ArgsState  (0, p, argc, "-T", "--train-model");
  P->modelFile  = ArgsFileGen(p, argc, "-M", "falcon_model", ".fcm"); // FCM = Falcon Compression Model

  if(P->loadModel){
    if(P->modelFile == NULL || strlen(P->modelFile) == 0){
      fprintf(stderr,
        "Error: Model loading enabled but no model file specified.\n"
              "Please provide a model file with -M option.\n");
      Free(P);
      exit(1);
    }
    TestReadFile(P->modelFile);
  }

  if(P->useMagnet) {
    if(P->magnetFilter == NULL || strlen(P->magnetFilter) == 0){
      fprintf(stderr,
        "Error: MAGNET filtering enabled but no filter file specified.\n"
              "Please provide a FASTA filter file with -mf option.\n");
      Free(P);
      exit(1);
    }

    TestReadFile(P->magnetFilter);

    if (P->verbose) {
      PrintMagnetVersion();
    }
  }

  P->nModels = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-m") == 0)
      P->nModels += 1;

  if(P->nModels == 0 && P->level == 0)
    P->level = DEFAULT_LEVEL;

  if(P->level != 0){
    xpl = GetLevels(P->level);
    xargc = StrToArgv(xpl, &xargv);
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-m") == 0)
        P->nModels += 1;
    }

  gamma = DEFAULT_GAMMA;
  for(n = 1 ; n < xargc ; ++n)
    if(strcmp(xargv[n], "-g") == 0)
      gamma = atof(xargv[n+1]);

  col = MAX_COLLISIONS;
  for(n = 1 ; n < xargc ; ++n)
    if(strcmp(xargv[n], "-c") == 0)
      col = atoi(xargv[n+1]);

  P->col       = ArgsNum    (col,   p, argc, "-c", 1, 253);
  P->gamma     = ArgsDouble (gamma, p, argc, "-g");
  P->gamma     = ((int) (P->gamma * 65536)) / 65536.0;
  P->output    = ArgsFileGen(p, argc, "-x", "top", ".csv");
  #ifdef LOCAL_SIMILARITY
  if(P->local == 1){
    P->outLoc  = ArgsFileGen(p, argc, "-y", "local", ".fal");
    }
  #endif

  FILE *OUTLOC = NULL;
  FILE *OUTPUT = NULL;

  TIME *Time = NULL;

  if(P->trainModel) {

    if(P->nModels == 0){
      fprintf(stderr, "Error: at least you need to use a context model!\n");
      exit(1);
    }

    // READ MODEL PARAMETERS FROM XARGS & ARGS
    T = (Threads *) Calloc(P->nThreads, sizeof(Threads));
    for(ref = 0 ; ref < P->nThreads ; ++ref){
      T[ref].model = (ModelPar *) Calloc(P->nModels, sizeof(ModelPar));
      T[ref].id    = ref;
      k = 0;
      for(n = 1 ; n < argc ; ++n)
        if(strcmp(argv[n], "-m") == 0)
          T[ref].model[k++] = ArgsUniqModel(argv[n+1], 0);
      if(P->level != 0){
        for(n = 1 ; n < xargc ; ++n)
          if(strcmp(xargv[n], "-m") == 0)
            T[ref].model[k++] = ArgsUniqModel(xargv[n+1], 0);
      }
    }

    P->nFiles     = ReadFNames (P, argv[argc-1], 0);
    fprintf(stderr, "\n");
    if(P->verbose) PrintArgsTrain(P, T[0], argv[argc-1]);

    fprintf(stderr, "==[ PROCESSING ]====================\n");
    Time = CreateClock(clock());

    CompressActionTraining(T, argv[argc-1]);

    fprintf(stderr, "  [+] Freeing compression models ... ");
    for(n = 0 ; n < P->nModels ; ++n)
#ifdef KMODELSUSAGE
      FreeKModel(KModels[n]);
#else
    FreeCModel(Models[n]);
#endif
#ifdef KMODELSUSAGE
    Free(KModels);
#else
    Free(Models);
#endif
    fprintf(stderr, "Done!\n");

    StopTimeNDRM(Time, clock());
    fprintf(stderr, "\n");

    fprintf(stderr, "==[ STATISTICS ]====================\n");
    StopCalcAll(Time, clock());
    fprintf(stderr, "\n");

    if (xargv) {
      if (xargv[0]) {
        Free(xargv[0]);  // Free the duplicated string
      }
      Free(xargv);         // Free the array itself
    }

    // Free P resources
    if (P->output) {
      Free(P->output);
    }

    // Free file names
    if (P->files) {
      Free(P->files);
    }

    // Free model file name
    if (P->modelFile) {
      Free(P->modelFile);
    }

    RemoveClock(Time);
    for(ref = 0 ; ref < P->nThreads ; ++ref){
      Free(T[ref].model);
    }
    Free(T);
    Free(P);

    return EXIT_SUCCESS;
  }

  if(!P->force)
    FAccessWPerm(P->output);
  OUTPUT = Fopen(P->output, "w");

#ifdef LOCAL_SIMILARITY
  if(P->local == 1){
    if(!P->force)
      FAccessWPerm(P->outLoc);
    OUTLOC = Fopen(P->outLoc, "w");
  }
#endif

  if(P->nModels == 0){
    fprintf(stderr, "Error: at least you need to use a context model!\n");
    exit(1);
  }

  // READ MODEL PARAMETERS FROM XARGS & ARGS
  T = (Threads *) Calloc(P->nThreads, sizeof(Threads));
  for(ref = 0 ; ref < P->nThreads ; ++ref){
    T[ref].model = (ModelPar *) Calloc(P->nModels, sizeof(ModelPar));
    T[ref].id    = ref;
    T[ref].top   = CreateTop(topSize);
    k = 0;
    for(n = 1 ; n < argc ; ++n)
      if(strcmp(argv[n], "-m") == 0)
        T[ref].model[k++] = ArgsUniqModel(argv[n+1], 0);
    if(P->level != 0){
      for(n = 1 ; n < xargc ; ++n)
        if(strcmp(xargv[n], "-m") == 0)
          T[ref].model[k++] = ArgsUniqModel(xargv[n+1], 0);
    }
  }

  P->base = argv[argc-1];
  P->nDatabases = ReadDBFNames (P, argv[argc-1], 0);
  P->nFiles     = ReadFNames (P, argv[argc-2], 0);
  fprintf(stderr, "\n");
  if(P->verbose) PrintArgs(P, T[0], argv[argc-2], argv[argc-1], topSize);

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  Time = CreateClock(clock());
  CompressAction(T, argv[argc-2], P->base);

  k = 0;
  P->top = CreateTop(topSize * P->nThreads);
  for(ref = 0 ; ref < P->nThreads ; ++ref){
    for(n = 0 ; n < T[ref].top->size-1 ; ++n){
      P->top->V[k].value = T[ref].top->V[n].value;
      P->top->V[k].size  = T[ref].top->V[n].size;
      P->top->V[k].dbIndex = T[ref].top->V[n].dbIndex;
      #ifdef LOCAL_SIMILARITY
      if(P->local == 1){
        P->top->V[k].iPos = T[ref].top->V[n].iPos;
        P->top->V[k].ePos = T[ref].top->V[n].ePos;
        }
      #endif
      CopyStringPart(P->top->V[k].name, T[ref].top->V[n].name);
      ++k;
      }
    }

  fprintf(stderr, "  [+] Sorting top .................. ");
  qsort(P->top->V, k, sizeof(VT), SortByValue);
  fprintf(stderr, "Done!\n");

  fprintf(stderr, "  [+] Printing to output file ...... ");
  #ifdef LOCAL_SIMILARITY
  if(P->local == 1)
    PrintTopWP(OUTPUT, P->top, topSize, P->dbFiles);
  else
    PrintTop(OUTPUT, P->top, topSize, P->dbFiles);
  #else
  PrintTop(OUTPUT, P->top, topSize);
  #endif
  fclose(OUTPUT);
  fprintf(stderr, "Done!\n");

  #ifdef LOCAL_SIMILARITY
  if(P->local == 1){
    fprintf(stderr, "  [+] Running local similarity:\n");
    #ifdef KMODELSUSAGE
    LocalComplexityWKM(T[0], P->top, topSize, OUTLOC);
    #else
    LocalComplexityVariousDBs(T[0], P->top, topSize, P->dbFiles, P->nDatabases, OUTLOC);
    fprintf(stderr, "Done!\n");
    #endif
    fclose(OUTLOC);
    }
  #endif

  fprintf(stderr, "  [+] Freeing compression models ... ");
  for(n = 0 ; n < P->nModels ; ++n)
    #ifdef KMODELSUSAGE
    FreeKModel(KModels[n]);
    #else
    FreeCModel(Models[n]);
    #endif
  #ifdef KMODELSUSAGE
  Free(KModels);
  #else
  Free(Models);
  #endif
  fprintf(stderr, "Done!\n");

  StopTimeNDRM(Time, clock());
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ RESULTS ]=======================\n");
  if(topSize <= 100){
    #ifdef LOCAL_SIMILARITY
    if(P->local == 1)
      PrintTopInfoWP(P->top, topSize, P->dbFiles);
    else
      PrintTopInfo(P->top, topSize, P->dbFiles);
    #else
    PrintTopInfo(P->top, topSize);
    #endif
    }
  else{
    fprintf(stderr, "Top results have been sent to file.\n");
    }
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ STATISTICS ]====================\n");
  StopCalcAll(Time, clock());
  fprintf(stderr, "\n");

  if (xargv) {
    if (xargv[0]) {
      Free(xargv[0]);  // Free the duplicated string
    }
    Free(xargv);         // Free the array itself
  }

  // Free P resources
  if (P->output) {
    Free(P->output);
  }

  // Free file names
  if (P->files) {
    Free(P->files);
  }

  // Free DB file names
  if (P->dbFiles) {
    Free(P->dbFiles);
  }

  // Free model file name
  if (P->modelFile) {
    Free(P->modelFile);
  }

  RemoveClock(Time);
  DeleteTop(P->top);
  for(ref = 0 ; ref < P->nThreads ; ++ref){
    DeleteTop(T[ref].top);
    Free(T[ref].model);
    }
  Free(T);
  Free(P);

  return EXIT_SUCCESS;
  }

int32_t P_Filter(char **argv, int argc){
  char **p = *&argv;
  FILE *OUTPUT = NULL, *INPUT = NULL;

  PEYE = (EYEPARAM *) Malloc(1 * sizeof(EYEPARAM));
  if((PEYE->help = ArgsState(DEFAULT_HELP, p, argc, "-h", "--help")) == 1 || argc < 2){
    PrintMenuFilter();
    Free(P);
    return EXIT_SUCCESS;
  }
  if(ArgsState(DEF_VERSION, p, argc, "-V", "--version")){
    PrintVersion();
    Free(P);
    return EXIT_SUCCESS;
  }

  PEYE->verbose    = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v", "--verbose");
  PEYE->force      = ArgsState  (DEFAULT_FORCE,   p, argc, "-F", "--force");
  PEYE->windowSize = ArgsNum    (100,             p, argc, "-s", 1, 999999);
  PEYE->windowType = ArgsNum    (1,               p, argc, "-w", 0, 3);
  PEYE->sampling   = ArgsNum    (10,              p, argc, "-x", 1, 999999);
  PEYE->threshold  = ArgsDouble (1.0,             p, argc, "-t");
  PEYE->lowerSimi  = ArgsDouble (0.00,            p, argc, "-sl");
  PEYE->upperSimi  = ArgsDouble (100.00,          p, argc, "-su");
  PEYE->lowerSize  = ArgsNum64  (1,               p, argc, "-dl", 1,
  9999999999);
  PEYE->upperSize  = ArgsNum64  (9999999999,      p, argc, "-du", 1,
  9999999999);
  PEYE->output     = ArgsFileGen(p, argc, "-o", "coords", ".fil");

  if(!PEYE->force)
    FAccessWPerm(PEYE->output);
  OUTPUT = Fopen(PEYE->output, "w");

  fprintf(stderr, "\n");
  if(PEYE->verbose){
    PrintArgsFilter(PEYE);
    }

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  TIME *Time = CreateClock(clock());

  int sym;
  char fname[MAX_NAME];
  double fvalue;
  uint64_t fsize;

  FILTER *FIL = CreateFilter(PEYE->windowSize, PEYE->sampling, PEYE->windowType,
  PEYE->threshold);
  INPUT = Fopen(argv[argc-1], "r");
  while((sym = fgetc(INPUT)) != EOF){
    if(sym == '#'){
      if(fscanf(INPUT, "\t%lf\t%"PRIu64"\t%s\n", &fvalue, &fsize, fname) != 3){
        fprintf(stderr, "  [x] Error: unknown type of file!\n");
        exit(1);
        }

      if(fsize > PEYE->upperSize ||  fsize < PEYE->lowerSize ||
        fvalue > PEYE->upperSimi || fvalue < PEYE->lowerSimi)
        continue;

      fprintf(stderr, "  [+] Filtering & segmenting %s ... ", fname);
      fprintf(OUTPUT, "$\t%lf\t%"PRIu64"\t%s\n", fvalue, fsize, fname);
      InitEntries(FIL, fsize, INPUT);
      FilterStream(FIL, OUTPUT);
      DeleteEntries(FIL);
      fprintf(stderr, "Done!\n");
      }
    }
  DeleteFilter(FIL);

  fclose(OUTPUT);
  fclose(INPUT);

  StopTimeNDRM(Time, clock());
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ STATISTICS ]====================\n");
  StopCalcAll(Time, clock());
  fprintf(stderr, "\n");

  RemoveClock(Time);

  return EXIT_SUCCESS;
  }

int32_t P_Filter_Visual(char **argv, int argc){
  char **p = *&argv, fname[MAX_NAME];
  FILE *OUTPUT = NULL, *INPUT = NULL;
  int sym;
  double fvalue;
  uint32_t n, tmp, maxName, extraLength, cmp;
  uint64_t maxSize = 0, fsize, iPos, ePos, nSeq, filtered;
  Painter *Paint;
  COLORS *CLR;

  PEYE = (EYEPARAM *) Malloc(1 * sizeof(EYEPARAM));
  if((PEYE->help = ArgsState(DEFAULT_HELP, p, argc, "-h", "--help")) == 1 || argc < 2){
    PrintMenuVisual();
    Free(P);
    return EXIT_SUCCESS;
  }

  if(ArgsState(DEF_VERSION, p, argc, "-V", "--version")){
    PrintVersion();
    Free(P);
    return EXIT_SUCCESS;
  }

  PEYE->verbose    = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v", "--verbose");
  PEYE->force      = ArgsState  (DEFAULT_FORCE,   p, argc, "-F", "--force");
  PEYE->width      = ArgsDouble (DEFAULT_WIDTH,   p, argc, "-w");
  PEYE->space      = ArgsDouble (DEFAULT_SPACE,   p, argc, "-s");
  PEYE->showScale  = ArgsState  (DEFAULT_SHOWS,   p, argc, "-ss", "-showScale");
  PEYE->showNames  = ArgsState  (DEFAULT_NAMES,   p, argc, "-sn", "-showNames");
  PEYE->sameScale  = ArgsState  (DEFAULT_RSCAL,   p, argc, "-rs", "-sameScale");
  PEYE->best       = ArgsState  (DEFAULT_GBEST,   p, argc, "-bg", "-best");
  PEYE->start      = ArgsDouble (0.35,            p, argc, "-i");
  PEYE->rotations  = ArgsDouble (1.50,            p, argc, "-r");
  PEYE->hue        = ArgsDouble (1.92,            p, argc, "-u");
  PEYE->gamma      = ArgsDouble (0.50,            p, argc, "-g");
  PEYE->proportion = ArgsDouble (500,             p, argc, "-p");
  PEYE->lowerSimi  = ArgsDouble (0.00,            p, argc, "-sl");
  PEYE->upperSimi  = ArgsDouble (100.00,          p, argc, "-su");
  PEYE->lowerSize  = ArgsNum64  (1,               p, argc, "-dl", 1,
  9999999999);
  PEYE->upperSize  = ArgsNum64  (9999999999,      p, argc, "-du", 1,
  9999999999);
  PEYE->enlarge    = ArgsNum64  (0,               p, argc, "-e",  0,
  9999999999);
  PEYE->output     = ArgsFileGen(p, argc, "-o", "femap", ".svg");

  if(!PEYE->force)
    FAccessWPerm(PEYE->output);
  OUTPUT = Fopen(PEYE->output, "w");

  fprintf(stderr, "\n");
  if(PEYE->verbose){
    PrintArgsEye(PEYE);
    }

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  TIME *Time = CreateClock(clock());

  // TODO: OPTION TO IGNORE SCALE AND SET THEM AT SAME SIZE

  CLR = (COLORS *) Calloc(1, sizeof(COLORS));
  CLR->start     = PEYE->start;
  CLR->rotations = PEYE->rotations;
  CLR->hue       = PEYE->hue;
  CLR->gamma     = PEYE->gamma;

  // BUILD UNIQUE ARRAY FOR NAMES USING REGULAR EXPRESSION:
  //
  SLABELS *SL = CreateSLabels();
  uint64_t unique = 0;
  // tested at: https://regex101.com/
  char *regexString = ".*\\|.*\\|.*\\|_([a-z A-Z]*_[a-z A-Z]*)";
  regex_t regexCompiled;
  regmatch_t groupArray[2];
  if(PEYE->best == 1){
    if(regcomp(&regexCompiled, regexString, REG_EXTENDED)){
      fprintf(stderr, "  [x] Error: regular expression compilation!\n");
      return 1;
      }
    }

  INPUT = Fopen(argv[argc-1], "r");
  nSeq = 0;
  maxName = 0;
  filtered = 0;
  while((sym = fgetc(INPUT)) != EOF){

    if(sym == '$' || sym == '#'){
      if(fscanf(INPUT, "\t%lf\t%"PRIu64"\t%s\n", &fvalue, &fsize, fname) != 3){
        fprintf(stderr, "  [x] Error: unknown type of file!\n");
        exit(1);
        }

      if(PEYE->best == 1){
        if(regexec(&regexCompiled, fname, 2, groupArray, 0) == 0){
          char sourceCopy[strlen(fname) + 1];
          strcpy(sourceCopy, fname);
          sourceCopy[groupArray[1].rm_eo] = 0;
          if(SearchSLabels(SL, sourceCopy + groupArray[1].rm_so) == 0){
            ++unique;
            AddSLabel(SL, sourceCopy + groupArray[1].rm_so);
            UpdateSLabels(SL);
            }
          else{
            ++filtered;
            continue;
            }
          }
        }

      if(fsize > PEYE->upperSize ||  fsize < PEYE->lowerSize ||
        fvalue > PEYE->upperSimi || fvalue < PEYE->lowerSimi){
        ++filtered;
        continue;
        }

      if(fsize > maxSize)
        maxSize = fsize;

      if(PEYE->best == 1){
        if((tmp = strlen(SL->names[SL->idx-1])) > maxName)
          maxName = tmp;
        }
      else{
        if((tmp = strlen(fname)) > maxName)
          maxName = tmp;
        }
      ++nSeq;
      }
    }
  rewind(INPUT);

  fprintf(stderr, "Skipping %"PRIu64" from %"PRIu64" entries.\n",
  filtered, filtered+nSeq);


  if(PEYE->best == 1){
    fprintf(stderr, "Number of unique existing species: %"PRIu64".\n", unique);
    fprintf(stderr, "Unique species:\n");
    for(n = 0 ; n < SL->idx ; ++n)
      fprintf(stderr, "  [+] %s\n", SL->names[n]);
    DeleteSLabels(SL);
    SL = CreateSLabels();
    }
  Paint = CreatePainter(maxSize, PEYE->width, PEYE->space, PEYE->proportion,
  "#ffffff");

  extraLength = 0;
  if(PEYE->showNames == 1){
    extraLength = 10 * maxName; // LETTER SIZE * MAXNAME
    Paint->cy += extraLength;
    }

  if(PEYE->showScale == 1)
    nSeq += 2;

  PrintHead(OUTPUT, (2 * DEFAULT_CX) + (((Paint->width + PEYE->space) * nSeq) -
  PEYE->space), Paint->size + EXTRA + Paint->width + extraLength);
  Rect(OUTPUT, (2 * DEFAULT_CX) + (((Paint->width + PEYE->space) * nSeq) -
  PEYE->space), Paint->size + EXTRA + Paint->width + extraLength, 0, 0,
  "#ffffff");

  if(PEYE->showScale == 1){
    nSeq -= 2;
    Paint->cx += (2 * Paint->width + PEYE->space);
    // PRINT HEATMAP SCALE
    uint32_t size = 4 * Paint->width;
    for(n = 0 ; n < size ; ++n){
      char color[12];
      Rect(OUTPUT, Paint->width, 1, Paint->cx - (Paint->width*2),
      Paint->cy + n, HeatMapColor(((double) n / size), color, CLR));
      }
    Text(OUTPUT, Paint->cx-(Paint->width*2 + 14), Paint->cy+18,     "+");
    Text(OUTPUT, Paint->cx-(Paint->width*2 + 12), Paint->cy+size-6, "-");

    // HIGH COMPLEX SCALE
    Rect(OUTPUT, Paint->width, Paint->width, Paint->cx - (Paint->width*2),
    Paint->cy + Paint->width + size, GetRgbColor(HIGH_COMPLEX));
    Text(OUTPUT, Paint->cx-(Paint->width*2 + 14), (Paint->cy+Paint->width+size)
    +(Paint->width/2)+5, "+");
    // MEDIUMH COMPLEX SCALE
    Rect(OUTPUT, Paint->width, Paint->width, Paint->cx - (Paint->width*2),
    Paint->cy + 2*Paint->width + size, GetRgbColor(MEDIUMH_COMPLEX));
    // MEDIUML COMPLEX SCALE
    Rect(OUTPUT, Paint->width, Paint->width, Paint->cx - (Paint->width*2),
    Paint->cy + 3*Paint->width + size, GetRgbColor(MEDIUML_COMPLEX));
    // LOW COMPLEX SCALE
    Rect(OUTPUT, Paint->width, Paint->width, Paint->cx - (Paint->width*2),
    Paint->cy + 4*Paint->width + size, GetRgbColor(LOW_COMPLEX));
    Text(OUTPUT, Paint->cx-(Paint->width*2+12), (Paint->cy+4*Paint->width+size)
    +(Paint->width/2)+4, "-");
    }

  if(nSeq > 0) fprintf(stderr, "Addressing regions individually:\n");
  while((sym = fgetc(INPUT)) != EOF){

    if(sym == '$' || sym == '#'){
      if(fscanf(INPUT, "\t%lf\t%"PRIu64"\t%s\n", &fvalue, &fsize, fname) != 3){
        fprintf(stderr, "  [x] Error: unknown type of file!\n");
        exit(1);
        }

      if(sym == '#'){
        if(PEYE->best == 1){
          if(regexec(&regexCompiled, fname, 2, groupArray, 0) == 0){
            char sourceCopy[strlen(fname) + 1];
            strcpy(sourceCopy, fname);
            sourceCopy[groupArray[1].rm_eo] = 0;
            if(SearchSLabels(SL, sourceCopy + groupArray[1].rm_so) == 0){
              AddSLabel(SL, sourceCopy + groupArray[1].rm_so);
              UpdateSLabels(SL);
              }
            else{ // SKIP
              continue;
              }
            }
          }

        // SKIPS: FILTERED ATTRIBUTES
        if(fsize > PEYE->upperSize ||  fsize < PEYE->lowerSize ||
          fvalue > PEYE->upperSimi || fvalue < PEYE->lowerSimi){
          continue;
          }

        if(PEYE->showNames == 1){  // PRINT NAMES 90D
          if(PEYE->best == 1){
            Text90d(OUTPUT, -(Paint->cy-32), Paint->cx+Paint->width-
            (Paint->width/2.0)+10, SL->names[SL->idx-1]);
            }
          else{
            Text90d(OUTPUT, -(Paint->cy-32), Paint->cx+Paint->width-
            (Paint->width/2.0)+10, fname);
            }
          }

        if(PEYE->best == 1)
          fprintf(stderr, "  [+] Painting %s (%s) ... ", SL->names[SL->idx-1],
          fname);
        else
          fprintf(stderr, "  [+] Painting %s ... ", fname);

        char tmpTxt[MAX_NAME], color[12];
        if(fvalue < 10){
          sprintf(tmpTxt, "%.1lf", fvalue);
          Text(OUTPUT, (Paint->cx+Paint->width/2)-12, Paint->cy-10, tmpTxt);
          }
        else{
          sprintf(tmpTxt, "%u", (unsigned) fvalue);
          Text(OUTPUT, (Paint->cx+Paint->width/2)-9, Paint->cy-10, tmpTxt);
          }
        RectWithBorder(OUTPUT, Paint->width, Paint->width, Paint->cx, Paint->cy,
        HeatMapColor(BoundDouble(0.0, 1-fvalue/100.0, 1.0), color, CLR));

        if(nSeq > 0)
          Paint->cx += Paint->width + Paint->space;
        fprintf(stderr, "Done!\n");
        continue; // CONTINUE WITHOUT WRITTING LOCAL COMPLEXITY
        }

      if(PEYE->best == 1){
        if(regexec(&regexCompiled, fname, 2, groupArray, 0) == 0){
          char sourceCopy[strlen(fname) + 1];
          strcpy(sourceCopy, fname);
          sourceCopy[groupArray[1].rm_eo] = 0;
          if(SearchSLabels(SL, sourceCopy + groupArray[1].rm_so) == 0){
            AddSLabel(SL, sourceCopy + groupArray[1].rm_so);
            UpdateSLabels(SL);
            }
          else{ // SKIP
            while(fscanf(INPUT, "%"PRIu64":%"PRIu64"\t%u\n", &iPos, &ePos, &cmp)
            == 3)
              ; // DO NOTHING
            continue;
            }
          }
        }

      // SKIPS: FILTERED ATTRIBUTES
      if(fsize > PEYE->upperSize ||  fsize < PEYE->lowerSize ||
        fvalue > PEYE->upperSimi || fvalue < PEYE->lowerSimi){
        while(fscanf(INPUT, "%"PRIu64":%"PRIu64"\t%u\n", &iPos, &ePos, &cmp)
        == 3)
          ; // DO NOTHING
        continue;
        }

      if(PEYE->showNames == 1){  // PRINT NAMES 90D
        if(PEYE->best == 1){
          Text90d(OUTPUT, -(Paint->cy-32), Paint->cx+Paint->width-
          (Paint->width/2.0)+10, SL->names[SL->idx-1]);
          }
        else{
          Text90d(OUTPUT, -(Paint->cy-32), Paint->cx+Paint->width-
          (Paint->width/2.0)+10, fname);
          }
        }

      char tmpTxt[MAX_NAME], color[12];
      if(fvalue < 10){
        sprintf(tmpTxt, "%.1lf", fvalue);
        Text(OUTPUT, (Paint->cx+Paint->width/2)-12, Paint->cy-10, tmpTxt);
        }
      else{
        sprintf(tmpTxt, "%u", (unsigned) fvalue);
        Text(OUTPUT, (Paint->cx+Paint->width/2)-9, Paint->cy-10, tmpTxt);
        }
      RectWithBorder(OUTPUT, Paint->width, Paint->width, Paint->cx, Paint->cy,
      HeatMapColor(BoundDouble(0.0, 1-fvalue/100.0, 1.0), color, CLR));

      Paint->cy += Paint->width + Paint->space;

      if(PEYE->best == 1)
        fprintf(stderr, "  [+] Painting %s (%s) ... ", SL->names[SL->idx-1],
        fname);
      else
        fprintf(stderr, "  [+] Painting %s ... ", fname);
      while(1){
        off_t beg = Ftello(INPUT);
        if(fscanf(INPUT, "%"PRIu64":%"PRIu64"\t%u\n", &iPos, &ePos, &cmp) != 3){
          Fseeko(INPUT, (off_t) beg, SEEK_SET);
          Chromosome(OUTPUT, Paint->width, GetPoint(Paint, fsize), Paint->cx,
          Paint->cy);
          if(nSeq > 0)
            Paint->cx += Paint->width + Paint->space;
          break;
          }
        else{
          int color = 0;
          switch(cmp){
            case 0: color = LOW_COMPLEX;     break;
            case 1: color = MEDIUML_COMPLEX; break;
            case 2: color = MEDIUMH_COMPLEX; break;
            case 3: color = HIGH_COMPLEX;    break;
            default: color = 0;
            }
          if(PEYE->enlarge == 0){
            Rect(OUTPUT, Paint->width, GetPoint(Paint, ePos-iPos+1), Paint->cx,
            Paint->cy + GetPoint(Paint, iPos), GetRgbColor(color));
            }
          else{
            Rect(OUTPUT, Paint->width, GetPoint(Paint, ePos-iPos+1+PEYE->enlarge),
            Paint->cx, Paint->cy + GetPoint(Paint, iPos), GetRgbColor(color));
            }
          }
        }

      Paint->cy -= Paint->width + Paint->space;
      fprintf(stderr, "Done!\n");
      }
    }

  PrintFinal(OUTPUT);

  fclose(OUTPUT);
  fclose(INPUT);
  Free(CLR);
  if(PEYE->best == 1)
    DeleteSLabels(SL);
  StopTimeNDRM(Time, clock());
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ STATISTICS ]====================\n");
  StopCalcAll(Time, clock());
  fprintf(stderr, "\n");

  RemoveClock(Time);
  return EXIT_SUCCESS;
  }

int32_t P_Inter(char **argv, int argc){
  char        **p = *&argv, **xargv, *xpl = NULL;
  int32_t     xargc = 0;
  uint32_t    n, k, col, ref;
  double      gamma;
  Threads     *T;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEFAULT_HELP, p, argc, "-h", "--help")) == 1 || argc < 2){
    PrintMenuInter();
    Free(P);
    return EXIT_SUCCESS;
  }

  if(ArgsState(DEF_VERSION, p, argc, "-V", "--version")){
    PrintVersion();
    Free(P);
    return EXIT_SUCCESS;
  }

  if(ArgsState(0, p, argc, "-s", "--show")){
    PrintLevels();
    Free(P);
    return EXIT_SUCCESS;
  }


  P->verbose  = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v", "--verbose");
  P->force    = ArgsState  (DEFAULT_FORCE,   p, argc, "-F", "--force");
  P->level    = ArgsNum    (0, p, argc, "-l", MIN_LEV, MAX_LEV);
  P->nThreads = ArgsNum    (DEFAULT_THREADS, p, argc, "-n", MIN_THREADS,
  MAX_THREADS);

  P->nModels = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-m") == 0)
      P->nModels += 1;

  if(P->nModels == 0 && P->level == 0)
    P->level = DEFAULT_LEVEL;

  if(P->level != 0){
    xpl = GetLevels(P->level);
    xargc = StrToArgv(xpl, &xargv);
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-m") == 0)
        P->nModels += 1;
    }

  gamma = DEFAULT_GAMMA;
  for(n = 1 ; n < xargc ; ++n)
    if(strcmp(xargv[n], "-g") == 0)
      gamma = atof(xargv[n+1]);

  col = MAX_COLLISIONS;
  for(n = 1 ; n < xargc ; ++n)
    if(strcmp(xargv[n], "-c") == 0)
      col = atoi(xargv[n+1]);

  P->col       = ArgsNum    (col,   p, argc, "-c", 1, 200);
  P->gamma     = ArgsDouble (gamma, p, argc, "-g");
  P->gamma     = ((int)(P->gamma * 65536)) / 65536.0;
  P->nFiles    = ReadFNames (P, argv[argc-1], 1);
  P->output    = ArgsFileGen(p, argc, "-x", "matrix", ".csv");
  P->labels    = ArgsFileGen(p, argc, "-o", "labels", ".csv");
  FILE *OUTPUT = Fopen(P->output, "w");
  FILE *LABELS = Fopen(P->labels, "w");

  if(P->nModels == 0){
    fprintf(stderr, "Error: at least you need to use a context model!\n");
    return EXIT_FAILURE;
    }

  // READ MODEL PARAMETERS FROM XARGS & ARGS
  T = (Threads *) Calloc(P->nFiles, sizeof(Threads));
  for(ref = 0 ; ref < P->nFiles ; ++ref){
    T[ref].model = (ModelPar *) Calloc(P->nModels, sizeof(ModelPar));
    T[ref].id = ref;
    k = 0;
    for(n = 1 ; n < argc ; ++n)
      if(strcmp(argv[n], "-m") == 0)
        T[ref].model[k++] = ArgsUniqModel(argv[n+1], 0);
    if(P->level != 0){
      for(n = 1 ; n < xargc ; ++n)
        if(strcmp(xargv[n], "-m") == 0)
          T[ref].model[k++] = ArgsUniqModel(xargv[n+1], 0);
      }
    }

  fprintf(stderr, "\n");
  if(P->verbose) PrintArgsInter(P, T[0]);

  P->matrix = (double  **) Calloc(P->nFiles, sizeof(double *));
  P->size   = (uint64_t *) Calloc(P->nFiles, sizeof(uint64_t));
  for(n = 0 ; n < P->nFiles ; ++n){
    P->size[n]   = FopenBytesInFile(P->files[n]);
    P->matrix[n] = (double *) Calloc(P->nFiles, sizeof(double));
    }

  if(P->nThreads > P->nFiles){
    fprintf(stderr, "Error: the number of threads must not be higher than the "
    "number of files\n");
    exit(1);
    }

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  TIME *Time = CreateClock(clock());
  for(n = 0 ; n < P->nFiles ; ++n){
    CompressActionInter(T, n);
    for(k = 0 ; k < P->nFiles ; ++k){
      fprintf(stderr, "%.4lf\t", P->matrix[n][k]);
      }
    fprintf(stderr, "\n");
    }
  StopTimeNDRM(Time, clock());
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ RESULTS ]=======================\n");
  fprintf(stderr, "Normalized Dissimilarity Rate matrix:\n");
  for(n = 0 ; n < P->nFiles ; ++n){
    for(k = 0 ; k < P->nFiles ; ++k){
      fprintf(stderr, "%.4lf\t", P->matrix[n][k]);
      fprintf(OUTPUT, "%.4lf\t", P->matrix[n][k]);
      }
    fprintf(stderr, "\n");
    fprintf(OUTPUT, "\n");
    fprintf(LABELS, "%s\t", P->files[n]);
    }
  fprintf(LABELS, "\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ STATISTICS ]====================\n");
  StopCalcAll(Time, clock());
  fprintf(stderr, "\n");

  // STOP &
  RemoveClock(Time);
  for(ref = 0 ; ref < P->nFiles ; ++ref)
    Free(T[ref].model);
  Free(T);
  if(!OUTPUT) fclose(OUTPUT);
  if(!LABELS) fclose(LABELS);

  return EXIT_SUCCESS;
  }

int32_t P_Inter_Visual(char **argv, int argc){
  char        **p = *&argv;
  uint32_t    n;
  double      start, rotations, hue, gamma, width, space;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEFAULT_HELP, p, argc, "-h", "--help")) == 1 || argc < 2){
    PrintMenuInterVisual();
    Free(P);
    return EXIT_SUCCESS;
  }

  if(ArgsState(DEF_VERSION, p, argc, "-V", "--version")){
    PrintVersion();
    Free(P);
    return EXIT_SUCCESS;
  }

  P->verbose  = ArgsState    (DEFAULT_VERBOSE, p, argc, "-v", "--verbose");
  P->force    = ArgsState    (DEFAULT_FORCE,   p, argc, "-F", "--force");
  P->image    = (uint8_t *)   ArgsFilesImg    (p, argc, "-x");
  P->labels   = (uint8_t *)   ArgsFile        (p, argc, "-l");
  width       = ArgsDouble   (DEFAULT_WIDTH,   p, argc, "-w");
  space       = ArgsDouble   (DEFAULT_SPACE,   p, argc, "-a");
  start       = ArgsDouble   (0.35,            p, argc, "-s");
  rotations   = ArgsDouble   (1.50,            p, argc, "-r");
  hue         = ArgsDouble   (1.92,            p, argc, "-u");
  gamma       = ArgsDouble   (0.50,            p, argc, "-g");

  fprintf(stderr, "\n");
  if(P->verbose){
    fprintf(stderr, "==[ CONFIGURATION ]=================\n");
    fprintf(stderr, "Verbose mode ....................... yes\n");
    fprintf(stderr, "Heatmap design characteristics:\n");
    fprintf(stderr, "  [+] Width ........................ %.3g\n", width);
    fprintf(stderr, "  [+] Space ........................ %.3g\n", space);
    fprintf(stderr, "  [+] Start ........................ %.3g\n", start);
    fprintf(stderr, "  [+] Rotations .................... %.3g\n", rotations);
    fprintf(stderr, "  [+] Hue .......................... %.3g\n", hue);
    fprintf(stderr, "  [+] Gamma ........................ %.3g\n", gamma);
    fprintf(stderr, "Input labels filename .............. %s\n",   P->labels);
    fprintf(stderr, "Output heatmap filename ............ %s\n",   P->image);
    fprintf(stderr, "\n");
    }

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  ReadMatrixSize(argv[argc-1]);
  P->matrix = (double  **) Calloc(P->nFiles, sizeof(double *));
  P->files  = (char    **) Calloc(P->nFiles, sizeof(char   *));
  for(n = 0 ; n < P->nFiles ; ++n){
    P->matrix[n] = (double *) Calloc(P->nFiles, sizeof(double));
    P->files [n] = (char   *) Calloc(MAX_LABEL, sizeof(char  ));
    P->files[n][0]  = '?';
    P->files[n][1]  = '\0';
    }
  fprintf(stderr, "  [+] Loading labels ... ");
  ReadLabels();
  fprintf(stderr, "Done!\n");
  fprintf(stderr, "  [+] Loading matrix ... ");
  ReadMatrix(argv[argc-1]);
  fprintf(stderr, "Done!\n");
  fprintf(stderr, "  [+] Painting heatmap ... ");
  PaintMatrix(start, rotations, hue, gamma, width, space);
  fprintf(stderr, "Done!\n\n");

  return EXIT_SUCCESS;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int argc, char *argv[])
{
  char **p = *&argv;

  if(ArgsState(0, p, argc, "-V", "--version"))
  {
    PrintVersion();
    return EXIT_SUCCESS;
  }

  if(ArgsState(DEFAULT_HELP, p, argc, "-h", "--help") && argc == 2 || argc < 2)
  {
    PrintMainMenu();
    return EXIT_SUCCESS;
  }

  switch(KeyString(argv[1]))
  {
    case K1: PrintMainMenu();                                   break;
    case K2: P_Falcon                        (argv+1, argc-1);  break;
    case K3: P_Filter                        (argv+1, argc-1);  break;
    case K4: P_Filter_Visual                 (argv+1, argc-1);  break;
    case K5: P_Inter                         (argv+1, argc-1);  break;
    case K6: P_Inter_Visual                  (argv+1, argc-1);  break;

    default:
      PrintWarning("unknown menu option!");
    PrintMainMenu();
  }

  return EXIT_SUCCESS;
}