
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "serialization.h"
#include "mem.h"
#include "common.h"

// Helper function to serialize a hashtable to file
static int SerializeHashTable(FILE *F, HashTable *HT) {
  // Write index array in one go
  if(fwrite(HT->index, sizeof(ENTMAX), HASH_SIZE, F) != HASH_SIZE)
    return -1;

  // Write entries more efficiently - single write per hash index
  for(uint32_t i = 0; i < HASH_SIZE; i++) {
    if(fwrite(HT->entries[i], sizeof(Entry), HT->maxC, F) != HT->maxC)
      return -2;
  }

  return 0;
}

// Helper function to deserialize a hashtable from file
static int DeserializeHashTable(FILE *F, HashTable *HT, uint32_t col) {
  // Initialize hash table
  HT->maxC = col;
  HT->index = (ENTMAX *) Calloc(HASH_SIZE, sizeof(ENTMAX));
  HT->entries = (Entry **) Calloc(HASH_SIZE, sizeof(Entry *));
  if(!HT->index || !HT->entries)
    return -1;

  // Read index array
  if(fread(HT->index, sizeof(ENTMAX), HASH_SIZE, F) != HASH_SIZE)
    return -2;

  // Allocate and read each hash entry block
  for(uint32_t i = 0; i < HASH_SIZE; i++) {
    HT->entries[i] = (Entry *) Calloc(HT->maxC, sizeof(Entry));
    if(!HT->entries[i])
      return -3;

    // Read entire entry array at once
    if(fread(HT->entries[i], sizeof(Entry), HT->maxC, F) != HT->maxC)
      return -4;
  }

  return 0;
}

// Helper function to serialize an array to file
static int SerializeArray(FILE *F, Array *AR, uint64_t nPModels) {
  uint64_t size = nPModels << 2; // * 4 for ACGT
  return fwrite(AR->counters, sizeof(ACC), size, F) != size ? -1 : 0;
}

// Helper function to deserialize an array from file
static int DeserializeArray(FILE *F, Array *AR, uint64_t nPModels) {
  uint64_t size = nPModels << 2; // * 4 for ACGT
  AR->counters = (ACC *) Calloc(size, sizeof(ACC));
  if(!AR->counters)
    return -1;
    
  return fread(AR->counters, sizeof(ACC), size, F) != size ? -2 : 0;
}

int SaveModels(const char *filename, CModel **Models, uint32_t nModels, uint32_t col) {
  // Validate input parameters
  if(!filename || !Models || nModels == 0) {
    fprintf(stderr, "Error: Invalid parameters for SaveModels\n");
    return -1;
  }

  // Use common file handler with error checking
  FILE *F = Fopen(filename, "wb");

  // Write file header
  ModelHeader header;
  memset(&header, 0, sizeof(ModelHeader)); // Ensure clean initialization
  header.magic = MODEL_MAGIC_NUMBER;
  header.version = MODEL_VERSION;
  header.nModels = nModels;
  header.alphabetSize = ALPHABET_SIZE;
  header.timestamp = (uint64_t)time(NULL);
  header.hashSize = HASH_SIZE;
  header.maxCollisions = col;

  if(fwrite(&header, sizeof(ModelHeader), 1, F) != 1) {
    fprintf(stderr, "Error writing model file header\n");
    Fclose(F);
    return -3;
  }

  // For each model
  for(uint32_t n = 0; n < nModels; n++) {
    CModel *M = Models[n];
    if(!M) {
      fprintf(stderr, "Error: NULL model at index %u\n", n);
      Fclose(F);
      return -4;
    }

    // Write model entry header
    ModelMeta entryHeader;
    memset(&entryHeader, 0, sizeof(ModelMeta)); // Ensure clean initialization
    entryHeader.ctx = M->ctx;
    entryHeader.alphaDen = M->alphaDen;
    entryHeader.ir = M->ir;
    entryHeader.edits = M->edits;
    entryHeader.eDen = M->edits != 0 ? M->SUBS.eDen : 0;
    entryHeader.mode = M->mode;
    entryHeader.nPModels = M->nPModels;
    entryHeader.maxCount = M->maxCount;
    entryHeader.multiplier = M->multiplier;

    if(fwrite(&entryHeader, sizeof(ModelMeta), 1, F) != 1) {
      fprintf(stderr, "Error writing model entry header for model %u\n", n);
      Fclose(F);
      return -5;
    }

    // Save model data
    int result = 0;
    switch(M->mode) {
      case HASH_TABLE_MODE:
        result = SerializeHashTable(F, &M->hTable);
        break;
      case ARRAY_MODE:
        result = SerializeArray(F, &M->array, M->nPModels);
        break;
      default:
        fprintf(stderr, "Unknown model mode: %u\n", M->mode);
        Fclose(F);
        return -6;
    }

    if(result != 0) {
      fprintf(stderr, "Error serializing model %u data: %d\n", n, result);
      Fclose(F);
      return -7;
    }
  }

  // Ensure data is committed to disk
  fflush(F);
  int fd = fileno(F);
  if(fd >= 0) {
    fsync(fd);
  }
  
  Fclose(F);
  return 0;
}

int LoadModels(const char *filename, CModel ***ModelsPtr, uint32_t *nModels, uint32_t *col) {
  // Validate input parameters
  if(!filename || !ModelsPtr || !nModels || !col) {
    fprintf(stderr, "Error: Invalid parameters for LoadModels\n");
    return -1;
  }
  
  *ModelsPtr = NULL;
  *nModels = 0;
  *col = 0;

  // Use common file handler with error checking
  FILE *F = Fopen(filename, "rb");

  // Read file header
  ModelHeader header;
  if(fread(&header, sizeof(ModelHeader), 1, F) != 1) {
    fprintf(stderr, "Error reading model file header\n");
    Fclose(F);
    return -3;
  }

  // Validate header
  if(header.magic != MODEL_MAGIC_NUMBER) {
    fprintf(stderr, "Error: Invalid model file format (wrong magic number)\n");
    Fclose(F);
    return -4;
  }

  if(header.version != MODEL_VERSION) {
    fprintf(stderr, "Error: Unsupported model file version: %u\n", header.version);
    Fclose(F);
    return -5;
  }

  if(header.alphabetSize != ALPHABET_SIZE) {
    fprintf(stderr, "Error: Model file has different alphabet size: %u (expected %u)\n",
            header.alphabetSize, ALPHABET_SIZE);
    Fclose(F);
    return -6;
  }

  // Allocate models array using common memory functions
  CModel **Models = (CModel **) Malloc(header.nModels * sizeof(CModel *));
  
  // Initialize to NULL for safe cleanup
  for(uint32_t i = 0; i < header.nModels; i++) {
    Models[i] = NULL;
  }

  // For each model
  for(uint32_t n = 0; n < header.nModels; n++) {
    // Read model entry header
    ModelMeta entryHeader;
    if(fread(&entryHeader, sizeof(ModelMeta), 1, F) != 1) {
      fprintf(stderr, "Error reading model entry header for model %u\n", n);
      // Free already loaded models
      FreeLoadedModels(Models, n);
      Fclose(F);
      return -8;
    }

    // Create model structure using common memory allocation
    Models[n] = (CModel *) Calloc(1, sizeof(CModel));
    CModel *M = Models[n];

    // Fill in basic model parameters
    M->ctx = entryHeader.ctx;
    M->alphaDen = entryHeader.alphaDen;
    M->ir = entryHeader.ir;
    M->edits = entryHeader.edits;
    M->mode = entryHeader.mode;
    M->nPModels = entryHeader.nPModels;
    M->maxCount = entryHeader.maxCount;
    M->multiplier = entryHeader.multiplier;
    M->pModelIdx = 0;
    M->pModelIdxIR = M->nPModels - 1;
    M->ref = 1; // This is a reference model

    // Initialize edits structure if needed
    if(M->edits != 0) {
      M->SUBS.seq = CreateCBuffer(BUFFER_SIZE, BGUARD);
      M->SUBS.in = 0;
      M->SUBS.idx = 0;
      M->SUBS.mask = (uint8_t *) Calloc(BGUARD, sizeof(uint8_t));
      M->SUBS.threshold = M->edits;
      M->SUBS.eDen = entryHeader.eDen;
    }

    // Load model data
    int result = 0;
    switch(M->mode) {
      case HASH_TABLE_MODE:
        result = DeserializeHashTable(F, &M->hTable, header.maxCollisions);
        break;
      case ARRAY_MODE:
        result = DeserializeArray(F, &M->array, M->nPModels);
        break;
      default:
        fprintf(stderr, "Unknown model mode: %u\n", M->mode);
        FreeLoadedModels(Models, n+1);
        Fclose(F);
        return -12;
    }

    if(result != 0) {
      fprintf(stderr, "Error deserializing model %u data: %d\n", n, result);
      FreeLoadedModels(Models, n+1);
      Fclose(F);
      return -13;
    }
  }

  *ModelsPtr = Models;
  *nModels = header.nModels;
  *col = header.maxCollisions;
  
  Fclose(F);
  return 0;
}

void FreeLoadedModels(CModel **Models, uint32_t nModels) {
  if(!Models) return;
  
  for(uint32_t n = 0; n < nModels; n++) {
    if(Models[n]) {
      FreeCModel(Models[n]);
    }
  }
  
  Free(Models);
}

void PrintModelInfo(const char *filename) {
  if(!filename) {
    fprintf(stderr, "Error: No filename provided\n");
    return;
  }
  
  FILE *F = Fopen(filename, "rb");

  // Read file header
  ModelHeader header;
  if(fread(&header, sizeof(ModelHeader), 1, F) != 1) {
    fprintf(stderr, "Error reading model file header\n");
    Fclose(F);
    return;
  }

  // Validate header
  if(header.magic != MODEL_MAGIC_NUMBER) {
    fprintf(stderr, "Error: Invalid model file format (wrong magic number)\n");
    Fclose(F);
    return;
  }

  // Print header info
  time_t timestamp = (time_t)header.timestamp;
  char timeStr[100];
  strftime(timeStr, sizeof(timeStr), "%Y-%m-%d %H:%M:%S", localtime(&timestamp));

  fprintf(stderr, "==[ MODEL FILE INFO ]=================\n");
  fprintf(stderr, "Format version ..................... %u\n", header.version);
  fprintf(stderr, "Number of models ................... %u\n", header.nModels);
  fprintf(stderr, "Alphabet size ...................... %u\n", header.alphabetSize);
  fprintf(stderr, "Max hash collisions ................ %u\n", header.maxCollisions);
  fprintf(stderr, "Created on ......................... %s\n", timeStr);
  fprintf(stderr, "\n");

  // Print model info
  fprintf(stderr, "==[ MODELS ]=======================\n");
  for(uint32_t n = 0; n < header.nModels; n++) {
    ModelMeta entryHeader;
    if(fread(&entryHeader, sizeof(ModelMeta), 1, F) != 1) {
      fprintf(stderr, "Error reading model entry header for model %u\n", n);
      Fclose(F);
      return;
    }

    fprintf(stderr, "[Model %u]\n", n+1);
    fprintf(stderr, "  [+] Context order ................ %u\n", entryHeader.ctx);
    fprintf(stderr, "  [+] Alpha denominator ............ %u\n", entryHeader.alphaDen);
    fprintf(stderr, "  [+] Inverted repeats ............. %s\n",
           entryHeader.ir == 0 ? "no" : "yes");
    fprintf(stderr, "  [+] Storage mode ................. %s\n",
           entryHeader.mode == ARRAY_MODE ? "array" : "hash table");
    fprintf(stderr, "  [+] Number of models ............. %lu\n", entryHeader.nPModels);

    if(entryHeader.edits != 0) {
      fprintf(stderr, "  [+] Allowable substitutions ...... %u\n", entryHeader.edits);
      fprintf(stderr, "  [+] Substitutions alpha den ...... %u\n", entryHeader.eDen);
    }

    // Skip model data for display purposes
    if(entryHeader.mode == HASH_TABLE_MODE) {
      // Skip index array
      Fseeko(F, HASH_SIZE * sizeof(ENTMAX), SEEK_CUR);

      // Skip hash entries (more efficiently)
      Fseeko(F, HASH_SIZE * header.maxCollisions * sizeof(Entry), SEEK_CUR);
    } else if(entryHeader.mode == ARRAY_MODE) {
      // Skip array
      Fseeko(F, (entryHeader.nPModels << 2) * sizeof(ACC), SEEK_CUR);
    }
  }
  fprintf(stderr, "\n");

  Fclose(F);
}