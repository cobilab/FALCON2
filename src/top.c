#include "top.h"
#include "mem.h"
#include <stdio.h>
#include <stdlib.h>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void CopyStringPart(uint8_t *a, uint8_t *b){
  uint32_t r;
  for(r = 0 ; r < MAX_NAME-1 ; ++r){
    a[r] = b[r];
    if(b[r] == '\0') 
      break;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int SortByValue(const void *a, const void *b){ 
  VT *ia = (VT *) a;
  VT *ib = (VT *) b;
  if     (ia->value < ib->value) return -1;
  else if(ia->value > ib->value) return 1;
  else                           return 0;
  } 

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

TOP *CreateTop(uint32_t size){
  uint32_t n;
  TOP *T  = (TOP *) Calloc(1, sizeof(TOP));
  T->size = size + 1;
  T->V    = (VT  *) Calloc(T->size, sizeof(VT));
  for(n = 0 ; n < T->size ; ++n){
    T->V[n].value = 1.0;
    T->V[n].name  = (uint8_t *) Calloc(MAX_NAME, sizeof(uint8_t));
    T->V[n].size  = 1;
    #ifdef LOCAL_SIMILARITY
    T->V[n].iPos  = 1;
    T->V[n].ePos  = 1;
    #endif
    }
  return T;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void AddElement(VT *Vt, double value, uint8_t *nm, uint64_t size){
  CopyStringPart(Vt->name, nm);
  Vt->value    = value;
  Vt->size     = size;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void AddElementWithDb(VT *Vt, double value, uint8_t *nm, uint64_t size, uint32_t dbIndex){
  CopyStringPart(Vt->name, nm);
  Vt->value    = value;
  Vt->size     = size;
  Vt->dbIndex = dbIndex;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef LOCAL_SIMILARITY
void AddElementWP(VT *Vt, double value, uint8_t *nm, uint64_t size, uint64_t 
iPos, uint64_t ePos){
  CopyStringPart(Vt->name, nm);
  Vt->value = value;
  Vt->size  = size;
  Vt->iPos  = iPos;
  Vt->ePos  = ePos;
  }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef LOCAL_SIMILARITY
void AddElementWPWithDb(VT *Vt, double value, uint8_t *nm, uint64_t size, uint64_t
iPos, uint64_t ePos, uint32_t dbIndex){
  CopyStringPart(Vt->name, nm);
  Vt->value = value;
  Vt->size  = size;
  Vt->iPos  = iPos;
  Vt->ePos  = ePos;
  Vt->dbIndex = dbIndex;
}
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateTop(double bits, uint8_t *nm, TOP *T, uint64_t size){
  uint32_t last = T->size - 1;
  if(T->id < last){
    AddElement(&T->V[T->id], bits, nm, size);
    qsort(T->V, T->id+1, sizeof(VT), SortByValue); 
    }
  else if(T->id == last){
    AddElement(&T->V[last], bits, nm, size);
    qsort(T->V, T->size, sizeof(VT), SortByValue);
    }
  else{ // real NRC = 1.0-bits
    if(T->V[last].value > bits){
      AddElement(&T->V[last], bits, nm, size);
      qsort(T->V, T->size, sizeof(VT), SortByValue);
      }
    }
  T->id++;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateTopWithDB(double bits, uint8_t *nm, TOP *T, uint64_t size, uint32_t dbIndex){
  uint32_t last = T->size - 1;
  if(T->id < last){
    AddElementWithDb(&T->V[T->id], bits, nm, size, dbIndex);
    qsort(T->V, T->id+1, sizeof(VT), SortByValue);
  }
  else if(T->id == last){
    AddElementWithDb(&T->V[last], bits, nm, size, dbIndex);
    qsort(T->V, T->size, sizeof(VT), SortByValue);
  }
  else{ // real NRC = 1.0-bits
    if(T->V[last].value > bits){
      AddElementWithDb(&T->V[last], bits, nm, size, dbIndex);
      qsort(T->V, T->size, sizeof(VT), SortByValue);
    }
  }
  T->id++;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef LOCAL_SIMILARITY
void UpdateTopWP(double bits, uint8_t *nm, TOP *T, uint64_t size, uint64_t 
iPos, uint64_t ePos){
  uint32_t last = T->size - 1;
  if(T->id < last){
    AddElementWP(&T->V[T->id], bits, nm, size, iPos, ePos);
    qsort(T->V, T->id+1, sizeof(VT), SortByValue);
    }
  else if(T->id == last){
    AddElementWP(&T->V[last], bits, nm, size, iPos, ePos);
    qsort(T->V, T->size, sizeof(VT), SortByValue);
    }
  else{ // real NRC = 1.0-bits
    if(T->V[last].value > bits){
      AddElementWP(&T->V[last], bits, nm, size, iPos, ePos);
      qsort(T->V, T->size, sizeof(VT), SortByValue);
      }
    }
  T->id++;
  }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef LOCAL_SIMILARITY
void UpdateTopWPWithDb(double bits, uint8_t *nm, TOP *T, uint64_t size, uint64_t
iPos, uint64_t ePos, uint32_t dbIndex){
  uint32_t last = T->size - 1;
  if(T->id < last){
    AddElementWPWithDb(&T->V[T->id], bits, nm, size, iPos, ePos, dbIndex);
    qsort(T->V, T->id+1, sizeof(VT), SortByValue);
  }
  else if(T->id == last){
    AddElementWPWithDb(&T->V[last], bits, nm, size, iPos, ePos, dbIndex);
    qsort(T->V, T->size, sizeof(VT), SortByValue);
  }
  else{ // real NRC = 1.0-bits
    if(T->V[last].value > bits){
      AddElementWPWithDb(&T->V[last], bits, nm, size, iPos, ePos, dbIndex);
      qsort(T->V, T->size, sizeof(VT), SortByValue);
    }
  }
  T->id++;
}
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintTop(FILE *F, TOP *Top, uint32_t size, char **dbFiles){
  uint32_t n;
  double pttmp;
  if(size > Top->size){
    fprintf(stderr, "  [x] Error: top is larger than size!\n");
    exit(1);
    }

  // Print header with added Database column
  fprintf(F, "Rank\tSize\tSimilarity\tSequence\tDatabase\n");

  for(n = 0 ; n < size ; ++n){
    pttmp = (1.0-Top->V[n].value) * 100.0;
    fprintf(F, "%u\t%"PRIu64"\t%6.3lf\t%s\t%s\n",
            n+1,
            Top->V[n].size,
            pttmp,
            Top->V[n].name,
            dbFiles ? dbFiles[Top->V[n].dbIndex] : "unknown");
    if(pttmp == 0.0)
      return;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef LOCAL_SIMILARITY
void PrintTopWP(FILE *F, TOP *Top, uint32_t size, char **dbFiles){
  uint32_t n;
  double pttmp;
  if(size > Top->size){
    fprintf(stderr, "  [x] Error: top is larger than size!\n");
    exit(1);
    }

  // Print header with added Database column
  fprintf(F, "Rank\tSize\tSimilarity\tSequence\tDatabase\n");

  for(n = 0 ; n < size ; ++n){
    pttmp = (1.0-Top->V[n].value) * 100.0;
    fprintf(F, "%u\t%"PRIu64"\t%6.3lf\t%s\t%"PRIu64"\t%"PRIu64"\t%s\n", n+1,
    Top->V[n].size, pttmp, Top->V[n].name, Top->V[n].iPos, Top->V[n].ePos, dbFiles ? dbFiles[Top->V[n].dbIndex] : "unknown");
    if(pttmp == 0.0)
      return;
    }
  }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintTopInfo(TOP *Top, uint32_t size, char **dbFiles){
  uint32_t n;
  if(size > Top->size){
    fprintf(stderr, "  [x] Error: top is larger than size!\n");
    exit(1);
    }
  fprintf(stderr, "  [*] Top %u:\n", size);

  // Print header with added Database column
  fprintf(stderr, "  [*] Rank\tSize\tSimilarity\tSequence\tDatabase\n");

  for(n = 0 ; n < size ; ++n)
    fprintf(stderr, "  [*] %u \t%"PRIu64"\t%7.4lf\t%s\t%s\n", n+1, Top->V[n].size,
    (1.0-Top->V[n].value)*100.0, Top->V[n].name, dbFiles ? dbFiles[Top->V[n].dbIndex] : "unknown");
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef LOCAL_SIMILARITY
void PrintTopInfoWP(TOP *Top, uint32_t size, char **dbFiles){
  uint32_t n;
  if(size > Top->size){
    fprintf(stderr, "  [x] Error: top is larger than size!\n");
    exit(1);
    }
  fprintf(stderr, "  [*] Top %u:\n", size);

  // Print header with added Database column
  fprintf(stderr, "  [*] Rank\tSize\tSimilarity\tSequence\tDatabase\n");

  for(n = 0 ; n < size ; ++n)
    fprintf(stderr, "  [*] %u \t%"PRIu64"\t%7.4lf\t%s\t%"PRIu64"\t%"PRIu64"\t%s\n",
    n+1, Top->V[n].size, (1.0-Top->V[n].value)*100.0, Top->V[n].name, 
    Top->V[n].iPos, Top->V[n].ePos, dbFiles ? dbFiles[Top->V[n].dbIndex] : "unknown");
  }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void DeleteTop(TOP *T){
  uint32_t n;
  for(n = 0 ; n < T->size ; ++n)
    Free(T->V[n].name);
  Free(T->V);
  Free(T);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
