
#ifndef SERIALIZATION_H
#define SERIALIZATION_H

#include "defs.h"
#include "models.h"
#include "param.h"

// Magic number to identify valid serialized model files
#define MODEL_MAGIC_NUMBER     0x46414C434F4E4D53 // "FALCONMS" in hex (FALCON Model Serialization)
#define MODEL_VERSION          1                  // Version of the serialization format

typedef struct {
    uint64_t magic;              // Magic number for validation
    uint32_t version;            // Serialization format version
    uint32_t nModels;            // Number of models
    uint32_t alphabetSize;       // Size of alphabet
    uint64_t timestamp;          // Time when serialization was performed
    uint32_t hashSize;           // Hash table size used
    uint32_t maxCollisions;      // Maximum collisions allowed
} ModelHeader;

typedef struct {
    uint32_t ctx;                // Context order
    uint32_t alphaDen;           // Alpha denominator
    uint8_t  ir;                 // Inverted repeats flag
    uint32_t edits;              // Number of substitution edits allowed
    uint32_t eDen;               // Edit denominator for substitution
    uint32_t mode;               // Storage mode (hash table or array)
    uint64_t nPModels;           // Number of probability models
    uint32_t maxCount;           // Max counter value
    uint64_t multiplier;         // Multiplier value
    uint64_t dataSize;           // Size of model data in bytes
} ModelMeta;

/**
 * Save compression models to a file
 * 
 * @param filename The name of the file to save the models to
 * @param Models Array of models to save
 * @param nModels Number of models
 * @param col Maximum allowed hash collisions
 * @return 0 on success, negative value on error
 */
int SaveModels(const char *filename, CModel **Models, uint32_t nModels, uint32_t col);

/**
 * Load compression models from a file
 * 
 * @param filename The name of the file to load the models from
 * @param Models Pointer to array of models to be allocated
 * @param nModels Pointer to store the number of models 
 * @param col Pointer to store the collisions parameter
 * @return 0 on success, negative value on error
 */
int LoadModels(const char *filename, CModel ***Models, uint32_t *nModels, uint32_t *col);

/**
 * Free all loaded models
 * 
 * @param Models Array of models to free
 * @param nModels Number of models
 */
void FreeLoadedModels(CModel **Models, uint32_t nModels);

/**
 * Print information about serialized models in a file
 * 
 * @param filename The name of the file to examine
 */
void PrintModelInfo(const char *filename);

#endif //SERIALIZATION_H
