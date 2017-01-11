////////////////////////////////////////////////////////////////////// 
// merlin/GenotypeInference.h 
// (c) 2000-2007 Goncalo Abecasis
// 
// This file is distributed as part of the MERLIN source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile MERLIN.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Tuesday December 18, 2007
// 
 
#ifndef __GENOTYPEINFERENCE_H__
#define __GENOTYPEINFERENCE_H__

#include "Mantra.h"
#include "TreeInfo.h"
#include "Tree.h"
#include "MathFloatVector.h"

class GenotypeInference
   {
   public:
      // These variables control exactly what should be infered
      static bool inferBest;          // infer the most likely genotype
      static bool inferExpected;      // infer the expected genotype
      static bool inferProbabilities; // infer genotype probabilities

      // This variable should be set when previously infered data is merged
      // into the pedigree
      static bool offline;

      // Initializes temporary memory allocation
      void     AllocateMemory(Pedigree & ped);
      int      SelectMarkers(IntArray & markers, Vector & markerPositions, Vector & analysisPositions);

      // These provide information on inferred genotypes for any individual
      double   GetExpectedGenotype(int individual, int marker);
      double   GetExpectedGenotype(Pedigree & ped, int individual, int marker);
      double   GetProbability(int individual, int marker, int genotype);

      // This function carries out genotype inference
      void     InferGenotypes(Mantra & mantra, TreeInfo & stats,
               Tree & withMarker, Tree & without, Tree & single,
               int marker);

      // Checks whether genotype inference was carried out for a particular marker
      bool     resultsAvailable(int marker);
      int      availableMarkers() { return markers; }

      // Generate pedigree and data files summarizing inferred genotypes
      void OutputGenotypes(Pedigree & ped, IntArray & markerList);

      // These functions are used to support pedigrees with precalculate
      // expected genotypes
      void CreateOfflineMarkers(Pedigree & ped);
      void CreateOfflineLookup();

      // This function checks that offline calculations for allele frequencies
      // and dosages are consistent
      void CheckFrequencies(Pedigree & ped);

      // This function updates allele frequencies based on estimated allele doses
      void UpdateFrequencies(Pedigree & ped, bool foundersOnly = true);

   private:
      // Key constants
      int markers;

      // Probability that an individual is homozygous for allele 1 or 2
      // Each vector has n_individuals * n_marker entries
      FloatVector probability[2];

      // This function maps an individual to a slot within the probability array
      int      GetSlot(int individual, int marker);

      // List of selected markers
      IntArray markerKey;

      // Frequencies for selected markers
      Vector   frequencies;

      // Track which individuals are males and females,
      // which is required for X chromosome inference
      String   isFemale;
   };

#endif

 
