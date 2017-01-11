////////////////////////////////////////////////////////////////////// 
// merlin/MerlinHaplotype.h 
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
 
#ifndef __MERLIN_HAPLOTYPE_H__
#define __MERLIN_HAPLOTYPE_H__

#include "ConquerHaplotyping.h"
#include "MerlinBitSet.h"
#include "TreeInfo.h"
#include "MathMatrix.h"
#include "MerlinCore.h"

#include <stdio.h>

// Haplotyping routines
//

class FamilyAnalysis;
class HaplotypeChain;

class MerlinHaplotype
   {
   public:
      // Constructor and destructor
      MerlinHaplotype();
      ~MerlinHaplotype();

      // User interfaces
      void MostLikely(FamilyAnalysis & f)
         { HaplotypeFamily(f, false); }
      void Sample(FamilyAnalysis & f)
         { HaplotypeFamily(f, true); }
      void All(FamilyAnalysis & f);

      // The work horse
      void HaplotypeFamily(FamilyAnalysis & which, bool sample);
      void UninformativeFamily(FamilyAnalysis & which, bool sample);

      // File management
      void OpenFile();
      void CloseFile();

      static bool outputFounders;
      static bool outputHorizontal;
      static bool outputGraph;

   private:
      // File for storing haplotypes
      FILE * output;
      FILE * foutput;
      FILE * goutput;

      // The output broker
      void OutputHaplotypes(StringArray * haplo, StringArray & recomb,
                            const char * header_format, ...);
      void HorizontalOutput(StringArray * haplo);

      // Founder haplotypes output broker
      void OutputFounders(StringArray * haplo, const char * label, int weight = 1);

      // Information for the current family
      FamilyAnalysis * family;

      // The most likely path through inheritance vector space
      // at each informative marker
      int        maximum_markers;
      IntArray * inheritanceVector;
      int        maximum_founders;

      // Memory allocation
      void AllocateVectors(int markers, int founders);

      // Conditional probabilities
      bool ScoreConditional(Tree * right, Vector & scale);

      // We will need to calculate likelihoods for each tree
      TreeInfo stats;

      // Listing all inheritance vectors
      void ListAllRecursively(Tree & tree, int node, int bit, double count,
                              HaplotypeChain & head);

      // Locating the best inheritance vector
      void FindBest(Tree & tree, IntArray & best_vector);
      void FindBestRecursively(Tree & tree, int node, int bit, IntArray & best);

      // Locating the best inheritance vector conditional on neighbour
      void FindBest(Tree & tree, IntArray & best_vector, IntArray & flanker);
      void FindBestRecursively(Tree & tree, int node, int bit, double lk,
                               IntArray & best_vector, IntArray & flanker);

      // Sampling among possible inheritance vectors
      void Sample(Tree & tree, IntArray & best_vector);
      void SampleRecursively(Tree & tree, int node, int bit, IntArray & best,
                             double weight = 1.0);

      // Sampling among possible inheritance vectors given a previous choice
      void Sample(Tree & tree, IntArray & best_vector, IntArray & flanker);
      void SampleRecursively(Tree & tree, int node, int bit, double lk,
                             IntArray & best, IntArray & flanker);

      // Label haplotypes according to inheritance vector
      int    LabelCluster(int marker, IntArray & vec, StringArray * states, bool sample);
      void   LabelChromosomes(int marker, IntArray & vec, StringArray & states);
      void   LabelDescent(IntArray & vec, StringArray & states);
      String NameAllele(MarkerInfo * markerInfo, longint binary_code);

      // Label recombinations according to neighbouring inheritance vectors
      void LabelRecombinants(IntArray& vec1, IntArray& vec2, String & recomb,
                             bool sample);

      // Used when searching for the best inheritance vector
      IntArray current;
      double   sum;

      // Distance between current inheritance vectors
      double   distance;
      double   theta[2];
      double   adjusted_theta;

      // These track recombination between a pair of inheritance vectors
      IntArray recombination;

      // This array tracks founder flips and symmetries
      IntArray flipKey;

      // MerlinBitSets separates inheritance patterns into independent
      // components to facilitate sampling
      //
      MerlinBitSets bitSet;

      // Utility function for managing output files
      FILE * OpenFile(const char * extension);
   };

class HaplotypeChain
   {
   public:
      StringArray    * haplo;
      HaplotypeChain * next;
      double count;
      int    index;

      HaplotypeChain(int serial = 1)
         {
         count = 0.0;
         next = NULL;
         haplo = NULL;
         index = serial;
         }

      ~HaplotypeChain()
         {
         if (haplo != NULL) delete [] haplo;
         if (next != NULL) delete next;
         }

      void Append(StringArray * states, double vectors, int markers);
      int  Length();
   };

#endif

 
