////////////////////////////////////////////////////////////////////// 
// merlin/Magic.h 
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
 
#ifndef __MAGIC_H__
#define __MAGIC_H__

#include "Pedigree.h"
#include "Tree.h"
#include "Mantra.h"

class InheritanceTree : public Tree
   {
   public:
      InheritanceTree();
      ~InheritanceTree();

      // Scores the likelihoods of all inheritance vectors
      void ScoreVectors(Mantra & m);

      // Scores the likelihood of a single inheritance vector
      double  EvaluateInheritanceVector(Mantra & m, int * vector,
              int * graph, longint * alleles, longint * alleles2, int * fixed);

      // When mergeStrategy is set to MERGE_ALL, branches with identical
      // likelihood outcomes are merged even if they imply different allelic
      // states for founders. This option is typically disabled when listing
      // all non-recombinant haplotypes.
      static int mergingStrategy;

      // These are the possible merging strategies ...
      #define MERGE_ALL      0
      #define MERGE_ZEROS    1
      #define MERGE_BOOLEAN  2

   private:
      // Stacks for the most commonly used array
      IntArray * graph_stack, * dirty_stack, * fixed_stack;
      LongArray * alleles_stack, * alleles2_stack;
      int stack_depth;

      void AllocateStack(int new_depth);

      // Temporary arrays for FindRedundancy and CalculateLikelihood
      IntArray fixed2, stack;

      // Setup the graph, alleles and fixed arrays
      bool PrepareForGraphing(Mantra & m, int * graph, longint * alleles,
                              int * fixed, longint * alleles2, int * dirty);

      // Recurse throught the pedigree and store the likelihood of all
      // possible inheritance vectors
      int    ScoreRecursive(Mantra & m, int current_bit, int last_pivot,
                            IntArray & graph, LongArray & alleles,
                            IntArray & fixed, LongArray & alleles2,
                            IntArray & dirty);

      // Check if branching is required (while recursing through the pedigree)
      int    FindRedundancy(Mantra & m, int pivot, int * graph, longint * alleles,
                            int * fixed, longint * alleles2, int * dirty);

      // Calculate the likelihood for a single inheritance vector
      double CalculateLikelihood(Mantra & m, int * graph, longint * alleles,
                                 int * fixed, longint * alleles2, int * dirty);
   };


#endif
 
