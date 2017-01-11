////////////////////////////////////////////////////////////////////// 
// clusters/SparseLikelihood.h 
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
 
#ifndef __SPARSELIKELIHOOD_H__
#define __SPARSELIKELIHOOD_H__

#include "HaploGraph.h"
#include "HaploSet.h"
#include "MathMatrix.h"
#include "HaploTree.h"

class SparseLikelihood
   {
   public:
      int    markers;
      int    haplos;
      Vector frequencies;

      SparseLikelihood()  { trees = NULL; index = NULL; tol = 1e-7; ztol = 1e-10; secondPass = false; }
      ~SparseLikelihood()
         { if (index != NULL) delete [] index;
           if (trees != NULL) delete [] trees; }

      void Initialize(const IntArray & alleleCounts);
      void Initialize(const IntArray & allelicCounts, const IntArray & haplos, const Vector & freqs);
      void RandomEM(HaplotypeSets * sets);
      void EM(HaplotypeSets * sets, bool quiet = false);

      void Print(double minfreq, HaplotypeSets * sets = NULL);
      void PrintExtras(double minfreq);

      void RetrieveHaplotypes(IntArray & haplos, Vector & freqs, double minFreq);

      double SetLikelihood(HaplotypeSet * set);
      double GraphLikelihood(HaplotypeGraph * graph);
      double GraphLikelihood(HaplotypeGraph * graph, int tree, int i);
      double CurrentLikelihood(HaplotypeGraph * graph, int tree);

      // Returns log of constant that has been used to scale likelihoods
      // for all component graphs in a haplotype set.
      double GetLogOffset(HaplotypeSet * set);

      void   ArgMax(HaplotypeGraph * graph, int i = 0);
      void   ExtractConfiguration(HaplotypeGraph * graph, IntArray & configuration, bool sample);
      void   ExtractConfiguration(HaplotypeGraph * graph, bool sample, int index = 0);
      void   SelectState(HaplotypeGraph * graph);
      void   SampleState(HaplotypeGraph * graph);

      double best_llk;

      double tol, ztol;

   private:
      IntArray * index;
      IntArray last;
      IntArray new_index;
      IntArray selectedState;
      Vector   marginals;

      bool     secondPass;
      double   setLikelihood;

      Vector   new_frequencies;

      bool UpdateIndex(HaplotypeGraph * graph, int tree, int i);
      bool ExpandIndex(HaplotypeGraph * graph, int tree, int i);
      bool SelectAllele(int tree, int i, int allele);

      void CoreEM(HaplotypeSets * sets, bool quiet = false);

      // Each tree has up to (2^treeBits) markers
      int treeBits;

      // Number of trees
      int treeCount;

      // Array of trees
      HaploTree * trees;

      // Likelihood factors within each graph
      Matrix   factorMatrix;
      Vector * factors;
   };

#endif

 
