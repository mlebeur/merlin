////////////////////////////////////////////////////////////////////// 
// clusters/Likelihood.h 
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
 
#ifndef __LIKELIHOOD_H__
#define __LIKELIHOOD_H__

#include "HaploGraph.h"
#include "HaploSet.h"
#include "MathVector.h"

#define PAIRWISE_COUPLING       1  /* The first two alleles are coupled */
#define PAIRWISE_EQUILIBRIUM    0  /* The markers are in equilibrium */
#define PAIRWISE_REPULSION     -1  /* The first two alleles are in repulsion */

class Likelihood
   {
   public:
      // Information about haplotypes
      int    haplos;
      Vector frequencies;

      // Information about individual markers
      int      markers;
      IntArray alleleCounts;

      Likelihood();
      ~Likelihood() { if (index != NULL) delete [] index; }

      void Initialize(const IntArray & allelicCounts);
      void Initialize(const IntArray & allelicCounts, const IntArray & haplos, const Vector & freqs);
      void EditedEM(HaplotypeSets * sets);
      void RandomEM(HaplotypeSets * sets);
      void EM(HaplotypeSets * sets);

      void Print(double minfreq, HaplotypeSets * sets = NULL);
      void PrintExtras();
      void SaveSummary();

      void RetrieveHaplotypes(IntArray & haplos, Vector & freqs, double minFreq);

      double SetLikelihood(HaplotypeSet * set);
      double GraphLikelihood(HaplotypeGraph * graph);
      double GraphLikelihood(HaplotypeGraph * graph, int i);
      double CurrentLikelihood(HaplotypeGraph * graph);

      // Calculate pairwise LD coefficients for bi-allelic markers
      bool   PairwiseLD(double & dprime, double & deltasq, int & coupling);
      bool   PairwiseLD(double & dprime, double & deltasq)
         { int coupling; return PairwiseLD(dprime, deltasq, coupling); }

      double best_llk;
      Vector best_frequencies;

      double tol, ztol;

      // The following variables are used for favoring similar haplotypes
      void   CheckSmoothing();
      double pseudoMutation;
      bool   gradual_smoothing;
      bool   two_pass_smoothing;

      const  char * failText;

   private:
      IntArray * index;
      IntArray last;
      Vector   marginals;

      bool     secondPass;
      double   setLikelihood;
      double   missingHaplotypes;

      Vector   new_frequencies;

      void UpdateIndex(HaplotypeGraph * graph, int i);
      void ExpandIndex(HaplotypeGraph * graph, int i);
      void ReverseIndex(HaplotypeGraph * graph, int i, int pos);

      void CoreEM(HaplotypeSets * sets);

      void PseudoCoalescent(int start, int stop);
   };

#endif

 
