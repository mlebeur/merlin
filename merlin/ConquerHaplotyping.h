////////////////////////////////////////////////////////////////////// 
// merlin/ConquerHaplotyping.h 
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
 
#ifndef __CONQUER_HAPLOTYPE_H__
#define __CONQUER_HAPLOTYPE_H__

#include "TreeFlips.h"
#include "Conquer.h"
#include "Mantra.h"

class MultipointHaplotyping : public TreeFlips
   {
   public:
      double theta;
      double oneminus;
      double x;
      double scale;

      void MoveAlong(Mantra & m, double distance, double rescale = 1.0);
      void MoveAlong(Mantra & m, double theta[2], double rescale = 1.0);

      // Conditions likelihoods at flanker on those at the current location
      void Condition(Tree & flanker, double distance, double rescale = 1.0);
      void Condition(Tree & flanker, double theta[2], double rescale = 1.0);

      // Allocates variables required for likelihood conditioning
      void SetupConditioning(Mantra & m);

   private:
      // These apply a variant of the Elston and Idury algorithm
      void DivideAndConquer(int node = 0);
      void DivideAndConquer(int * isMale, int node = 0);
      void ReUnite(int left, int right);
      void SpecialReUnite(int left, int right, const int * flips);
      void DoubleReUnite(int left, int right, const int * flips);

      // These apply an approximation allowing for up to n recombinants
      void Approximation(Mantra & m, double distance, int order, double rescale);
      void Approximation(Mantra & m, double theta[2], int order, double rescale);
      void DivideAndConquer(Tree & tree, int node = 0);
      void DivideAndConquer(Tree & tree, int * isMale, int node = 0);
      void ReUnite(Tree & original, int left, int right);
      void SpecialReUnite(Tree & tree, int left, int right, const int * flips);
      void DoubleReUnite(TreeFlips & original, int left, int right, const int *flips);

      // These are used for conditioning
      double ConditionalLikelihood(int bit, double lk = 1.0, int node = 0);
      void Condition(Tree & flanker, int bit, double scale = 1.0, int node = 0);

      // Auxiliary variables for conditioning
      MerlinBitSets bitSet;
      IntArray      current;
      IntArray      previous;

      // Quantities for handling differences in male-female recombination
      double        recombination[2];
      double        complement[2];
      double        ratio[2];

      void SetMeiosisSex(bool isMale)
         {
         theta = recombination[isMale];
         oneminus = complement[isMale];
         x = ratio[isMale];
         }
   };

#endif
 
