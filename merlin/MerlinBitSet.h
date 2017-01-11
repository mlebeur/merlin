////////////////////////////////////////////////////////////////////// 
// merlin/MerlinBitSet.h 
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
 
#ifndef __MERLINBITSET_H__
#define __MERLINBITSET_H__

#include "Mantra.h"

#define BS_DEFAULT     1
#define BS_FOUNDER     2
#define BS_COUPLE      8

class MerlinBitSet;
class MerlinBitSets;

class MerlinBitSet
   {
   public:
      // Identification for this bit_set
      int type;
      int index;
      int last_bit;

      // Sex of parent where recombination occurs
      bool maleMeiosis;

      // Index for bits in this set
      IntArray bits;

      // Calculate likelihood this bit set conditional on previous
      double MaxLikelihood(IntArray & current, IntArray & previous, double theta);
      double Likelihood(IntArray & current, IntArray & previous, double theta);

      // Count recombinants in this bit set
      void  CountRecombinants(IntArray & current, IntArray & previous);

      // Set state of ambiguous meiosis conditional on previous selection
      void   SelectMostLikely(IntArray & current, IntArray & previous);
      void   SelectBySampling(IntArray & current, IntArray & previous, double theta);

      // Setup auxiliary data structures for founder couples
      void   SetupCouple(Mantra & m, int couple);

      // Mark recombinants and couple symmetries for specific transition
      void   MapTransition(Mantra & m, int path,
                      IntArray & current, IntArray & previous,
                      IntArray & crossover, IntArray & key);

   private:
      // For symmetric founder couples only
      IntArray p1, p2, p3;
      IntArray couple_flips;
      int couple_founder1, couple_founder2;

      // Likelihoods for deconvoluted inheritance vectors
      double lk[8];

      // Recombinant counts for each deconvoluted vector
      int rec[8];

      // Maximum likelihood vector
      int best;

      // Number of ambiguous outcomes
      int ambiguous;

      // Update individual likelihood components
      void Update(IntArray & current, IntArray & previous, double theta);
      void Rescale(double theta);

      friend class MerlinBitSets;

      static void Swap(int & a, int & b)
         {
         int tmp = a;
         a = b;
         b = tmp;
         }
   };

class MerlinBitSets
   {
   public:
      // Each BitSet contains a set of grouped bits
      int            size, count;
      MerlinBitSet * sets;

      // Index maps each bit to a corresponding set
      IntArray       Index;

      MerlinBitSets()
         {
         size = count = 0;
         sets = NULL;
         };

      ~MerlinBitSets()
         {
         if (sets) delete [] sets;
         }

      // Setup non-overlapping bit sets of founder, founder couple and other bits
      void BuildSets(Mantra & m);

      // Complete likelihood update
      void UpdateComponents(IntArray & current, IntArray & previous, double theta[2]);
      void CountRecombinants(IntArray & current, IntArray & previous);

      // Update recombinant counts within each bit set
      void UpdateRecombinationCounts(IntArray & current, IntArray & previous);

      // Partial likelihood updates for navigating down a tree
      void UpdateMaxLikelihood(int bit, double theta[2], double & likelihood,
                               IntArray & current, IntArray & previous);
      void UpdateLikelihood(int bit, double theta[2], double & likelihood,
                            IntArray & current, IntArray & previous);

      // Ambiguous meioses get resolved by these two functions
      void SelectMostLikely(IntArray & current, IntArray & previous);
      void SelectBySampling(IntArray & current, IntArray & previous, double theta[2]);

      // Short-hand access to individual sets
      MerlinBitSet & operator [] (int index)
         { return sets[index]; }

      // Label recombinants
      void BestTransition(Mantra & m,
           IntArray & current, IntArray & previous, double theta[2],
           IntArray & crossover, IntArray & key);
      void SampleTransition(Mantra & m,
           IntArray & current, IntArray & previous, double theta[2],
           IntArray & crossover, IntArray & key);
      void SampleTransition(Mantra & m,
           IntArray & current, IntArray & previous, double theta[2],
           IntArray & crossover, IntArray & key, int MaxRecombination);

   private:
      // Allocate bit sets as required
      void Dimension(int no_of_sets);

      bool RecursiveSampler(Mantra & m, IntArray & current, IntArray & previous,
                            int recombinants, double lk = 1.0, int level = 0);

      double sampled_so_far;

   };


#endif


 
