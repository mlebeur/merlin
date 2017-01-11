////////////////////////////////////////////////////////////////////// 
// clusters/HaploGraph.h 
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
 
#ifndef __HAPLOGRAPH_H__
#define __HAPLOGRAPH_H__

#include "HaploFamily.h"
#include "Unknown.h"

#include <stdio.h>

class HaplotypeGraph
   {
   public:
      // Number of haplotypes in this graph
      int count;

      // Allele state for each marker
      IntArray * haplotypes;

      // Index of missing genotypes and other ambiguities
      SetOfUnknowns unknowns;

      // Index of last genotyped marker for each haplotype
      IntArray lastGenotype;

      // Relative weight for this graph
      double weight;

      // Probability of observing this haplotype set
      double likelihood;

      // Scale factor for the likelihood, in 2^-128 units
      int    scale;

      HaplotypeGraph()
         { haplotypes = NULL; count = 0;}

      ~HaplotypeGraph()
         { if (haplotypes != NULL) delete [] haplotypes; }

      void LoadFromFile(FILE * file, int haplos);
      void LoadFromMemory(FounderGraph * founderGraph, int haplos, int markers);

      // Number of alternative solutions for this haplotype graph
      double Complexity(IntArray & alleleCounts);

   private:
      void Dimension(int size, int markers);
   };

#endif


 
