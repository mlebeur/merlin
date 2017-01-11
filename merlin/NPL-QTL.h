////////////////////////////////////////////////////////////////////// 
// merlin/NPL-QTL.h 
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
 
#ifndef __NPL_QTL_H__
#define __NPL_QTL_H__

#include "Pedigree.h"
#include "Mantra.h"
#include "Tree.h"

class NPL_QTL : public Tree
   {
   public:
      void   ScoreWithMeanZero(Mantra & m, int trait);
      void   ScoreWithSampleMean(Mantra & m, int trait);

   private:
      void   ScoreNPL(Mantra & m, int trait, double mean);
      void   ScoreMeans(Pedigree & ped);

      int    RecursivelyScoreNPL(Mantra & m, int bit);
      bool   isRedundantNPL(Mantra & m, int pivot);
      double CalculateNPL(Mantra & m);

      Vector avg;

      // Phenotypic means calculated once per run
      static Vector traitMeans;
   };

#endif
 
