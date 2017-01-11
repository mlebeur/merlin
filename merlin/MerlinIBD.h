////////////////////////////////////////////////////////////////////// 
// merlin/MerlinIBD.h 
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
 
#ifndef __MERLIN_IBD_H__
#define __MERLIN_IBD_H__

#include "MathMatrix.h"
#include "Pedigree.h"
#include "Mantra.h"
#include "Tree.h"
#include "MerlinCore.h"

#include <stdio.h>

class MerlinIBD
   {
   public:
      MerlinIBD(Mantra & m) : mantra(m) {}

      double Calculate(Tree & tree, const char * label);
      double Score(Tree & tree, int node, int bit, int start = 0);
      void   Output(const char * label, double scale);

      void   OpenFile();
      void   CloseFile();

      void   SelectFamily(Pedigree * p, Family * f);

   private:
      Pedigree * ped;
      Family   * family;
      Mantra   & mantra;

      // Specialized routines for trees with symmetries
      double WeightedScore(Tree & tree, int node, int bit, int start,
                           double weight = 2.0);
      double PairWiseIBD1(int a, int b, int c, int d);
      double PairWiseIBD2(int a, int b, int c, int d);

      // Output files
      FILE     * ibdfile;

      // Temporary storage
      Matrix   IBD1, IBD2;

      // Uninformative meiosis in trees with symmetries
      IntArray uninformative;
   };

#endif

 
 
