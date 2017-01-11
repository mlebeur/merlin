////////////////////////////////////////////////////////////////////// 
// merlin/MerlinKinship15.h 
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
 
#ifndef __MERLIN_KINSHIP15_H__
#define __MERLIN_KINSHIP15_H__

#include "MathMatrix.h"
#include "Pedigree.h"
#include "Mantra.h"
#include "Tree.h"
#include "MerlinCore.h"

#include <stdio.h>

class MerlinKinship15
   {
   public:
      MerlinKinship15(Mantra & m);
      ~MerlinKinship15();

      double Calculate(Tree & tree, const char * label = NULL);
      double Score(Tree & tree, int node, int bit, int start = 0);
      double WScore(Tree & tree, int node, int bit, int start, double weight);
      void   Output(const char * label, double scale);

      void   OpenFile();
      void   CloseFile();

      void   SelectFamily(Pedigree * p, Family * f);

   private:
      Pedigree * ped;
      Family   * family;
      Mantra   & mantra;

      // Useful information about the current pedigree
      //
      bool isTwoGeneration;
      bool isInbred;

      // Routines for scoring generalized kinships in small sets of alleles
      //

      // Single blocks of 2, 3 and 4 alleles IBD
      double   Kinship2(int allele1, int allele2);
      double   Kinship3(int allele1, int allele2, int allele3);
      double   Kinship4(int allele1, int allele2, int allele3, int allele4);

      // Two blocks with (3,1), (2,2), (2,1) and (1,1) alleles IBD
      double   Kinship31(int allele1, int allele2, int allele3, int allele4);
      double   Kinship22(int allele1, int allele2, int allele3, int allele4);
      double   Kinship21(int allele1, int allele2, int allele3);
      double   Kinship11(int allele1, int allele2)
         { return 1.0 - Kinship2(allele1, allele2); }

      // Three blocks with (2,1,1) or (1,1,1) alleles IBD
      double   Kinship211(int allele1, int allele2, int allele3, int allele4);
      double   Kinship111(int allele1, int allele2, int allele3);

      // Routine for ordering pairs of alleles
      void Order(int &allele1, int &allele2)
         {
         if (allele1 > allele2)
            {
            int swap = allele1;
            allele1 = allele2;
            allele2 = swap;
            }
         }

      void     Order(int &allele1, int &state1, int &allele2, int &state2)
         {
         int swap;

         if (state1 > state2)
            {
            swap = state1; state1 = state2; state2 = swap;
            swap = allele1; allele1 = allele2; allele2 = swap;
            }
         }

      // Output files
      FILE     * kinfile;

      // Temporary storage
      Matrix   kinship[14];

      // Routine for initialize kinship coeffiecients
      void InitMatrix(int coefficient)
         {
         kinship[coefficient].Dimension(mantra.n, mantra.n);
         kinship[coefficient].Zero();
         }
   };

#endif


 
 
