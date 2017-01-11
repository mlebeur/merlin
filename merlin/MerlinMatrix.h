////////////////////////////////////////////////////////////////////// 
// merlin/MerlinMatrix.h 
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
 
#ifndef __MERLIN_MATRIX_H__
#define __MERLIN_MATRIX_H__

#include "Pedigree.h"
#include "Mantra.h"
#include "Tree.h"
#include "MerlinCore.h"
#include "StringMap.h"

#include <stdio.h>

class MerlinMatrix
   {
   public:
      MerlinMatrix(Mantra & m);
      ~MerlinMatrix();

      double Calculate(Tree & tree, const char * label);
      double Score(Tree & tree, int node, int bit, int start = 0, int pos = 0);
      void   Output(const char * label, double scale);

      void   OpenFile();
      void   CloseFile();

      void   SelectFamily(Pedigree * p, Family * f);

   private:
      Pedigree * ped;
      Family   * family;
      Mantra   & mantra;

      // Output files
      FILE       * matrixfile;

      // The current kinship matrix, in linear form
      String       matrix;

      // Kinship matrices for the current pedigree
      StringIntMap matrices;

      // Matrix probabilites for the current pedigree
      Vector       probabilities;

      // These are used to track couple symmetries
      IntArray * symmetries;
      int        size;
      int        couples;

      // These routines unravel equivalent kinship matrices
      void ResetSymmetries();
      const char * Unravel(const String & in, String & out, int index);
   };

#endif

 
 
