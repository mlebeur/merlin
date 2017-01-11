////////////////////////////////////////////////////////////////////// 
// regress/RegressKinship.h 
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
 
#ifndef __MERLIN_KINSHIP_H__
#define __MERLIN_KINSHIP_H__

#include "MathMatrix.h"
#include "Pedigree.h"
#include "Mantra.h"
#include "Tree.h"

#include <stdio.h>

class RegressKinship
   {
   public:
      RegressKinship(Mantra & m);
      ~RegressKinship();

      double Calculate(Tree & tree);
      double Score(Tree & tree, int node, int bit, int start = 0, int pos = 0);

      void   SelectFamily(Pedigree * p, Family * f);

      double Retrieve(int person1, int person2);
      double Retrieve(int person1, int person2, int person3, int person4);

      void   ResetSymmetries();
      void   Unravel(int state);

      void   UpdateScores(IntArray & k, int start_pos, int end_pos, double sum);

      Pedigree * GetPedigree() { return ped; }
      Family   * GetFamily()   { return family; }

   private:
      Pedigree * ped;
      Family   * family;
      Mantra   & mantra;

      int count;
      int founders;

      // Array of kinship coefficients for the current gene flow pattern
      IntArray kin;

      // Vector for storing observed kinship coefficients
      Vector   kinship;

      // Vector for storing cross-products of kinship coefficients
      Vector   kinship2;

      // Extra storage for unravelling founder couple symmetries
      IntArray * symmetries;
      IntArray   unraveled;
      int        sym_size;

      // How many equivalent gene flow states must be unraveled?
      int        alternate_states;
   };

#endif

 
