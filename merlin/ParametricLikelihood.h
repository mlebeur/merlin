////////////////////////////////////////////////////////////////////// 
// merlin/ParametricLikelihood.h 
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
 
#ifndef _PARAMETRICLIKELIHOOD_H_
#define _PARAMETRICLIKELIHOOD_H_

#include "Houdini.h"
#include "DiseaseModel.h"

class ParametricLikelihood : public FuzzyInheritanceTree
   {
   public:
      ParametricLikelihood();
      ~ParametricLikelihood();

      DisModel * model;

   protected:
      virtual int    CountAlleles(Mantra & m);
      virtual bool   ShortScoreVectors(Mantra & m);
      virtual void   FillPenetranceMatrix(Mantra & m);
      virtual bool   isPhenotyped(Mantra & m, int i);

      virtual double GetAlleleFrequency(int founder, int allele);
      virtual void   UpdateAlleles(int founder, int carrier);
      virtual void   UndoAlleles(int founder, int carrier);

      virtual double JointProbability(Mantra & m);

      virtual void   SetupAlleleList(Mantra & m);
      virtual void   CleanUpAlleleList(Mantra & m);

      virtual bool IdenticalPenetrances(int founder1, int founder2);
      virtual void SwapPenetrances(int founder1, int founder2);

   private:
       IntArray compatibleAlleles;
   };

#endif

 
