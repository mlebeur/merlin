////////////////////////////////////////////////////////////////////// 
// merlin/MerlinError.h 
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
 
#ifndef __MERLINERROR_H__
#define __MERLINERROR_H__

#include "MathGold.h"
#include "Pedigree.h"

class ErrorRateEstimator : public ScalarMinimizer
   {
   public:
      ErrorRateEstimator(Pedigree & p) : ped(p)
         { alleleError = trace = false; }

      virtual double f (double error_rate);
      void Estimate();
      void EstimateModel();

   private:
      Pedigree & ped;

      bool alleleError;
      bool trace;

      double CalculateLikelihood();
   };

#endif

 
