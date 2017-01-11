////////////////////////////////////////////////////////////////////// 
// regress/AutoFit.h 
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
 
#ifndef __AUTOFIT_H__
#define __AUTOFIT_H__

#include "MathNormal.h"
#include "Pedigree.h"

class NormalEquations;

class AutoFit
   {
   public:
      AutoFit();

      void FitPolygenicModels(Pedigree & ped);

      static bool   useCovariates;
      static bool   useProbands;

      Vector means, variances, heritabilities;
   };

#endif

  
