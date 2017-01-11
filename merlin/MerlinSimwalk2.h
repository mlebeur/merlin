////////////////////////////////////////////////////////////////////// 
// merlin/MerlinSimwalk2.h 
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
 
#ifndef __MERLIN_SIMWALK2__
#define __MERLIN_SIMWALK2__

#include "MerlinCore.h"
#include "KongAndCox.h"

#include <stdio.h>

class FamilyAnalysis;

// This class outputs non-parametric statistics in a format
// that Simwalk2 can use

class NPLforSimwalk2
   {
   public:
      NPLforSimwalk2(FamilyAnalysis & f) : family(f)
         { }

      void OpenFiles();
      void CloseFiles();
      void Output();
      void Skip();

   private:
      FamilyAnalysis & family;
      FILE * output;

      void OutputValue(double value);
   };

#endif

 
 
