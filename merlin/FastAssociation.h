////////////////////////////////////////////////////////////////////// 
// merlin/FastAssociation.h 
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
 
#ifndef __FASTASSOCIATION_H__
#define __FASTASSOCIATION_H__

#include "AssociationAnalysis.h"

class FastAssociationAnalysis : protected AssociationAnalysis
   {
   public:
      void AnalyseAssociation(FamilyAnalysis & engine);

      virtual ~FastAssociationAnalysis() {}

      static double fastFilter;

      static void OutputFastModels(const String & prefix, Pedigree & ped);

   private:
      // This variable is static so it can be updated across chromosomes
      static RefinedQtlModel refinedFastModels;
   };

#endif

 
