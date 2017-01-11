////////////////////////////////////////////////////////////////////// 
// merlin/AssociationAnalysis.h 
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
 
#ifndef __ASSOCIATION_ANALYSIS_H__
#define __ASSOCIATION_ANALYSIS_H__

#include "VarianceComponents.h"
#include "GenotypeInference.h"
#include "MerlinKinship.h"

class FamilyAnalysis;

class AssociationAnalysis : protected VarianceComponents
   {
   public:
      virtual ~AssociationAnalysis() {}

      void AnalyseAssociation(FamilyAnalysis & engine);
      void PolygenicModel(Pedigree & ped, int t, NormalSet & mvn, IntArray * pheno);

      static void OutputRefinedModels(const String & prefix, Pedigree & ped);

      void OpenPerFamilyFile(int chromosome);
      void ClosePerFamilyFile(int chromosome);

   protected:
      virtual void SetupPDF(MerlinPDF & pdf, String & trait);

   private:
      // This variable is static so it can be updated across chromosomes
      static RefinedQtlModel refinedAssociationModels;

   };

#endif


 
