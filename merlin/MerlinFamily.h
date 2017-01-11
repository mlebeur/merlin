////////////////////////////////////////////////////////////////////// 
// merlin/MerlinFamily.h 
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
 
#ifndef __FAMILY_ANALYSIS_H__
#define __FAMILY_ANALYSIS_H__

// Task specific modules
#include "MerlinHaplotype.h"
#include "MerlinSimwalk2.h"
#include "MerlinKinship.h"
#include "MerlinKinship15.h"
#include "MerlinMatrix.h"
#include "MerlinIBD.h"
#include "VarianceComponents.h"
#include "GenotypeInference.h"
#include "AnalysisTask.h"
#include "MerlinPDF.h"
#include "TreeIndex.h"

// The MerlinCore class manages basic likelihood calculations
#include "MerlinCore.h"

class FamilyAnalysis : public MerlinCore
   {
   public:
      FamilyAnalysis(Pedigree & pedigree);
      virtual ~FamilyAnalysis();

      virtual void SetupGlobals();
      virtual int  SetupMap(int chromosome = -1);

      void SetupFiles();
      void CloseFiles();

      virtual bool SelectFamily(Family * f, bool warnOnSkip = true);
      void Analyse();

      void ShowLODs();
      void ShowLikelihood();

      void FreeMemory();
      void AbortAnalysis();
      void SkipFamily();
      void ScoreNPL();
      void ScoreNPLBounds(int index);

      int      nplCount;
      Vector   nplMean, nplInv;
      Tree      *npl;
      IndexTree *nplIndex;

      // Flags for selecting analyses to perform
      static bool estimateIBD;
      static bool estimateKinship;
      static bool estimateKinship15;
      static bool estimateInformation;
      static bool estimateMatrices;
      static bool selectCases;
      static bool storeKinshipForVc;
      static bool storeKinshipForAssoc;
      static bool fastAssociationAnalysis;
      static bool calcLikelihood;

      // Error detection and genotype inference flags
      static bool outputInferredGenotypes;
      static bool inferGenotypes;
      static bool findErrors;

      // Flags for controlling haplotyping
      static bool bestHaplotype;
      static int  sampledHaplotypes;
      static bool allHaplotypes;

      // Flags for controlling linkage analyses
      static bool nplAll;
      static bool nplPairs;
      static bool nplQtl;
      static bool nplDeviates;
      static bool nplExponential;

      // Flags for controlling output
      static bool simwalk2;
      static bool perFamily;
      static bool writePDF;

      // Specialized engines
      MerlinMatrix       matrix;
      MerlinKinship      kinship;
      MerlinKinship15    kinship15;
      MerlinIBD          ibd;
      MerlinHaplotype    haplo;
      NPLforSimwalk2     hybrid;
      VarianceComponents vc;
      GenotypeInference  imputationEngine;

      // PDF output file
      MerlinPDF pdf;

      // Pointer to handler for NPL information
      KongAndCox * kac;

   protected:
      // This function actually carries out most analyses
      virtual void AnalyseLocation(int pos, Tree & inheritance);

      // Error checking and genotype inference routines
      void ErrorCheck(Tree & withMarker, Tree & without, int marker);
      void InferGenotypes(Tree & withMarker, Tree & without, int marker);

      // Interface between error checking routines and MerlinCore
      virtual void GenotypeAnalysisHook(Tree & withMarker, Tree & without, int marker);
      virtual bool GenotypeAnalysisEnabled(int informativeMarker, bool global);

   private:
      // Generic task list ... new functionality should go here
      AnalysisTask * taskList;
      AnalysisInfo   taskInfo;

      // Output file listing possible errors
      void OpenErrorFile();
      FILE * errorfile;

      // Output file for listing information for each family
      void OpenPerFamilyFiles();

      // Overall likelihood for each chromosome
      double lkSum;
      int lkCount;

      // Functions for managing task list
      void NewTask(AnalysisTask * task);
   };

#endif

 
