////////////////////////////////////////////////////////////////////// 
// regress/RegressAnalysis.h 
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
 
#ifndef __REGRESSANALYSIS_H__
#define __REGRESSANALYSIS_H__

#include "MerlinCore.h"
#include "FancyRegression.h"
#include "MerlinPDF.h"

class RegressionAnalysis : public MerlinCore
   {
   public:
      // Constructor and destructor
      RegressionAnalysis(Pedigree & ped);
      ~RegressionAnalysis();

      // Trait model parameters
      static double mean;
      static double variance;
      static double heritability;
      static double testRetestCorrel;

      // Trait modelling options
      static String modelsFile;
      static bool   autoFit;

      static bool unrestricted;

      static bool rankFamilies;

      // Optional flag for outputing PDF files
      static bool writePDF;

      // Setup global variables
      virtual void SetupGlobals();

      // Setup the current map
      virtual int SetupMap(int chromosome);

      // Setup analysis for the current family
      virtual bool SelectFamily(Family * f, bool warnOnSkip = true);

      // Carry out the core analyses
      virtual void AnalyseLocation(int location, Tree & inheritance);

      // Analyse the currently selected family
      void Analyse();

      // Print LOD scores at each location
      void PrintScores();

      // Routines for setting up trait models
      void AutoFitModels();
      bool ReadModelsFromFile();
      void SetupDefaultModels();

      // Manages pdf output
      MerlinPDF pdf;

   private:
      void RankFamilies();

      int             modelCount;
      FancyRegression * regress;
      RegressKinship  kinship;
      Tree            tree;
      String          label;
   };

#endif


 
