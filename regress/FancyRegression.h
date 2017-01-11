////////////////////////////////////////////////////////////////////// 
// regress/FancyRegression.h 
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
 
#ifndef __FANCYREGRESSION_H__
#define __FANCYREGRESSION_H__

#include "RegressKinship.h"
#include "MathCholesky.h"
#include "MerlinPDF.h"

class RegressionAnalysis;

#define COVAR_SEX       1
#define COVAR_USER      2
#define COVAR_TRAIT     3
#define COVAR_AFFECTION 4

class FancyRegression
   {
   public:
      String label, shortLabel;

      int    trait;
      double heritability;
      double mean;
      double variance;
      double testRetestCorrel;

      double maxInfo;
      double sumMaxInfo;

      FancyRegression();

      void  SetPositionCount(int count);
      void  SetupFamily(RegressKinship & kin);
      void  Analyse(RegressKinship & kin, int position);
      void  PrintScores(FILE * tablefile, int chromosome, MerlinPDF & pdf, StringArray & labels);

      void  ResetMaxInfo() { sumMaxInfo = 0; }
      int   CountPairs()   { return pairs1.Length(); }

      friend class RegressionAnalysis;

      bool  bounded;

   private:
      Cholesky chol;

      Vector scores, information;
      Vector sum, diff, y;
      Vector weights, pi, local_pi, E, intermediate;
      Matrix sigma_ss, sigma_dd, sigma_sd;
      Matrix sigma_ds, sigma_yy, B, H, temp;
      Matrix sigma_pi, local_sigma_pi;
      IntArray repeatCounts;

      IntArray pairs1, pairs2;

      double square(double x) { return x * x; }

      double calc_correl(RegressKinship & kin, int person1, int person2);

      // Used for tests of heterogeneity based on a binary covariate
      bool    covariate;
      int     covarType, covarId;
      double  covarThreshold;

      // Helper functions related to covariates
      bool    CovariateObserved(Person & p);
      int     Classify(Person & p);
   };

#endif

 
