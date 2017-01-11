////////////////////////////////////////////////////////////////////// 
// merlin/KongAndCox.h 
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
 
#ifndef __KONG_AND_COX_H__
#define __KONG_AND_COX_H__

#include "AnalysisTask.h"
#include "TreeIndex.h"
#include "TreeInfo.h"
#include "MathGold.h"

class LinearModel : public ScalarMinimizer
   {
   public:
      // The first dimension in these arrays corresponds to
      // a particular scoring statistic, the second dimension
      // corresponds to families.
      Vector *  min;
      Vector *  max;
      Vector ** scores;

      int statistics, families, positions;
      int stat, pos;

      double minDelta, maxDelta, scale;
      int    informativeFamilies;

      LinearModel();
      ~LinearModel();

      void Dimension(int statistics, int families, int positions);

      double Evaluate(double x) { return f(x); }
      virtual double f(double x);

      void UpdateDeltaBounds();
      void CalculateScore(double & Z, double & delta, double & chisq);
      void CalculateMinScore(double & Z, double & delta, double & chisq);
      void CalculateMaxScore(double & Z, double & delta, double & chisq);

   private:
      void FreeMemory();
      void AllocateMemory();
   };

class ExponentialModel : public ScalarMinimizer
   {
   public:
      Vector ** scores;
      Vector ** probs;

      int statistics, families, positions;
      int stat, pos;

      double minDelta, maxDelta;

      ExponentialModel();
      ~ExponentialModel();

      void Dimension(int statistics, int families, int positions);

      void NullDistribution(int statistic, int family, Vector & z, Vector & p);
      void UpdateDistribution(int statistic, int family, int position, Vector & p);

      void UpdateDeltaBounds();
      void CalculateScore(double & delta, double & chisq);
      void CalculateMinScore(double & delta, double & chisq);
      void CalculateMaxScore(double & delta, double & chisq);

      double Evaluate(double x) { return f(x); }
      virtual double f(double x);

   private:
      void FreeMemory();
      void AllocateMemory();
   };

class KongAndCox : public AnalysisTask
   {
   public:
      LinearModel        lm;
      ExponentialModel   em;
      IndexTree *        index;
      Tree *             npl;
      TreeInfo           stats;

      StringArray pheno, phenoTag;

      int    statistics, families, positions;
      Vector means, scales;

      static bool nplAll;
      static bool nplPairs;
      static bool nplQtl;
      static bool nplDeviates;
      static bool nplExtras;
      static bool nplExponential;
      static bool rawOutput;

      KongAndCox(Pedigree & ped);
      ~KongAndCox();

      static bool CalculateNPL(Pedigree & ped);

      virtual void Setup(AnalysisInfo & info);
      virtual void SetupFamily(AnalysisInfo & info);
      virtual void SkipFamily(AnalysisInfo & info);
      virtual void Analyse(AnalysisInfo & info, Tree & tree, int position);
      virtual void Report(AnalysisInfo & info);

      virtual void OpenFiles(AnalysisInfo & info, String & prefix);
      virtual void CloseFiles();

      virtual const char * TaskDescription();

   private:
      void AllocateMemory();
      void LabelAnalyses(Pedigree & ped);
      void ScoreNPLBounds(int statistics, int family);

      void OutputLabel(int chromosome, double position,
                       const char * what, const char * poslabel, double Z);
      void OutputLOD(double delta, double chisq);

      void OutputPerFamily(Pedigree & ped, const char * label, double delta);
      void OutputRawScores(Pedigree & ped, int chromosome, StringArray & labels);

      String filename, tablename, rawOutputFilename;
      FILE * file, *tablefile, * rawOutputFile;
   };

#endif
 
