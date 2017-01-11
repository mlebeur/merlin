////////////////////////////////////////////////////////////////////// 
// merlin/MerlinModel.h 
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
 
#ifndef __MERLINMODEL_H__
#define __MERLINMODEL_H__

#include "ParametricLikelihood.h"
#include "DiseaseModel.h"
#include "AnalysisTask.h"
#include "MathGold.h"

// This class is used for managing parametric trait models
//

class MerlinModel
   {
   public:
      MerlinModel();
      ~MerlinModel();

      // Load parametric models from file
      void LoadModels();

      // Return number of parametric models to be considered
      int  CountModels();

      // Static variables
      static String modelsFile;

      // Creates a task object capable of carrying out parametric analyses
      AnalysisTask * CreateTask();

   private:
      StringArray    labels;
      DisModel * models;
      int modelCount;
   };

class ParametricAnalysis : public AnalysisTask, private ScalarMinimizer
   {
   public:
      ParametricAnalysis(DisModel * models, StringArray & labels);
      virtual ~ParametricAnalysis();

      virtual void Setup(AnalysisInfo & info);
      virtual void SetupFamily(AnalysisInfo & info);
      virtual void Analyse(AnalysisInfo & info, Tree & tree, int position);
      virtual void AnalyseUninformative(AnalysisInfo & info, Tree & tree);
      virtual void ReportFamily();
      virtual void Report(AnalysisInfo & info);

      virtual void OpenFiles(AnalysisInfo & info, String & prefix);
      virtual void CloseFiles();

      virtual const char * TaskDescription();

      // Estimates proportion of admixed families to maximize LOD score
      void   EstimateAlpha();

   private:
      // Disease model information
      int            modelCount;
      DisModel     * models;
      StringArray    labels;

      // Analysis results
      IntArray * impossible;
      Vector * scores;
      Matrix * rawScores;
      Vector baseline;

      // Tree with model information
      Tree * likelihoods;

      // File for storing results for individual families
      FILE * scoreFile, * scoreTable;
      String filename, tablename;
      static double scale;

      // Calculate likelihood ratio at a specific location and assuming
      // a specific admixture proportion alpha
      virtual double f(double alpha);

      int    pos, model;
      double alpha, hlod;
   };

#endif
 
