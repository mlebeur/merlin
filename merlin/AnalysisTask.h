////////////////////////////////////////////////////////////////////// 
// merlin/AnalysisTask.h 
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
 
#ifndef __ANALYSISTASK_H__
#define __ANALYSISTASK_H__

#include "Tree.h"
#include "Mantra.h"

class MerlinPDF;

class AnalysisInfo
   {
   public:
      // Generic information about task
      int chromosome;
      int positions;
      int families;

      // Likelihood for current set of inheritance vectors
      double lk;

      // Serial number for the family currently being analysed
      int famno;

      // Label for the family currently being analyzed
      String * famid;

      // Other information
      Mantra * m;
      StringArray * labels;
      Vector * analysisPositions;

      // PDF information
      bool drawPDF;
      MerlinPDF * pdf;

      // Information on requested analyses
      bool perFamily;
   };

class AnalysisTask
   {
   public:
      AnalysisTask();
      virtual ~AnalysisTask();

      virtual void Setup(AnalysisInfo & info);
      virtual void SetupFamily(AnalysisInfo & info);
      virtual void Analyse(AnalysisInfo & info, Tree & tree, int position);
      virtual void AnalyseUninformative(AnalysisInfo & info, Tree & dummy);
      virtual void ReportFamily();
      virtual void Report(AnalysisInfo & info);
      virtual void SkipFamily(AnalysisInfo & info);
      virtual void * ProcessMessage(const char * message);

      virtual void OpenFiles(AnalysisInfo & info, String & prefix);
      virtual void CloseFiles();

      virtual const char * TaskDescription();

      AnalysisTask * next;

      AnalysisTask(AnalysisTask * nextTask)
         { nextTask = next; }
   };

#endif
 
