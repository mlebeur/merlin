////////////////////////////////////////////////////////////////////// 
// merlin/AnalysisTask.cpp 
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
 
#include "AnalysisTask.h"

AnalysisTask::AnalysisTask()
   {
   }

AnalysisTask::~AnalysisTask()
   {
   }

void AnalysisTask::Setup(AnalysisInfo &)
   {
   }

void AnalysisTask::SetupFamily(AnalysisInfo &)
   {
   }

void AnalysisTask::Analyse(AnalysisInfo &, Tree &, int)
   {
   }

void AnalysisTask::AnalyseUninformative(AnalysisInfo &, Tree &)
   {
   }

void AnalysisTask::ReportFamily()
   {
   }

void AnalysisTask::Report(AnalysisInfo & info)
   {
   }

void AnalysisTask::OpenFiles(AnalysisInfo & info, String & prefix)
   {
   }

void AnalysisTask::CloseFiles()
   {
   }

void AnalysisTask::SkipFamily(AnalysisInfo & info)
   {
   }

const char * AnalysisTask::TaskDescription()
   {
   return NULL;
   }

void * AnalysisTask::ProcessMessage(const char * message)
   {
   if (next == NULL)
      return NULL;

   return next->ProcessMessage(message);
   }
 
