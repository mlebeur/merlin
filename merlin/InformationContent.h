////////////////////////////////////////////////////////////////////// 
// merlin/InformationContent.h 
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
 
#ifndef __INFORMATIONCONTENT_H__
#define __INFORMATIONCONTENT_H__

#include "AnalysisTask.h"
#include "MathVector.h"
#include "TreeInfo.h"

class InformationContent : public AnalysisTask
   {
   public:
      InformationContent();
      virtual ~InformationContent();

      virtual void Setup(AnalysisInfo & info);
      virtual void Analyse(AnalysisInfo & info, Tree & tree, int position);
      virtual void AnalyseUninformative(AnalysisInfo & info, Tree & tree);
      virtual void Report(AnalysisInfo & info);

      virtual void OpenFiles(AnalysisInfo & info, String & prefix);
      virtual void CloseFiles();

   private:
      Vector   infoBits, totalBits;
      TreeInfo stats;
      FILE *   infoFile, * infoTable;
      String   filename, tablename;
   };

#endif
  
