////////////////////////////////////////////////////////////////////// 
// merlin/QtlModel.h 
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
 
#ifndef __QTLMODEL_H__
#define __QTLMODEL_H__

#include "StringArray.h"
#include "MathVector.h"
#include "IntArray.h"

class Pedigree;

class QtlModel
   {
   public:
      int modelCount;
      IntArray traitList;
      IntArray * covariateList;
      IntArray * candidateList;

      QtlModel();
      virtual ~QtlModel();

      void LoadFromFile(const char * filename);
      void LoadFromFile(IFILE & f);

      bool HaveModels()
         { return modelCount; }

      void PrintSummary();

      bool SkipMarker(int model, int marker);

   private:
   };

class RefinedQtlModel
   {
   public:
      int modelCount;
      StringArray labels;
      StringArray comments;
      Vector      scores;
      Vector *    values;

      RefinedQtlModel();
      virtual ~RefinedQtlModel();

      void Dimension(int count);

      void WriteFiles(const String & prefix, Pedigree & ped, QtlModel & base);

      void PrintDataFile(FILE * f);
      void PrintPedigreeFile(FILE * f, Pedigree & ped);
      void PrintModelFile(FILE * f, QtlModel & base);

      void PrintDataFile(const char * filename);
      void PrintPedigreeFile(const char * filename, Pedigree & ped);
      void PrintModelFile(const char * filename, QtlModel & base);

      bool CheckForImprovement(int model, double score);

   private:
   };

#endif


 
