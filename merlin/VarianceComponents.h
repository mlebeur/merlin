////////////////////////////////////////////////////////////////////// 
// merlin/VarianceComponents.h 
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
 
#ifndef __VARIANCECOMPONENTS_H__
#define __VARIANCECOMPONENTS_H__

#include "MathNormal.h"
#include "Pedigree.h"
#include "QtlModel.h"

class MerlinPDF;
class FamilyAnalysis;
class NormalEquations;

class VarianceComponents
   {
   public:
      VarianceComponents();

      void Analyse(FamilyAnalysis & engine);

      static bool   useCovariates;
      static bool   useProbands;
      static double unlinkedFraction;

      static QtlModel customModels;

      void OpenPerFamilyFile();
      void ClosePerFamilyFile();

   protected:
      void ListCovariates();
      void PrintCovariates();

      void PrintPvalue(double pvalue);

      bool CheckCovariates(Person & p);
      bool CheckRepeatedMeasures(Pedigree & ped, IntArray & traitIndex, IntArray * pheno);
      bool CheckAveragedMeasures(Pedigree & ped, IntArray * pheno, int repeatedMeasures);

      void SetupPDF(MerlinPDF & pdf, String & trait, int probandCount);

      void DuplicateModel(NormalEquations & dest, NormalEquations & src,
                          int vc_count, int beta_count);
      double TotalVariance(Vector & variances, int vc_count);

      FILE * perFamily;

      void WritePerFamilyLOD(Pedigree & ped, IntArray * pheno,
           const char * position, Vector & nullLLK, Vector & fullLLK);

      IntArray covariates;

      // Utility functions for dealing with repeated measurements,
      // when they have not been averaged
      bool   SkipTrait(int trait);
      bool   IsPhenotyped(Person & person, IntArray & measurements, int trait);
      bool   IndexMeasurements(IntArray & index, String & label, int trait);
      double RetrieveMeasurement(Person & person, IntArray & measurements, int & measurement);
      bool   ExtractMeasurementNumber(const String & label, int & prefix_len, int & number);
      int    FindMeasurement(const String & label, int measurement);
      void   SummarizeRepeatedMeasures(Pedigree & ped, String & label, IntArray * pheno);
   };

#endif

 
