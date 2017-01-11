////////////////////////////////////////////////////////////////////// 
// merlin/MerlinKinship.h 
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
 
#ifndef __MERLIN_KINSHIP_H__
#define __MERLIN_KINSHIP_H__

#include "MathMatrix.h"
#include "Pedigree.h"
#include "Mantra.h"
#include "Tree.h"
#include "MerlinCore.h"
#include "AnalysisTask.h"

#include <stdio.h>

class MerlinKinship : public AnalysisTask
   {
   public:
      MerlinKinship(Mantra & m);
      ~MerlinKinship();

      virtual void Analyse(AnalysisInfo & info, Tree & tree, int position);
      virtual void AnalyseUninformative(AnalysisInfo & info, Tree & dummy);
      double Score(Tree & tree, int node, int bit, int start = 0);

      virtual void OpenFiles(AnalysisInfo & info, String & prefix);
      virtual void CloseFiles();

      void   Store(int position, double scale);
      void   Output(const char * label, double scale);
      void   SelectCases(const char * label, double scale);

      virtual void Setup(AnalysisInfo & info);
      virtual void SetupFamily(AnalysisInfo & info);

      bool   CheckFamily(int serial);
      double Retrieve(int family, int position, int index1, int index2);
      double RetrieveFromScratch(int i, int j, double scale);

      // Discards kinship information, for families
      // were analysis could not be completed
      virtual void SkipFamily(AnalysisInfo & info);

      bool   write, store, selectCases;

   private:
      Pedigree * ped;
      Family   * family;
      Mantra   & mantra;

      // Routines for scoring kinship in trees with symmetries
      double   AlleleKinship(int allele1, int allele2);
      double   WeightedScore(Tree & tree, int node, int bit, int start,
                             double weight = 2.0);

      // Output files
      FILE     * kinfile;
      FILE     * casefile;

      // Output file names
      String   kinfilename;
      String   casefilename;

      // Temporary storage
      Matrix   kinship;

      // Track uninformative meiosis in the pedigree
      IntArray uninformative;

      // Store kinship matrices for all families
      // prior to variance components analysis
      Vector * matrices;

      // Array for storing sharing scores [one per individual]
      Vector   scores;

      // Array for storing expect sum of sharing scores [one per trait]
      Vector   prior;

      // Scale constant for converting sum of kinships to NPL scores
      Vector   priorToZ;
   };

#endif


 
