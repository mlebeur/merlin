////////////////////////////////////////////////////////////////////// 
// clusters/HaploSet.h 
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
 
#ifndef __HAPLOSET_H__
#define __HAPLOSET_H__

#include "HaploGraph.h"
#include "HaploFamily.h"

class HaplotypeSet
   {
   public:
      String            setid;
      HaplotypeGraph ** graphs;
      int               graphCount;
      int               haploCount;
      double            weight;
      int               scale;
      int               copies;

      HaplotypeSet();
      ~HaplotypeSet();

      void LoadFromMemory(FounderGraph * input, int markers);
      void LoadFromMemory(FamilyHaplos * input, int markers);

      void StandardizeWeights();

      double Complexity(IntArray & alleleCounts);

      static int  maxHaploCount;

   private:
      int  size;
      void Grow();
   };

class HaplotypeSets : public StringMap
   {
   public:
      HaplotypeSet * FindSet(const ::String & tag)
         {
         int set = Find(tag, create_set);
         return (HaplotypeSet *) objects[set];
         }

      void LoadFromFile(FILE * file);
      HaplotypeSet * LoadFromMemory(FamilyHaplos * input, int markers);

      void ClearContents();

      ~HaplotypeSets();
      HaplotypeSets();

      void Filter(IntArray & alleleCounts, int maxBits);
      void ShowCounts(IntArray & alleleCounts);

      void AllocateAlleleLabels(int markers);
      void SetAlleleLabels(int marker, StringArray & labels);
      ::String GetAlleleLabel(int marker, int allele);

      ::StringArray * alleleLabels;

       bool quiet;

      StringArray discarded;

   private:
      static void * create_set();
   };

#endif
 
