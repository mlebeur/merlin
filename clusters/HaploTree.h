////////////////////////////////////////////////////////////////////// 
// clusters/HaploTree.h 
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
 
#ifndef __HAPLOTREE_H__
#define __HAPLOTREE_H__

#include "IntArray.h"
#include "StringArray.h"
#include "MathVector.h"

class HaploTree
   {
   public:
      int size;
      int count;
      int depth;

      IntArray **branches;
      IntArray   entries;
      IntArray   levels;

      Vector new_freqs;
      Vector freqs;

      HaploTree();
      ~HaploTree();

      void Clear();

      int  AddHaplotype(const IntArray & state);
      void MakeShrub(int alleles);
      void MakeEmptyTree(const IntArray & alleleCounts);

      int  GetBranch(int branch, int allele);
      int  PeekBranch(int branch, int allele) const;

      void SetupTraversal(IntArray & pointer, IntArray & state) const;
      int  Traverse(IntArray & pointer, IntArray & state) const;
      int  Traverse(IntArray & pointer, IntArray & state, double minfreq) const;

      void Copy(const HaploTree & source);
      void Merge(const HaploTree & top, const HaploTree & bottom);

      void Print(double minfreq) const;
      void Print(StringArray * labels, double minfreq) const;

      void SetupFrequencies();
      void SetAlleleCounts(const IntArray & counts)
         { alleleCounts = counts; }

   private:
      void Grow(int minsize);
      int  AllocateBranch(int level);
      int  AllocateFrequency();

      int  leafs;

      IntArray   alleleCounts;
   };

#endif

 
