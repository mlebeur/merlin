////////////////////////////////////////////////////////////////////// 
// merlin/NPL-ASP.h 
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
 
#ifndef __NPL_ASP_H__
#define __NPL_ASP_H__

#include "Pedigree.h"
#include "Mantra.h"
#include "Tree.h"

class NPL_ASP_Tree : public Tree
   {
   public:
      void    ScoreNPL(Mantra & f, int affection);

      // Dummy virtual destructor
      virtual ~NPL_ASP_Tree() { }

   private:
      int    RecursivelyScoreNPL(Mantra & m, int bit);
      bool   isRedundantNPL(Mantra & m, int pivot);
      virtual double CalculateNPL(Mantra & m) = 0;
      virtual void   PrepareForScoring(Mantra & m);

   };

class NPL_ALL : public NPL_ASP_Tree
   {
   public:
      // Dummy virtual destructor
      virtual ~NPL_ALL() { }

   protected:
      virtual double CalculateNPL(Mantra & m);
      virtual void   PrepareForScoring(Mantra & m);

      static const double factors[20];
      static const double reciprocals[20];

      IntArray counts;
      int      permutations;
      double   permutationWeight;

   private:
      double   ScoreSubset(Mantra & m, IntArray & affecteds);
      IntArray subset, graph;
   };

class NPL_Pairs : public NPL_ASP_Tree
   {
   public:
      // Dummy virtual destructor
      virtual ~NPL_Pairs() { }

   protected:
      virtual double CalculateNPL(Mantra & m);
   };

class NPL_Pairs_For_Simwalk2 : public NPL_ASP_Tree
   {
   public:
      // Dummy virtual destructor
      virtual ~NPL_Pairs_For_Simwalk2() { }

   protected:
      virtual double CalculateNPL(Mantra & m);
   };

class NPL_Pairs_Maternal : public NPL_ASP_Tree
   {
   public:
      // Dummy virtual destructor
      virtual ~NPL_Pairs_Maternal() { }

   protected:
      virtual double CalculateNPL(Mantra & m);
   };

class NPL_Pairs_Paternal : public NPL_ASP_Tree
   {
   public:
      // Dummy virtual destructor
      virtual ~NPL_Pairs_Paternal() { }

   protected:
      virtual double CalculateNPL(Mantra & m);
   };

#endif 
