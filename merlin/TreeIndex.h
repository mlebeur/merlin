////////////////////////////////////////////////////////////////////// 
// merlin/TreeIndex.h 
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
 
#ifndef __TREEINDEX_H__
#define __TREEINDEX_H__

#include "Tree.h"
#include "MathVector.h"

class IndexTree : public Tree
   {
   private:
      void UpdateFrequencies(const Tree & freqs, int n1, int n2, double w);
      void UpdateFrequencies(int node, double weigth);

      void EvenFrequencies(int node, double weigth);

   public:
      void MakeIndex(const Tree & source);
      void UpdateFrequencies(const Tree & freqs, double scale = 1.0);
      void EvenFrequencies();

      Vector uniqValues;
      Vector freqs;
   };

#endif
 
