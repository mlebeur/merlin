////////////////////////////////////////////////////////////////////// 
// merlin/TreeInfo.h 
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
 
#ifndef __TREEINFO_H__
#define __TREEINFO_H__

#include "TreeBasics.h"
#include "MathVector.h"
#include "IntArray.h"

class TreeInfo
   {
   public:
      int    nodeCount;
      int    zeroCount;
      int    leafCount;

      double min;
      double max;
      double mean;
      double var;
      double information;

      TreeInfo();

      double GetMean(BasicTree & tree, int node = 0);
      void   GetMeanVar(BasicTree & tree, int node = 0);
      void   GetBounds(BasicTree & tree, int node = 0);

      void   GetTreeInfo(BasicTree & tree, int node = 0);
      void   GetDetailedInfo(BasicTree & tree, int node = 0);
      double GetMinInfo(BasicTree & tree, int node = 0);

      // Useful for avoiding recursion
      IntArray stack;
      Vector   weights;

   private:
      double  sumInfo;
   };

#endif

 
