////////////////////////////////////////////////////////////////////// 
// merlin/Tree.h 
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
 
#ifndef __TREE_H__
#define __TREE_H__

#include "TreeBasics.h"
#include "IntArray.h"

class Tree : public BasicTree
   {
   public:
      void Graft(int target, int source);

      void Floor(int node, double floor);
      void Add(int node, double value);
      void Add(Tree & tree, int node1 = 0, int node2 = 0);
      void Multiply(int node, double scale);
      void Multiply(Tree & tree, int node1 = 0, int node2 = 0);

      double MeanProduct(Tree & t, int node1=0, int node2=0, double weight=1.0);
      double Mean(int node1, double weight=1.0) const;
      double WeightedSum() const { return Mean(0, 1.0); }

      void ScaleSum(int node, double scale, double sum);

      int  CopyAndScale(int node, double scale);
      int  CopyAndScale(Tree & tree, int node, double scale);
      int  CopyAndAdd(Tree & tree, int node, double summand);
      void GraftAndScale(int target, int source, double scale);
      void GraftAndScale(int target, Tree & tree, int source, double scale);

      int  CopyAndScaleAndAdd(int node,
           double scale, double sum_target, double sum_source);
      void GraftAndScaleAndAdd(int target, int source,
           double scale, double sum_target, double sum_source);

      int  CopyAndScaleAndAdd(Tree & source, int node, double mul, double sum);
      void GraftAndScaleAndAdd(int target, Tree & source, int node,
           double mul, double sum);

      int  CopyAndScaleAndChoose(Tree & source, int node, double mul, double min);
      void GraftAndScaleAndChoose(int target, Tree & source, int node,
           double mul, double min);

      int  CopyAndScaleAndChoose(int node,
           double scale, double lo_target, double lo_source);
      void GraftAndScaleAndChoose(int target, int source,
           double scale, double lo_target, double lo_source);

      Tree & operator *= (Tree & rhs)
         { Multiply(rhs, 0, 0); return *this; }

      void MakeTree(const int * bits, int depth, double value);
   };

#endif
 
