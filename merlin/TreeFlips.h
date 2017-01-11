////////////////////////////////////////////////////////////////////// 
// merlin/TreeFlips.h 
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
 
#ifndef __TREEFLIPS_H__
#define __TREEFLIPS_H__

#include "Magic.h"

class TreeFlips : public InheritanceTree
   {
   public:
      // Routines for using standard flips reflecting symmetry of founder
      // alleles

      int  FlipCopyAndScale(const int * flips, int node, double scale);
      void FlipGraftAndScale(const int * flips, int target, int source,
                             double scale);
      int  FlipCopyAndScale(const int * flips, Tree & tree, int node, double scale);
      void FlipGraftAndScale(const int * flips, int target, Tree & tree,
                             int source, double scale);

      int  FlipCopyAndScaleAndAdd(const int * flips, int node,
           double scale, double sum_target, double sum_source);
      void FlipGraftAndScaleAndAdd(const int * flips, int target, int source,
           double scale, double sum_target, double sum_source);
      int  FlipCopyAndScaleAndAdd(const int * flips, Tree & tree, int node,
           double scale, double sum);
      void FlipGraftAndScaleAndAdd(const int * flips, int target,
           Tree & tree, int source, double scale, double sum);

      int  FlipCopyAndScaleAndChoose(const int * flips, int node,
           double scale, double sum_target, double sum_source);
      void FlipGraftAndScaleAndChoose(const int * flips, int target, int source,
           double scale, double sum_target, double sum_source);
      int  FlipCopyAndScaleAndChoose(const int * flips, Tree & tree, int node,
           double scale, double sum);
      void FlipGraftAndScaleAndChoose(const int * flips, int target,
           Tree & tree, int source, double scale, double sum);

      // Routines for handling double flips used in founder
      // couple reduction

      void DoubleFlipBranch(const int * flips, int node);
      void UpgradeNode(int node, int new_type);      
      void PseudoUpgradeNode(int node, int new_type);      


   };

#endif

 
