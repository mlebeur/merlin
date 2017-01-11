////////////////////////////////////////////////////////////////////// 
// merlin/TreeNode.h 
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
 
#ifndef __TREENODE_H__
#define __TREENODE_H__

#define TREE_NODE_ZERO      0
#define TREE_NODE_LEAF      1
#define TREE_NODE_ONE       2
#define TREE_NODE_TWO       3
#define TREE_NODE_INT_VALUE 4

#define PAIR_NODES(a, b)   (((a) << 2) + (b))

struct TreeNode
   {
   int type;
   union
      {
      double value;
      int    child[2];
      int    integer;
      void * info;
      };
   };

#endif


 
