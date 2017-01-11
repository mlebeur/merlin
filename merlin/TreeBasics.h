////////////////////////////////////////////////////////////////////// 
// merlin/TreeBasics.h 
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
 
#ifndef __TREEBASE_H__
#define __TREEBASE_H__

#include "TreeNode.h"
#include "TreeManager.h"
#include "IntArray.h"

#include <stdlib.h>
#include <stdio.h>

class TreeManager;

class BasicTree
   {
   public:
     // Maximum number of nodes to store in memory
     static int maxNodes;
     static int totalNodes;

     // Log of constant that has been used to scale all the values
     // stored in the tree -- we usually assume this is zero
     double logOffset;

     // Tree is stored as an array of nodes
     TreeNode * nodes;
     int        nextFree;

     // Depth of tree
     int        bit_count;

     // Total number of nodes allocated
     int        count;

     BasicTree()
       { Initialize(); }

     BasicTree(BasicTree & tree)
       { Initialize(); Copy(tree); }

     ~BasicTree()
       {
       if (nodes != NULL)
         {
         totalNodes -= count;
         free(nodes);
         }
       }

     void Clear()
       { nextFree = 0; if (nodes == NULL) count = 0; logOffset = 0; }

     void Free()
       { nextFree = count = 0; if (nodes != NULL) free(nodes); nodes = NULL; }

     int NewNode()
       {
       if (nextFree == count)
         Grow();
       return nextFree++;
       }

     void Grow();

     void Trim(int node = 0);
     void TrimNoMerge(int node = 0);
     void MakeBooleanTree(int node = 0);

     int  Copy(int node);
     int  Copy(const BasicTree & source, int node);
     void CarbonCopy(const BasicTree & source);
     void Copy(const BasicTree & source);

     BasicTree & operator = (BasicTree & rhs);

     void Print(int node = 0, int pad = 0);

     void MakeMinimalTree(double value, int bits);
     void MakeEmptyTree(int bits);
     bool IsEmpty()
       { return nextFree == 0 || nodes[0].type == TREE_NODE_ZERO; }

     // Routines for managing swap files
     static void SetupSwap();
     static void FreeSwap();
     static void CloseSwap(bool quiet = false);
     static double SwapFileSize() { return tmpfileInfo.GetFileSize(); }

     // Discards the contents of the current tree, to save memory
     void Discard();

     // Routines for updating thresholds for swapping
     static void UpdateSwapThreshold(int nTrees);

     // Routines for swapping trees in and out of memory
     void Pack();             // Compresses tree and frees memory
     void PackOnly();         // Compresses tree but does not free memory
     void RePack();           // Frees memory for previously packed tree
     void UnPack();           // Decompresses tree

     // Routines for swapping trees to a specific file
     void Pack(TreeManager & output);
     void PackOnly(TreeManager & output);
     void PackOrDiscard(TreeManager & output);
     void UnPack(TreeManager & input);

     // Routines for writing and reading trees to user files
     void WriteToFile(FILE * output);
     void ReadFromFile(FILE * input);

     // Routines to allocate input / output buffers
     static void AllocateBuffers();
     static void FreeBuffers();

     // Routine to check whether tree has non-zero nodes
     bool FindNonZero(int node = 0);

     // Routine to check whether tree is candidate for swapping
     bool IsSwappable()
         { return count > SWAP_MIN; }

     bool IsInMemory()
         { return nodes != NULL; }

     bool IsInSwapFile()
         { return nodes == NULL && count > 0; }

     static int LargeTreeSize()
         { return SWAP_MIN; }

   protected:
     // Used for avoiding recursion
     IntArray stack;

   private:
     void Initialize()
       {
       bit_count = nextFree = count = 0;
       nodes = NULL;
       logOffset = 0.0;
       }

     // Routines for managing memory allocation failures
     void MemoryCeiling();
     void OutOfMemory();

     // Information used to swap tree in and out of memory
     int           leafCount;
     TreePosition  tmpfileOffset;

     // Global file output buffers
     //
     static char        * skeleton;
     static double      * leaves;
     static MiniDeflate   zip;

     static TreeManager  tmpfileInfo;

     // Minimum size of inheritance trees that are swapped out of memory
     //
     static int          SWAP_MIN;
   };

class TreesTooBig
   {
   public:
     int memory_request;

     TreesTooBig(const TreesTooBig & source)
       { memory_request = source.memory_request; }
     TreesTooBig(int request)
       { memory_request = request; }

   };

///////////////////////////////////////////////////////////////
// The SAFE_SET macro is crucial for setting node values
// when MERLIN is compiled by gcc-3.0 and above. SAFE_SET
// ensures that the left side (variable) of an assignment
// expression is  evaluated after the right side (value)
//
// Without it, the following piece of code can raise
// segmentation faults whenever NewNode() changes the
// nodes pointer!
//
//   nodes[i].child[0] = NewNode();
//
// Note that NewNode() and many other tree functions can
// change the value of the nodes pointer!
//

#define SAFE_SET(a,b) {int tmp = b; a = tmp; }

#endif


 
