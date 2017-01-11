////////////////////////////////////////////////////////////////////// 
// merlin/TreeIndex.cpp 
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
 
#include "TreeIndex.h"
#include "IntArray.h"

void IndexTree::MakeIndex(const Tree & source)
   {
   Copy(source);

   int leafCount = 0;

   // Count number of leaf nodes
   for (int i = 0; i < nextFree; i++)
     if (nodes[i].type == TREE_NODE_LEAF || nodes[i].type == TREE_NODE_ZERO)
       leafCount++;

   // List all leaf node values in values array
   uniqValues.Dimension(leafCount);

   for (int i = 0, next = 0; i < nextFree; i++)
     if (nodes[i].type == TREE_NODE_LEAF || nodes[i].type == TREE_NODE_ZERO)
       uniqValues[next++] = nodes[i].value;

   // Sort values and delete duplicates
   uniqValues.Sort();
   uniqValues.RemoveDuplicates();

   // Recode leaf values to ranks
   for (int i = 0; i < nextFree; i++)
     if (nodes[i].type == TREE_NODE_LEAF || nodes[i].type == TREE_NODE_ZERO)
       {
       nodes[i].type = TREE_NODE_INT_VALUE;
       nodes[i].integer = uniqValues.BinarySearch(nodes[i].value);
       }
   }

void IndexTree::UpdateFrequencies(const Tree & source, double scale)
   {
   freqs.Dimension(uniqValues.Length());
   freqs.Zero();

   UpdateFrequencies(source, 0, 0, scale);
   }

void IndexTree::UpdateFrequencies(int node, double weight)
   {
   switch (nodes[node].type)
    {
    case TREE_NODE_INT_VALUE :
      freqs[nodes[node].integer] += weight;
      break;
    case TREE_NODE_ONE :
      UpdateFrequencies(nodes[node].child[0], weight);
      break;
    case TREE_NODE_TWO :
      UpdateFrequencies(nodes[node].child[0], weight * 0.5);
      UpdateFrequencies(nodes[node].child[1], weight * 0.5);
      break;
    }
   }

void IndexTree::UpdateFrequencies(const Tree & tree, int node1, int node2, double weight)
   {
   switch (PAIR_NODES(nodes[node1].type, tree.nodes[node2].type))
     {
     case PAIR_NODES(TREE_NODE_INT_VALUE, TREE_NODE_ZERO) :
     case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ZERO)  :
     case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ZERO)  :
       return;
     case PAIR_NODES(TREE_NODE_INT_VALUE, TREE_NODE_LEAF) :
       freqs[nodes[node1].integer] += tree.nodes[node2].value * weight;
       return;
     case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_LEAF) :
     case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_LEAF) :
       UpdateFrequencies(node1, tree.nodes[node2].value * weight);
       return;
     case PAIR_NODES(TREE_NODE_INT_VALUE, TREE_NODE_ONE) :
     case PAIR_NODES(TREE_NODE_INT_VALUE, TREE_NODE_TWO) :
       freqs[nodes[node1].integer] += tree.Mean(node2, weight);
       return;
     case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
       UpdateFrequencies(tree, nodes[node1].child[0],
                     tree.nodes[node2].child[0], weight);
       return;
     case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_TWO) :
       UpdateFrequencies(tree, nodes[node1].child[0],
                     tree.nodes[node2].child[0], weight * 0.5);
       UpdateFrequencies(tree, nodes[node1].child[0],
                     tree.nodes[node2].child[1], weight * 0.5);
       return;
     case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ONE) :
       UpdateFrequencies(tree, nodes[node1].child[0],
                           tree.nodes[node2].child[0], weight * 0.5);
         UpdateFrequencies(tree, nodes[node1].child[1],
                           tree.nodes[node2].child[0], weight * 0.5);
         return;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_TWO) :
         UpdateFrequencies(tree, nodes[node1].child[0],
                           tree.nodes[node2].child[0], weight * 0.5);
         UpdateFrequencies(tree, nodes[node1].child[1],
                           tree.nodes[node2].child[1], weight * 0.5);
         return;
      }
   }

void IndexTree::EvenFrequencies()
   {
   freqs.Dimension(uniqValues.Length());
   freqs.Zero();

   EvenFrequencies(0, 1.0);
   }

void IndexTree::EvenFrequencies(int node, double weight)
   {
   switch (nodes[node].type)
    {
    case TREE_NODE_INT_VALUE :
      freqs[nodes[node].integer] += weight;
      break;
    case TREE_NODE_ONE :
      EvenFrequencies(nodes[node].child[0], weight);
      break;
    case TREE_NODE_TWO :
      EvenFrequencies(nodes[node].child[0], weight * 0.5);
      EvenFrequencies(nodes[node].child[1], weight * 0.5);
      break;
    }
   }



 
