////////////////////////////////////////////////////////////////////// 
// merlin/TreeFlips.cpp 
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
 
#include "TreeFlips.h"
#include "MathConstant.h"

int TreeFlips::FlipCopyAndScale(const int * flips, int node, double scale)
   {
   int new_node = NewNode();

   switch (nodes[new_node].type = nodes[node].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         nodes[new_node].value  = nodes[node].value * scale;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0],
            FlipCopyAndScale(flips + 1, nodes[node].child[*flips], scale));
         SAFE_SET(nodes[new_node].child[1],
            FlipCopyAndScale(flips + 1, nodes[node].child[!*flips], scale));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0],
            FlipCopyAndScale(flips + 1, nodes[node].child[0], scale));
         break;
      }
   return new_node;
   }

int TreeFlips::FlipCopyAndScale(const int * flips, Tree & tree, int node, double scale)
   {
   int new_node = NewNode();

   switch (nodes[new_node].type = tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         nodes[new_node].value  = tree.nodes[node].value * scale;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0],
            FlipCopyAndScale(flips + 1, tree, tree.nodes[node].child[*flips], scale));
         SAFE_SET(nodes[new_node].child[1],
            FlipCopyAndScale(flips + 1, tree, tree.nodes[node].child[!*flips], scale));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0],
            FlipCopyAndScale(flips + 1, tree, tree.nodes[node].child[0], scale));
         break;
      }
   return new_node;
   }

int TreeFlips::FlipCopyAndScaleAndAdd(const int * flips, int node,
    double scale, double sum_target, double sum_source)
   {
   int new_node = NewNode();

   switch (nodes[new_node].type = nodes[node].type)
      {
      case TREE_NODE_ZERO :
         nodes[new_node].type = TREE_NODE_LEAF;
         nodes[new_node].value = sum_target;
         nodes[node].type = TREE_NODE_LEAF;
         nodes[node].value = sum_source;
         break;
      case TREE_NODE_LEAF :
         nodes[new_node].value  = sum_target + nodes[node].value * scale;
         nodes[node].value     += sum_source;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0],
            FlipCopyAndScaleAndAdd(flips + 1, nodes[node].child[*flips], scale, sum_target, sum_source));
         SAFE_SET(nodes[new_node].child[1],
            FlipCopyAndScaleAndAdd(flips + 1, nodes[node].child[!*flips], scale, sum_target, sum_source));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0],
            FlipCopyAndScaleAndAdd(flips + 1, nodes[node].child[0], scale, sum_target, sum_source));
         break;
      }

   return new_node;
   }

int TreeFlips::FlipCopyAndScaleAndAdd(const int * flips, Tree & tree, int node,
    double scale, double sum)
   {
   int new_node = NewNode();

   switch (nodes[new_node].type = tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         nodes[new_node].type = TREE_NODE_LEAF;
         nodes[new_node].value = sum;
         break;
      case TREE_NODE_LEAF :
         nodes[new_node].value = sum + tree.nodes[node].value * scale;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0],
            FlipCopyAndScaleAndAdd(flips + 1, tree, tree.nodes[node].child[*flips], scale, sum));
         SAFE_SET(nodes[new_node].child[1],
            FlipCopyAndScaleAndAdd(flips + 1, tree, tree.nodes[node].child[!*flips], scale, sum));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0],
            FlipCopyAndScaleAndAdd(flips + 1, tree, tree.nodes[node].child[0], scale, sum));
         break;
      }

   return new_node;
   }

void TreeFlips::FlipGraftAndScale(const int * flips, int target, int source,
     double scale)
   {
   switch (nodes[target].type = nodes[source].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         nodes[target].value = nodes[source].value * scale;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[target].child[0],
            FlipCopyAndScale(flips + 1, nodes[source].child[*flips], scale));
         SAFE_SET(nodes[target].child[1],
            FlipCopyAndScale(flips + 1, nodes[source].child[!*flips], scale));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0],
            FlipCopyAndScale(flips + 1, nodes[source].child[0], scale));
         break;
      }
   }

void TreeFlips::FlipGraftAndScale(const int * flips, int target,
     Tree & tree, int source, double scale)
   {
   switch (nodes[target].type = tree.nodes[source].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         nodes[target].value = tree.nodes[source].value * scale;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[target].child[0],
            FlipCopyAndScale(flips + 1, tree, tree.nodes[source].child[*flips], scale));
         SAFE_SET(nodes[target].child[1],
            FlipCopyAndScale(flips + 1, tree, tree.nodes[source].child[!*flips], scale));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0],
            FlipCopyAndScale(flips + 1, tree, tree.nodes[source].child[0], scale));
         break;
      }
   }

void TreeFlips::FlipGraftAndScaleAndAdd(const int * flips, int target, int source,
     double scale, double sum_target, double sum_source)
   {
   switch (nodes[target].type = nodes[source].type)
      {
      case TREE_NODE_ZERO :
         nodes[target].type = TREE_NODE_LEAF;
         nodes[target].value = sum_target;
         nodes[source].type = TREE_NODE_LEAF;
         nodes[source].value = sum_source;
         break;
      case TREE_NODE_LEAF :
         nodes[target].value  = sum_target + nodes[source].value * scale;
         nodes[source].value += sum_source;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[target].child[0],
            FlipCopyAndScaleAndAdd(flips + 1, nodes[source].child[*flips], scale, sum_target, sum_source));
         SAFE_SET(nodes[target].child[1],
            FlipCopyAndScaleAndAdd(flips + 1, nodes[source].child[!*flips], scale, sum_target, sum_source));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0],
            FlipCopyAndScaleAndAdd(flips + 1, nodes[source].child[0], scale, sum_target, sum_source));
         break;
      }
   }

void TreeFlips::FlipGraftAndScaleAndAdd(const int * flips, int target,
   Tree & tree, int source, double scale, double sum)
   {
   switch (nodes[target].type = tree.nodes[source].type)
      {
      case TREE_NODE_ZERO :
         nodes[target].type = TREE_NODE_LEAF;
         nodes[target].value = sum;
         break;
      case TREE_NODE_LEAF :
         nodes[target].value = sum + nodes[source].value * scale;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[target].child[0],
            FlipCopyAndScaleAndAdd(flips + 1, tree, tree.nodes[source].child[*flips], scale, sum));
         SAFE_SET(nodes[target].child[1],
            FlipCopyAndScaleAndAdd(flips + 1, tree, tree.nodes[source].child[!*flips], scale, sum));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0],
            FlipCopyAndScaleAndAdd(flips + 1, tree, tree.nodes[source].child[0], scale, sum));
         break;
      }
   }

// Routines used by the haplotyping engine
//
//    When haplotyping the objective is to find the single most likely path
//    through inheritance vector space, so the likelihood at any point is
//    the likelihood for the most likely route there rather than the sum
//    over all possible routes -- a call to max replaces straight add --.
//

void TreeFlips::FlipGraftAndScaleAndChoose(const int * flips, int target,
                int source, double scale, double sum_target, double sum_source)
   {
   switch (nodes[target].type = nodes[source].type)
      {
      case TREE_NODE_ZERO :
         nodes[target].type = TREE_NODE_LEAF;
         nodes[target].value = sum_target;
         nodes[source].type = TREE_NODE_LEAF;
         nodes[source].value = sum_source;
         break;
      case TREE_NODE_LEAF :
         nodes[target].value = max(nodes[source].value * scale, sum_target);
         nodes[source].value = max(nodes[source].value, sum_source);
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[target].child[0],
            FlipCopyAndScaleAndChoose(flips + 1, nodes[source].child[*flips], scale, sum_target, sum_source));
         SAFE_SET(nodes[target].child[1],
            FlipCopyAndScaleAndChoose(flips + 1, nodes[source].child[!*flips], scale, sum_target, sum_source));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0],
            FlipCopyAndScaleAndChoose(flips + 1, nodes[source].child[0], scale, sum_target, sum_source));
         break;
      }
   }

int TreeFlips::FlipCopyAndScaleAndChoose(const int * flips, int node,
    double scale, double sum_target, double sum_source)
   {
   int new_node = NewNode();

   switch (nodes[new_node].type = nodes[node].type)
      {
      case TREE_NODE_ZERO :
         nodes[new_node].type = TREE_NODE_LEAF;
         nodes[new_node].value = sum_target;
         nodes[node].type = TREE_NODE_LEAF;
         nodes[node].value = sum_source;
         break;
      case TREE_NODE_LEAF :
         nodes[new_node].value = max(sum_target, nodes[node].value * scale);
         nodes[node].value     = max(sum_source, nodes[node].value);
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0],
            FlipCopyAndScaleAndChoose(flips + 1, nodes[node].child[*flips], scale, sum_target, sum_source));
         SAFE_SET(nodes[new_node].child[1],
            FlipCopyAndScaleAndChoose(flips + 1, nodes[node].child[!*flips], scale, sum_target, sum_source));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0],
            FlipCopyAndScaleAndChoose(flips + 1, nodes[node].child[0], scale, sum_target, sum_source));
         break;
      }

   return new_node;
   }

int TreeFlips::FlipCopyAndScaleAndChoose(const int * flips, Tree & tree, int node,
    double scale, double floor)
   {
   int new_node = NewNode();

   switch (nodes[new_node].type = tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         nodes[new_node].type = TREE_NODE_LEAF;
         nodes[new_node].value = floor;
         break;
      case TREE_NODE_LEAF :
         nodes[new_node].value = max(floor, tree.nodes[node].value * scale);
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0],
            FlipCopyAndScaleAndChoose(flips + 1, tree, tree.nodes[node].child[*flips], scale, floor));
         SAFE_SET(nodes[new_node].child[1],
            FlipCopyAndScaleAndChoose(flips + 1, tree, tree.nodes[node].child[!*flips], scale, floor));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0],
            FlipCopyAndScaleAndChoose(flips + 1, tree, tree.nodes[node].child[0], scale, floor));
         break;
      }

   return new_node;
   }

void TreeFlips::FlipGraftAndScaleAndChoose(const int * flips, int target,
   Tree & tree, int source, double scale, double floor)
   {
   switch (nodes[target].type = tree.nodes[source].type)
      {
      case TREE_NODE_ZERO :
         nodes[target].type = TREE_NODE_LEAF;
         nodes[target].value = floor;
         break;
      case TREE_NODE_LEAF :
         nodes[target].value = max(floor, nodes[source].value * scale);
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[target].child[0],
            FlipCopyAndScaleAndChoose(flips + 1, tree, tree.nodes[source].child[*flips], scale, floor));
         SAFE_SET(nodes[target].child[1],
            FlipCopyAndScaleAndChoose(flips + 1, tree, tree.nodes[source].child[!*flips], scale, floor));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0],
            FlipCopyAndScaleAndChoose(flips + 1, tree, tree.nodes[source].child[0], scale, floor));
         break;
      }
   }

void TreeFlips::DoubleFlipBranch(const int * flips, int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
      case TREE_NODE_LEAF :
         return;

      case TREE_NODE_ONE  :
         {
         int child = nodes[node].child[0];
         if ( *flips < 2 )
            DoubleFlipBranch(flips + 1, child);
         else
            switch (nodes[child].type)
               {
               case TREE_NODE_ZERO :
               case TREE_NODE_LEAF :
                  return;
               case TREE_NODE_ONE  :
                  DoubleFlipBranch(flips + 2, nodes[child].child[0]);
                  break;
               case TREE_NODE_TWO :
                  {
                  int new_node = NewNode();
                  nodes[node].type = TREE_NODE_TWO;
                  nodes[node].child[1] = new_node;
                  nodes[child].type = TREE_NODE_ONE;
                  nodes[new_node].type = TREE_NODE_ONE;
                  nodes[new_node].child[0] = nodes[child].child[1];
                  DoubleFlipBranch(flips + 2, nodes[child].child[0]);
                  DoubleFlipBranch(flips + 2, nodes[child].child[1]);
                  }
               }
         return;
         }

      case TREE_NODE_TWO :
         switch (*flips)
            {
            case 1 :
               {
               int temp = nodes[node].child[0];
               nodes[node].child[0] = nodes[node].child[1];
               nodes[node].child[1] = temp;
               }
            case 0 :
               DoubleFlipBranch(flips + 1, nodes[node].child[0]);
               DoubleFlipBranch(flips + 1, nodes[node].child[1]);
               return;
            case 2 :
               {
               int child = nodes[node].child[0];
               UpgradeNode(child, nodes[nodes[node].child[1]].type);
               UpgradeNode(nodes[node].child[1], nodes[child].type);

               if (nodes[child].type <= TREE_NODE_LEAF)
                  {
                  int new_node = NewNode();
                  nodes[new_node].type = TREE_NODE_TWO;
                  nodes[new_node].child[0] = nodes[node].child[0];
                  nodes[new_node].child[1] = nodes[node].child[1];
                  nodes[node].type = TREE_NODE_ONE;
                  nodes[node].child[0] = new_node;
                  }
               else if (nodes[child].type == TREE_NODE_ONE)
                  {
                  nodes[node].type = TREE_NODE_ONE;
                  nodes[child].type = TREE_NODE_TWO;
                  nodes[child].child[1] = nodes[nodes[node].child[1]].child[0];
                  DoubleFlipBranch(flips + 2, nodes[child].child[0]);
                  DoubleFlipBranch(flips + 2, nodes[child].child[1]);
                  }
               else
                  {
                  int other = nodes[node].child[1];
                  int temp  = nodes[child].child[1];
                  nodes[child].child[1] = nodes[other].child[0];
                  nodes[other].child[0] = temp;
                  DoubleFlipBranch(flips + 2, nodes[child].child[0]);
                  DoubleFlipBranch(flips + 2, nodes[child].child[1]);
                  DoubleFlipBranch(flips + 2, nodes[other].child[0]);
                  DoubleFlipBranch(flips + 2, nodes[other].child[1]);
                  }
               }
            }
      }
   }

void TreeFlips::UpgradeNode(int node, int new_type)
   {
   int new_node;

   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
      case TREE_NODE_LEAF :
         if (new_type < TREE_NODE_ONE) return;
         new_node = NewNode();
         if ((nodes[new_node].type = nodes[node].type) == TREE_NODE_LEAF)
            nodes[new_node].value = nodes[node].value;
         nodes[node].type = TREE_NODE_ONE;
         nodes[node].child[0] = new_node;
      case TREE_NODE_ONE :
         if (new_type < TREE_NODE_TWO) return;
         nodes[node].type = TREE_NODE_TWO;
         SAFE_SET(nodes[node].child[1], Copy(nodes[node].child[0]));
      case TREE_NODE_TWO :
         return;
      }
   }

void TreeFlips::PseudoUpgradeNode(int node, int new_type)
   {
   int new_node;

   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
      case TREE_NODE_LEAF :
         if (new_type < TREE_NODE_ONE) return;
         new_node = NewNode();
         if ((nodes[new_node].type = nodes[node].type) == TREE_NODE_LEAF)
            nodes[new_node].value = nodes[node].value;
         nodes[node].type = TREE_NODE_ONE;
         nodes[node].child[0] = new_node;
      case TREE_NODE_ONE :
         if (new_type < TREE_NODE_TWO) return;
         nodes[node].child[1] = nodes[node].child[0];
      case TREE_NODE_TWO :
         return;
      }
   }

 
