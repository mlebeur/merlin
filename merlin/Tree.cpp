////////////////////////////////////////////////////////////////////// 
// merlin/Tree.cpp 
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
 
#include "Tree.h"
#include "MathConstant.h"

void Tree::Multiply(int node, double scale)
   {
   stack.Dimension(bit_count + 1);
   int ptr = 1;

   stack[0] = node;

   while (ptr)
      {
      node = stack[--ptr];

      switch (nodes[node].type)
         {
         case TREE_NODE_ZERO :
            break;
         case TREE_NODE_LEAF :
            nodes[node].value *= scale;
            break;
         case TREE_NODE_TWO :
            stack[ptr++] = nodes[node].child[1];
         case TREE_NODE_ONE :
            stack[ptr++] = nodes[node].child[0];
            break;
         }
      }
   }

void Tree::Add(int node, double constant)
   {
   stack.Dimension(bit_count + 1);
   int ptr = 1;

   stack[0] = node;

   while (ptr)
      {
      node = stack[--ptr];

      switch (nodes[node].type)
         {
         case TREE_NODE_ZERO :
            nodes[node].type = TREE_NODE_LEAF;
            nodes[node].value = constant;
            break;
         case TREE_NODE_LEAF :
            nodes[node].value += constant;
            break;
         case TREE_NODE_ONE :
            stack[ptr++] = nodes[node].child[0];
            break;
         case TREE_NODE_TWO :
            stack[ptr++] = nodes[node].child[0];
            stack[ptr++] = nodes[node].child[1];
            break;
         }
      }
   }

void Tree::Floor(int node, double floor)
   {
   stack.Dimension(bit_count + 1);
   int ptr = 1;

   stack[0] = node;

   while (ptr)
      {
      node = stack[--ptr];

      switch (nodes[node].type)
         {
         case TREE_NODE_ZERO :
            nodes[node].type = TREE_NODE_LEAF;
            nodes[node].value = floor;
            break;
         case TREE_NODE_LEAF :
            nodes[node].value = max(nodes[node].value, floor);
            break;
         case TREE_NODE_ONE :
            stack[ptr++] = nodes[node].child[0];
            break;
         case TREE_NODE_TWO :
            stack[ptr++] = nodes[node].child[0];
            stack[ptr++] = nodes[node].child[1];
            break;
         }
      }
   }

void Tree::ScaleSum(int node, double scale, double sum)
   {
   stack.Dimension(bit_count + 1);
   int ptr = 1;

   stack[0] = node;

   while (ptr)
      {
      node = stack[--ptr];

      switch (nodes[node].type)
         {
         case TREE_NODE_ZERO :
            nodes[node].type = TREE_NODE_LEAF;
            nodes[node].value = sum;
            break;
         case TREE_NODE_LEAF :
            nodes[node].value = nodes[node].value * scale + sum;
            break;
         case TREE_NODE_TWO :
            stack[ptr++] = nodes[node].child[1];
         case TREE_NODE_ONE :
            stack[ptr++] = nodes[node].child[0];
            break;
         }
      }
   }

void Tree::Graft(int target, int source)
   {
   switch (nodes[target].type = nodes[source].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         nodes[target].value = nodes[source].value;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[target].child[0], Copy(nodes[source].child[0]));
         SAFE_SET(nodes[target].child[1], Copy(nodes[source].child[1]));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0], Copy(nodes[source].child[0]));
         break;
      }
   }

void Tree::GraftAndScale(int target, int source, double scale)
   {
   switch (nodes[target].type = nodes[source].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         nodes[target].value  = nodes[source].value * scale;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[target].child[0], CopyAndScale(nodes[source].child[0], scale));
         SAFE_SET(nodes[target].child[1], CopyAndScale(nodes[source].child[1], scale));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0], CopyAndScale(nodes[source].child[0], scale));
         break;
      }
   }

void Tree::GraftAndScale(int target, Tree & tree, int source, double scale)
   {
   switch (nodes[target].type = tree.nodes[source].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         nodes[target].value = tree.nodes[source].value * scale;
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0],
            CopyAndScale(tree, tree.nodes[source].child[0], scale));
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[target].child[0],
            CopyAndScale(tree, tree.nodes[source].child[0], scale));
         SAFE_SET(nodes[target].child[1],
            CopyAndScale(tree, tree.nodes[source].child[1], scale));
         break;
      }
   }

void Tree::GraftAndScaleAndAdd
   (int target, Tree & tree, int source, double scale, double add)
   {
   switch (nodes[target].type = tree.nodes[source].type)
      {
      case TREE_NODE_ZERO :
         nodes[target].type = TREE_NODE_LEAF;
         nodes[target].value = add;
         break;
      case TREE_NODE_LEAF :
         nodes[target].value = tree.nodes[source].value * scale + add;
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0],
            CopyAndScaleAndAdd(tree, tree.nodes[source].child[0], scale, add));
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[target].child[0],
            CopyAndScaleAndAdd(tree, tree.nodes[source].child[0], scale, add));
         SAFE_SET(nodes[target].child[1],
            CopyAndScaleAndAdd(tree, tree.nodes[source].child[1], scale, add));
         break;
      }
   }

void Tree::GraftAndScaleAndChoose
   (int target, Tree & tree, int source, double scale, double floor)
   {
   switch (nodes[target].type = tree.nodes[source].type)
      {
      case TREE_NODE_ZERO :
         nodes[target].type = TREE_NODE_LEAF;
         nodes[target].value = floor;
         break;
      case TREE_NODE_LEAF :
         nodes[target].value = max(tree.nodes[source].value * scale, floor);
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0],
            CopyAndScaleAndChoose(tree, tree.nodes[source].child[0], scale, floor));
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[target].child[0],
            CopyAndScaleAndChoose(tree, tree.nodes[source].child[0], scale, floor));
         SAFE_SET(nodes[target].child[1],
            CopyAndScaleAndChoose(tree, tree.nodes[source].child[1], scale, floor));
         break;
      }
   }

void Tree::GraftAndScaleAndAdd(int target, int source, double scale,
                               double sum_target, double sum_source)
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
            CopyAndScaleAndAdd(nodes[source].child[0], scale, sum_target, sum_source));
         SAFE_SET(nodes[target].child[1],
            CopyAndScaleAndAdd(nodes[source].child[1], scale, sum_target, sum_source));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0],
            CopyAndScaleAndAdd(nodes[source].child[0], scale, sum_target, sum_source));
         break;
      }
   }

int Tree::CopyAndScale(int node, double scale)
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
         SAFE_SET(nodes[new_node].child[0], CopyAndScale(nodes[node].child[0], scale));
         SAFE_SET(nodes[new_node].child[1], CopyAndScale(nodes[node].child[1], scale));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0], CopyAndScale(nodes[node].child[0], scale));
         break;
      }

   return new_node;
   }

int Tree::CopyAndScale(Tree & tree, int node, double scale)
   {
   int new_node = NewNode();

   switch (nodes[new_node].type = tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         nodes[new_node].value = tree.nodes[node].value * scale;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0],
            CopyAndScale(tree, tree.nodes[node].child[0], scale));
         SAFE_SET(nodes[new_node].child[1],
            CopyAndScale(tree, tree.nodes[node].child[1], scale));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0],
            CopyAndScale(tree, tree.nodes[node].child[0], scale));
         break;
      }

   return new_node;
   }

int Tree::CopyAndScaleAndAdd(Tree & tree, int node, double scale, double add)
   {
   int new_node = NewNode();

   switch (nodes[new_node].type = tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         nodes[new_node].type = TREE_NODE_LEAF;
         nodes[new_node].value = add;
         break;
      case TREE_NODE_LEAF :
         nodes[new_node].value = tree.nodes[node].value * scale + add;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0],
            CopyAndScaleAndAdd(tree, tree.nodes[node].child[0], scale, add));
         SAFE_SET(nodes[new_node].child[1],
            CopyAndScaleAndAdd(tree, tree.nodes[node].child[1], scale, add));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0],
            CopyAndScaleAndAdd(tree, tree.nodes[node].child[0], scale, add));
         break;
      }

   return new_node;
   }

int Tree::CopyAndScaleAndChoose(Tree & tree, int node, double scale, double floor)
   {
   int new_node = NewNode();

   switch (nodes[new_node].type = tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         nodes[new_node].type = TREE_NODE_LEAF;
         nodes[new_node].value = floor;
         break;
      case TREE_NODE_LEAF :
         nodes[new_node].value = max(tree.nodes[node].value * scale, floor);
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0],
            CopyAndScaleAndChoose(tree, tree.nodes[node].child[0], scale, floor));
         SAFE_SET(nodes[new_node].child[1],
            CopyAndScaleAndChoose(tree, tree.nodes[node].child[1], scale, floor));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0],
            CopyAndScaleAndChoose(tree, tree.nodes[node].child[0], scale, floor));
         break;
      }

   return new_node;
   }

int Tree::CopyAndScaleAndAdd(int node, double scale, double sum_target,
                             double sum_source)
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
            CopyAndScaleAndAdd(nodes[node].child[0], scale, sum_target, sum_source));
         SAFE_SET(nodes[new_node].child[1],
            CopyAndScaleAndAdd(nodes[node].child[1], scale, sum_target, sum_source));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0],
            CopyAndScaleAndAdd(nodes[node].child[0], scale, sum_target, sum_source));
         break;
      }

   return new_node;
   }

void Tree::Multiply(Tree & tree, int node1, int node2)
   {
   switch (PAIR_NODES(nodes[node1].type, tree.nodes[node2].type))
      {
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ONE)  :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_TWO)  :
         break;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ZERO)  :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ZERO)  :
         nodes[node1].type = TREE_NODE_ZERO;
         break;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_LEAF) :
         Multiply(nodes[node1].child[0], tree.nodes[node2].value );
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_LEAF) :
         Multiply(nodes[node1].child[0], tree.nodes[node2].value );
         Multiply(nodes[node1].child[1], tree.nodes[node2].value );
         break;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_LEAF) :
         nodes[node1].value *= tree.nodes[node2].value;
         break;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
         nodes[node1].type = TREE_NODE_ONE;
         SAFE_SET(nodes[node1].child[0],
            CopyAndScale(tree, tree.nodes[node2].child[0], nodes[node1].value));
         break;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         {
         double old_value = nodes[node1].value;
         nodes[node1].type = TREE_NODE_TWO;
         SAFE_SET(nodes[node1].child[0],
            CopyAndScale(tree, tree.nodes[node2].child[0], old_value));
         SAFE_SET(nodes[node1].child[1],
            CopyAndScale(tree, tree.nodes[node2].child[1], old_value));
         }
         break;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
         Multiply(tree, nodes[node1].child[0], tree.nodes[node2].child[0]);
         break;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_TWO) :
         nodes[node1].type = TREE_NODE_TWO;
         SAFE_SET(nodes[node1].child[1], Copy(nodes[node1].child[0]));
         Multiply(tree, nodes[node1].child[0], tree.nodes[node2].child[0]);
         Multiply(tree, nodes[node1].child[1], tree.nodes[node2].child[1]);
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ONE) :
         Multiply(tree, nodes[node1].child[0], tree.nodes[node2].child[0]);
         Multiply(tree, nodes[node1].child[1], tree.nodes[node2].child[0]);
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_TWO) :
         Multiply(tree, nodes[node1].child[0], tree.nodes[node2].child[0]);
         Multiply(tree, nodes[node1].child[1], tree.nodes[node2].child[1]);
         break;
      }
   }

double Tree::Mean(int node, double weight) const
   {
   switch (nodes[node].type)
     {
     case TREE_NODE_ZERO :
        return 0.0;
     case TREE_NODE_LEAF :
        return nodes[node].value * weight;
     case TREE_NODE_ONE :
        return Mean(nodes[node].child[0], weight);
     case TREE_NODE_TWO :
        return Mean(nodes[node].child[0], weight * 0.5) +
               Mean(nodes[node].child[1], weight * 0.5);
     }
   return 0.0;
   }

double Tree::MeanProduct(Tree & tree, int node1, int node2, double weight)
   {
   switch (PAIR_NODES(nodes[node1].type, tree.nodes[node2].type))
      {
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ONE)  :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_TWO)  :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ZERO)  :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ZERO)  :
         return 0.0;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_LEAF) :
         return nodes[node1].value * tree.nodes[node2].value * weight;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_LEAF) :
         return Mean(node1, weight) * tree.nodes[node2].value;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         return tree.Mean(node2, weight) * nodes[node1].value;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
         return MeanProduct(tree, nodes[node1].child[0],
                            tree.nodes[node2].child[0], weight);
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_TWO) :
         return MeanProduct(tree, nodes[node1].child[0],
                            tree.nodes[node2].child[0], weight * 0.5) +
                MeanProduct(tree, nodes[node1].child[0],
                            tree.nodes[node2].child[1], weight * 0.5);
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ONE) :
         return MeanProduct(tree, nodes[node1].child[0],
                            tree.nodes[node2].child[0], weight * 0.5) +
                MeanProduct(tree, nodes[node1].child[1],
                            tree.nodes[node2].child[0], weight * 0.5);
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_TWO) :
         return MeanProduct(tree, nodes[node1].child[0],
                            tree.nodes[node2].child[0], weight * 0.5) +
                MeanProduct(tree, nodes[node1].child[1],
                            tree.nodes[node2].child[1], weight * 0.5);
      }
   return 0.0;
   }

// Branch grafting and copying routines used by the haplotyping engine
//
//    When haplotyping the objective is to find the single most likely path
//    through inheritance vector space, so the likelihood at any point is
//    the likelihood for the most likely route there rather than the sum
//    over all possible routes -- a call to max replaces straight add --.
//

int Tree::CopyAndScaleAndChoose(int node, double scale, double sum_target,
                                double sum_source)
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
         nodes[new_node].value  = max(sum_target, nodes[node].value * scale);
         nodes[node].value = max(nodes[node].value, sum_source);
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0],
            CopyAndScaleAndChoose(nodes[node].child[0], scale, sum_target, sum_source));
         SAFE_SET(nodes[new_node].child[1],
            CopyAndScaleAndChoose(nodes[node].child[1], scale, sum_target, sum_source));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0],
            CopyAndScaleAndChoose(nodes[node].child[0], scale, sum_target, sum_source));
         break;
      }

   return new_node;
   }

void Tree::GraftAndScaleAndChoose(int target, int source,
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
         nodes[target].value  = max(nodes[source].value * scale, sum_target);
         nodes[source].value  = max(nodes[source].value, sum_source);
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[target].child[0],
            CopyAndScaleAndChoose(nodes[source].child[0], scale, sum_target, sum_source));
         SAFE_SET(nodes[target].child[1],
            CopyAndScaleAndChoose(nodes[source].child[1], scale, sum_target, sum_source));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[target].child[0],
            CopyAndScaleAndChoose(nodes[source].child[0], scale, sum_target, sum_source));
         break;
      }
   }

void Tree::MakeTree(const int * bits, int depth, double value)
   {
   Clear();

   bit_count = depth;

   for (int i = depth - 1; i >= 0; i--)
      {
      int node = NewNode();
      nodes[node].type = TREE_NODE_TWO;
      nodes[node].child[bits[i] != 0] = (depth - i) * 2;
      nodes[node].child[bits[i] == 0] = (depth - i) * 2 - 1;
      node = NewNode();
      nodes[node].type = TREE_NODE_ZERO;
      }

   int node = NewNode();
   nodes[node].type = TREE_NODE_LEAF;
   nodes[node].value = value;
   }

void Tree::Add(Tree & tree, int node1, int node2)
   {
   switch (PAIR_NODES(nodes[node1].type, tree.nodes[node2].type))
      {
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_LEAF) :
         nodes[node1].type = TREE_NODE_LEAF;
         nodes[node1].value = tree.nodes[node2].value;
         break;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ONE)  :
         nodes[node1].type = TREE_NODE_ONE;
         SAFE_SET(nodes[node1].child[0], Copy(tree, tree.nodes[node2].child[0]));
         break;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_TWO)  :
         nodes[node1].type = TREE_NODE_ONE;
         SAFE_SET(nodes[node1].child[0], Copy(tree, tree.nodes[node2].child[0]));
         SAFE_SET(nodes[node1].child[1], Copy(tree, tree.nodes[node2].child[1]));
         break;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ZERO)  :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ZERO)  :
         break;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_LEAF) :
         Add(node1, tree.nodes[node2].value);
         break;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_LEAF) :
         nodes[node1].value += tree.nodes[node2].value;
         break;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
         nodes[node1].type = TREE_NODE_ONE;
         SAFE_SET(nodes[node1].child[0],
            CopyAndAdd(tree, tree.nodes[node2].child[0], nodes[node1].value));
         break;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         {
         double old_value = nodes[node1].value;
         nodes[node1].type = TREE_NODE_TWO;
         SAFE_SET(nodes[node1].child[0],
            CopyAndAdd(tree, tree.nodes[node2].child[0], old_value));
         SAFE_SET(nodes[node1].child[1],
            CopyAndAdd(tree, tree.nodes[node2].child[1], old_value));
         }
         break;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
         Add(tree, nodes[node1].child[0], tree.nodes[node2].child[0]);
         break;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_TWO) :
         nodes[node1].type = TREE_NODE_TWO;
         SAFE_SET(nodes[node1].child[1], Copy(nodes[node1].child[0]));
         Add(tree, nodes[node1].child[0], tree.nodes[node2].child[0]);
         Add(tree, nodes[node1].child[1], tree.nodes[node2].child[1]);
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ONE) :
         Add(tree, nodes[node1].child[0], tree.nodes[node2].child[0]);
         Add(tree, nodes[node1].child[1], tree.nodes[node2].child[0]);
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_TWO) :
         Add(tree, nodes[node1].child[0], tree.nodes[node2].child[0]);
         Add(tree, nodes[node1].child[1], tree.nodes[node2].child[1]);
         break;
      }
   }

int Tree::CopyAndAdd(Tree & tree, int node, double summand)
   {
   int new_node = NewNode();

   switch (nodes[new_node].type = tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         nodes[new_node].type  = TREE_NODE_LEAF;
         nodes[new_node].value = summand;
         break;
      case TREE_NODE_LEAF :
         nodes[new_node].value = tree.nodes[node].value + summand;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0],
            CopyAndAdd(tree, tree.nodes[node].child[0], summand));
         SAFE_SET(nodes[new_node].child[1],
            CopyAndAdd(tree, tree.nodes[node].child[1], summand));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0],
            CopyAndAdd(tree, tree.nodes[node].child[0], summand));
         break;
      }

   return new_node;
   }






 
