////////////////////////////////////////////////////////////////////// 
// merlin/TreeInfo.cpp 
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
 
#include "TreeInfo.h"
#include "MathConstant.h"

#include <math.h>

TreeInfo::TreeInfo()
   {
   nodeCount = zeroCount = leafCount = 0;
   sumInfo = information = 0.0;
   min = max = mean = var = 0.0;
   }

double TreeInfo::GetMean(BasicTree & tree, int node)
   {
   min =  1e100;
   max = -1e100;
   mean = var = 0.0;

   stack.Dimension(tree.bit_count + 1);
   weights.Dimension(tree.bit_count + 1);

   int    ptr = 1;
   double weight;

   stack[0] = node;
   weights[0] = 1.0;

   while (ptr--)
      {
      weight = weights[ptr];
      node   = stack[ptr];

      switch (tree.nodes[node].type)
         {
         case TREE_NODE_ZERO :
            break;
         case TREE_NODE_LEAF :
            mean += tree.nodes[node].value * weight;
            break;
         case TREE_NODE_ONE :
            weights[ptr] = weight;
            stack[ptr++] = tree.nodes[node].child[0];
            break;
         case TREE_NODE_TWO :
            weights[ptr] = weight * 0.5;
            stack[ptr++] = tree.nodes[node].child[1];
            weights[ptr] = weight * 0.5;
            stack[ptr++] = tree.nodes[node].child[0];
            break;
         }
      }

   return mean;
   }

void TreeInfo::GetMeanVar(BasicTree & tree, int node)
   {
   min =  1e100;
   max = -1e100;
   mean = var = 0.0;

   stack.Dimension(tree.bit_count + 1);
   int ptr = 1;

   weights.Dimension(tree.bit_count + 1);
   double weight;

   stack[0] = node;
   weights[0] = pow(2.0, tree.bit_count);

   while (ptr--)
      {
      weight = weights[ptr];
      node = stack[ptr];

      switch (tree.nodes[node].type)
         {
         case TREE_NODE_ZERO :
            break;
         case TREE_NODE_LEAF :
            mean += tree.nodes[node].value * weight;
            var  += tree.nodes[node].value * tree.nodes[node].value * weight;
            break;
         case TREE_NODE_ONE :
            weights[ptr] = weight;
            stack[ptr++] = tree.nodes[node].child[0];
            break;
         case TREE_NODE_TWO :
            weights[ptr] = weight * 0.5;
            stack[ptr++] = tree.nodes[node].child[1];
            weights[ptr] = weight * 0.5;
            stack[ptr++] = tree.nodes[node].child[0];
            break;
         }
      }

   if (tree.bit_count)
      var = (var - mean * mean * pow( 2.0, -tree.bit_count))
            * (pow (2.0, -tree.bit_count)) + TINY;
   else
      var = TINY;
   mean *= pow(2.0, -tree.bit_count);
   }

void TreeInfo::GetBounds(BasicTree & tree, int node)
   {
   min =  1e100;
   max = -1e100;
   mean = var = 0.0;

   stack.Dimension(tree.bit_count + 1);
   int ptr = 1;

   weights.Dimension(tree.bit_count + 1);
   double weight;
   bool   have_zero = false;

   stack[0] = node;
   weights[0] = pow(2.0, tree.bit_count);

   while (ptr--)
      {
      weight = weights[ptr];
      node = stack[ptr];

      switch (tree.nodes[node].type)
         {
         case TREE_NODE_ZERO :
            have_zero = true;
            break;
         case TREE_NODE_LEAF :
            if (tree.nodes[node].value < min) min = tree.nodes[node].value;
            if (tree.nodes[node].value > max) max = tree.nodes[node].value;
            mean += tree.nodes[node].value * weight;
            var  += tree.nodes[node].value * tree.nodes[node].value * weight;
            break;
         case TREE_NODE_ONE :
            weights[ptr] = weight;
            stack[ptr++] = tree.nodes[node].child[0];
            break;
         case TREE_NODE_TWO :
            weights[ptr] = weight * 0.5;
            stack[ptr++] = tree.nodes[node].child[1];
            weights[ptr] = weight * 0.5;
            stack[ptr++] = tree.nodes[node].child[0];
            break;
         }
      }

   if (tree.bit_count)
      var = (var - mean * mean * pow( 2.0, -tree.bit_count))
            * (pow (2.0, -tree.bit_count)) + TINY;
   else
      var = TINY;

   if (have_zero && min > 0.0)
      min = 0.0;
   if (have_zero && max < 0.0)
      max = 0.0;

   mean *= pow(2.0, -tree.bit_count);
   }

void TreeInfo::GetTreeInfo(BasicTree & tree, int node)
   {
   nodeCount = zeroCount = leafCount = 0;
   mean = sumInfo = 0.0;
   min =  1e100;
   max = -1e100;

   weights.Dimension(tree.bit_count + 1);
   double weight;

   stack.Dimension(tree.bit_count + 1);
   int ptr = 1;

   stack[0] = node;
   weights[0] = pow(2.0, tree.bit_count);

   while (ptr--)
      {
      weight = weights[ptr];
      node = stack[ptr];

      switch (tree.nodes[node].type)
         {
         case TREE_NODE_ZERO :
            break;
         case TREE_NODE_LEAF :
            if (tree.nodes[node].value == 0.0) break;
            if (weight > 1)
               {
               mean += tree.nodes[node].value * weight;
               sumInfo += tree.nodes[node].value *
                          log(tree.nodes[node].value) * weight;
               break;
               }
            mean += tree.nodes[node].value;
            sumInfo += tree.nodes[node].value * log(tree.nodes[node].value);
            break;
         case TREE_NODE_ONE :
            weights[ptr] = weight;
            stack[ptr++] = tree.nodes[node].child[0];
            break;
         case TREE_NODE_TWO :
            weights[ptr] = weight * 0.5;
            stack[ptr++] = tree.nodes[node].child[1];
            weights[ptr] = weight * 0.5;
            stack[ptr++] = tree.nodes[node].child[0];
            break;
         }
      }

   if (mean && tree.bit_count)
      {
      // As usual, consider entropy to be
      //   Sum p * ln (p) for all p
      //   All p sum to 1.0

      // The prior probalities are assumed to be 1 / N
      //    E = SUM OVER N of 1 / N * LOG(1/N)
      //      = LOG(1/N)
      double prior = log(pow(2.0, -tree.bit_count)) + TINY;

      // The posterior probabilities are the value in each leaf
      // divide by the sum of all leafs, so...
      //    E = SUM (leaf / TOTAL) * log (leaf / TOTAL) =
      //      = 1 / TOTAL * SUM leaf * { log(leaf) - LOG(TOTAL) } =
      //      = 1 / TOTAL * {SUM leaf * log(leaf) - LOG(TOTAL) * TOTAL }
      //      = 1 / TOTAL * SUM leaf * log(leaf) - LOG(TOTAL)
      //
      // Obviously, TOTAL = Sum (leaf)

      // Stop negative information due to roundoff!
      double posterior = sumInfo / mean - log(mean);
      information = ::max(1.0 - posterior / prior, 0.0);
      }
   else
      information = 0.0;
   }

double TreeInfo::GetMinInfo(BasicTree & tree, int node)
   {
   double minInfo = 0.0;

   weights.Dimension(tree.bit_count + 1);
   double weight;

   stack.Dimension(tree.bit_count + 1);
   int ptr = 1;

   stack[0] = node;
   weights[0] = pow(2.0, tree.bit_count);

   while (ptr--)
      {
      weight = weights[ptr];
      node = stack[ptr];

      switch (tree.nodes[node].type)
         {
         case TREE_NODE_ZERO :
            break;
         case TREE_NODE_LEAF :
            if (tree.nodes[node].value == 0.0) break;
            minInfo += weight;
            break;
         case TREE_NODE_ONE :
            weights[ptr] = weight;
            stack[ptr++] = tree.nodes[node].child[0];
            break;
         case TREE_NODE_TWO :
            weights[ptr] = weight * 0.5;
            stack[ptr++] = tree.nodes[node].child[1];
            weights[ptr] = weight * 0.5;
            stack[ptr++] = tree.nodes[node].child[0];
            break;
         }
      }

   if (minInfo && tree.bit_count)
      return 1.0 - log (1.0 / minInfo) / (log(pow(2.0, -tree.bit_count)) + TINY);
   else
      return 0.0;
   }

void TreeInfo::GetDetailedInfo(BasicTree & tree, int node)
   {
   nodeCount = zeroCount = leafCount = 0;
   mean = sumInfo = 0.0;
   min =  1e100;
   max = -1e100;

   weights.Dimension(tree.bit_count + 1);
   double weight;

   stack.Dimension(tree.bit_count + 1);
   int ptr = 1;

   stack[0] = node;
   weights[0] = pow(2.0, tree.bit_count);

   while (ptr--)
      {
      nodeCount++;

      weight = weights[ptr];
      node = stack[ptr];

      switch (tree.nodes[node].type)
         {
         case TREE_NODE_ZERO :
            zeroCount++;
            break;
         case TREE_NODE_LEAF :
            leafCount++;
            if (tree.nodes[node].value == 0) break;
            if (tree.nodes[node].value > max)
               max = tree.nodes[node].value;
            if (tree.nodes[node].value < min)
               min = tree.nodes[node].value;
            if (weight > 1)
               {
               mean += tree.nodes[node].value * weight;
               sumInfo += tree.nodes[node].value *
                          log(tree.nodes[node].value) * weight;
               break;
               }
            mean += tree.nodes[node].value;
            sumInfo += tree.nodes[node].value * log(tree.nodes[node].value);
            break;
         case TREE_NODE_ONE :
            weights[ptr] = weight;
            stack[ptr++] = tree.nodes[node].child[0];
            break;
         case TREE_NODE_TWO :
            weights[ptr] = weight * 0.5;
            stack[ptr++] = tree.nodes[node].child[1];
            weights[ptr] = weight * 0.5;
            stack[ptr++] = tree.nodes[node].child[0];
            break;
         }
      }

   if (mean && tree.bit_count)
      {
      // As usual, consider entropy to be
      //   Sum p * ln (p) for all p
      //   All p sum to 1.0

      // The prior probalities are assumed to be 1 / N
      //    E = SUM OVER N of 1 / N * LOG(1/N)
      //      = LOG(1/N)
      double prior = log(pow(2.0, -tree.bit_count)) + TINY;

      // The posterior probabilities are the value in each leaf
      // divide by the sum of all leafs, so...
      //    E = SUM (leaf / TOTAL) * log (leaf / TOTAL) =
      //      = 1 / TOTAL * SUM leaf * { log(leaf) - LOG(TOTAL) } =
      //      = 1 / TOTAL * {SUM leaf * log(leaf) - LOG(TOTAL) * TOTAL }
      //      = 1 / TOTAL * SUM leaf * log(leaf) - LOG(TOTAL)
      //
      // Obviously, TOTAL = Sum (leaf)

      // Stop negative information due to roundoff!
      double posterior = sumInfo / mean - log(mean);
      information = ::max(1.0 - posterior / prior, 0.0);
      }
   else
      information = 0.0;
   }

 
