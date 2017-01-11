////////////////////////////////////////////////////////////////////// 
// regress/RegressKinship.cpp 
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
 
#include "RegressKinship.h"

#include <math.h>

// The routines in this section perform Kinship calculations
//    Calculate(...) calculates and outputs kinship coefficients
//                   returning the overall likelihood of all vectors in the tree
//    Score(...)  traverses an inheritance tree and fills kinship matrix
//                   with unscaled kinship coefficients

// Constructor and destructor
//

RegressKinship::RegressKinship(Mantra & m)
   : mantra(m)
   {
   ped = NULL;
   family = NULL;
   sym_size = 0;
   }

RegressKinship::~RegressKinship()
   {
   if (sym_size) delete [] symmetries;
   }

// Call this routine to calculate pairwise kinship coefficients and their
// cross products conditional on some gene flow tree
//

double RegressKinship::Calculate(Tree & tree)
   {
   int pairs = (count - founders) * (count + founders - 1) / 2;

   // Dimension Kinship matrices
   kin.Dimension(pairs);
   kin.Zero();

   kinship.Dimension(pairs);
   kinship.Zero();

   kinship2.Dimension(pairs * (pairs + 1) / 2);
   kinship2.Zero();

   // Calculate IBDs as well as the overall likelihood
   double likelihood = Score(tree, 0, tree.bit_count - 1, mantra.two_f);

   kinship.Multiply(0.25 / (likelihood * alternate_states));
   kinship2.Multiply(0.0625 / (likelihood * alternate_states));

   return likelihood * pow(2.0, -mantra.bit_count);
   }

// This is routine is the actual work-horse. It traverse the gene flow tree
// and updates kinship coefficients at each possible leaf node
//

double RegressKinship::Score(Tree & tree, int node, int bit, int start, int pos)
   {
   int pivot = (bit >= 0) ? mantra.bits[bit] : 0;
   int end = (bit >= 0) ? pivot & ~1 : mantra.two_n;

   if (bit != mantra.bit_count - 1)
      {
      int pivot = bit == -1 ? mantra.two_n : mantra.bits[bit];
      int last_pivot = mantra.bits[bit+1] + 1;

      for (int i = last_pivot; i < pivot; i++)
         mantra.state[i] = mantra.state[mantra.vector[i]];
      }

   int start_pos = pos;

   for ( int i = start, ii = start >> 1 ; i < end; i++, i++, ii++)
      for ( int j = 0, jj = 0; j < i; j++, j++, jj++, pos++)
         if (mantra.ibd[ii][jj] == MANTRA_IBD_UNKNOWN)
            kin[pos] =
               ((mantra.state[i]   == mantra.state[j]   ) +
                (mantra.state[i+1] == mantra.state[j+1] ) +
                (mantra.state[i]   == mantra.state[j+1] ) +
                (mantra.state[i+1] == mantra.state[j]   ));

   // This variable is initialized to prevent bogus compiler warnings
   double sum = 0.0;

   if (bit >= 0)
      {
      switch (tree.nodes[node].type)
         {
         case TREE_NODE_ZERO :
            return 0.0;
         case TREE_NODE_ONE  :
            node = tree.nodes[node].child[0];
         case TREE_NODE_LEAF :
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] |= 1];
            sum = Score(tree, node, bit - 1, end, pos);
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] &= ~1];
            sum += Score(tree, node, bit - 1, end, pos);
            break;
         case TREE_NODE_TWO :
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] |= 1];
            sum = Score(tree, tree.nodes[node].child[1],bit - 1, end, pos);
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] &= ~1];
            sum += Score(tree, tree.nodes[node].child[0], bit - 1, end, pos);
            break;
         }
      }
   else
      {
      // At the base of the tree there are only zero nodes
      // and leaf nodes ...

      // Zero nodes we ignore
      if (tree.nodes[node].type == TREE_NODE_ZERO)
         return 0.0;

      // So this is a leaf node
      sum = tree.nodes[node].value;
      }

   if (start_pos != pos)
      {
      UpdateScores(kin, start_pos, pos, sum);

      for (int state = 1; state < alternate_states; state++)
         {
         Unravel(state);
         UpdateScores(unraveled, start_pos, pos, sum);
         }
      }

   return sum;
   }

// Select a family for analysis
//

void RegressKinship::SelectFamily(Pedigree * p, Family * f)
   {
   ped = p;
   family = f;

   founders = f->founders;
   count = f->count;

   alternate_states = 1 << mantra.couples ;

   if (alternate_states == 1) return;

   ResetSymmetries();

   for (int i = family->founders, position = 0; i < family->count; i++)
      for (int j = 0; j < i; j++)
         {
         if (j < family->founders && mantra.couple_index[j] != -1)
            symmetries[mantra.couple_index[j]].Push(position);
         position++;
         }
   }

// Routines for retrieving kinship coefficients and their cross-products
//

double RegressKinship::Retrieve(int person1, int person2)
   {
   if (person2 > person1)
      { int swap = person1; person1 = person2; person2 = swap; }

   switch (mantra.ibd[person1][person2])
      {
      case MANTRA_IBD_ZERO :
         return 0.0;
      case MANTRA_IBD_HALF :
         return 0.25;
      case MANTRA_IBD_ONE :
         return 0.50;
      default :
         {
         int index = (person1 - founders) * (founders + person1 - 1) / 2 + person2;

         return kinship[index];
         }
      }
   }

double RegressKinship::Retrieve(int person1, int person2, int person3, int person4)
   {
   if (person2 > person1)
      { int swap = person1; person1 = person2; person2 = swap; }

   if (person4 > person3)
      { int swap = person3; person3 = person4; person4 = swap; }

   if (mantra.ibd[person1][person2] != MANTRA_IBD_UNKNOWN ||
       mantra.ibd[person3][person4] != MANTRA_IBD_UNKNOWN)
      return Retrieve(person1, person2) * Retrieve(person3, person4);

   int index1 = (person1 - founders) * (founders + person1 - 1) / 2 + person2;
   int index2 = (person3 - founders) * (founders + person3 - 1) / 2 + person4;

   if (index1 > index2)
      return kinship2[index1 * (index1 + 1) / 2 + index2];
   else
      return kinship2[index2 * (index2 + 1) / 2 + index1];
   }

// Additional routines for handling symmetries

void RegressKinship::ResetSymmetries()
   {
   int couples = mantra.couples;

   if (couples > sym_size)
      {
      if (sym_size) delete [] symmetries;

      sym_size = (couples + 3) & (~3);
      symmetries = new IntArray[sym_size];
      }

   for (int i = 0; i < couples; i++)
      symmetries[i].Clear();
   }

void RegressKinship::Unravel(int idx)
   {
   unraveled = kin;

   int couple = 0;

   while (idx > 0)
      {
      if (idx & 1)
         for (int i = 0; i < symmetries[couple].Length(); i+=2)
            {
            unraveled[symmetries[couple][i]] = kin[symmetries[couple][i+1]];
            unraveled[symmetries[couple][i+1]] = kin[symmetries[couple][i]];
            }

      idx >>= 1;
      couple++;
      }
   }

void RegressKinship::UpdateScores(IntArray & k, int start_pos, int end_pos, double sum)
   {
   int pair_index = start_pos * (start_pos + 1) / 2;

   // Update cross-products
   for (int i = start_pos; i < end_pos; i++)
      if (kin[i])
         for (int j = 0; j <= i; j++)
            kinship2[pair_index++] += k[i] * k[j] * sum;
      else
         pair_index += i + 1;

   // Update kinship coefficients
   for (int i = start_pos; i < end_pos; i++)
      kinship[i] += k[i] * sum;
   }
 
