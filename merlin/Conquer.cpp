////////////////////////////////////////////////////////////////////// 
// merlin/Conquer.cpp 
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
 
#include "Conquer.h"
#include "MapFunction.h"

#include <math.h>

int Multipoint::maximum_recombinants = 0;

////////////////////////////////////////////////////////////////////////////
// This unit implements the Elston and Idury divide-and-conquer algorithm
// for general multipoint calculations
//

void Multipoint::MoveAlong(Mantra & m, double distance, double rescale)
   {
   if (distance == 0.0)
      return;

   if (maximum_recombinants)
      {
      Approximation(m, distance, maximum_recombinants, rescale);
      return;
      }

   theta = DistanceToRecombination(distance);
   oneminus = 1.0 - theta;
   x = theta / oneminus;
   onex = 1.0 + x;

   // The rescale constant must take hidden bits into accounts
   // if likelihoods are going to be used for something other
   // than comparing alternative gene flow patterns.
   pow_onex.Dimension(bit_count + 1);
   pow_onex[0] = pow(oneminus, bit_count + m.hidden_bit_count) * rescale;
   for (int i = 1; i <= bit_count; i++)
      pow_onex[i] = pow_onex[i - 1] * onex;

   for (int i = 0; i < m.two_f; i+=2)
      if (m.founder_bit_count[i])
         SpecialReUnite(0, 0, m.founder_bits[i]);

   DivideAndConquer(bit_count);
   for (int i = 0; i < m.couples; i++)
      if (m.couple_bits[i][m.two_n - 1])
         DoubleReUnite(0, 0, m.couple_bits[i]);
   }

////////////////////////////////////////////////////////////////////////////
// The version below allows for different male and female recombination
// fractions.
//

void Multipoint::MoveAlong(Mantra & m, double theta[2], double rescale)
   {
   if (theta[0] == 0.0 && theta[1] == 0.0)
      {
      if (rescale != 1.0)
         Multiply(0, rescale);
      return;
      }

   if (maximum_recombinants)
      {
      Approximation(m, theta, maximum_recombinants, rescale);
      return;
      }

   recombination[0] = theta[0];
   recombination[1] = theta[1];

   complement[0] = 1.0 - theta[0];
   complement[1] = 1.0 - theta[1];

   ratio[0] = theta[0] / complement[0];
   ratio[1] = theta[1] / complement[1];

   factor[0] = 1.0 + ratio[0];
   factor[1] = 1.0 + ratio[1];

   // Find out if some founders have only two children, in
   // which case additional speed-ups apply
   if (m.quick_bit_count)
      {
      QuickMoveAlong(m, rescale);
      return;
      }

   // The rescale constant must take hidden bits into accounts
   // if likelihoods are going to be used for something other
   // than comparing alternative gene flow patterns.
   male_scale.Dimension(m.male_bit_count + 1);
   male_scale[0] = pow(complement[1], m.male_bit_count + m.hidden_male_bit_count);
   for (int i = 1; i <= m.male_bit_count; i++)
      male_scale[i] = male_scale[i - 1] * factor[1];

   female_scale.Dimension(m.bit_count - m.male_bit_count + 1);
   female_scale[0] = pow(complement[0], m.bit_count + m.hidden_bit_count -
                                        m.male_bit_count - m.hidden_male_bit_count)
                     * rescale;
   for (int i = 1; i < female_scale.Length(); i++)
      female_scale[i] = female_scale[i - 1] * factor[0];

   for (int i = 0; i < m.two_f; i+=2)
      if (m.founder_bit_count[i])
         {
         // SetMeiosisSex(m.founder_male[i]);
         x = ratio[m.founder_male[i]];
         onex = factor[m.founder_male[i]];

         SpecialReUnite(0, 0, m.founder_bits[i]);
         }

   DivideAndConquer(&female_scale[m.bit_count - m.male_bit_count],
                    &male_scale[m.male_bit_count],
                    &m.bit_sex[m.bit_count]);

   SetMeiosisSex(0);

   for (int i = 0; i < m.couples; i++)
      if (m.couple_bits[i][m.two_n - 1])
         DoubleReUnite(0, 0, m.couple_bits[i]);
   }

void Multipoint::QuickMoveAlong(Mantra & m, double rescale)
   {
   // Increase apparent recombination rate for offspring of founders
   // with two children, after fixing one meiosis
   recombination[2] = 2 * complement[0] * recombination[0];
   recombination[3] = 2 * complement[1] * recombination[1];

   complement[2] = 1.0 - recombination[2];
   complement[3] = 1.0 - recombination[3];

   ratio[2] = recombination[2] / complement[2];
   ratio[3] = recombination[3] / complement[3];

   factor[2] = 1.0 + ratio[2];
   factor[3] = 1.0 + ratio[3];

   // Calculate rescaling constants ...
   double scale = pow(complement[1], m.male_bit_count + m.hidden_male_bit_count)
                * pow(complement[0], m.bit_count + m.hidden_bit_count -
                                     m.male_bit_count - m.hidden_male_bit_count)
                * rescale;

   // Process complicated founder symmetries and update meiosis types
   meiosis_type = m.bit_sex;

   for (int i = 0; i < m.two_f; i+=2)
      switch (m.founder_bit_count[i])
         {
         case 0 :
            break;
         case 1 :
            for (int j = 0; j < m.bit_count; j++)
               if (m.founder_bits[i][j])
                  {
                  j = m.bit_count - j;
                  if (meiosis_type[j] == 0)
                     {
                     meiosis_type[j] = 2;
                     scale *= complement[2] / (complement[0] * complement[0]);
                     }
                  else
                     {
                     meiosis_type[j] = 3;
                     scale *= complement[3] / (complement[1] * complement[1]);
                     }
                  break;
                  }
            break;
         default :
            x = ratio[m.founder_male[i]];
            onex = factor[m.founder_male[i]];

            SpecialReUnite(0, 0, m.founder_bits[i]);
         }

   // Handle regular bits, using increased recombination rates for
   // some hidden bits (as defined in the meiosis_type array).
   QuickDivideAndConquer(scale, &meiosis_type[m.bit_count], m.bit_count);

   SetMeiosisSex(0);

   for (int i = 0; i < m.couples; i++)
      if (m.couple_bits[i][m.two_n - 1])
         DoubleReUnite(0, 0, m.couple_bits[i]);
   }

///////////////////////////////////////////////////////////////////////
// This first set of routines allows for exact solutions and is the
// most general. DivideAndConquer() successively bissects inheritance
// trees, while ReUnite() joins transformed halves and SpecialReUnite()
// and DoubleReUnite() handle founder and founder couple symmetries
// respectively.
//

void Multipoint::DivideAndConquer(int splits_to_go, int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return;
      case TREE_NODE_LEAF :
         nodes[node].value *= pow_onex[splits_to_go];
         return;
      case TREE_NODE_ONE :
         DivideAndConquer(splits_to_go, nodes[node].child[0]);
         return;
      case TREE_NODE_TWO :
         ReUnite(nodes[node].child[0], nodes[node].child[1]);
         DivideAndConquer(splits_to_go - 1, nodes[node].child[0]);
         DivideAndConquer(splits_to_go - 1, nodes[node].child[1]);
         return;
      }
   }

void Multipoint::DivideAndConquer(double * femaleScale, double * maleScale,
                                  int * isMale, int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return;
      case TREE_NODE_LEAF :
         nodes[node].value *= (*femaleScale) * (*maleScale);
         return;
      case TREE_NODE_ONE :
         DivideAndConquer(femaleScale, maleScale, isMale - 1, nodes[node].child[0]);
         return;
      case TREE_NODE_TWO :
         // SetMeiosisSex(*isMale);
         x = ratio[*isMale];
         ReUnite(nodes[node].child[0], nodes[node].child[1]);
         if (*isMale)
            {
            DivideAndConquer(femaleScale, maleScale - 1, isMale - 1, nodes[node].child[0]);
            DivideAndConquer(femaleScale, maleScale - 1, isMale - 1, nodes[node].child[1]);
            }
         else
            {
            DivideAndConquer(femaleScale - 1, maleScale, isMale - 1, nodes[node].child[0]);
            DivideAndConquer(femaleScale - 1, maleScale, isMale - 1, nodes[node].child[1]);
            }
         return;
      }
   }

void Multipoint::QuickDivideAndConquer(double scale, int * meiosis_kind, int bit, int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return;
      case TREE_NODE_LEAF :
         while (bit)
            {
            scale *= factor[*meiosis_kind];
            meiosis_kind --;
            bit --;
            }
         nodes[node].value *= scale;
         return;
      case TREE_NODE_ONE :
         scale *= factor[*meiosis_kind];
         QuickDivideAndConquer(scale, meiosis_kind - 1, bit - 1, nodes[node].child[0]);
         return;
      case TREE_NODE_TWO :
         x = ratio[*meiosis_kind];
         ReUnite(nodes[node].child[0], nodes[node].child[1]);

         QuickDivideAndConquer(scale, meiosis_kind - 1, bit - 1, nodes[node].child[0]);
         QuickDivideAndConquer(scale, meiosis_kind - 1, bit - 1, nodes[node].child[1]);
         return;
      }
   }

void Multipoint::ReUnite(int left, int right)
   {
   int type_left  = nodes[left].type;
   int type_right = nodes[right].type;

   switch (PAIR_NODES(type_left, type_right))
      {
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ZERO) :
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_LEAF) :
         {
         double value_left = nodes[left].value;

         nodes[left].value  += nodes[right].value * x;
         nodes[right].value += value_left * x;
         }
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_LEAF) :
         nodes[left].type = TREE_NODE_LEAF;
         nodes[left].value = nodes[right].value * x;
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ZERO) :
         nodes[right].type = TREE_NODE_LEAF;
         nodes[right].value  = nodes[left].value * x;
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
         ReUnite(nodes[left].child[0], nodes[right].child[0]);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ZERO) :
         GraftAndScale(right, left, x);
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_TWO) :
         GraftAndScale(left, right, x);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_LEAF) :
         GraftAndScaleAndAdd(right, left, x,
               nodes[right].value, nodes[right].value * x);
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         GraftAndScaleAndAdd(left, right, x,
               nodes[left].value, nodes[left].value * x);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_TWO) :
         nodes[left].type = TREE_NODE_TWO;
         SAFE_SET(nodes[left].child[1], Copy(nodes[left].child[0]));
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ONE) :
         nodes[right].type = TREE_NODE_TWO;
         SAFE_SET(nodes[right].child[1], Copy(nodes[right].child[0]));
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_TWO) :
         break;
      }

   ReUnite(nodes[left].child[0], nodes[right].child[0]);
   ReUnite(nodes[left].child[1], nodes[right].child[1]);
   }

void Multipoint::SpecialReUnite(int left, int right, const int * flips)
   {
   int type_left  = nodes[left].type;
   int type_right = nodes[right].type;

   switch (PAIR_NODES(type_left, type_right))
      {
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ZERO) :
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_LEAF) :
         {
         double value_left = nodes[left].value;

         if (left != right)
            {
            nodes[left].value  += nodes[right].value * x;
            nodes[right].value += value_left * x;
            }
         else
            nodes[left].value  *= onex;
         }
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_LEAF) :
         nodes[left].type = TREE_NODE_LEAF;
         nodes[left].value = nodes[right].value * x;
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ZERO) :
         nodes[right].type = TREE_NODE_LEAF;
         nodes[right].value = nodes[left].value * x;
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
         SpecialReUnite(nodes[left].child[0], nodes[right].child[0], flips + 1);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ZERO) :
         FlipGraftAndScale(flips, right, left, x);
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_TWO) :
         FlipGraftAndScale(flips, left, right, x);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_LEAF) :
         FlipGraftAndScaleAndAdd(flips, right, left, x,
               nodes[right].value, nodes[right].value * x);
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         FlipGraftAndScaleAndAdd(flips, left, right, x,
               nodes[left].value, nodes[left].value * x);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_TWO) :
         nodes[left].type = TREE_NODE_TWO;
         SAFE_SET(nodes[left].child[1], Copy(nodes[left].child[0]));
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ONE) :
         nodes[right].type = TREE_NODE_TWO;
         SAFE_SET(nodes[right].child[1], Copy(nodes[right].child[0]));
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_TWO) :
         break;
      }

   if (*flips)
      {
      SpecialReUnite(nodes[left].child[1], nodes[right].child[0], flips + 1);
      if (right == left) return;
      SpecialReUnite(nodes[left].child[0], nodes[right].child[1], flips + 1);
      }
   else
      {
      SpecialReUnite(nodes[left].child[0], nodes[right].child[0], flips + 1);
      SpecialReUnite(nodes[left].child[1], nodes[right].child[1], flips + 1);
      }
   }

void Multipoint::DoubleReUnite(int left, int right, const int * flips)
   {
   int type_left  = nodes[left].type;
   int type_right = nodes[right].type;

   // If either side is a trivial node (eg. ZERO or LEAF) then ...
   switch (PAIR_NODES(type_left, type_right))
      {
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ZERO) :
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_LEAF) :
         {
         double value_left = nodes[left].value;

         if (left != right)
            {
            nodes[left].value  += nodes[right].value * x;
            nodes[right].value += value_left * x;
            }
         else
            nodes[left].value *= onex;
         }
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_LEAF) :
         nodes[left].type = TREE_NODE_LEAF;
         nodes[left].value = nodes[right].value * x;
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ZERO) :
         nodes[right].type = TREE_NODE_LEAF;
         nodes[right].value = nodes[left].value * x;
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ZERO) :
         GraftAndScale(right, left, x);
         DoubleFlipBranch(flips, right);
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_TWO) :
         GraftAndScale(left, right, x);
         DoubleFlipBranch(flips, left);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_LEAF) :
         GraftAndScaleAndAdd(right, left, x,
         nodes[right].value, nodes[right].value * x);
         DoubleFlipBranch(flips, right);
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         GraftAndScaleAndAdd(left, right, x,
         nodes[left].value, nodes[left].value * x);
         DoubleFlipBranch(flips, left);
         return;
      }

   // Otherwise, flips may need special handling
   if (*flips < 2)
      {
      switch (PAIR_NODES(type_left, type_right))
         {
         case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
            DoubleReUnite(nodes[left].child[0], nodes[right].child[0], flips + 1);
            return;
         case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_TWO) :
            nodes[left].type = TREE_NODE_TWO;
            SAFE_SET(nodes[left].child[1], Copy(nodes[left].child[0]));
            break;
         case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ONE) :
            nodes[right].type = TREE_NODE_TWO;
            SAFE_SET(nodes[left].child[1], Copy(nodes[right].child[0]));
            break;
         case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_TWO) :
            break;
         }

      DoubleReUnite(nodes[left].child[0], nodes[right].child[*flips], flips+1);
      if (*flips && right == left) return;
      DoubleReUnite(nodes[left].child[1], nodes[right].child[!*flips], flips+1);
      return;
      }

   // This is the hard part, since we must simultaneously consider 2 meiosis
   int right_child = nodes[right].child[0];
   int right_other = type_right == TREE_NODE_TWO ? nodes[right].child[1] : 0;
   int left_child  = nodes[left].child[0];
   int left_other  = type_left == TREE_NODE_TWO ? nodes[left].child[1] : 0;

   int left_level = type_left == TREE_NODE_TWO ||
                    nodes[right_child].type == TREE_NODE_TWO ||
                    right_other && nodes[right_other].type == TREE_NODE_TWO ?
                    TREE_NODE_TWO : TREE_NODE_ONE;

   int right_level = type_right == TREE_NODE_TWO ||
                     nodes[left_child].type == TREE_NODE_TWO ||
                     left_other && nodes[left_other].type == TREE_NODE_TWO ?
                     TREE_NODE_TWO : TREE_NODE_ONE;

   UpgradeNode(right_child, left_level);
   if (right_other) UpgradeNode(right_other, left_level);
   UpgradeNode(left_child, right_level);
   if (left_other) UpgradeNode(left_other, right_level);
   UpgradeNode(left, left_level);
   UpgradeNode(right, right_level);

   switch (PAIR_NODES(left_level, right_level))
      {
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
         DoubleReUnite(nodes[nodes[left].child[0]].child[0],
                       nodes[nodes[right].child[0]].child[0],
                       flips + 2);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_TWO) :
         DoubleReUnite(nodes[nodes[left].child[0]].child[0],
                       nodes[nodes[right].child[0]].child[0],
                       flips + 2);
         DoubleReUnite(nodes[nodes[left].child[0]].child[1],
                       nodes[nodes[right].child[1]].child[0],
                       flips + 2);
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ONE) :
         DoubleReUnite(nodes[nodes[left].child[0]].child[0],
                       nodes[nodes[right].child[0]].child[0],
                       flips + 2);
         DoubleReUnite(nodes[nodes[left].child[1]].child[0],
                       nodes[nodes[right].child[0]].child[1],
                       flips + 2);
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_TWO) :
         DoubleReUnite(nodes[nodes[left].child[0]].child[0],
                       nodes[nodes[right].child[0]].child[0],
                       flips + 2);
         DoubleReUnite(nodes[nodes[left].child[1]].child[1],
                       nodes[nodes[right].child[1]].child[1],
                       flips + 2);
         DoubleReUnite(nodes[nodes[left].child[0]].child[1],
                       nodes[nodes[right].child[1]].child[0],
                       flips + 2);
         if (left == right) return;
         DoubleReUnite(nodes[nodes[left].child[1]].child[0],
                       nodes[nodes[right].child[0]].child[1],
                       flips + 2);
         break;
      }
   }

///////////////////////////////////////////////////////////////////
// This set of routines allows for approximate solutions for dense
// maps (i.e. small theta). It provides equivalent functionality to
// the routines in the previous section.
//

void Multipoint::Approximation
    (Mantra & m, double distance, int order, double rescale)
   {
   theta = DistanceToRecombination(distance/order);
   oneminus = 1.0 - theta;
   x = theta / oneminus;
   onex = 1.0 + x;

   pow_onex.Dimension(bit_count + 1);
   pow_onex[0] = 1.0;
   for (int i = 1; i <= bit_count; i++)
      pow_onex[i] = pow_onex[i - 1] + x;

   TreeFlips original;

   Multiply(0, rescale);

   while (order--)
      {
      original.CarbonCopy(*this);

      DivideAndConquer(original, bit_count);

      for (int i = 0; i < m.two_f; i+=2)
         if (m.founder_bit_count[i])
            SpecialReUnite(original, 0, 0, m.founder_bits[i]);

      for (int i = 0; i < m.couples; i++)
         if (m.couple_bits[i][m.two_n - 1])
            DoubleReUnite(original, 0, 0, m.couple_bits[i]);
      }
   }

void Multipoint::Approximation
    (Mantra & m, double theta[2], int order, double rescale)
   {
   if (order > 1)
      {
      theta[0] = 0.5 - 0.5 * pow(1.0 - 2.0 * theta[0], 1.0 / order);
      theta[1] = 0.5 - 0.5 * pow(1.0 - 2.0 * theta[1], 1.0 / order);
      }

   // Setup sex specific recombination fraction and related quantities
   recombination[0] = theta[0];
   recombination[1] = theta[1];

   complement[0] = 1.0 - theta[0];
   complement[1] = 1.0 - theta[1];

   ratio[0] = theta[0] / complement[0];
   ratio[1] = theta[1] / complement[1];

   factor[0] = 1.0 + ratio[0];
   factor[1] = 1.0 + ratio[1];

   // The rescale constant must take hidden bits into accounts
   // if likelihoods are going to be used for something other
   // than comparing alternative gene flow patterns.
   male_scale.Dimension(m.male_bit_count + 1);
   male_scale[0] = 1.0;
   for (int i = 1; i <= m.male_bit_count; i++)
      male_scale[i] = male_scale[i - 1] + ratio[1];

   female_scale.Dimension(m.bit_count - m.male_bit_count + 1);
   female_scale[0] = 0.0;
   for (int i = 1; i < female_scale.Length(); i++)
      female_scale[i] = female_scale[i - 1] + ratio[0];

   TreeFlips original;

   Multiply(0, rescale);

   while (order--)
      {
      original.CarbonCopy(*this);

      DivideAndConquer(original,
                       &female_scale[m.bit_count - m.male_bit_count],
                       &male_scale[m.male_bit_count],
                       &m.bit_sex[m.bit_count]);

      for (int i = 0; i < m.two_f; i+=2)
         if (m.founder_bit_count[i])
            {
            SetMeiosisSex(m.founder_male[i]);
            SpecialReUnite(original, 0, 0, m.founder_bits[i]);
            }

      SetMeiosisSex(0);
      for (int i = 0; i < m.couples; i++)
         if (m.couple_bits[i][m.two_n - 1])
            DoubleReUnite(original, 0, 0, m.couple_bits[i]);
      }
   }

void Multipoint::DivideAndConquer(Tree & original, int splits_to_go, int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return;
      case TREE_NODE_LEAF :
         nodes[node].value *= pow_onex[splits_to_go];
         return;
      case TREE_NODE_ONE :
         DivideAndConquer(original, splits_to_go, nodes[node].child[0]);
         return;
      case TREE_NODE_TWO :
         DivideAndConquer(original, splits_to_go - 1, nodes[node].child[0]);
         DivideAndConquer(original, splits_to_go - 1, nodes[node].child[1]);
         ReUnite(original, nodes[node].child[0], nodes[node].child[1]);
         ReUnite(original, nodes[node].child[1], nodes[node].child[0]);
         return;
      }
   }

void Multipoint::DivideAndConquer(Tree & tree,
    double * femaleScale, double * maleScale, int * isMale, int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return;
      case TREE_NODE_LEAF :
         nodes[node].value *= (*femaleScale) + (*maleScale);
         return;
      case TREE_NODE_ONE :
         DivideAndConquer(tree, femaleScale, maleScale, isMale - 1, nodes[node].child[0]);
         return;
      case TREE_NODE_TWO :
         if (*isMale)
            {
            DivideAndConquer(tree, femaleScale, maleScale - 1, isMale - 1, nodes[node].child[0]);
            DivideAndConquer(tree, femaleScale, maleScale - 1, isMale - 1, nodes[node].child[1]);
            }
         else
            {
            DivideAndConquer(tree, femaleScale - 1, maleScale, isMale - 1, nodes[node].child[0]);
            DivideAndConquer(tree, femaleScale - 1, maleScale, isMale - 1, nodes[node].child[1]);
            }
         SetMeiosisSex(*isMale);
         ReUnite(tree, nodes[node].child[0], nodes[node].child[1]);
         ReUnite(tree, nodes[node].child[1], nodes[node].child[0]);
         return;
      }
   }

void Multipoint::ReUnite(Tree & original, int left, int right)
   {
   int type_left  = nodes[left].type;
   int type_right = original.nodes[right].type;

   switch (PAIR_NODES(type_left, type_right))
      {
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ZERO) :
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_LEAF) :
         Add(left, original.nodes[right].value * x);
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_TWO) :
         GraftAndScale(left, original, right, x);
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         GraftAndScaleAndAdd(left, original, right, x, nodes[left].value);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
         ReUnite(original, nodes[left].child[0], original.nodes[right].child[0]);
         return;
      case PAIR_NODES(TREE_NODE_TWO , TREE_NODE_ONE) :
         ReUnite(original, nodes[left].child[0], original.nodes[right].child[0]);
         ReUnite(original, nodes[left].child[1], original.nodes[right].child[0]);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_TWO) :
         nodes[left].type = TREE_NODE_TWO;
         SAFE_SET(nodes[left].child[1], Copy(nodes[left].child[0]));
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_TWO) :
         ReUnite(original, nodes[left].child[0], original.nodes[right].child[0]);
         ReUnite(original, nodes[left].child[1], original.nodes[right].child[1]);
      }
   }

void Multipoint::SpecialReUnite
   (Tree & original, int left, int right, const int * flips)
   {
   int type_left  = nodes[left].type;
   int type_right = original.nodes[right].type;

   switch (PAIR_NODES(type_left, type_right))
      {
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ZERO) :
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_LEAF) :
         Add(left, original.nodes[right].value * x);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
         SpecialReUnite(original, nodes[left].child[0], original.nodes[right].child[0], flips + 1);
         return;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ONE) :
         SpecialReUnite(original, nodes[left].child[0], original.nodes[right].child[0], flips + 1);
         SpecialReUnite(original, nodes[left].child[1], original.nodes[right].child[0], flips + 1);
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_TWO) :
         FlipGraftAndScale(flips, left, original, right, x);
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         FlipGraftAndScaleAndAdd(flips, left, original, right, x, nodes[left].value);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_TWO) :
         nodes[left].type = TREE_NODE_TWO;
         SAFE_SET(nodes[left].child[1], Copy(nodes[left].child[0]));
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_TWO) :
         ;
      }

   SpecialReUnite(original, nodes[left].child[0], original.nodes[right].child[*flips], flips + 1);
   SpecialReUnite(original, nodes[left].child[1], original.nodes[right].child[!*flips], flips + 1);
   }

void Multipoint::DoubleReUnite(TreeFlips & original, int left, int right, const int *flips)
   {
   int type_left  = nodes[left].type;
   int type_right = original.nodes[right].type;

   switch (PAIR_NODES(type_left, type_right))
      {
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ZERO) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ZERO) :
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_LEAF) :
         Add(left, original.nodes[right].value * x);
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_TWO) :
         GraftAndScale(left, original, right, x);
         DoubleFlipBranch(flips, left);
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         GraftAndScaleAndAdd(left, original, right, x, nodes[left].value);
         DoubleFlipBranch(flips, left);
         return;
      }

   if (*flips < 2)
      switch (PAIR_NODES(type_left, type_right))
         {
         case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
            DoubleReUnite(original, nodes[left].child[0],
                          original.nodes[right].child[0], flips + 1);
            return;
         case PAIR_NODES(TREE_NODE_TWO , TREE_NODE_ONE) :
            DoubleReUnite(original, nodes[left].child[0],
                          original.nodes[right].child[0], flips + 1);
            DoubleReUnite(original, nodes[left].child[1],
                          original.nodes[right].child[0], flips + 1);
            return;
         case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_TWO) :
            nodes[left].type = TREE_NODE_TWO;
            SAFE_SET(nodes[left].child[1], Copy(nodes[left].child[0]));
         case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_TWO) :
            DoubleReUnite(original, nodes[left].child[0],
                          original.nodes[right].child[*flips], flips + 1);
            DoubleReUnite(original, nodes[left].child[1],
                          original.nodes[right].child[!*flips], flips + 1);
            return;
         }

   // This is the hard part ... need to consider two meiosis at once!
   int right_child = original.nodes[right].child[0];
   int right_other = type_right == TREE_NODE_TWO ? original.nodes[right].child[1] : 0;
   int left_child  = nodes[left].child[0];
   int left_other  = type_left == TREE_NODE_TWO ? nodes[left].child[1] : 0;

   int left_level = type_left == TREE_NODE_TWO ||
                    original.nodes[right_child].type == TREE_NODE_TWO ||
                    right_other && original.nodes[right_other].type == TREE_NODE_TWO ?
                    TREE_NODE_TWO : TREE_NODE_ONE;

   int right_level = type_right == TREE_NODE_TWO ||
                     nodes[left_child].type == TREE_NODE_TWO ||
                     left_other && nodes[left_other].type == TREE_NODE_TWO ?
                     TREE_NODE_TWO : TREE_NODE_ONE;

   UpgradeNode(left_child, right_level);
   if (left_other) UpgradeNode(left_other, right_level);
   UpgradeNode(left, left_level);

   // As far as possible, leave original tree unchanged...
   original.PseudoUpgradeNode(right_child, left_level);
   if (right_other) original.PseudoUpgradeNode(right_other, left_level);
   original.PseudoUpgradeNode(right, right_level);

   switch (PAIR_NODES(left_level, right_level))
      {
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
         DoubleReUnite(original, nodes[nodes[left].child[0]].child[0],
                       original.nodes[original.nodes[right].child[0]].child[0],
                       flips + 2);
         break;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_TWO) :
         DoubleReUnite(original, nodes[nodes[left].child[0]].child[0],
                       original.nodes[original.nodes[right].child[0]].child[0],
                       flips + 2);
         DoubleReUnite(original, nodes[nodes[left].child[0]].child[1],
                       original.nodes[original.nodes[right].child[1]].child[0],
                       flips + 2);
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ONE) :
         DoubleReUnite(original, nodes[nodes[left].child[0]].child[0],
                       original.nodes[original.nodes[right].child[0]].child[0],
                       flips + 2);
         DoubleReUnite(original, nodes[nodes[left].child[1]].child[0],
                       original.nodes[original.nodes[right].child[0]].child[1],
                       flips + 2);
         break;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_TWO) :
         DoubleReUnite(original, nodes[nodes[left].child[0]].child[0],
                       original.nodes[original.nodes[right].child[0]].child[0],
                       flips + 2);
         DoubleReUnite(original, nodes[nodes[left].child[1]].child[1],
                       original.nodes[original.nodes[right].child[1]].child[1],
                       flips + 2);
         DoubleReUnite(original, nodes[nodes[left].child[0]].child[1],
                       original.nodes[original.nodes[right].child[1]].child[0],
                       flips + 2);
         DoubleReUnite(original, nodes[nodes[left].child[1]].child[0],
                       original.nodes[original.nodes[right].child[0]].child[1],
                       flips + 2);
         break;
      }
   }

///////////////////////////////////////////////////////////////////
// This final set does not provide additional functionality but
// provides for faster calculation when Markov chains are walking
// through a set of highly informative markers.
//

void Multipoint::SetupConditioning(Mantra & m)
   {
   bitSet.BuildSets(m);

   current.Dimension(m.bit_count);
   previous.Dimension(m.bit_count);
   }

void Multipoint::Condition(Tree & flanker, double distance, double rescale)
   {
   if (distance == 0.0 || bit_count == 0)
      {
      flanker.Multiply(*this);
      flanker.Multiply(0, rescale);
      return;
      }

   theta = DistanceToRecombination(distance);
   recombination[0] = recombination[1] = theta;

   Condition(flanker, bit_count - 1, rescale);
   }

void Multipoint::Condition(Tree & flanker, double theta[2], double rescale)
   {
   if (theta[0] == 0.0 && theta[1] == 0.0 || bit_count == 0)
      {
      flanker.Multiply(*this);
      flanker.Multiply(0, rescale);
      return;
      }

   recombination[0] = theta[0];
   recombination[1] = theta[1];

   Condition(flanker, bit_count - 1, rescale);
   }

void Multipoint::Condition(Tree & tree, int bit, double scale, int node)
   {
   switch (tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return;
      case TREE_NODE_LEAF :
         if (bit >= 0)
            {
            double value = tree.nodes[node].value;

            tree.nodes[node].type = TREE_NODE_TWO;
            SAFE_SET(tree.nodes[node].child[0], tree.NewNode());
            SAFE_SET(tree.nodes[node].child[1], tree.NewNode());

            tree.nodes[tree.nodes[node].child[0]].type = TREE_NODE_LEAF;
            tree.nodes[tree.nodes[node].child[1]].type = TREE_NODE_LEAF;
            tree.nodes[tree.nodes[node].child[0]].value = value;
            tree.nodes[tree.nodes[node].child[1]].value = value;
            }
         else
            {
            tree.nodes[node].value *= ConditionalLikelihood(bit_count - 1, scale);
            return;
            }
         break;
      case TREE_NODE_ONE :
         tree.nodes[node].type = TREE_NODE_TWO;
         SAFE_SET(tree.nodes[node].child[1], tree.Copy(tree.nodes[node].child[0]));
         break;
      case TREE_NODE_TWO :
         break;
      }

   // For internal nodes, just need to condition subtrees
   previous[bit] = 0;
   Condition(tree, bit - 1, scale, tree.nodes[node].child[0]);

   previous[bit] = 1;
   Condition(tree, bit - 1, scale, tree.nodes[node].child[1]);
   }

double Multipoint::ConditionalLikelihood(int bit, double lk, int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return 0.0;
      case TREE_NODE_LEAF :
         while (bit >= 0)
            {
            current[bit] = 2;

            bitSet.UpdateLikelihood(bit, recombination, lk, current, previous);
            bit--;
            }
         return lk * nodes[node].value;
      case TREE_NODE_ONE :
         {
         current[bit] = 2;

         bitSet.UpdateLikelihood(bit, recombination, lk, current, previous);
         return ConditionalLikelihood(bit - 1, lk, nodes[node].child[0]);
         }
      case TREE_NODE_TWO :
         {
         double lk2 = lk;

         // First search the more likely path (with no recombinants)
         current[bit] = previous[bit];
         bitSet.UpdateLikelihood(bit, recombination, lk, current, previous);
         lk = ConditionalLikelihood(bit - 1, lk, nodes[node].child[previous[bit]]);

         // Now search the other path
         current[bit] = !previous[bit];
         bitSet.UpdateLikelihood(bit, recombination, lk2, current, previous);
         lk2= ConditionalLikelihood(bit - 1, lk2, nodes[node].child[!previous[bit]]);

         return lk + lk2;
         }
      }
   // Should never get here, but keeps compiler happy!
   return 0.0;
   }


 
