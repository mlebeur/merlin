////////////////////////////////////////////////////////////////////// 
// merlin/ConquerHaplotyping.cpp 
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
 
#include "ConquerHaplotyping.h"
#include "MapFunction.h"
#include "MathConstant.h"

void MultipointHaplotyping::MoveAlong(Mantra & m, double distance, double rescale)
   {
   if (distance == 0.0)
      return;

   if (Multipoint::maximum_recombinants)
      {
      Approximation(m, distance, Multipoint::maximum_recombinants, rescale);
      return;
      }

   theta = DistanceToRecombination(distance);
   oneminus = 1.0 - theta;
   x = theta / oneminus;
   scale = rescale;

   for (int i = 0; i < m.two_f; i+=2)
      if (m.founder_bit_count[i])
         SpecialReUnite(0, 0, m.founder_bits[i]);

   DivideAndConquer();

   for (int i = 0; i < m.couples; i++)
      if (m.couple_bits[i][m.two_n - 1])
         DoubleReUnite(0, 0, m.couple_bits[i]);
   }

void MultipointHaplotyping::MoveAlong(Mantra & m, double theta[2], double rescale)
   {
   if (theta[0] == 0.0 && theta[1] == 0)
      return;

   if (Multipoint::maximum_recombinants)
      {
      Approximation(m, theta, Multipoint::maximum_recombinants, rescale);
      return;
      }

   recombination[0] = theta[0];
   recombination[1] = theta[1];

   complement[0] = 1.0 - theta[0];
   complement[1] = 1.0 - theta[1];

   ratio[0] = theta[0] / complement[0];
   ratio[1] = theta[1] / complement[1];

   scale = rescale;

   for (int i = 0; i < m.two_f; i+=2)
      if (m.founder_bit_count[i])
         {
         SetMeiosisSex(m.founder_male[i]);
         SpecialReUnite(0, 0, m.founder_bits[i]);
         }

   DivideAndConquer(&m.bit_sex[m.bit_count]);

   SetMeiosisSex(0);
   for (int i = 0; i < m.couples; i++)
      if (m.couple_bits[i][m.two_n - 1])
         DoubleReUnite(0, 0, m.couple_bits[i]);
   }

void MultipointHaplotyping::DivideAndConquer(int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return;
      case TREE_NODE_LEAF :
         nodes[node].value *= scale;
         return;
      case TREE_NODE_ONE :
         DivideAndConquer(nodes[node].child[0]);
         return;
      case TREE_NODE_TWO :
         ReUnite(nodes[node].child[0], nodes[node].child[1]);
         DivideAndConquer(nodes[node].child[0]);
         DivideAndConquer(nodes[node].child[1]);
         return;
      }
   }

void MultipointHaplotyping::DivideAndConquer(int * isMale, int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return;
      case TREE_NODE_LEAF :
         nodes[node].value *= scale;
         return;
      case TREE_NODE_ONE :
         DivideAndConquer(isMale - 1, nodes[node].child[0]);
         return;
      case TREE_NODE_TWO :
         SetMeiosisSex(*isMale);
         ReUnite(nodes[node].child[0], nodes[node].child[1]);
         DivideAndConquer(isMale - 1, nodes[node].child[0]);
         DivideAndConquer(isMale - 1, nodes[node].child[1]);
         return;
      }
   }

void MultipointHaplotyping::ReUnite(int left, int right)
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

         nodes[left].value  = max(value_left, nodes[right].value * x);
         nodes[right].value = max(nodes[right].value, value_left * x);
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
         GraftAndScaleAndChoose(right, left, x,
               nodes[right].value, nodes[right].value * x);
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         GraftAndScaleAndChoose(left, right, x,
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

void MultipointHaplotyping::SpecialReUnite(int left, int right, const int * flips)
   {
   int type_left  = nodes[left].type;
   int type_right = nodes[right].type;

   switch (PAIR_NODES(type_left, type_right))
      {
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ZERO) :
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_LEAF) :
         if (left != right)
            {
            double value_left = nodes[left].value;

            nodes[left].value  = max(value_left, nodes[right].value * x);
            nodes[right].value = max(nodes[right].value, value_left * x);
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
         FlipGraftAndScaleAndChoose(flips, right, left, x,
               nodes[right].value, nodes[right].value * x);
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         FlipGraftAndScaleAndChoose(flips, left, right, x,
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

void MultipointHaplotyping::DoubleReUnite(int left, int right, const int * flips)
   {
   int type_left  = nodes[left].type;
   int type_right = nodes[right].type;

   // If either side is a trivial node (eg. ZERO or LEAF) then ...
   switch (PAIR_NODES(type_left, type_right))
      {
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ZERO) :
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_LEAF) :
         if (left != right)
            {
            double value_left = nodes[left].value;

            nodes[left].value  = max(value_left, nodes[right].value * x);
            nodes[right].value = max(nodes[right].value, value_left * x);
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
         GraftAndScaleAndChoose(right, left, x,
         nodes[right].value, nodes[right].value * x);
         DoubleFlipBranch(flips, right);
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         GraftAndScaleAndChoose(left, right, x,
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
            SAFE_SET(nodes[right].child[1], Copy(nodes[right].child[0]));
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

void MultipointHaplotyping::Approximation
    (Mantra & m, double distance, int order, double rescale)
   {
   theta = DistanceToRecombination(distance/order);
   oneminus = 1.0 - theta;
   x = theta / oneminus;

   TreeFlips original;

   Multiply(0, rescale);

   while (order--)
      {
      original.CarbonCopy(*this);

      DivideAndConquer(original);

      for (int i = 0; i < m.two_f; i+=2)
         if (m.founder_bit_count[i])
            SpecialReUnite(original, 0, 0, m.founder_bits[i]);

      for (int i = 0; i < m.couples; i++)
         if (m.couple_bits[i][m.two_n - 1])
            DoubleReUnite(original, 0, 0, m.couple_bits[i]);
      }
   }

void MultipointHaplotyping::Approximation
    (Mantra & m, double theta[2], int order, double rescale)
   {
   if (order > 1)
      {
      theta[0] = 0.5 - 0.5 * pow(1.0 - 2.0 * theta[0], 1.0 / order);
      theta[1] = 0.5 - 0.5 * pow(1.0 - 2.0 * theta[1], 1.0 / order);
      }

   recombination[0] = theta[0];
   recombination[1] = theta[1];

   complement[0] = 1.0 - theta[0];
   complement[1] = 1.0 - theta[1];

   ratio[0] = theta[0] / complement[0];
   ratio[1] = theta[1] / complement[1];

   TreeFlips original;

   Multiply(0, rescale);

   while (order--)
      {
      original.CarbonCopy(*this);

      DivideAndConquer(original, &m.bit_sex[m.bit_count]);

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

void MultipointHaplotyping::DivideAndConquer(Tree & original, int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
      case TREE_NODE_LEAF :
         return;
      case TREE_NODE_ONE :
         DivideAndConquer(original, nodes[node].child[0]);
         return;
      case TREE_NODE_TWO :
         DivideAndConquer(original, nodes[node].child[0]);
         DivideAndConquer(original, nodes[node].child[1]);
         ReUnite(original, nodes[node].child[0], nodes[node].child[1]);
         ReUnite(original, nodes[node].child[1], nodes[node].child[0]);
         return;
      }
   }

void MultipointHaplotyping::DivideAndConquer(Tree & original, int * isMale, int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
      case TREE_NODE_LEAF :
         return;
      case TREE_NODE_ONE :
         DivideAndConquer(original, isMale - 1, nodes[node].child[0]);
         return;
      case TREE_NODE_TWO :
         DivideAndConquer(original, isMale - 1, nodes[node].child[0]);
         DivideAndConquer(original, isMale - 1, nodes[node].child[1]);
         SetMeiosisSex(*isMale);
         ReUnite(original, nodes[node].child[0], nodes[node].child[1]);
         ReUnite(original, nodes[node].child[1], nodes[node].child[0]);
         return;
      }
   }

void MultipointHaplotyping::ReUnite(Tree & original, int left, int right)
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
         Floor(left, original.nodes[right].value * x);
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_LEAF) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_TWO) :
         GraftAndScale(left, original, right, x);
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         GraftAndScaleAndChoose(left, original, right, x, nodes[left].value);
         return;
      case PAIR_NODES(TREE_NODE_ONE, TREE_NODE_ONE) :
         ReUnite(original, nodes[left].child[0], original.nodes[right].child[0]);
         return;
      case PAIR_NODES(TREE_NODE_TWO, TREE_NODE_ONE) :
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

void MultipointHaplotyping::SpecialReUnite
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
         Floor(left, original.nodes[right].value * x);
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
         FlipGraftAndScaleAndChoose(flips, left, original, right, x, nodes[left].value);
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

void MultipointHaplotyping::DoubleReUnite(TreeFlips & original, int left, int right, const int *flips)
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
         Floor(left, original.nodes[right].value * x);
         return;
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_ZERO, TREE_NODE_TWO) :
         GraftAndScale(left, original, right, x);
         DoubleFlipBranch(flips, left);
         return;
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_ONE) :
      case PAIR_NODES(TREE_NODE_LEAF, TREE_NODE_TWO) :
         GraftAndScaleAndChoose(left, original, right, x, nodes[left].value);
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
         return;
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

void MultipointHaplotyping::SetupConditioning(Mantra & m)
   {
   bitSet.BuildSets(m);

   current.Dimension(m.bit_count);
   previous.Dimension(m.bit_count);
   }


void MultipointHaplotyping::Condition(Tree & flanker, double distance, double rescale)
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

void MultipointHaplotyping::Condition(Tree & flanker, double theta[2], double rescale)
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

void MultipointHaplotyping::Condition(Tree & tree, int bit, double scaling, int node)
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
            tree.nodes[node].value *= ConditionalLikelihood(bit_count - 1, scaling);
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
   Condition(tree, bit - 1, scaling, tree.nodes[node].child[0]);

   previous[bit] = 1;
   Condition(tree, bit - 1, scaling, tree.nodes[node].child[1]);
   }

double MultipointHaplotyping::ConditionalLikelihood(int bit, double lk, int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return 0.0;
      case TREE_NODE_LEAF :
         while (bit >= 0)
            {
            current[bit] = 2;

            bitSet.UpdateMaxLikelihood(bit, recombination, lk, current, previous);
            bit--;
            }
         return lk * nodes[node].value;
      case TREE_NODE_ONE :
         {
         current[bit] = 2;

         bitSet.UpdateMaxLikelihood(bit, recombination, lk, current, previous);
         return ConditionalLikelihood(bit - 1, lk, nodes[node].child[0]);
         }
      case TREE_NODE_TWO :
         {
         double lk2 = lk;

         // First search the more likely path (with no recombinants)
         current[bit] = previous[bit];
         bitSet.UpdateMaxLikelihood(bit, recombination, lk, current, previous);
         lk = ConditionalLikelihood(bit - 1, lk, nodes[node].child[previous[bit]]);

         // Now search the other path
         current[bit] = !previous[bit];
         bitSet.UpdateMaxLikelihood(bit, recombination, lk2, current, previous);
         lk2= ConditionalLikelihood(bit - 1, lk2, nodes[node].child[!previous[bit]]);

         return lk > lk2 ? lk : lk2;
         }
      }
   // Should never get here, but keeps compiler happy!
   return 0.0;
   }


 
