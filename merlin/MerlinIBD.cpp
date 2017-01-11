////////////////////////////////////////////////////////////////////// 
// merlin/MerlinIBD.cpp 
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
 
#include "MerlinIBD.h"
#include "MerlinCore.h"

#include <math.h>

// The routines in this section perform IBD calculations
//    Calculate(...) calculates and outputs IBD state probabilities
//                   returning the overall likelihood of all vectors in the tree
//    Score(...)  traverses an inheritance tree and fills IBD1 and IBD2
//                   with relative sharing probabities
//    Output(...) rescales sharing probabilities to the 0.0 .. 1.0 range
//                   and stores them in a file
//    Open(...) and Close(...) manage the IBD file
//

double MerlinIBD::Calculate(Tree & tree, const char * label)
   {
   // Dimension IBD matrices
   IBD1.Dimension(mantra.n, mantra.n);
   IBD2.Dimension(mantra.n, mantra.n);
   IBD1.Zero();
   IBD2.Zero();

   // Calculate IBDs as well as the overall likelihood
   double likelihood  = Score(tree, 0, tree.bit_count - 1);

   Output(label, likelihood);

   return likelihood * pow(2.0, -mantra.bit_count);
   }

double MerlinIBD::Score(Tree & tree, int node, int bit, int start)
   {
   double sum = 0.0; // Initialization avoids compiler warning
   int pivot, end;

   if (bit != mantra.bit_count - 1)
      {
      int pivot = bit == -1 ? mantra.two_n : mantra.bits[bit];
      int last_pivot = mantra.bits[bit+1] + 1;

      for (int i = last_pivot; i < pivot; i++)
         mantra.state[i] = mantra.state[mantra.vector[i]];
      }

   if (bit >= 0)
      {
      pivot = mantra.bits[bit];
      end = pivot & ~1;

      switch (tree.nodes[node].type)
         {
         case TREE_NODE_ZERO :
            return 0.0;
         case TREE_NODE_LEAF :
            uninformative.Dimension(mantra.two_n);
            uninformative.Zero();
            uninformative[pivot] = true;

            mantra.state[pivot] = mantra.state[mantra.vector[pivot] &= ~1];
            sum = WeightedScore(tree, node, bit - 1, end);
            break;
         case TREE_NODE_ONE :
            uninformative.Dimension(mantra.two_n);
            uninformative.Zero();
            uninformative[pivot] = true;

            mantra.state[pivot] = mantra.state[mantra.vector[pivot] &= ~1];
            sum = WeightedScore(tree, tree.nodes[node].child[0], bit - 1, end);
            break;
         case TREE_NODE_TWO :
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] &= ~1];
            sum = Score(tree, tree.nodes[node].child[0], bit - 1, end);
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] |= 1];
            sum += Score(tree, tree.nodes[node].child[1],bit - 1, end);
            break;
         }
      }
   else
      {
      // At the base of the tree there are only zero nodes
      // and leaf nodes ...

      if (tree.nodes[node].type == TREE_NODE_ZERO)
         return 0.0;

      end = mantra.two_n;
      sum = tree.nodes[node].value;
      }


   for ( int i = start, ii = start >> 1 ; i < end; i++, i++, ii++)
      for ( int j = 0, jj = 0; j <= i; j++, j++, jj++)
         if (mantra.ibd[ii][jj] == MANTRA_IBD_UNKNOWN)
            {
            if (mantra.state[i]   == mantra.state[j]   &&
                mantra.state[i+1] == mantra.state[j+1] ||
                mantra.state[i]   == mantra.state[j+1] &&
                mantra.state[i+1] == mantra.state[j]   )
                IBD2[ii][jj] += sum;
            else
            if (mantra.state[i]   == mantra.state[j]   ||
                mantra.state[i+1] == mantra.state[j+1] ||
                mantra.state[i]   == mantra.state[j+1] ||
                mantra.state[i+1] == mantra.state[j])
               IBD1[ii][jj] += sum;
            }

   return sum;
   }

void MerlinIBD::Output(const char * label, double scale)
   {
   scale = 1.0 / scale;

   for (int i = 0; i < mantra.n; i++)
      for (int j = 0; j <= i; j++)
         {
         fprintf(ibdfile, "%s %s %s %s ",
                (const char *) family->famid,
                (const char *) (*ped)[family->path[i]].pid,
                (const char *) (*ped)[family->path[j]].pid,
                (const char *) label);

         switch (mantra.ibd[i][j])
            {
            case MANTRA_IBD_ZERO :
               fprintf(ibdfile, " 1.0 0.0 0.0\n");
               break;
            case MANTRA_IBD_HALF :
            case MANTRA_IBD_HALF_MALE :
               fprintf(ibdfile, " 0.0 1.0 0.0\n");
               break;
            case MANTRA_IBD_ONE :
            case MANTRA_IBD_ONE_MALE :
               fprintf(ibdfile, " 0.0 0.0 1.0\n");
               break;
            default :
               {
               double p2 = IBD2[i][j] * scale;
               double p1 = IBD1[i][j] * scale;

               if (mantra.couples && j < mantra.f && mantra.couple_index[j]!=-1)
                  {
                  p2 = 0.5 * (p2 + IBD2[i][mantra.partners[j]] * scale);
                  p1 = 0.5 * (p1 + IBD1[i][mantra.partners[j]] * scale);
                  }

               double p0 = 1.0 - p1 - p2;

               fprintf(ibdfile, " %.5f %.5f %.5f\n", p0, p1, p2);
               }
            }
         }
   }

void MerlinIBD::OpenFile()
   {
   String filename(MerlinCore::filePrefix);
   filename += ".ibd";

   ibdfile = fopen(filename, "wt");

   if (ibdfile == NULL)
      error("Error opening file [%s]", (const char *) filename);

   fprintf(ibdfile, "FAMILY ID1 ID2 MARKER P0 P1 P2\n");
   }

void MerlinIBD::CloseFile()
   {
   fclose(ibdfile);
   printf("IBD probabilities stored in file [%s.ibd]\n", (const char *) MerlinCore::filePrefix);
   }

void MerlinIBD::SelectFamily(Pedigree * p, Family * f)
   {
   ped = p;
   family = f;
   }

double MerlinIBD::PairWiseIBD1(int a, int b, int c, int d)
   {
   // First sort alleles so a is the bottom-most one
   int temp;

   if (a < b) { temp = a; a = b; b = temp; }
   if (c < d) { temp = c; c = d; d = temp; }
   if (a < c) { temp = a; a = c; c = temp;
                temp = b; b = d; d = temp; }

   // If all alleles are founder alleles
   if (a < mantra.two_f)
      return (a == c) ? (b != d) : (a == d || b == c || b == d);

   int parent = mantra.vector[a];

   if (uninformative[a])
      return 0.5 * ( PairWiseIBD1(parent, b == a ? parent : b,
                     c == a ? parent : c, d == a ? parent : d) +
                     PairWiseIBD1(parent + 1, b == a ? parent + 1 : b,
                     c == a ? parent + 1 : c, d == a ? parent + 1 : d));
   else
      return PairWiseIBD1(parent, b == a ? parent : b,
             c == a ? parent : c, d == a ? parent : d);
   }

double MerlinIBD::PairWiseIBD2(int a, int b, int c, int d)
   {
   // First sort alleles so a is the bottom-most one
   int temp;

   if (a < b) { temp = a; a = b; b = temp; }
   if (c < d) { temp = c; c = d; d = temp; }
   if (a < c) { temp = a; a = c; c = temp;
                temp = b; b = d; d = temp; }

   // Check for IBD2
   if (a == c && b == d)
      return 1.0;

   // Else, if all alleles are founder alleles
   if (a < mantra.two_f)
      return 0.0;

   int parent = mantra.vector[a];

   if (uninformative[a])
      return 0.5 * ( PairWiseIBD2(parent, b == a ? parent : b,
                     c == a ? parent : c, d == a ? parent : d) +
                     PairWiseIBD2(parent + 1, b == a ? parent + 1 : b,
                     c == a ? parent + 1 : c, d == a ? parent + 1 : d));
   else
      return PairWiseIBD2(parent, b == a ? parent : b,
             c == a ? parent : c, d == a ? parent : d);
   }

double MerlinIBD::WeightedScore(Tree & tree, int node, int bit, int start,
                        double weight)
   {
   double sum = 0.0; // Initialization avoids compiler warning
   int pivot, end;

   if (bit != mantra.bit_count - 1)
      {
      int pivot = bit == -1 ? mantra.two_n : mantra.bits[bit];
      int last_pivot = mantra.bits[bit+1] + 1;

      for (int i = last_pivot; i < pivot; i++)
         mantra.state[i] = mantra.state[mantra.vector[i]];
      }

   if (bit >= 0)
      {
      pivot = mantra.bits[bit];
      end = pivot & ~1;

      switch (tree.nodes[node].type)
         {
         case TREE_NODE_ZERO :
            return 0.0;
         case TREE_NODE_LEAF :
            uninformative[pivot] = true;
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] &= ~1];
            sum = WeightedScore(tree, node, bit - 1, end, weight * 2.0);
            break;
         case TREE_NODE_ONE :
            uninformative[pivot] = true;
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] &= ~1];
            sum = WeightedScore(tree, tree.nodes[node].child[0], bit - 1, end,
                  weight * 2.0);
            break;
         case TREE_NODE_TWO :
            uninformative[pivot] = false;
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] &= ~1];
            sum = WeightedScore(tree, tree.nodes[node].child[0], bit - 1, end,
                  weight);
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] |= 1];
            sum += WeightedScore(tree, tree.nodes[node].child[1],bit - 1, end,
                   weight);
            break;
         }
      }
   else
      {
      // At the base of the tree there are only zero nodes
      // and leaf nodes ...

      if (tree.nodes[node].type == TREE_NODE_ZERO)
         return 0.0;

      end = mantra.two_n;
      sum = tree.nodes[node].value * weight;
      }

   for ( int i = start, ii = start >> 1 ; i < end; i++, i++, ii++)
      for ( int j = 0, jj = 0; j <= i; j++, j++, jj++)
         if (mantra.ibd[ii][jj] == MANTRA_IBD_UNKNOWN)
            {
#ifndef __CHROMOSOME_X__
            IBD2[ii][jj] += PairWiseIBD2(i, i + 1, j, j + 1) * sum;
            IBD1[ii][jj] += PairWiseIBD1(i, i + 1, j, j + 1) * sum;
#else
            bool isMale = (*ped)[family->path[jj]].sex == SEX_MALE;

            IBD2[ii][jj] += PairWiseIBD2(i, i + 1, j, j + (isMale ? 0 : 1)) * sum;
            IBD1[ii][jj] += PairWiseIBD1(i, i + 1, j, j + (isMale ? 0 : 1)) * sum;
#endif
            }

   return sum;
   }


 
