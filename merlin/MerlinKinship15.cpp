////////////////////////////////////////////////////////////////////// 
// merlin/MerlinKinship15.cpp 
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
 
#include "MerlinKinship15.h"
#include "MerlinCore.h"

#include <math.h>

// The routines in this section perform Kinship calculations
//    Calculate(...) calculates and outputs kinship coefficients
//                   returning the overall likelihood of all vectors in the tree
//    Score(...)  traverses an inheritance tree and fills kinship matrix
//                   with unscaled kinship coefficients
//    Output(...) rescales coefficients to the 0.0 .. 1.0 range
//                   and stores them in a file
//    Open(...) and Close(...) manage the kinship file
//

// Flag for marking meiosis with ambiguous outcome
//

MerlinKinship15::MerlinKinship15(Mantra & m)
   : mantra(m)
   {
   ped = NULL;
   family = NULL;
   kinfile = NULL;

   isInbred = false;
   }

MerlinKinship15::~MerlinKinship15()
   {
   }

double MerlinKinship15::Calculate(Tree & tree, const char * label)
   {
   // Dimension Kinship matrices
   for (int i = (isInbred) ? 0 : 8; i < 14; i++)
      InitMatrix(i);

   // Calculate IBDs as well as the overall likelihood
   double likelihood  = Score(tree, 0, tree.bit_count - 1);

   if (label != NULL) Output(label, likelihood);

   return likelihood * pow(2.0, -mantra.bit_count);
   }

double MerlinKinship15::WScore(Tree & tree, int node, int bit, int start, double weight)
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
            mantra.state[pivot] = pivot;
            sum = WScore(tree, node, bit - 1, end, weight * 2.0);
            break;
         case TREE_NODE_ONE :
            mantra.state[pivot] = pivot;
            sum = WScore(tree, tree.nodes[node].child[0], bit - 1, end, weight * 2.0);
            break;
         case TREE_NODE_TWO :
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] & ~1];
            sum = WScore(tree, tree.nodes[node].child[0], bit - 1, end, weight);
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] | 1];
            sum += WScore(tree, tree.nodes[node].child[1],bit - 1, end, weight);
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
         if (mantra.ibd[ii][jj] != MANTRA_IBD_ONE &&
             mantra.ibd[ii][jj] != MANTRA_IBD_ZERO)
             {
             int si = mantra.state[i], si1 = mantra.state[i + 1];
             int sj = mantra.state[j], sj1 = mantra.state[j + 1];

             if (isInbred)
               {
               kinship[0][ii][jj] += Kinship4(si, si1, sj, sj1) * sum;
               kinship[1][ii][jj] += Kinship31(si, si1, sj, sj1) * sum;
               kinship[2][ii][jj] += Kinship31(si, si1, sj1, sj) * sum;
               kinship[3][ii][jj] += Kinship31(si, sj, sj1, si1) * sum;
               kinship[4][ii][jj] += Kinship31(si1, sj, sj1, si) * sum;
               kinship[5][ii][jj] += Kinship22(si, si1, sj, sj1) * sum;
               kinship[6][ii][jj] += Kinship211(si, si1, sj, sj1) * sum;
               kinship[7][ii][jj] += Kinship211(sj, sj1, si, si1) * sum;
               }

             kinship[8][ii][jj] += Kinship22(si, sj, si1, sj1) * sum;
             kinship[9][ii][jj] += Kinship211(si, sj, si1, sj1) * sum;
             kinship[10][ii][jj] += Kinship211(si1, sj1, si, sj) * sum;
             kinship[11][ii][jj] += Kinship22(si, sj1, si1, sj) * sum;
             kinship[12][ii][jj] += Kinship211(si, sj1, si1, sj) * sum;
             kinship[13][ii][jj] += Kinship211(si1, sj, si, sj1) * sum;
             }

   return sum;
   }

double MerlinKinship15::Score(Tree & tree, int node, int bit, int start)
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
            mantra.state[pivot] = pivot;
            sum = WScore(tree, node, bit - 1, end, 2.0);
            break;
         case TREE_NODE_ONE :
            mantra.state[pivot] = pivot;
            sum = WScore(tree, tree.nodes[node].child[0], bit - 1, end, 2.0);
            break;
         case TREE_NODE_TWO :
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] & ~1];
            sum = Score(tree, tree.nodes[node].child[0], bit - 1, end);
            mantra.state[pivot] = mantra.state[mantra.vector[pivot] | 1];
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
         if (mantra.ibd[ii][jj] != MANTRA_IBD_ONE &&
             mantra.ibd[ii][jj] != MANTRA_IBD_ZERO)
             {
             int a = mantra.state[i], b = mantra.state[i + 1];
             int c = mantra.state[j], d = mantra.state[j + 1];

             if (a == c)      // States 1, 2, 4, 9 or 10
               {
               if (b == d)    // States 1 or 9
                  if (a == b)
                     kinship[0][ii][jj] += sum;
                  else
                     kinship[8][ii][jj] += sum;
               else if (a == b)
                  kinship[1][ii][jj] += sum;
               else if (c == d)
                  kinship[3][ii][jj] += sum;
               else
                  kinship[9][ii][jj] += sum;
               }
             else if (b == d) // States 3, 5, 11
               {
               if (a == b)
                  kinship[2][ii][jj] += sum;
               else if (c == d)
                  kinship[4][ii][jj] += sum;
               else
                  kinship[10][ii][jj] += sum;
               }
             else if (a == d) // States 12 or 13
               {
               if (b == c)
                  kinship[11][ii][jj] += sum;
               else
                  kinship[12][ii][jj] += sum;
               }
             else if (b == c)
               kinship[13][ii][jj] += sum;
             else if (a == b)
               if (c == d)
                  kinship[5][ii][jj] += sum;
               else
                  kinship[6][ii][jj] += sum;
             else if (c == d)
               kinship[7][ii][jj] += sum;
             }

   return sum;
   }

// File I/O routines
//

void MerlinKinship15::Output(const char * label, double scale)
   {
   scale = 1.0 / scale;

   // Before output, we may need to apply founder couple symmetry to
   // estimated coefficients ...
   if (mantra.couples)
      for (int i = 0; i < mantra.f; i++)
         if (mantra.couple_index[i] != -1 && mantra.partners[i] < i)
            {
            int partner = mantra.partners[i];

            // First we average all estimates with those of partner
            // Estimates for sharing against first generation offspring are unchanged
            for (int k = (isInbred ? 0 : 8); k < 14; k++)
               for (int j = mantra.f; j < mantra.n; j++)
                  if ((mantra.vector[j * 2] >> 1) != i && (mantra.vector[j * 2] >> 1) != partner)
                     kinship[k][j][i] = kinship[k][j][partner] =
                        (kinship[k][j][i] + kinship[k][j][partner]) * 0.5;

            // Then we need to consider what happends to offspring when
            // we flip phase -- note that these individuals can never be
            // inbred
            for (int j = mantra.f; j < mantra.n; j++)
               if ((mantra.vector[j * 2] >> 1) == i || (mantra.vector[j * 2] >> 1) == partner)
                  {
                  for (int k = mantra.f; k < j; k++)
                     if ((mantra.vector[k * 2] >> 1) == i || (mantra.vector[k * 2] >> 1) == partner)
                        kinship[9][j][k] = kinship[10][j][k] =
                           (kinship[9][j][k] + kinship[10][j][k]) * 0.5;
                     else
                        {
                        kinship[8][j][k] = kinship[11][j][k] =
                           (kinship[8][j][k] + kinship[11][j][k]) * 0.5;
                        kinship[9][j][k] = kinship[13][j][k] =
                           (kinship[9][j][k] + kinship[13][j][k]) * 0.5;
                        kinship[10][j][k] = kinship[12][j][k] =
                           (kinship[10][j][k] + kinship[12][j][k]) * 0.5;
                        if (!isInbred) continue;
                        kinship[3][j][k] = kinship[4][j][k] =
                           (kinship[3][j][k] + kinship[4][j][k]) * 0.5;
                        }

                  for (int k = j + 1; k < mantra.n; k++)
                     if ((mantra.vector[k * 2] >> 1) == i || (mantra.vector[k * 2] >> 1) == partner)
                        kinship[9][k][j] = kinship[10][k][j] =
                           (kinship[9][k][j] + kinship[10][k][j]) * 0.5;
                     else
                        {
                        kinship[8][k][j] = kinship[11][k][j] =
                           (kinship[8][k][j] + kinship[11][k][j]) * 0.5;
                        kinship[9][k][j] = kinship[12][k][j] =
                           (kinship[9][k][j] + kinship[12][k][j]) * 0.5;
                        kinship[10][k][j] = kinship[13][k][j] =
                           (kinship[10][k][j] + kinship[13][k][j]) * 0.5;
                        if (!isInbred) continue;
                        kinship[1][k][j] = kinship[2][k][j] =
                           (kinship[1][k][j] + kinship[2][k][j]) * 0.5;
                        }
                  }
            }

   for (int i = 0; i < mantra.n; i++)
      for (int j = 0; j <= i; j++)
         {
         fprintf(kinfile, "%s %s %s %s ",
                (const char *) family->famid,
                (const char *) (*ped)[family->path[i]].pid,
                (const char *) (*ped)[family->path[j]].pid,
                (const char *) label);

         switch (mantra.ibd[i][j])
            {
            case MANTRA_IBD_ZERO :
               fprintf(kinfile, " 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000"
                                " 0.000 0.000 0.000 0.000 0.000 0.000 1.000\n");
               break;
            case MANTRA_IBD_ONE :
               fprintf(kinfile, " 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000"
                                " 1.000 0.000 0.000 0.000 0.000 0.000 0.000\n");
               break;
            default :
               {
               double kin[14];

               for (int k = 0; k < (isInbred ? 0 : 8); k++)
                  kin[k] = 0.0;

               for (int k = (isInbred ? 0 : 8); k < 14; k++)
                  kin[k] = kinship[k][i][j] * scale;

               if (j < mantra.f)
                  {
                  kin[1] = kin[2] = (kin[1] + kin[2]) * 0.5;
                  kin[8] = kin[11] = (kin[8] + kin[11]) * 0.5;
                  kin[9] = kin[12] = (kin[9] + kin[12]) * 0.5;
                  kin[10] = kin[13] = (kin[10] + kin[13]) * 0.5;
                  }

               for (int k = 0; k < 14; k++)
                  fprintf(kinfile, " %5.3f", kin[k]);

               double kin15 = 0.0;

               for (int k = (isInbred ? 0 : 8); k < 14; k++)
                  kin15 += kin[k];

               kin15 = kin15 > 1.0 ? 0.0 : 1.0 - kin15;

               fprintf(kinfile, " %5.3f\n", kin15);
               }
            }
         }
   }

void MerlinKinship15::OpenFile()
   {
   String filename(MerlinCore::filePrefix);
   filename += ".s15";

   kinfile = fopen(filename, "wt");

   if (kinfile == NULL)
      error("Error opening file [%s]", (const char *) filename);

   fprintf(kinfile, "FAMILY ID1 ID2 POSITION S01 S02 S03 S04 S05 S06 S07 S08 "
                    "S09 S10 S11 S12 S13 S14 S15\n");
   }

void MerlinKinship15::CloseFile()
   {
   fclose(kinfile);
   printf("Extended identity state probabilities stored in file [%s.s15]\n", (const char *) MerlinCore::filePrefix);
   }

void MerlinKinship15::SelectFamily(Pedigree * p, Family * f)
   {
   ped = p;
   family = f;

   isInbred = false;

   // Find out if there are any inbred individuals
   for (int i = mantra.f; i < mantra.n; i++)
      if (mantra.ibd[i][i] != MANTRA_IBD_ONE)
         {
         isInbred = true;
         break;
         }
   }

double MerlinKinship15::Kinship2(int allele1, int allele2)
   {
   if (allele1 == allele2)
      return 1.0;

   Order(allele1, allele2);

   if (allele2 < mantra.two_f)
      return 0.0;

   int stateA = mantra.state[mantra.vector[allele2] & ~1];
   int stateB = mantra.state[mantra.vector[allele2] | 1];

   return 0.5 * (Kinship2(allele1, stateA) + Kinship2(allele1, stateB));
   }

double MerlinKinship15::Kinship3(int allele1, int allele2, int allele3)
   {
   Order(allele1, allele2);
   Order(allele2, allele3);
   Order(allele1, allele2);

   if (allele1 == allele2)
      return Kinship2(allele1, allele3);

   if (allele2 < mantra.two_f)
      return 0.0;

   if (allele2 == allele3)
      return Kinship2(allele1, allele2);

   int stateA = mantra.state[mantra.vector[allele3] & ~1];
   int stateB = mantra.state[mantra.vector[allele3] | 1];

   return 0.5 * (Kinship3(allele1, allele2, stateA) +
                 Kinship3(allele1, allele2, stateB));
   }

double MerlinKinship15::Kinship4(int allele1, int allele2, int allele3, int allele4)
   {
   Order(allele1, allele2);
   Order(allele2, allele3);
   Order(allele3, allele4);
   Order(allele1, allele2);
   Order(allele2, allele3);

   if (allele1 == allele2)
      return Kinship3(allele1, allele3, allele4);

   if (allele2 < mantra.two_f)
      return 0.0;

   if (allele1 == allele3 || allele2 == allele3)
      return Kinship3(allele1, allele2, allele4);

   if (allele3 == allele4)
      return Kinship3(allele1, allele2, allele3);

   int stateA = mantra.state[mantra.vector[allele4] & ~1];
   int stateB = mantra.state[mantra.vector[allele4] | 1];

   return 0.5 * (Kinship4(allele1, allele2, allele3, stateA) +
                 Kinship4(allele1, allele2, allele3, stateB));
   }

double MerlinKinship15::Kinship31(int allele1, int allele2, int allele3, int allele4)
   {
   if (allele1 == allele4 || allele2 == allele4 || allele3 == allele4)
      return 0.0;

   Order(allele1, allele2);
   Order(allele2, allele3);

   if (allele1 == allele2)
      return Kinship21(allele1, allele3, allele4);

   Order(allele1, allele2);

   if (allele2 < mantra.two_f)
      return 0.0;

   if (allele2 == allele3)
      return Kinship21(allele1, allele2, allele4);

   if (allele4 > allele3)
      {
      int stateA = mantra.state[mantra.vector[allele4] & ~1];
      int stateB = mantra.state[mantra.vector[allele4] | 1];

      return 0.5 * (Kinship31(allele1, allele2, allele3, stateA) +
                    Kinship31(allele1, allele2, allele3, stateB));
      }
   else
      {
      int stateA = mantra.state[mantra.vector[allele3] & ~1];
      int stateB = mantra.state[mantra.vector[allele3] | 1];

      return 0.5 * (Kinship31(allele1, allele2, stateA, allele4) +
                    Kinship31(allele1, allele2, stateB, allele4));
      }
   }

double MerlinKinship15::Kinship22(int allele1, int allele2, int allele3, int allele4)
   {
   if (allele1 == allele2)
      return Kinship21(allele3, allele4, allele1);

   if (allele3 == allele4)
      return Kinship21(allele1, allele2, allele3);

   if (allele1 == allele3 || allele1 == allele4 || allele2 == allele3 || allele2 == allele4)
      return 0.0;

   Order(allele1, allele2);
   Order(allele3, allele4);

   if (allele2 < mantra.two_f || allele4 < mantra.two_f)
      return 0.0;

   if (allele4 > allele2)
      {
      int stateA = mantra.state[mantra.vector[allele4] & ~1];
      int stateB = mantra.state[mantra.vector[allele4] | 1];

      return 0.5 * (Kinship22(allele1, allele2, allele3, stateA) +
                    Kinship22(allele1, allele2, allele3, stateB));
      }
   else
      {
      int stateA = mantra.state[mantra.vector[allele2] & ~1];
      int stateB = mantra.state[mantra.vector[allele2] | 1];

      return 0.5 * (Kinship22(allele1, stateA, allele3, allele4) +
                    Kinship22(allele1, stateB, allele3, allele4));
      }
   }

double MerlinKinship15::Kinship21(int allele1, int allele2, int allele3)
   {
   if (allele1 == allele3 || allele2 == allele3)
      return 0.0;

   if (allele1 == allele2)
      return Kinship11(allele1, allele3);

   Order(allele1, allele2);

   if (allele2 < mantra.two_f)
      return 0.0;

   if (allele3 > allele2)
      {
      int stateA = mantra.state[mantra.vector[allele3] & ~1];
      int stateB = mantra.state[mantra.vector[allele3] | 1];

      return 0.5 * (Kinship21(allele1, allele2, stateA) +
                    Kinship21(allele1, allele2, stateB));
      }
   else
      {
      int stateA = mantra.state[mantra.vector[allele2] & ~1];
      int stateB = mantra.state[mantra.vector[allele2] |= 1];

      return 0.5 * (Kinship21(allele1, stateA, allele3) +
                    Kinship21(allele1, stateB, allele3));
      }
   }

double MerlinKinship15::Kinship211(int allele1, int allele2, int allele3, int allele4)
   {
   if (allele1 == allele2)
      return Kinship111(allele1, allele3, allele4);

   if (allele1 == allele3 || allele1 == allele4 || allele3 == allele4 ||
       allele2 == allele3 || allele2 == allele4)
      return 0.0;

   Order(allele1, allele2);
   Order(allele3, allele4);

   if (allele2 < mantra.two_f)
      return 0.0;

   if (allele4 > allele2)
      {
      int stateA = mantra.state[mantra.vector[allele4] & ~1];
      int stateB = mantra.state[mantra.vector[allele4] | 1];

      return 0.5 * (Kinship211(allele1, allele2, allele3, stateA) +
                    Kinship211(allele1, allele2, allele3, stateB));
      }
   else
      {
      int stateA = mantra.state[mantra.vector[allele2] & ~1];
      int stateB = mantra.state[mantra.vector[allele2] | 1];

      return 0.5 * (Kinship211(allele1, stateA, allele3, allele4) +
                    Kinship211(allele1, stateB, allele3, allele4));
      }
   }

double MerlinKinship15::Kinship111(int allele1, int allele2, int allele3)
   {
   if (allele1 == allele2 || allele2 == allele3 || allele1 == allele3)
      return 0.0;

   Order(allele1, allele2);
   Order(allele2, allele3);

   if (allele3 < mantra.two_f)
      return 1.0;

   int stateA = mantra.state[mantra.vector[allele3] & ~1];
   int stateB = mantra.state[mantra.vector[allele3] | 1];

   return 0.5 * (Kinship111(allele1, allele2, stateA) +
                 Kinship111(allele1, allele2, stateB));
   }


 
 
