////////////////////////////////////////////////////////////////////// 
// merlin/MerlinMatrix.cpp 
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
 
#include "MerlinMatrix.h"
#include "MerlinCore.h"
#include "Error.h"

#include <math.h>

// The routines in this section perform Kinship matrix calculations
//    Calculate(...) calculates and outputs kinship matrices
//                   returning the overall likelihood of all vectors in the tree
//    Score(...)  traverses an inheritance tree and stores individual matrices
//                   in the matrices stringmap
//    Output(...) rescales likelihoods to the 0.0 .. 1.0 range
//                   and stores them in a file
//    Open(...) and Close(...) manage the kinship matrix file
//

MerlinMatrix::MerlinMatrix(Mantra & m) : mantra(m)
   {
   matrixfile = NULL;
   symmetries = NULL;
   size = 0;
   }

MerlinMatrix::~MerlinMatrix()
   {
   if (symmetries != NULL) delete [] symmetries;
   }

double MerlinMatrix::Calculate(Tree & tree, const char * label)
   {
   // Clear probability vector
   probabilities.Clear();
   matrices.Clear();

   // Calculate IBDs as well as the overall likelihood
   double likelihood  = Score(tree, 0, tree.bit_count - 1, mantra.two_f);

   Output(label, likelihood);

   return likelihood * pow(2.0, -mantra.bit_count);
   }

double MerlinMatrix::Score(Tree & tree, int node, int bit, int start, int pos)
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

   for ( int i = start, ii = start >> 1 ; i < end; i++, i++, ii++)
      for ( int j = 0, jj = 0; j <= i; j++, j++, jj++)
         if (mantra.ibd[ii][jj] == MANTRA_IBD_UNKNOWN)
            matrix[pos++] = '0' +
               ((mantra.state[i]   == mantra.state[j]   ) +
                (mantra.state[i+1] == mantra.state[j+1] ) +
                (mantra.state[i]   == mantra.state[j+1] ) +
                (mantra.state[i+1] == mantra.state[j]   ));

   double sum = 0.0; // Initialization avoids compiler warning

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

      // And for leaf nodes, we either ...
      int index = matrices.Integer(matrix);

      // Update matrix likelihood ...
      if (index != -1)
         probabilities[index] += sum;
      // Or add a new entry for matrices we haven't seen before ...
      else
         {
         probabilities.Push(sum);
         matrices.SetInteger(matrix, probabilities.dim - 1);
         }
      }

   return sum;
   }

// File I/O routines
//

void MerlinMatrix::Output(const char * label, double scale)
   {
   int unraveled = 1 << couples;

   scale = 1.0 / (scale * unraveled);

   fprintf(matrixfile, "Position %s\n", (const char *) label);

   for (int i = 0; i < matrices.Length(); i++)
      {
      double prob = probabilities[matrices.Integer(i)] * scale;
      fprintf(matrixfile, "%s %g\n", (const char *) matrices[i], prob);

      for (int j = 1; j < unraveled; j++)
         fprintf(matrixfile, "%s %g\n", Unravel(matrices[i], matrix, j), prob);
      }
   }

void MerlinMatrix::OpenFile()
   {
   String filename(MerlinCore::filePrefix);
   filename += ".kmx";

   matrixfile = fopen(filename, "wt");

   if (matrixfile == NULL)
      error("Error opening file [%s]", (const char *) filename);
   }

void MerlinMatrix::CloseFile()
   {
   fclose(matrixfile);
   printf("Kinship matrices stored in file [%s.kmx]\n", (const char *) MerlinCore::filePrefix);
   }

void MerlinMatrix::SelectFamily(Pedigree * p, Family * f)
   {
   if (matrixfile == NULL) return;

   ped = p;
   family = f;

   fprintf(matrixfile, "Family %s\n", (const char *) f->famid);
   fprintf(matrixfile, "Pairs  ");

   matrix.Clear();
   for (int i = family->founders; i < family->count; i++)
      for (int j = 0; j <= i; j++)
         if (mantra.ibd[i][j] == MANTRA_IBD_UNKNOWN)
            {
            matrix += 'X';
            fprintf(matrixfile, "%s-%s ",
                   (const char *) (*ped)[family->path[i]].pid,
                   (const char *) (*ped)[family->path[j]].pid);
            }
   fprintf(matrixfile, "\n");

   if ((couples = mantra.couples) == 0) return;

   ResetSymmetries();

   int position = 0;

   for (int i = family->founders; i < family->count; i++)
      for (int j = 0; j <= i; j++)
         if (mantra.ibd[i][j] == MANTRA_IBD_UNKNOWN)
            {
            if (j < family->founders && mantra.couple_index[j] != -1)
               symmetries[mantra.couple_index[j]].Push(position);
            position++;
            }
   }

// These routines unravel kinship matrices based on pedigree symmetries
//

void MerlinMatrix::ResetSymmetries()
   {
   if (couples > size)
      {
      if (symmetries) delete [] symmetries;

      size = (couples + 3) & 4;
      symmetries = new IntArray[size];
      }

   for (int i = 0; i < couples; i++)
      symmetries[i].Clear();
   }

const char * MerlinMatrix::Unravel(const String & original, String & unraveled, int idx)
   {
   int couple = 0;

   unraveled = original;

   while (idx > 0)
      {
      if (idx & 1)
         for (int i = 0; i < symmetries[couple].Length(); i+=2)
            {
            unraveled[symmetries[couple][i]] = original[symmetries[couple][i+1]];
            unraveled[symmetries[couple][i+1]] = original[symmetries[couple][i]];
            }

      idx >>= 1;
      couple++;
      }

   return unraveled;
   }
 
 
