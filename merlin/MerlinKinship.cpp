////////////////////////////////////////////////////////////////////// 
// merlin/MerlinKinship.cpp 
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
 
#include "MerlinKinship.h"
#include "MerlinCore.h"
#include "NPL-ASP.h"
#include "Kinship.h"

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

MerlinKinship::MerlinKinship(Mantra & m)
   : mantra(m)
   {
   ped = NULL;
   family = NULL;
   casefile = kinfile = NULL;
   matrices = NULL;

   write = store = selectCases = false;
   }

MerlinKinship::~MerlinKinship()
   {
   if (matrices != NULL)
      delete [] matrices;
   }

void MerlinKinship::Analyse(AnalysisInfo & info, Tree & tree, int position)
   {
   // Dimension Kinship matrices
   kinship.Dimension(mantra.n, mantra.n);
   kinship.Zero();

   // Calculate IBDs as well as the overall likelihood
   double likelihood  = Score(tree, 0, tree.bit_count - 1);
   double scale = 1.0 / likelihood * 0.25;

   if (write) Output((*info.labels)[position], scale);
   if (selectCases) SelectCases((*info.labels)[position], scale);
   if (store) Store(position, scale);

   info.lk = likelihood * pow(2.0, -mantra.bit_count);
   }

void MerlinKinship::AnalyseUninformative(AnalysisInfo & info, Tree & tree)
   {
   // Allocates memory, etc...
   SetupFamily(info);

   // Dimension Kinship matrices
   kinship.Dimension(mantra.n, mantra.n);
   kinship.Zero();

   // Calculate IBDs as well as the overall likelihood
   double likelihood  = Score(tree, 0, tree.bit_count - 1);
   double scale = 1.0 / likelihood * 0.25;

   for (int position = 0; position < info.positions; position++)
      {
      if (write) Output((*info.labels)[position], scale);
      if (selectCases) SelectCases((*info.labels)[position], scale);
      if (store) Store(position, scale);
      }

   info.lk = likelihood * pow(2.0, -mantra.bit_count);
   }

double MerlinKinship::Score(Tree & tree, int node, int bit, int start)
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
            kinship[ii][jj] += sum *
               ((mantra.state[i]   == mantra.state[j]   ) +
                (mantra.state[i+1] == mantra.state[j+1] ) +
                (mantra.state[i]   == mantra.state[j+1] ) +
                (mantra.state[i+1] == mantra.state[j]   ));

   return sum;
   }

// Routines for scoring kinships in trees with symmetries
//

double MerlinKinship::AlleleKinship(int a, int b)
   {
   int parent_a = a < mantra.two_f ? a : mantra.vector[a];
   int parent_b = b < mantra.two_f ? b : mantra.vector[b];

   while (true)
      {
      // If the two alleles are the same
      if (a == b) return 1.0;

      // If both alleles descent from the same parent
      if ((parent_a & ~1) == (parent_b & ~1))
         // And the pedigree is not inbred
         if (mantra.ibd[parent_a >> 1][parent_a >> 1] == MANTRA_IBD_ONE)
            if (uninformative[a] || uninformative[b])
               return 0.5;
            else
               return parent_a == parent_b;
#ifdef __CHROMOSOME_X__
         else
            if (mantra.ibd[parent_a >> 1][parent_a >> 1] == MANTRA_IBD_ONE_MALE)
               return 1.0;
#endif
         // If the parent is inbred
         else
            if (uninformative[a] || uninformative[b])
               return 0.5 + AlleleKinship(parent_a & ~1, parent_a | 1) * 0.5;
            else
               return parent_a == parent_b ?
                      1.0 : AlleleKinship(parent_a & ~1, parent_a | 1);

      // If parents are different founders, return
      if (parent_a < mantra.two_f && parent_b < mantra.two_f)
         return 0.0;

      // Otherwise, move up the pedigree
      if (parent_a > parent_b)
         if (!uninformative[a])
            {
            a = parent_a;
            parent_a = mantra.vector[a];
            }
         else
            return 0.5 * (AlleleKinship(b, parent_a) +
                          AlleleKinship(b, parent_a + 1));
      else
         if (!uninformative[b])
            {
            b = parent_b;
            parent_b = mantra.vector[b];
            }
         else
            return 0.5 * (AlleleKinship(a, parent_b) +
                          AlleleKinship(a, parent_b + 1));
      }
   }

double MerlinKinship::WeightedScore(Tree & tree, int node, int bit, int start,
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
            kinship[ii][jj] += sum *
               (AlleleKinship(i, j) + AlleleKinship(i + 1, j + 1) +
                AlleleKinship(i + 1, j) + AlleleKinship(i, j + 1));

   return sum;
   }

// File I/O routines
//

void MerlinKinship::Output(const char * label, double scale)
   {
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
               fprintf(kinfile, " 0.000\n");
               break;
            case MANTRA_IBD_HALF :
               fprintf(kinfile, " 0.250\n");
               break;
            case MANTRA_IBD_ONE :
            case MANTRA_IBD_HALF_MALE :
               fprintf(kinfile, " 0.500\n");
               break;
            case MANTRA_IBD_ONE_MALE :
               fprintf(kinfile, " 1.000\n");
               break;
            default :
               if (!mantra.couples || j>=mantra.f || mantra.couple_index[j]==-1)
                  fprintf(kinfile, " %.5f\n", kinship[i][j] * scale);
               else
                  fprintf(kinfile, " %.5f\n", 0.5 *
                     (kinship[i][j] + kinship[i][mantra.partners[j]]) * scale);
            }
         }
   }

void MerlinKinship::OpenFiles(AnalysisInfo &, String & prefix)
   {
   if (write)
      {
      kinfilename = prefix + ".kin";
      kinfile = fopen(kinfilename, "wt");

      if (kinfile == NULL)
         error("Error opening file [%s]", (const char *) kinfilename);

      fprintf(kinfile, "FAMILY ID1 ID2 MARKER KINSHIP\n");
      }

   if (selectCases)
      {
      casefilename = prefix + ".sel";
      casefile = fopen(casefilename, "wt");

      if (casefile == NULL)
         error("Error opening file [%s]", (const char *) casefilename);

      fprintf(casefile, "FAMILY ID POSITION TRAIT NPL-SCORE CASE-SCORE\n");
      }
   }

void MerlinKinship::CloseFiles()
   {
   if (kinfile != NULL)
      {
      fclose(kinfile);
      printf("Kinship coefficients stored in file [%s]\n",
             (const char *) kinfilename);
      }

   if (casefile != NULL)
      {
      fclose(casefile);
      printf("Selected case information recorded in file [%s]\n",
             (const char *) casefilename);
      }

   kinfile = casefile = NULL;
   }

void MerlinKinship::SetupFamily(AnalysisInfo & info)
   {
   ped = mantra.pedigree;
   family = mantra.family;

   if (store)
      {
      int pairs = (family->count - family->founders) *
                  (family->count + family->founders + 1) / 2;

      // Always allocate at least one position, so we can distinguish
      // uninformative families (for VC analyses), from other groups
      // of unrelated individuals, which can help estimate overall trait
      // mean and variance
      matrices[family->serial].Dimension(pairs ? info.positions * pairs : 1);
      }

   if (selectCases)
      {
      // Calculate prior IBD sharing scores for this family (one per trait)
      Tree dummy;
      dummy.MakeMinimalTree(1.0, mantra.bit_count);

      // For calculating the variance of the NPL Pairs statistic
      NPL_Pairs npl;
      TreeInfo  stats;

      // Dimension Kinship matrices
      kinship.Dimension(mantra.n, mantra.n);
      kinship.Zero();

      // Calculate IBDs as well as the overall likelihood
      double likelihood = Score(dummy, 0, dummy.bit_count - 1);
      double scale = 1.0 / likelihood * 0.25;

      prior.Dimension(ped->affectionCount);
      priorToZ.Dimension(ped->affectionCount);
      for (int a = 0; a < ped->affectionCount; a++)
         {
         prior[a] = 0.0;

         // Calculate expected sharing score for each individual
         for (int i = 0; i < mantra.n; i++)
            if ((*ped)[family->path[i]].affections[a] == 2)
               {
               double expected = 0.0;

               for (int j = 0; j < mantra.n; j++)
                  if ((*ped)[family->path[j]].affections[a] == 2)
                     {
                     double pair = RetrieveFromScratch(i, j, scale);

                     expected += pair;
                     if (j <= i) prior[a] += pair;
                     }

               fprintf(casefile, "%s %s %s %s %s %.3f\n",
                       (const char *) family->famid,
                       (const char *) (*ped)[family->path[i]].pid,
                       (const char *) "EXPECTED",
                       (const char *) ped->affectionNames[a],
                       "0.000", expected);
               }

         // Score the NPL pairs statistic for this pedigree
         npl.ScoreNPL(mantra, a);

         // Calculate distribution summaries
         stats.GetMeanVar(npl);

         // Scale for converting NPL pairs scores into Z scores
         priorToZ[a] = (stats.var > 1E-10) ? sqrt(16 / stats.var) : 0.0;
         }
      }
   }

void MerlinKinship::Setup(AnalysisInfo & info)
   {
   if (matrices == NULL)
      matrices = new Vector [info.families];

   for (int i = 0; i < info.families; i++)
      matrices[i].Dimension(0);
   }

void MerlinKinship::SkipFamily(AnalysisInfo &)
   {
   if (store) matrices[family->serial].Dimension(0);
   }

void MerlinKinship::Store(int position, double scale)
   {
   int index = position * (family->count - family->founders) *
               (family->count + family->founders + 1) / 2;

   for (int i = mantra.f; i < mantra.n; i++)
      for (int j = 0; j <= i; j++)
         switch (mantra.ibd[i][j])
            {
            case MANTRA_IBD_ZERO :
               matrices[family->serial][index++] = 0.00;
               break;
            case MANTRA_IBD_HALF :
               matrices[family->serial][index++] = 0.25;
               break;
            case MANTRA_IBD_HALF_MALE :
            case MANTRA_IBD_ONE :
               matrices[family->serial][index++] = 0.50;
               break;
            case MANTRA_IBD_ONE_MALE :
               matrices[family->serial][index++] = 1.00;
               break;
            default :
               matrices[family->serial][index++] =
                 (!mantra.couples || j>=mantra.f || mantra.couple_index[j]==-1) ?
                  kinship[i][j] * scale :
                 (kinship[i][j] + kinship[i][mantra.partners[j]]) * scale * 0.5;
            }
   }

bool MerlinKinship::CheckFamily(int index)
   {
   return (matrices[index].Length() != 0);
   }

double MerlinKinship::Retrieve(int fam, int pos, int serial1, int serial2)
   {
   int founders = ped->families[fam]->founders;

   if (serial1 < founders && serial2 < founders)
#ifndef __CHROMOSOME_X__
       return (serial1 == serial2) ? 0.5 : 0.0;
#else
       {
       bool isMale = (*ped)[ped->families[fam]->path[serial1]].sex == SEX_MALE;

       return (serial1 == serial2) ? (isMale ? 1.0 : 0.5) : 0.0;
       }
#endif

   int count = ped->families[fam]->count;

   int index = pos * (count - founders) * (count + founders + 1) / 2;

   if (serial1 > serial2)
      index += (serial1 - founders) * (serial1 + founders + 1) / 2 + serial2;
   else
      index += (serial2 - founders) * (serial2 + founders + 1) / 2 + serial1;

   return matrices[fam][index];
   }

double MerlinKinship::RetrieveFromScratch(int i, int j, double scale)
   {
   if (j > i)
      return RetrieveFromScratch(j, i, scale);

   switch (mantra.ibd[i][j])
      {
      case MANTRA_IBD_ZERO :
         return 0.000;
      case MANTRA_IBD_HALF :
         return 0.250;
      case MANTRA_IBD_ONE :
      case MANTRA_IBD_HALF_MALE :
         return 0.500;
      case MANTRA_IBD_ONE_MALE :
         return 1.000;
      default :
         if (!mantra.couples || j>=mantra.f || mantra.couple_index[j]==-1)
            return kinship[i][j] * scale;
         else
            return 0.5 * (kinship[i][j] + kinship[i][mantra.partners[j]]) * scale;
      }
   }

void MerlinKinship::SelectCases(const char * label, double scale)
   {
   for (int a = 0; a < ped->affectionCount; a++)
      {
      // Total observed sharing
      double observed = 0.0;

      // Observed sharing per individual
      scores.Dimension(family->count);
      scores.Zero();

      for (int i = 0; i < mantra.n; i++)
         if ((*ped)[family->path[i]].affections[a] == 2)
            for (int j = 0; j < mantra.n; j++)
               if ((*ped)[family->path[j]].affections[a] == 2)
                  {
                  double pair = RetrieveFromScratch(i, j, scale);

                  scores[i] += pair;
                  if (i <= j) observed += pair;
                  }

      double max   = scores.Max();
      double npl   = (observed - prior[a]) * priorToZ[a];

      const char * tag = npl > -1E-7 ? "LINKED BEST" : "BEST";

      for (int i = 0; i < mantra.n; i++)
         if ((*ped)[family->path[i]].affections[a] == 2)
            fprintf(casefile, "%s %s %s %s %.3f %.3f %s\n",
                (const char *) family->famid,
                (const char *) (*ped)[family->path[i]].pid,
                (const char *) label,
                (const char *) ped->affectionNames[a],
                npl, scores[i],
                scores[i] == max ? (max++, tag) : "");
      }
   }

 
