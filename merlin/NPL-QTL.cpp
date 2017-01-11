////////////////////////////////////////////////////////////////////// 
// merlin/NPL-QTL.cpp 
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
 
#include "NPL-QTL.h"
#include "Error.h"

Vector NPL_QTL::traitMeans;

void NPL_QTL::ScoreMeans(Pedigree & ped)
   {
   // Unweighted means for each phenotype
   traitMeans.Dimension(ped.traitCount);
   traitMeans.Zero();

   // Number of observations for each phenotype
   IntArray counts(ped.traitCount);
   counts.Zero();

   // Simple loop through pedigree
   for (int i = 0; i < ped.count; i++)
      for (int t = 0; t < ped.traitCount; t++)
         if (ped[i].traits[t] != _NAN_)
            {
            traitMeans[t] += ped[i].traits[t];
            counts[t]++;
            }

   // Divide totals by sum to get average
   for (int t = 0; t < ped.traitCount; t++)
      if (counts[t] !=  0)
         traitMeans[t] /= counts[t];
   }

void NPL_QTL::ScoreWithMeanZero(Mantra & m, int trait)
   { ScoreNPL(m, trait, 0.0); }

void NPL_QTL::ScoreWithSampleMean(Mantra & m, int trait)
   {
   if (m.pedigree->traitCount != traitMeans.dim)
      ScoreMeans(*m.pedigree);
   ScoreNPL(m, trait, traitMeans[trait]);
   }

void NPL_QTL::ScoreNPL(Mantra & m, int trait, double mean)
   {
   // Setup a framework for traversing inheritance vectors
   // in this family
   m.SelectTrait(trait, mean);

   // Clear the tree
   Clear();

   // The bit_count refers to the depth of the tree
   bit_count = m.bit_count;

   // The average vector stores the average deviation for each allele
   avg.Dimension(m.two_f);

   // The work horse...
   RecursivelyScoreNPL(m, bit_count - 1);

   // Save some memory?
   Trim();
   }

double NPL_QTL::CalculateNPL(Mantra & m)
   {
   // Calculate the average displacement for each founder allele
   avg.Zero();

   // Score the NPL
//   for (int i = m.two_f; i < m.two_n; i+=2)
   for (int i = 0; i < m.two_n; i++)
     avg[m.state[i]] += m.pheno[i];

   double sum = 0.0;
   for (int i = 0; i < m.two_f; i++)
      sum += avg[i] * avg[i];

   // Done!
   return sum;
   }

int NPL_QTL::RecursivelyScoreNPL(Mantra & m, int bit)
   {
   int node  = NewNode();

   if (bit != m.bit_count - 1)
      {
      int pivot = bit == -1 ? m.two_n : m.bits[bit];
      int last_pivot = m.bits[bit+1] + 1;

      for (int i = last_pivot; i < pivot; i++)
         m.state[i] = m.state[m.vector[i]];
      }

   // Bottom of the pedigree
   if (bit < 0)
      {
      nodes[node].type = TREE_NODE_LEAF;
      nodes[node].value = CalculateNPL(m);
      return node;
      }

   int pivot = bit >= 0 ? m.bits[bit] : 0;

   // Find out if we can save time by using symmetries in the data
   if (isRedundantNPL(m, pivot))
      {
      m.state[pivot] = m.state[m.vector[pivot] &= ~1];
      nodes[node].type = TREE_NODE_ONE;
      SAFE_SET(nodes[node].child[0], RecursivelyScoreNPL(m, bit - 1));
      }
   else
      {
      nodes[node].type = TREE_NODE_TWO;

      m.state[pivot] = m.state[m.vector[pivot] &= ~1];
      SAFE_SET(nodes[node].child[0], RecursivelyScoreNPL(m, bit - 1));

      m.state[pivot] = m.state[m.vector[pivot] |= 1];
      SAFE_SET(nodes[node].child[1], RecursivelyScoreNPL(m, bit - 1));
      }

   return node;
   }

bool NPL_QTL::isRedundantNPL(Mantra & m, int pivot)
   {
   // No affected descendants, no effect on NPL
   if (m.branch_phenotyped[pivot] == 0)
      return true;

   int fatallele = m.state[m.vector[pivot] & ~1];
   int motallele = m.state[m.vector[pivot] |  1];

   // No real choice, since parent is inbred
   if (motallele == fatallele)
      return true;

   // Bottleneck before reaching affected descendants
   if (m.branch_phenotyped[fatallele] == m.branch_phenotyped[pivot] &&
       m.branch_phenotyped[motallele] == m.branch_phenotyped[pivot] ||
       m.branch_phenotyped[fatallele] == m.branch_phenotyped[motallele] &&
       m.branch_phenotyped[fatallele] == m.branch_phenotyped[pivot] - 1 &&
       m.pheno[fatallele] != _NAN_ && m.pheno[fatallele] == m.pheno[motallele])
      return true;

   return false;
   }


 
