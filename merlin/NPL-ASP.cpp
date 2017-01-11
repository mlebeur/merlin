////////////////////////////////////////////////////////////////////// 
// merlin/NPL-ASP.cpp 
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
 
#include "NPL-ASP.h"
#include "Error.h"

const double NPL_ALL::factors[20] =
   {  1.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,
      8.0,  9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0,
     16.0, 17.0, 18.0, 19.0 };

const double NPL_ALL::reciprocals[20] =
   {  1.0 / 1.0,  1.0 / 1.0,  1.0 / 2.0,  1.0 / 3.0,
      1.0 / 4.0,  1.0 / 5.0,  1.0 / 6.0,  1.0 / 7.0,
      1.0 / 8.0,  1.0 / 9.0,  1.0 / 10.0, 1.0 / 11.0,
      1.0 / 12.0, 1.0 / 13.0, 1.0 / 14.0, 1.0 / 15.0,
      1.0 / 16.0, 1.0 / 17.0, 1.0 / 18.0, 1.0 / 19.0 };

void NPL_ASP_Tree::ScoreNPL(Mantra & m, int affection)
   {
   // Setup a framework for traversing inheritance vectors
   // in this family
   m.SelectAffection(affection);

   // Clear the tree
   Clear();

   // The bit_count refers to the depth of the tree
   bit_count = m.bit_count;

   // If there is less than one affected in the family
   if (m.aff_count < 1)
      {
      int node = NewNode();
      nodes[node].type = TREE_NODE_LEAF;
      nodes[node].value = 1.0;
      return;
      }

   // Do specific initialization for this statistic
   PrepareForScoring(m);

   // The work horse...
   RecursivelyScoreNPL(m, bit_count - 1);

   // Save some memory?
   Trim();
   }

int NPL_ASP_Tree::RecursivelyScoreNPL(Mantra & m, int bit)
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

bool NPL_ASP_Tree::isRedundantNPL(Mantra & m, int pivot)
   {
   // No affected descendants, no effect on NPL
   if (m.branch_affected[pivot] == 0)
      return true;

   int fatallele = m.state[m.vector[pivot] & ~1];
   int motallele = m.state[m.vector[pivot] |  1];

   // No real choice, since parent is inbred
   if (motallele == fatallele)
      return true;

   // Bottleneck before reaching affected descendants
   if (m.branch_affected[fatallele] == m.branch_affected[pivot] &&
       m.branch_affected[motallele] == m.branch_affected[pivot] ||
       m.branch_affected[fatallele] == m.branch_affected[motallele] &&
       m.branch_affected[fatallele] == m.branch_affected[pivot] - 1 &&
       m.aff[fatallele] == 1 && m.aff[motallele] == 1)
      return true;

   return false;
   }

void NPL_ASP_Tree::PrepareForScoring(Mantra & ) { }

double NPL_ALL::CalculateNPL(Mantra & m)
   {
   if (m.aff_count < 8)
      // Evaluate statistic directly
      return ScoreSubset(m, m.affecteds);
   else
      {
      // Calculate statistic for connected subsets of affected individuals

      // First identify subsets of connected individuals
      graph.SetSequence(0, 1);

      // Quick Union Algorithm used to build graph
      for (int i = 0; i < m.aff_count; i ++)
         {
         int p  = m.affecteds[i];
         int g1 = m.state[p]; while (graph[g1] != g1) g1 = graph[g1];
         int g2 = m.state[p + 1]; while (graph[g2] != g2) g2 = graph[g2];

         graph[g1] = g2;
         }

      // Next we flatten the graph
      for (int i = 0; i < m.two_f; i++)
         {
         int g1 = graph[i]; while (graph[g1] != g1) g1 = graph[g1];

         graph[i] = g1;
         }

      // Finally we cycle through graph components
      double score = 1.0;
      for (int i = 0; i < m.two_f; i++)
         if (graph[i] == i)
            {
            subset.Clear();

            // List affected individuals in graph
            for (int j = 0; j < m.aff_count; j++)
               if (graph[m.state[m.affecteds[j]]] == i)
                  subset.Push(m.affecteds[j]);

            // Finally, calculate statistic for each subset of affecteds
            if (subset.Length())
               {
               permutations = 1 << subset.Length();
               permutationWeight = 1.0 / permutations;
               score *= ScoreSubset(m, subset);
               }
            }

      return score;
      }
   }

double NPL_ALL::ScoreSubset(Mantra & m, IntArray & affecteds)
   {
   int count = affecteds.Length();

   // Count the number of times each founder allele appears
   counts.Zero();

   // Score the NPL for the first permutation
   double score = permutationWeight;
   for (int i = 0; i < count; i++)
     score *= factors[++counts[m.state[affecteds[i]]]];

   double sum = score;

   // To get the new score from the previous one,
   // multiply by 1.0 / count[outgoing_allele]
   // and by count[incoming_allele]
   for (int i = 0; i < permutations - 1; i++)
      {
      int j = i, bit = 0;
      while (j & 1) { j >>= 1; bit++; }
      score *= reciprocals[counts[m.state[affecteds[bit]]]--];
      score *= factors[++counts[m.state[affecteds[bit]^=1]]];
      sum += score;
      }

   // Done!
   return sum;
   }

void NPL_ALL::PrepareForScoring(Mantra & m)
   {
   // Check that we can handle this pedigree
   if (m.aff_count >= 20)
      for (int i = 0; i < m.two_f; i++)
         if (m.branch_affected[i] >= 20 || m.aff_count > 30)
            error("Calculating Whittemore and Halpern Statistic\n\n"
                  "There are too many affecteds in family %s\n",
                  (const char *) m.family->famid);

   // The counts array stores the number of times each allele appears
   counts.Dimension(m.two_f);

   // How many different sets of alleles can we pick amongst affecteds?
   permutations = 1 << m.aff_count;
   permutationWeight = 1.0 / permutations;

   // Auxiliary structures for handling pedigrees with several affecteds
   graph.Dimension(m.two_f);
   }

double NPL_Pairs::CalculateNPL(Mantra & m)
   {
   double score = 0.0;

   // This statistic is simply the sum of kinship coefficients
   // among affected individuals in the pedigree
   for (int i = 0; i < m.aff_count; i++)
      {
      int maternal = m.state[m.affecteds[i]];
      int paternal = m.state[m.affecteds[i] + 1];

      for (int j = 0; j <= i; j++)
         {
         if (m.state[m.affecteds[j]] == maternal)
            score++;
         if (m.state[m.affecteds[j]] == paternal)
            score++;
         if (m.state[m.affecteds[j]+1] == maternal)
            score++;
         if (m.state[m.affecteds[j]+1] == paternal)
            score++;
         }
      }

   return score;
   }

double NPL_Pairs_For_Simwalk2::CalculateNPL(Mantra & m)
   {
   double score = 0.0;

   // This statistic is simply the sum of kinship coefficients
   // among affected individuals in the pedigree
   for (int i = 1; i < m.aff_count; i++)
      {
      int maternal = m.state[m.affecteds[i]];
      int paternal = m.state[m.affecteds[i] + 1];

      for (int j = 0; j < i; j++)
         {
         if (m.state[m.affecteds[j]] == maternal)
            score++;
         if (m.state[m.affecteds[j]] == paternal)
            score++;
         if (m.state[m.affecteds[j]+1] == maternal)
            score++;
         if (m.state[m.affecteds[j]+1] == paternal)
            score++;
         }
      }

   return score;
   }

double NPL_Pairs_Maternal::CalculateNPL(Mantra & m)
   {
   double score = 0.0;

   // This statistic is simply the sum of kinship coefficients
   // among affected individuals in the pedigree
   for (int i = 1; i < m.aff_count; i++)
      {
      bool founder  = m.affecteds[i] < m.two_f;
      int  maternal = m.state[m.affecteds[i]];
      int  paternal = m.state[m.affecteds[i] + 1];

      for (int j = 0; j < i; j++)
         {
         bool founder2  = m.affecteds[j] < m.two_f;
         int  maternal2 = m.state[m.affecteds[j]];
         int  paternal2 = m.state[m.affecteds[j] + 1];
         double weight  = (founder || founder2) ? 0.5 : 1.0;

         // Founders can't share alleles IBD
         if (founder && founder2) continue;

         // Score IBD, allowing for founder symmetries
         if ((maternal == maternal2) ||
             (paternal == maternal2 && founder) ||
             (maternal == paternal2 && founder2))
            score += weight;
         }
      }

   return score;
   }

double NPL_Pairs_Paternal::CalculateNPL(Mantra & m)
   {
   double score = 0.0;

   // This statistic is simply the sum of kinship coefficients
   // among affected individuals in the pedigree
   for (int i = 1; i < m.aff_count; i++)
      {
      bool founder  = m.affecteds[i] < m.two_f;
      int  maternal = m.state[m.affecteds[i]];
      int  paternal = m.state[m.affecteds[i] + 1];

      for (int j = 0; j < i; j++)
         {
         bool founder2  = m.affecteds[j] < m.two_f;
         int  maternal2 = m.state[m.affecteds[j]];
         int  paternal2 = m.state[m.affecteds[j] + 1];
         double weight  = (founder || founder2) ? 0.5 : 1.0;

         // Founders can't share alleles IBD
         if (founder && founder2) continue;

         // Score IBD, allowing for founder symmetries
         if ((paternal == paternal2) ||
             (maternal == paternal2 && founder) ||
             (paternal == maternal2 && founder2))
            score += weight;
         }
      }

   return score;
   }

 
