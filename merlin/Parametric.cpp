////////////////////////////////////////////////////////////////////// 
// merlin/Parametric.cpp 
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
 
#include "Parametric.h"
#include "MathStats.h"

#include <math.h>

LikelihoodTree::LikelihoodTree()
   {
   stack_depth = 0;
   }

LikelihoodTree::~LikelihoodTree()
   {
   FreeStack();
   }

void LikelihoodTree::ScoreModel(Mantra & m)
   {
   // Setup a framework for traversing inheritance vectors
   m.SelectBinaryTrait(model.affection);

   // Clear the tree and define its depth
   Clear();
   bit_count = m.bit_count;

   // Setup up basic graphs and allele states, using founder genotypes
   AllocateStack(m);

   // The work horse
   ScoreRecursive(m, m.bit_count - 1, m.two_f, graphs, graph_phenotypes, m.two_f);

   // Generate more compact form of tree, if possible
   Trim();
   }

int LikelihoodTree::ScoreRecursive
    (Mantra & m, int bit, int last_pivot,
    IntArray * graph, IntArray * graph_phenos, int freec)
   {
   // Allocate a new tree node
   int node = NewNode();

   // Where in the inheritance vector are we now?
   int pivot = bit >= 0 ? m.bits[bit] : m.two_n;

   // Ensure graph consistency and track whether additional genotyped descendants
   // are possible for each graph
   graph_phenos[1] = graph_phenos[0];
   graph_phenos++;

   for (int i = last_pivot; i < pivot; i++)
      {
      m.state[i] = m.state[m.vector[i]];

      // Adjust number of genotyped descendants, taking care to avoid twins
      if (m.isBit[i])
         // This line updates the number of potential genotyped descendants
         (*graph_phenos)[(*graph)[m.state[m.vector[i] ^ 1]]] -= m.branch_affected[i];
      }

   // Update the graph of connected components in pedigree
   bool save_graph = false;
   for ( ; last_pivot < pivot; last_pivot++)
      {
      int phenotype = m.aff[last_pivot];

      if (phenotype == 0) continue;

      if (last_pivot & 1)
         {
         int componentA = (*graph)[m.state[last_pivot - 1]];
         int componentB = (*graph)[m.state[last_pivot]];

         if (componentA != componentB)
            {
            if (!save_graph)
               {
               save_graph = true;
               graph[1] = graph[0];
               graph++;
               }

            (*graph_phenos)[freec] = (*graph_phenos)[componentA] +
                                    (*graph_phenos)[componentB] - 2;

            components[freec] = components[componentA];
            components[freec].Append(components[componentB]);

            for (int i = 0; i < components[freec].Length(); i++)
               (*graph)[components[freec][i]] = freec;

            freec++;
            }
         else
            (*graph_phenos)[componentA] -= 2;
         }
      }

   // Update likelihood components
   UpdateLikelihood(m, *graph, *graph_phenos, last_pivot);

   // Is this a leaf node?
   if (bit < 0)
      {
      nodes[node].type = TREE_NODE_LEAF;
      nodes[node].value = CalculateLikelihood(m, *graph);

      return node;
      }

   // Find out if we can save time by using symmetries in the data
   switch (nodes[node].type = FindRedundancy(m, pivot))
      {
      case TREE_NODE_ZERO :
         nodes[node].value = 0.0;
         break;
      case TREE_NODE_ONE :
         m.state[pivot] = m.state[m.vector[pivot] &= ~1];
         SAFE_SET(nodes[node].child[0],
                  ScoreRecursive(m, bit - 1, last_pivot, graph, graph_phenos,
                                 freec));
         break;
      case TREE_NODE_TWO :
         m.state[pivot] = m.state[m.vector[pivot] &= ~1];
         SAFE_SET(nodes[node].child[0],
                  ScoreRecursive(m, bit - 1, last_pivot, graph, graph_phenos,
                                 freec));

         m.state[pivot] = m.state[m.vector[pivot] |= 1];
         SAFE_SET(nodes[node].child[1],
                  ScoreRecursive(m, bit - 1, last_pivot, graph, graph_phenos,
                                 freec));
         break;
      }

   return node;
   }

int LikelihoodTree::FindRedundancy(Mantra & m, int pivot)
   {
   // Are there are any genotyped descendants?
   if (m.branch_affected[pivot] == 0)
      return TREE_NODE_ONE;

   // Are parents homozygous?
   int parent = m.vector[pivot];
   int matallele = m.state[parent & ~1];
   int fatallele = m.state[parent | 1];

   if (matallele == fatallele)
      return TREE_NODE_ONE;

   // If all the phenotyped descendants of the two alternative
   // founder alleles are on this branch, it doesn't matter which
   // we pick -- one will go down, the other will be fixed (or not)
   if (m.branch_affected[matallele] == m.branch_affected[pivot] &&
       m.branch_affected[fatallele] == m.branch_affected[pivot])
       return TREE_NODE_ONE;

   return TREE_NODE_TWO;
   }

double LikelihoodTree::CalculateLikelihood(Mantra & m, IntArray & graph)
   {
   double likelihood = 1.0;

   /* Go through each component */
   for (int i = 0; i < m.two_f; i ++)
      if (components[graph[i]][0] == i)
         likelihood *= graph_likelihood[graph[i]];

   return likelihood;
   }

void LikelihoodTree::UpdateLikelihood(Mantra & m,
      IntArray & graph, IntArray & graph_phenos, int pivot)
   {
   for (int i = 0; i < m.two_f; i ++)
      {
      int c = graph[i];

      if ((graph_phenos[c] == 0) && (m.vector[i] == -1))
         {
         graph_phenos[c] = -1;

         PrepareComponent(m, graph, c, pivot);
         graph_likelihood[c] = ComponentLikelihood(m, graph, c);

         #ifdef __DEBUG_ERROR_GRAPH__
         DebugGraph(m, graph, c, i);
         #endif
         }
      }
   }

double LikelihoodTree::ComponentLikelihood(Mantra & m, IntArray & graph, int component)
   {
   IntArray & original = activeComponent;

   int alleles = original.Length();

   // Initialize allele states
   allele.Dimension(m.two_n);
   for (int i = 0; i < alleles; i++)
      allele[original[i]] = 0;

   double lk = 0.0;

   while (true)
      {
      // Evaluate the likelihood for the current allele state
      lk += PartialLikelihood(m, graph, component);

      // Cycle through to the next possible allelic state
      int i = 0;
      while (i < alleles && allele[original[i]] == 1)
         i++;

      // If there are no more possible allelic states, we are finished
      if (i == alleles)
         break;

      // Otherwise keep iterating
      for (int j = 0; j < i; j++)
         allele[original[j]] = 0;

      allele[original[i]] = 1;
      }

   return lk;
   }

double LikelihoodTree::PartialLikelihood(Mantra & m, IntArray & graph, int component)
   {
   double  lk = 1.0;

   // Prior probability of founder allele states
   for (int i = 0; i < activeComponent.Length(); i++)
      lk *= (allele[activeComponent[i]] == 0) ? 1.0 - model.p : model.p;

   // Probability of observed genotypes conditional on founder allele states
   for (int j = 0; j < cPheno.Length(); j++)
      {
      int pheno = cPheno[j];
      int geno = allele[cState[j * 2]] + allele[cState[j * 2 + 1]];

      if (pheno == 1)
         lk *= 1.0 - model.pen[geno];
      else
         lk *= model.pen[geno];
      }

   return lk;
   }

void LikelihoodTree::PrepareComponent(Mantra & m, IntArray & graph,
      int c, int pivot)
   {
   // Setup temporary list of component descent graphs
   activeComponent = components[c];

   // List phenotyped descendants
   cPheno.Clear(); cState.Clear();
   for (int i = 0; i < pivot; i += 2)
      if ((m.aff[i]) && (graph[m.state[i]] == c))
         {
         cPheno.Push(m.aff[i]);
         cState.Push(m.state[i]);
         cState.Push(m.state[i+1]);
         }
   }

void LikelihoodTree::AllocateStack(Mantra & m)
   {
   // Work out a reasonable stack depth
   int new_depth = m.bit_count;

   if (m.two_n > new_depth)
      new_depth = m.two_n;

   // new_depth = (new_depth & ~31) + 32;

   // Allocate stack if necessary
   if (new_depth > stack_depth)
      {
      FreeStack();

      stack_depth = new_depth;

      components = new IntArray [stack_depth];
      graphs = new IntArray [stack_depth];
      graph_phenotypes = new IntArray [stack_depth];

      graph_likelihood.Dimension(stack_depth);
      }

   // Initialize graph
   graphs[0].Dimension(m.two_f);
   graphs[0].SetSequence(0, 1);

   // Initialize components for each singleton graph
   for (int i = 0; i < m.two_f; i++)
      {
      components[i].Dimension(1);
      components[i][0] = m.vector[i] >= 0 ? -1 : i;
      }

   graph_phenotypes[0] = m.branch_affected;

   // Setup graphs and founder components for genotyped founders
   for (int i = 0; i < m.two_f; i+=2)
      if (m.aff[i])
         {
         if (m.vector[i + 1] == -1)
            {
            graphs[0][i + 1] = graphs[0][i];
            components[i].Append(components[i + 1]);
            graph_phenotypes[0][i] += graph_phenotypes[0][i + 1] - 2;
            }
         else
            graph_phenotypes[0][i] -= 2;

         }
   }

void LikelihoodTree::FreeStack()
   {
   if (stack_depth)
      {
      // Delete temporary variables from stack
      delete [] components;
      delete [] graphs;
      delete [] graph_phenotypes;
      }
   }

void DiseaseModel::Copy(DiseaseModel & source)
   {
   p = source.p;
   pen[0] = source.pen[0];
   pen[1] = source.pen[1];
   pen[2] = source.pen[2];
   affection = source.affection;
   }

bool DiseaseModel::CheckModel()
   {
   return  !(affection < 0 || affection > Pedigree::affectionCount ||
             p <= 0.0 || p >= 1.0 || pen[0] < 0.0 || pen[0] > 1.0 ||
             pen[1] < 0.0 || pen[1] > 1.0 || pen[2] < 0.0 || pen[2] > 1.0);
   }

 
