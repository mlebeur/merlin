////////////////////////////////////////////////////////////////////// 
// merlin/Magic.cpp 
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
 
#include "Magic.h"
#include "Mantra.h"

#include <math.h>

int InheritanceTree::mergingStrategy = MERGE_ALL;

InheritanceTree::InheritanceTree()
   {
   stack_depth = 0;
   }

InheritanceTree::~InheritanceTree()
   {
   if (stack_depth)
      {
      delete [] graph_stack;
      delete [] dirty_stack;
      delete [] alleles_stack;
      delete [] alleles2_stack;
      delete [] fixed_stack;
      }
   }

void InheritanceTree::ScoreVectors(Mantra & m)
   {
   // Clear the tree and define its depth
   Clear();
   bit_count = m.bit_count;

   // Setup up basic graphs and allele states, using founder genotypes
   AllocateStack(m.bit_count);

   LongArray & alleles = alleles_stack[m.bit_count];
   LongArray & alleles2 = alleles2_stack[m.bit_count];
   IntArray & graph = graph_stack[m.bit_count];
   IntArray & fixed = fixed_stack[m.bit_count];
   IntArray & dirty = dirty_stack[m.bit_count];

   alleles.Dimension(m.two_f);
   alleles2.Dimension(m.two_f);
   graph.Dimension(m.two_f);
   fixed.Dimension(m.two_f);
   dirty.Dimension(m.two_f);

   if (!PrepareForGraphing(m, graph, alleles, fixed, alleles2, dirty))
      {
      // Problem detected among founder genotypes [e.g. male heterozygous for X]
      MakeEmptyTree(bit_count);
      return;
      }

   // The work horse
   ScoreRecursive(m, m.bit_count - 1, m.two_f,
                  graph, alleles, fixed, alleles2, dirty);

#ifndef __PRINTOUT_STATS__
   // Compact form for calculations
   switch (mergingStrategy)
      {
      case MERGE_ALL :
         Trim(); break;
      case MERGE_ZEROS :
         TrimNoMerge(); break;
      case MERGE_BOOLEAN :
         MakeBooleanTree(); break;
      }
#endif
   }

bool InheritanceTree::PrepareForGraphing
     (Mantra & m, int * graph, longint * alleles,
      int * fixed, longint * alleles2, int * dirty)
   {
   // Each allele is a bitset listing possible states for each founder allele
   // and corresponds to the founder genotype by default

   // Each graph includes alleles that pass through the same individual,
   // and includes a single element (for ungenotyped founders) or the
   // two founder alleles (for genotyped founders) by default

   // Graphs that descent from homozygous founders are fixed by default,
   //

   for (int i = 0; i < m.two_f; i+=2)
      {
      alleles[i] = alleles[i + 1] = m.genotype[i];
      dirty[i]   = dirty[i+1]     = false;
      graph[i]   = i;
      graph[i+1] = m.genotype[i] != NOTZERO ? i : i + 1;
      fixed[i]   =
      fixed[i+1] = !m.hetero[i] && (m.genotype[i] != NOTZERO);
      // For heterozygous individuals, initialize the alleles2 array
      if (m.hetero[i])
         // Make sure that we this isn't a male "heterozygous" for X
         if (m.vector[i + 1] == -1)
            {
            alleles2[i] = alleles[i] & (alleles[i] - 1);
            alleles2[i + 1] = alleles[i] ^ alleles2[i];
            }
         // If it is, we have probably encountered a genotyping error...
         else
            return false;
      }

   // Temporary storage for FindRedundancy and CalculateLikelihood functions
   fixed2.Dimension(m.two_f);
   stack.Dimension(m.two_f);

   return true;
   }

int InheritanceTree::ScoreRecursive
    (Mantra & m, int bit, int last_pivot, IntArray & graph, LongArray & alleles,
     IntArray & fixed, LongArray & alleles2, IntArray & dirty)
   {
   // Allocate a new tree node
   int node = NewNode();

   // Where in the inheritance vector are we now?
   int pivot = bit >= 0 ? m.bits[bit] : m.two_n;

   // Has this bit been masked?
   if (bit >= 0 && m.bit_mask[bit])
      {
      m.state[pivot] = m.state[m.vector[pivot] &= ~1];

      nodes[node].type = TREE_NODE_ONE;
      SAFE_SET(nodes[node].child[0],
               ScoreRecursive(m, bit - 1, last_pivot, graph,
                              alleles, fixed, alleles2, dirty));
      return node;
      }

   if (bit != m.bit_count - 1)
      for (int i = last_pivot; i < pivot; i++)
         m.state[i] = m.state[m.vector[i]];

   // Update the graph up to current pivot
   for ( ; last_pivot < pivot; last_pivot++)
      {
      longint genotype = m.genotype[last_pivot];
      int state = m.state[last_pivot];

      // Nothing to do for individuals who are not genotyped
      if (genotype == NOTZERO) continue;

      // Work out what graph this allele belongs to
      int graph_left = state;

      while (graph[graph_left] != graph_left)
         graph_left = graph[graph_left];

      // And if it needs recalculating
      dirty[graph_left] |= m.hetero[last_pivot] || alleles[state] != genotype;

      // list possible alleles at each locus
      if ((alleles[state] &= genotype) == 0) break;

      fixed[state] |= !m.hetero[last_pivot] || alleles[state] != genotype;

      // update graph -- only need to do this once per genotyped person
      if (last_pivot & 1)
         {
         int graph_right = m.state[last_pivot - 1];

         while (graph[graph_right] != graph_right)
            graph_right = graph[graph_right];

         graph[graph_right] = graph_left;
         dirty[graph_left] |= dirty[graph_right];
         }
      }

   // Is this an impossible state?
   if (last_pivot < pivot)
      {
      nodes[node].type  = TREE_NODE_ZERO;
      nodes[node].value = 0.0;

      return node;
      }

   // Is this a leaf node?
   if (bit < 0)
      {
      nodes[node].type = TREE_NODE_LEAF;
      nodes[node].value =
         CalculateLikelihood(m, graph, alleles, fixed, alleles2, dirty);

#ifdef __FUGUE__
      // Fugue estimates haplotype frequencies in the presence of LD
      // and does not require the probability of genotypes given gene
      // flow to be calculated at a single marker
      if (nodes[node].value)
         nodes[node].value = 1.0;
#endif

      return node;
      }

   // Find out if we can save time by using symmetries in the data
   switch (nodes[node].type =
           FindRedundancy(m, pivot, graph, alleles, fixed, alleles2, dirty))
      {
      case TREE_NODE_ZERO :
         nodes[node].value = 0.0;
         break;
      case TREE_NODE_ONE :
         m.state[pivot] = m.state[m.vector[pivot] &= ~1];
         SAFE_SET(nodes[node].child[0],
                  ScoreRecursive(m, bit - 1, last_pivot, graph,
                                 alleles, fixed, alleles2, dirty));
         break;
      case TREE_NODE_TWO :
         graph_stack[bit] = graph; dirty_stack[bit] = dirty;
         alleles_stack[bit] = alleles; fixed_stack[bit] = fixed;
         alleles2_stack[bit] = alleles2;

         m.state[pivot] = m.state[m.vector[pivot] &= ~1];
         SAFE_SET(nodes[node].child[0],
                  ScoreRecursive(m, bit - 1, last_pivot,
                                 graph, alleles, fixed, alleles2, dirty));

         m.state[pivot] = m.state[m.vector[pivot] |= 1];
         SAFE_SET(nodes[node].child[1],
                  ScoreRecursive(m, bit - 1, last_pivot,
                                graph_stack[bit], alleles_stack[bit],
                                fixed_stack[bit], alleles2_stack[bit],
                                dirty_stack[bit]));
         break;
      }

   return node;
   }

int InheritanceTree::FindRedundancy
    (Mantra & m, int pivot, int * graph, longint * alleles,
     int * fixed, longint * alleles2, int * dirty)
   {
   int limit = pivot & ~1;

   // Squash the graph
   for (int i = 0; i < m.two_f; i++)
      {
      int head = i;
      while (graph[head] != head) head = graph[head];
      graph[i] = head;
      }

   // Do the graphing...
   // Loop over all non-trivial graphs (eg. with multiple alleles)
   for (int i = 0, j; i < m.two_f; i++)
      if (graph[i] == i && alleles[i] != NOTZERO && dirty[i])
         {
         // Don't need to retrace this graph until it changes
         dirty[i] = false;

         // sp is the stack pointer for our stack of unresolved alleles
         int sp = 0;

         // If at least one allele is known we should be able
         // to fix all alleles in the graph
         for (j = 0; j < m.two_f; j++)
            if (graph[j] == i)
               {
               stack[sp++] = j;

               while (sp && fixed[stack[sp - 1]])
                  {
                  int state = stack[--sp];
                  longint allele = alleles[state];

                  // Since k is known, all alleles in heterozygotes
                  // with k are also known...
                  for (int k = 0, other; k < limit; k += 2)
                     if (m.hetero[k])
                        {
                        if (m.state[k] == state)
                           other = m.state[k + 1];
                        else if (m.state[k + 1] == state)
                           other = m.state[k];
                        else continue;

                        if ((alleles[other] &= ~allele) == 0)
                           return TREE_NODE_ZERO;

                        fixed[other] = true;
                        }

                  // Bring alleles that are fixed so far to the fore
                  int catch22 = sp;
                  while ( --catch22 >= 0 )
                     if (fixed[stack[catch22]])
                        {
                        int tmp = stack[catch22];
                        stack[catch22] = stack[sp - 1];
                        stack[sp - 1] = tmp;
                        break;
                        }
                  }
               }

         if (sp)
            // if we get here, we have no homozygotes
            {
            // We won't be able to really call any alleles, so
            // make a scratch copy of the allele states tables

            for (int j = 0; j < m.two_f; j++)
               if (graph[j] == i)
                  {
                  alleles2[j] = alleles[i];
                  fixed2[j] = false;
                  }

            // Choose a random allele at one marker and see what
            // happens
            int random = stack[sp - 1];
            alleles2[random] = alleles[random] & (alleles[random] - 1);

            // This should fix all other alleles
            // The last allele is guaranteed to be compatible with all
            // others, so consider sp - 1 graph components only
            while (--sp)
               {
               int state = stack[sp];
               longint allele = alleles2[state];

               // After fixing an allele, the state of all alleles
               // all alleles in heterozygotes who carry it is also
               // also known...
               for (int k = 0, other; k < limit; k += 2)
                  if (m.hetero[k])
                     {
                     if (m.state[k] == state)
                        other = m.state[k + 1];
                     else if (m.state[k + 1] == state)
                        other = m.state[k];
                     else continue;

                     if ((alleles2[other] &= ~allele) == 0)
                        return TREE_NODE_ZERO;

                     fixed2[other] = true;
                     }

               // Bring alleles that are fixed so far to the fore
               int catch22 = sp;
               while ( --catch22 >= 0 )
                  if (fixed2[stack[catch22]])
                     {
                     int tmp = stack[catch22];
                     stack[catch22] = stack[sp - 1];
                     stack[sp - 1] = tmp;
                     break;
                     }
               }
            }
         }

   // Are parents homozygous?
   int parent = m.vector[pivot];
   int matallele = m.state[parent & ~1];
   int fatallele = m.state[parent | 1];

   // Are the possible states of the parental alleles identical
   if ( alleles[matallele] == alleles[fatallele] )
      {
      // Then if one is fixed, the parent is homozygous
      if ( fixed[matallele] )
        return TREE_NODE_ONE;

      // Or if the parents are in the same graph, and
      // in the same state, they are also homozygous
      if ( graph[matallele] == graph[fatallele] &&
           alleles2[matallele] == alleles2[fatallele])
         return TREE_NODE_ONE;

      // Or if the two alleles are the same!
      if ( matallele == fatallele )
         return TREE_NODE_ONE;
      }

   // If all the genotyped descendants of the two alternative
   // founder alleles are on this branch, it doesn't matter which
   // we pick -- one will go down, the other will be fixed (or not)
   if (m.branch_genotyped[matallele] == m.branch_genotyped[pivot] &&
       m.branch_genotyped[fatallele] == m.branch_genotyped[pivot] &&
       mergingStrategy == MERGE_ALL)
       return TREE_NODE_ONE;

   return TREE_NODE_TWO;
   }

double InheritanceTree::CalculateLikelihood(Mantra & m, int * graph,
       longint * alleles, int * fixed, longint * alleles2, int * dirty )
   {
   // Squash the graph
   for (int i = 0; i < m.two_f; i++)
      {
      int head = i;
      while (graph[head] != head) head = graph[head];
      graph[i] = head;
      }

   double likelihood = 1.0;

   // Do the graphing...
   // Loop over all non-trivial graphs (eg. with multiple alleles)
   for (int i = 0, j; i < m.two_f; i++)
      if (graph[i] == i && alleles[i] != NOTZERO && (dirty[i] || !fixed[i]))
         {
         // sp is the stack pointer for stack of unresolved alleles
         int sp = 0;

         // If at least one allele is known we should be able
         // to fix all alleles in the graph
         for (j = 0; j < m.two_f; j++)
            // Find all founder alleles in this graph
            if (graph[j] == i && m.vector[j] == -1)
               {
               stack[sp++] = j;

               while (sp && fixed[stack[sp - 1]])
                  {
                  int state = stack[--sp];
                  longint allele = alleles[state];

                  // Since k is known, all alleles in heterozygotes
                  // with k are also known...
                  for (int k = 0, other; k < m.two_n; k += 2)
                     if (m.hetero[k])
                        {
                        if (m.state[k] == state)
                           other = m.state[k + 1];
                        else if (m.state[k + 1] == state)
                           other = m.state[k];
                        else continue;

                        if ((alleles[other] &= ~allele) == 0)
                           return 0.0;

                        fixed[other] = true;
                        }

                  // Bring alleles that are fixed so far to the fore
                  int catch22 = sp;
                  while ( --catch22 >= 0 )
                     if (fixed[stack[catch22]])
                        {
                        int tmp = stack[catch22];
                        stack[catch22] = stack[sp - 1];
                        stack[sp - 1] = tmp;
                        break;
                        }
                  }
               }

         if (sp)
            // if we get here, we have no homozygotes
            {
            // We won't be able to really call any alleles, so
            // make a scratch copy of the allele state tables

            for (int j = 0; j < m.two_f; j++)
               if (graph[j] == i)
                  {
                  alleles2[j] = alleles[i];
                  fixed2[j] = false;
                  }

            // There will either be zero or two possible states
            double likelihoods[2] = { 1.0, 1.0 };

            // Choose a random allele at one marker and see what
            // happens
            int random = stack[sp - 1];
            alleles2[random] = alleles[random] & (alleles[random] - 1);

            // This should fix all other alleles
            while (sp--)
               {
               int state = stack[sp], freq;
               longint allele = alleles2[state], mask;

               // Update likelihood based on current choice
               mask = allele; freq = 1;
               while ((bool)(mask & NOTONE)) { mask >>= 1; freq++; }
               likelihoods[0] *= m.frequencies[freq];

               // Update likelihood based on alternative choice
               mask = alleles[state] ^ allele; freq = 1;
               while ((bool)(mask & NOTONE)) { mask >>= 1; freq++; }
               likelihoods[1] *= m.frequencies[freq];

               // Since k is known, all alleles in heterozygotes
               // with k are also known...
               for (int k = 0, other; k < m.two_n; k += 2)
                  if (m.hetero[k])
                     {
                     if (m.state[k] == state)
                        other = m.state[k + 1];
                     else if (m.state[k + 1] == state)
                        other = m.state[k];
                     else continue;

                     if ((alleles2[other] &= ~allele) == 0)
                        return 0.0;

                     fixed2[other] = true;
                     }

               // Bring alleles that are fixed so far to the fore
               int catch22 = sp;
               while ( --catch22 >= 0 )
                  if (fixed2[stack[catch22]])
                     {
                     int tmp = stack[catch22];
                     stack[catch22] = stack[sp - 1];
                     stack[sp - 1] = tmp;
                     break;
                     }
               }
            likelihood *= likelihoods[0] + likelihoods[1];
            }
         }

   // Used fixed vectors to update likelihood
   for (int i = 0; i < m.two_f; i++)
      // Find all fixed alleles
      if (fixed[i] && m.vector[i] == -1)
         {
         int freq = 1;
         longint allele = alleles[i];
         while (bool(allele & NOTONE)) { allele >>= 1; freq++; }
         likelihood *= m.frequencies[freq];
         }

   return likelihood;
   }

double InheritanceTree::EvaluateInheritanceVector(Mantra & m, int * vector,
       int * graph, longint * alleles, longint * alleles2, int * fixed)
   {
   // Setup storage
   IntArray dirty(m.two_f);

   // Setup the inheritance vector
   for (int i = 0; i < m.bit_count; i++)
      {
      int bit = m.bits[i];

      m.vector[bit] = vector[i] == 0 ? m.vector[bit] & ~1 : m.vector[bit] | 1;
      }

   // Identify founder alleles for each individual
   for (int i = m.two_f; i < m.two_n; i++)
      m.state[i] = m.state[m.vector[i]];

   if (!PrepareForGraphing(m, graph, alleles, fixed, alleles2, dirty))
      // Problem in founder genotypes
      return 0.0;

   // Update the graph up to current pivot
   for (int counter = m.two_f; counter < m.two_n; counter++)
      {
      longint genotype = m.genotype[counter];
      int state = m.state[counter];

      if (genotype == NOTZERO) continue;

      // list possible alleles at each locus
      if ((alleles[state] &= genotype) == 0) return 0.0;

      fixed[state] |= !m.hetero[counter] || alleles[state] != genotype;

      // update graph -- only need to do this once per genotyped person
      // TODO -- Use squashing when updating graphs!
      if (counter & 1)
         {
         int graph_left  = state;
         int graph_right = m.state[counter - 1];

         while (graph[graph_left] != graph_left)
            graph_left = graph[graph_left];

         while (graph[graph_right] != graph_right)
            graph_right = graph[graph_right];

         graph[graph_right] = graph_left;
         }
      }

   dirty.Set(1);

   return CalculateLikelihood(m, graph, alleles, fixed, alleles2, dirty);
   }

void InheritanceTree::AllocateStack(int new_depth)
   {
   if (new_depth >= stack_depth)
      {
      // Free old storage if required
      if (stack_depth)
         {
         delete [] graph_stack;
         delete [] dirty_stack;
         delete [] alleles_stack;
         delete [] alleles2_stack;
         delete [] fixed_stack;
         }

      // Round up at 8 bits to limit reallocations
      stack_depth = (new_depth + 8) & ~7;

      // Allocate stacks
      graph_stack = new IntArray[stack_depth];
      dirty_stack = new IntArray[stack_depth];
      alleles_stack = new LongArray[stack_depth];
      alleles2_stack = new LongArray[stack_depth];
      fixed_stack = new IntArray[stack_depth];
      }
   }
 
