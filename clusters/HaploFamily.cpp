////////////////////////////////////////////////////////////////////// 
// clusters/HaploFamily.cpp 
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
 
#include "HaploFamily.h"
#include "MerlinCore.h"
#include "Error.h"
#include "Magic.h"

// Static variable definitions
//

int FamilyHaplos::maxBits = 24;

// Family Haplo routines
//

FamilyHaplos::FamilyHaplos() : table(8)
   {
   graphs = NULL;
   graph = fixed = NULL;
   alleles = alleles2 = NULL;
   allocatedMarkers = 0;
   }

FamilyHaplos::~FamilyHaplos()
   {
   DeleteGraphs();

   if (graph != NULL) delete [] graph;
   if (fixed != NULL) delete [] fixed;
   if (alleles != NULL) delete [] alleles;
   if (alleles2 != NULL) delete [] alleles2;
   }

void FamilyHaplos::DeleteGraphs()
   {
   // First delete chain of founder allele graphs
   // We don't do this recursively to avoid stack
   // overflows
   FounderGraph * current = graphs;

   while (current != NULL)
      {
      FounderGraph * next = current->next;
      delete current;
      current = next;
      }
   }

String & FamilyHaplos::RetrieveGenotypes(String & genotypes, Person & p)
   {
   genotypes.Clear();

#ifdef __CHROMOSOME_X__
   genotypes.catprintf("%c\n", p.sex == 1 ? 'M' : 'F');
#endif

   for (int i = 0; i < markers.Length(); i++)
      {
      int m = markers[i];

      genotypes.catprintf("%d,", p.markers[m].SequenceCoded());
      }

   return genotypes;
   }

String & FamilyHaplos::RetrieveGenotypes(String & genotypes, Family & f)
   {
   genotypes.Clear();

   for (int j = 0; j < f.count; j++)
      {
      Person & p = f.ped[f.path[j]];

      #ifdef __CHROMOSOME_X__
         genotypes += p.sex == 1 ? 'M' : 'F';
      #endif

      if (!p.isFounder())
         genotypes.catprintf("(%d,%d,%d),", p.serial, p.father->serial, p.mother->serial);
      else
         genotypes.catprintf("(%d)", p.serial);

      for (int i = 0; i < markers.Length(); i++)
         {
         int m = markers[i];

         genotypes.catprintf("%d,", p.markers[m].SequenceCoded());
         }

      genotypes += '\n';
      }

   return genotypes;
   }

bool FamilyHaplos::Haplotype(Pedigree & ped, Person & p)
   {
   Tree joint;

   return Haplotype(joint, ped, p);
   }

bool FamilyHaplos::Haplotype(Tree & joint, Pedigree & ped, Person & p)
   {
#ifdef __CHROMOSOME_X__
   if (p.sex == SEX_MALE)
      return HaplotypeMaleX(joint, ped, p);
#endif

   label = p.famid + "->" + p.pid;

   AllocateMemory(2, 0);

   bool first_heterozygote = false;
   for (int i = 0; i < markers.Length(); i++)
      {
      int m = markers[i];

      if (p.markers[m].isKnown())
         {
         if (p.markers[m].isHeterozygous())
            {
            if (first_heterozygote)
               {
               graph[i].SetSequence(0, 1);
               alleles[i][0] = 1 << (p.markers[m].Lo() - 1);
               alleles[i][1] = 1 << (p.markers[m].Hi() - 1);
               alleles2[i].Zero();
               fixed[i].Set(1);
               first_heterozygote = false;
               }
            else /* multiple heterozygous genotypes */
               {
               graph[i].Zero();
               alleles[i][0] = alleles[i][1] = p.markers[m].BinaryCoded();
               alleles2[i][0] = 1 << (p.markers[m].Hi() - 1);
               alleles2[i][1] = alleles[i][0] ^ alleles2[i][0];
               fixed[i].Zero();
               }
            }
         else /* Homozygous genotype */
            {
            graph[i].SetSequence(0, 1);
            alleles[i].Set(p.markers[m].BinaryCoded());
            alleles2[i].Zero();
            fixed[i].Set(1);
            }
         }
      else /* Missing genotype */
         {
         graph[i].SetSequence(0, 1);
         alleles[i].Set(NOTZERO);
         alleles2[i].Set(NOTZERO);
         fixed[i].Zero();
         }
      }

   joint.MakeMinimalTree(1.0, 0);
   joint.nodes[0].info = (void *)
      graphs->Append(graph, fixed, alleles, alleles2, first_heterozygote ? 1.0 : 2.0);

   return true;
   }

bool FamilyHaplos::HaplotypeMaleX(Tree & joint, Pedigree & ped, Person & p)
   {
   label = p.famid + "->" + p.pid;

   AllocateMemory(1, 0);

   for (int i = 0; i < markers.Length(); i++)
      {
      int m = markers[i];

      if (p.markers[m].isKnown())
         {
         if (p.markers[m].isHomozygous())
            {
            /* Known haplotype */
            graph[i][0] = 0;
            alleles[i].Set(p.markers[m].BinaryCoded());
            alleles2[i][0] = 0;
            fixed[i][0] = 1;
            }
         else
            {
            /* Discard heterozygous males */
            mendel_errors.Push(m);
            return false;
            }
         }
      else
         {
         /* Missing genotype */
         graph[i][0] = 0;
         alleles[i][0] = NOTZERO;
         alleles2[i][0] = NOTZERO;
         fixed[i][0] = 0;
         }
      }

   joint.MakeMinimalTree(1.0, 0);
   joint.nodes[0].info = (void *)
      graphs->Append(graph, fixed, alleles, alleles2, 1.0);

   return true;
   }

bool FamilyHaplos::Haplotype(Pedigree & ped, Family & f)
   {
   Tree joint;

   return Haplotype(joint, ped, f);
   }

bool FamilyHaplos::Haplotype(Tree & joint, Pedigree & ped, Family & f)
   {
   mendel_errors.Clear();

   // Special purpose routines for haplotyping single individuals
   if (f.count == 1)
      return Haplotype(joint, ped, ped[f.first]);

   label = f.famid;

   m.Prepare(ped, f);

   if (m.bit_count > maxBits)
      return false;

   InheritanceTree tree;

   int oldStrategy = tree.mergingStrategy;
   tree.mergingStrategy = MERGE_BOOLEAN;

   m.SelectMarker(markers[0]);
   tree.ScoreVectors(m);
   joint.Copy(tree);

   if (!tree.FindNonZero())
      mendel_errors.Push(markers[0]);

   for (int i = 1; i < markers.Length(); i++)
      {
      m.SelectMarker(markers[i]);
      tree.ScoreVectors(m);
      joint.Multiply(tree);

      if (!tree.FindNonZero())
         mendel_errors.Push(markers[0]);
      }

   obligate_recombinant = !joint.FindNonZero();

   AllocateMemory(m.two_f, m.bit_count);

   ListAllRecursively(joint, 0, m.bit_count - 1, 1);
   table.Clear();

#if __DEBUG__
   if (graphs->weight)
      {
      int    list = 0;
      double sum  = 0;
      for (FounderGraph * gr = graphs; gr != NULL; gr = gr->next)
         list++, sum += gr->weight * (2 << (m.bit_count - 1));

      printf("%10s %2d %2d %5.0f\n",
             (const char *) f.famid, m.bit_count, list, sum);

      return list;
      }
#endif

   tree.mergingStrategy = oldStrategy;

   return graphs->weight > 0;
   }

bool FamilyHaplos::Haplotype(Pedigree & ped, Family & f, IntArray & inheritance)
   {
   // Select the current family
   m.Prepare(ped, f);

   // Allocate auxiliary storage
   AllocateMemory(m.two_f, m.bit_count);

   // Select the appropriate segregation pattern
   vector = inheritance;

   // Retrieve haplotype graph
   if (!RetrieveGraph(true))
      return false;
   graphs->Append(graph, fixed, alleles, alleles2, 1.0);

   // Check if we were successful
   return graphs->weight > 0;
   }

void FamilyHaplos::FreeMemory()
   {
   if (graph != NULL) delete [] graph;
   if (fixed != NULL) delete [] fixed;
   if (alleles != NULL) delete [] alleles;
   if (alleles2 != NULL) delete [] alleles2;
   }

void FamilyHaplos::AllocateMemory(int haplos, int bits)
   {
   int mrk = markers.Length();

   if (mrk > allocatedMarkers)
      {
      FreeMemory();

      graph = new IntArray [mrk];
      fixed = new IntArray [mrk];
      alleles = new LongArray [mrk];
      alleles2 = new LongArray [mrk];

      allocatedMarkers = mrk;
      }

   for (int i = 0; i < mrk; i++)
      graph[i].Dimension(haplos),
      fixed[i].Dimension(haplos),
      alleles[i].Dimension(haplos),
      alleles2[i].Dimension(haplos);

   vector.Dimension(bits);
   vector.Zero();

   if (graphs != NULL)
      DeleteGraphs();

   graphs = new FounderGraph(markers.Length());
   }

void FamilyHaplos::RetrieveLikelihoods(Tree & tree)
   {
   RetrieveLikelihoods(tree, 0, tree.bit_count - 1);
   }

void FamilyHaplos::ListAllRecursively(Tree & tree, int node, int bit, double prob)
   {
   switch (tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         RetrieveGraph();
         tree.nodes[node].info = (void *)
            graphs->Append(graph, fixed, alleles, alleles2, prob, &table);
         break;
      case TREE_NODE_ONE :
         vector[bit] = 0;
         ListAllRecursively(tree, tree.nodes[node].child[0], bit-1, prob);
         break;
      case TREE_NODE_TWO :
         vector[bit] = 0;
         ListAllRecursively(tree, tree.nodes[node].child[0], bit-1, prob * .5);
         vector[bit] = 1;
         ListAllRecursively(tree, tree.nodes[node].child[1], bit-1, prob * .5);
         break;
      }
   }

void FamilyHaplos::RetrieveLikelihoods(Tree & tree, int node, int bit)
   {
   switch (tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         tree.nodes[node].value = *((FounderGraph *) tree.nodes[node].info)->likelihood;
         break;
      case TREE_NODE_ONE :
         vector[bit] = 0;
         RetrieveLikelihoods(tree, tree.nodes[node].child[0], bit-1);
         break;
      case TREE_NODE_TWO :
         vector[bit] = 0;
         RetrieveLikelihoods(tree, tree.nodes[node].child[0], bit-1);
         vector[bit] = 1;
         RetrieveLikelihoods(tree, tree.nodes[node].child[1], bit-1);
         break;
      }
   }

bool FamilyHaplos::RetrieveGraph(bool canFail)
   {
   for (int mrk = 0; mrk < markers.Length(); mrk++)
      {
      m.SelectMarker(markers[mrk]);

      double likelihood = engine.EvaluateInheritanceVector(
             m, vector, graph[mrk], alleles[mrk], alleles2[mrk], fixed[mrk]);

      // The selected inheritance vector is incompatible with this marker
      // (either we made a mistake, or there is bad inheritance)
      if (likelihood == 0.0)
         if (canFail)
            return false;
         else
            error("Invalid state encountered haplotyping marker %s",
                  (const char *) Pedigree::markerNames[markers[mrk]]);

      #ifdef __CHROMOSOME_X__
      // Edit out the second male allele
      int out = 0;
      for (int in = 0; in < m.two_f; in++)
         if (m.vector[in] == -1)
            {
            graph[mrk][out] = graph[mrk][in];
            alleles[mrk][out] = alleles[mrk][in];
            alleles2[mrk][out] = alleles2[mrk][in];
            fixed[mrk][out] = fixed[mrk][in];
            out++;
            }

      graph[mrk].Dimension(out);
      alleles[mrk].Dimension(out);
      alleles2[mrk].Dimension(out);
      fixed[mrk].Dimension(out);
      #endif
      }

   return true;
   }

// Allele Graph Routines
//

FounderGraph::FounderGraph(int mrk)
   {
   // Key Member Variables
   markers = mrk;
   weight  = 0;

   // Dynamic Arrays
   alleles = alleles2 = NULL;
   graph   = NULL;
   next    = NULL;

   // Pointer to likelihood
   likelihood = NULL;
   }

FounderGraph::~FounderGraph()
   {
   if (alleles != NULL) delete [] alleles;
   if (alleles2 != NULL) delete [] alleles2;
   if (graph != NULL) delete [] graph;
   }

int FounderGraph::Length()
   {
   FounderGraph * pointer = this;
   int length = 0;

   while (pointer->next != NULL)
      {
      pointer = pointer->next;
      length++;
      }

   return length;
   }

FounderGraph * FounderGraph::Append(
                  IntArray * gr, IntArray * fix,
                  LongArray * al, LongArray * al2, double prior,
                  BasicHash * table)
   {
   // First we create a standardized graph, for easy comparison
   // A standard graph has alleles of known state in singleton components
   // and other components named after the lowest numbered founder allele
   // they contain
   IntArray * standard_graph = new IntArray[markers];
#ifdef __CHROMOSOME_X__
   IntArray components(gr[0].Length() * 2);
#else
   IntArray components(gr[0].Length());
#endif

   for (int mrk = 0; mrk < markers; mrk++)
      {
      standard_graph[mrk].Dimension(gr[mrk].Length());
      components.SetSequence(0, 1);

      for (int i = 0; i < gr[mrk].Length(); i++)
         if (fix[mrk][i] || al[mrk][i] == NOTZERO)
            {
            standard_graph[mrk][i] = -1;
            al2[mrk][i] = 0;
            }
         else
            {
            al[mrk][i] ^= al2[mrk][i];

            if (components[gr[mrk][i]] <= i)
               standard_graph[mrk][i] = components[gr[mrk][i]];
            else
               standard_graph[mrk][i] = components[gr[mrk][i]] = i;
            }
      }

   FounderGraph * current = this;

   int h = Hash(standard_graph, al, al2);

   if (current->graph != NULL)
      if (table == NULL)
         while (true)
            {
            if (current->Compare(standard_graph, al, al2))
               {
               delete [] standard_graph;
               current->weight += prior;
               return current;
               }

            if (current->next == NULL)
               {
               current->next = new FounderGraph(markers);
               current = current->next;
               break;
               }

            current = current->next;
            }
      else
         {
         int pos = table->Find(h);

         while (pos >= 0)
            {
            current = (FounderGraph *) (*table)[pos];

            if (current->Compare(standard_graph, al, al2))
               {
               delete [] standard_graph;
               current->weight += prior;
               return current;
               }

            pos = table->Rehash(h, pos);
            }

         current = new FounderGraph(markers);
         current->next = this->next;
         this->next = current;

         table->Add(h, current);
         }

   if (current == this && table != NULL)
      table->Add(h, current);

   current->graph  = standard_graph;
   current->weight = prior;

   current->alleles = new LongArray [markers];
   current->alleles2 = new LongArray [markers];
   for (int mrk = 0; mrk < markers; mrk++)
      {
      current->alleles[mrk] = al[mrk];
      current->alleles2[mrk] = al2[mrk];
      }

   return current;
   }

bool FounderGraph::Compare(IntArray * gr, LongArray * al, LongArray * al2)
   {
   for (int mrk = 0; mrk < markers; mrk++)
      {
      if (graph[mrk] != gr[mrk]) return false;
      if (alleles[mrk] != al[mrk]) return false;
      if (alleles2[mrk] != al2[mrk]) return false;
      }
   return true;
   }

int FounderGraph::Hash(IntArray * gr, LongArray * al, LongArray * al2)
   {
   int h = 0;

   for (int mrk = 0; mrk < markers; mrk++)
      {
      h = gr[mrk].Hash(h);
      h = al[mrk].Hash(h);
      h = al2[mrk].Hash(h);
      }
   return h;
   }

int FounderGraph::GetAllele(int mrk, int founder)
   {
   longint binary = alleles[mrk][founder];
   longint mask = 1;
   int allele = 1;

   while ((mask & binary) == 0)
      {
      mask <<= 1;
      allele ++;
      }

   return allele;
   }

int FounderGraph::GetAllele2(int mrk, int founder)
   {
   longint binary = alleles2[mrk][founder];
   longint mask = 1;
   int allele = 1;

   while ((mask & binary) == 0)
      {
      mask <<= 1;
      allele ++;
      }

   return allele;
   }


 
