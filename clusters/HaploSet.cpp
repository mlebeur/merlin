////////////////////////////////////////////////////////////////////// 
// clusters/HaploSet.cpp 
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
 
#include "HaploSet.h"
#include "Error.h"
#include "StringArray.h"

#include <math.h>

int HaplotypeSet::maxHaploCount = 0;

HaplotypeSet::HaplotypeSet()
   {
   size = haploCount = graphCount = 0;
   graphs = NULL;
   weight = 0.0;
   copies = 1;
   }

HaplotypeSet::~HaplotypeSet()
   {
   if (graphs != NULL)
      {
      for (int i = 0; i < graphCount; i++)
         delete graphs[i];
      delete [] graphs;
      }
   }

void HaplotypeSet::Grow()
   {
   size = size == 0 ? 1 : size * 2;

   HaplotypeGraph ** new_graphs = new HaplotypeGraph * [size];

   for (int i = 0; i < graphCount; i++)
      new_graphs[i] = graphs[i];

   if (graphs != NULL) delete [] graphs;

   graphs = new_graphs;
   }

void HaplotypeSet::LoadFromMemory(FounderGraph * founderGraph, int markers)
   {
   HaplotypeGraph * graph = new HaplotypeGraph;

   int haplos = founderGraph->graph[0].Length();

   if (haploCount != 0 && haploCount != haplos)
      error("Expecting %d Haplotypes in Family %s\n",
            haploCount, (const char *) setid);

   haploCount = haplos;
   if (haplos > maxHaploCount)
      maxHaploCount = haploCount;

   graph->LoadFromMemory(founderGraph, haplos, markers);
   weight += graph->weight = founderGraph->weight;

   if (graphCount >= size) Grow();

   graphs[graphCount++] = graph;
   }

void HaplotypeSet::LoadFromMemory(FamilyHaplos * input, int markers)
   {
   setid = input->label;

   FounderGraph * graph = input->graphs;

   do {
      LoadFromMemory(graph, markers);
      graph = graph->next;
   } while (graph != NULL);
   }

void HaplotypeSet::StandardizeWeights()
   {
   if (weight != 1.0)
      {
      for (int i = 0; i < graphCount; i++)
         graphs[i]->weight /= weight;
      weight = 1.0;
      }
   }


double HaplotypeSet::Complexity(IntArray & alleleCounts)
   {
   double c = 0;

   for (int i = 0; i < graphCount; i++)
      c += graphs[i]->Complexity(alleleCounts);

   return c;
   }

HaplotypeSet * HaplotypeSets::LoadFromMemory(FamilyHaplos * family, int markers)
   {
   ::String input;

   HaplotypeSet * set = FindSet(family->label);

   set->setid = family->label;

   FounderGraph * graph = family->graphs;

   do {
      set->LoadFromMemory(graph, markers);
      graph = graph->next;
   } while (graph != NULL);

   return set;
   }

void HaplotypeSets::AllocateAlleleLabels(int markers)
   {
   if (alleleLabels != NULL)
      delete [] alleleLabels;

   alleleLabels = new StringArray[markers];
   }

void HaplotypeSets::SetAlleleLabels(int marker, StringArray & labels)
   {
   alleleLabels[marker].Dimension(labels.Length());

   for (int i = 0; i < labels.Length(); i++)
      alleleLabels[marker][i] = labels[i];
   }

String HaplotypeSets::GetAlleleLabel(int marker, int allele)
   {
   if (alleleLabels == NULL || alleleLabels[marker].Length() == 0)
      {
      ::String label;
      return label = allele;
      }

   return alleleLabels[marker][allele];
   }

void * HaplotypeSets::create_set()
   {
   return new HaplotypeSet;
   }

HaplotypeSets::HaplotypeSets()
   {
   alleleLabels = NULL;
   quiet = false;
   }

HaplotypeSets::~HaplotypeSets()
   {
   for (int i = 0; i < Length(); i++)
      delete (HaplotypeSet *) objects[i];

   if (alleleLabels != NULL)
      delete [] alleleLabels;
   }

void HaplotypeSets::ShowCounts(IntArray & alleleCounts)
   {
   int totalHaplotypes = 0;
   int totalSets = 0;
   int knownHaplotypes = 0;
   int knownSets = 0;

   for (int i = Length() - 1; i >= 0; i--)
      {
      HaplotypeSet * set = (HaplotypeSet *) objects[i];

      if (set->Complexity(alleleCounts) == 1.0)
         {
         knownHaplotypes += set->haploCount;
         knownSets++;
         }

      totalHaplotypes += set->haploCount;
      totalSets++;
      }

   printf("Total: %d Haplotypes in %d Sets\n"
          "Known: %d Haplotypes in %d Sets [UNDERESTIMATE]\n\n",
          totalHaplotypes, totalSets,
          knownHaplotypes, knownSets);
   }

void HaplotypeSets::Filter(IntArray & alleleCounts, int maxBits)
   {
   int totalHaplotypes = 0;
   int totalSets = 0;
   int knownHaplotypes = 0;
   int knownSets = 0;

//   printf("Filtering data...\n");
   for (int i = Length() - 1; i >= 0; i--)
      {
      HaplotypeSet * set = (HaplotypeSet *) objects[i];

      double c = set->Complexity(alleleCounts);
      double bits = log(c) / log(2.0);

      if (bits > maxBits)
         {
         if (!quiet)
            printf("DISCARDING family %s (%.1f bits)\n",
                   (const char *) set->setid, bits);
         else
            discarded.Add(set->setid);

         Delete(i);
         delete set;
         continue;
         }

      if (bits == 0.0)
         {
         knownHaplotypes += set->haploCount;
         knownSets++;
         }

      totalHaplotypes += set->haploCount;
      totalSets++;
      }

/*   printf("\nTotal: %d Haplotypes in %d Sets\n"
          "Known: %d Haplotypes in %d Sets [UNDERESTIMATE]\n\n",
          totalHaplotypes, totalSets,
          knownHaplotypes, knownSets); */
   }

void HaplotypeSets::ClearContents()
   {
   for (int i = 0; i < Length(); i++)
      delete (HaplotypeSet *) objects[i];

   if (alleleLabels != NULL)
      {
      delete [] alleleLabels;
      alleleLabels = NULL;
      }

   discarded.Clear();
   Clear();
   }
 
