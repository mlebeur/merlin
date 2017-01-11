////////////////////////////////////////////////////////////////////// 
// clusters/SparseLikelihood.cpp 
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
 
#include "SparseLikelihood.h"
#include "MathConstant.h"
#include "Random.h"
#include "Error.h"

#include <math.h>

#define    SCALE_UP     1e7  //pow(2.0, 128)
#define    SCALE_DOWN   1e-7 //pow(2.0, -128)

void SparseLikelihood::Initialize(const IntArray & alleleCounts)
   {
   if (trees != NULL)
      delete [] trees;
            
   markers = alleleCounts.Length();

   treeCount = markers;
   treeBits = 0;
   trees = new HaploTree [treeCount];

   for (int i = 0; i < markers; i++)
      trees[i].MakeShrub(alleleCounts[i]);

   if (index != NULL) delete [] index;
   index = new IntArray[HaplotypeSet::maxHaploCount];

   best_llk = 0;
   }

void SparseLikelihood::Initialize(const IntArray & counts, const IntArray  & haplos, const Vector & freqs)
   {
   if (trees != NULL)
      delete [] trees;

   markers = counts.Length();

   treeCount = 1;
   treeBits = 31;

   trees = new HaploTree [1];
   trees[0].MakeEmptyTree(counts);

   for (int base = 0, f = 0; base < haplos.Length(); base += markers)
      {
      int branch = 0;

      for (int i = 0; i < markers; i++)
         branch = trees[0].GetBranch(branch, haplos[base + i] - 1);

      if (branch != trees[0].freqs.Length())
         error("Internal error while organizing haplotype frequencies\n");

      trees[0].freqs.Push(freqs[f++]);
      }

   if (index != NULL) delete [] index;
   index = new IntArray[HaplotypeSet::maxHaploCount];

   best_llk = 0;
   }

void SparseLikelihood::Print(double minfreq, HaplotypeSets * sets)
   {
   String label, allele;

   printf("Haplotypes with estimated frequency > %g\n", minfreq);

   if (sets != NULL)
      trees[0].Print(sets->alleleLabels, minfreq);
   else
      trees[0].Print(minfreq);

   printf("The logLikelihood of the data is %.4f\n\n", best_llk);
   }

void SparseLikelihood::EM(HaplotypeSets * sets, bool quiet)
   {
   CoreEM(sets, quiet);

   if (!quiet) printf("\n");
   }

void SparseLikelihood::CoreEM(HaplotypeSets * sets, bool quiet)
   {
   int count = 0;
   for (int i = 0; i < sets->Length(); i++)
      {
      HaplotypeSet * set = (HaplotypeSet *) sets->Object(i);
      count += set->haploCount * set->copies;
      set->StandardizeWeights();
      }

   if (count == 0) return;

   double llk0, llk1;

//   double families_llk0, unrel_llk0;
   double families_llk1, unrel_llk1;

   while (true)
      {
      if (!quiet)
         printf("Using %d set%s%s %d marker%s...\n",
             treeCount,
             treeCount == 1 ? "" : "s",
             (trees[0].depth * treeCount == markers) ? " of" : " of up to",
             trees[0].depth, trees[0].depth == 1 ? "" : "s");

//      for (int i = 0; i < treeCount; i++)
//         trees[i].freqs.Set(1.0 / trees[i].freqs.dim);

      for (int i = 0; i < treeCount; i++)
         {
         trees[i].freqs.Set(1.0 / trees[i].freqs.Length());
         // for (int j = 0; j < trees[i].freqs.Length(); j++)
         //   trees[i].freqs[j] = globalRandom.Next();
         // trees[i].freqs.Multiply(1.0 / trees[i].freqs.Sum());
         }

      llk0 = 0.0;
      int pass = 0;

      while (true)
         {
         unrel_llk1 = families_llk1 = 0.0;

         for (int i = 0; i < treeCount; i++)
            trees[i].new_freqs.Zero();

         for (int i = 0; i < sets->Length(); i++)
            {
            HaplotypeSet * set = (HaplotypeSet *) sets->Object(i);

            secondPass = false;
            setLikelihood = SetLikelihood(set) * count / set->copies;

            if (setLikelihood == 0.0) continue;

            secondPass = true;
            SetLikelihood(set);

            // printf("%d - %g [%d]\n", i, setLikelihood, set->scale);

            double llk = log(setLikelihood * set->copies / count) * set->copies;

            // printf("XXX\n");

            if (set->scale)
               llk += log(SCALE_DOWN) * set->scale;

            // printf("XXX\n");

            if (set->haploCount <= 2)
               unrel_llk1 += llk;
            else
               families_llk1 += llk;
            }

         llk1 = unrel_llk1 + families_llk1;

         // if (treeBits == 0)
         //   {
         //   unrel_llk0 = unrel_llk1;
         //   families_llk0 = families_llk1;
         //   }

         for (int i = 0; i < treeCount; i++)
            trees[i].freqs = trees[i].new_freqs;

         if (!quiet)
            {
            printf("  Pass %2d, log(lk) = %.3f    \r", ++pass, llk1);
            fflush(stdout);
            }

         double precision = max(fabs(llk1 + llk0) * tol, ztol);

         if (fabs(llk1 - llk0) < precision)
            break;

         llk0 = llk1;
         }

      // Turn second pass flag off for safety
      secondPass = false;

      if (!quiet) printf("\n");
      // for (int i = 0; i < treeCount; i++)
      //   trees[i].Print(0.0001);

      if (treeCount == 1)
         break;

      int newTreeCount = (treeCount / 2) + (treeCount & 1);

      HaploTree * newTrees = new HaploTree [newTreeCount];

      for (int i = 0; i < newTreeCount; i++)
         if ((i * 2 + 1) < treeCount)
            newTrees[i].Merge(trees[i * 2], trees[i * 2 + 1]);
         else
            newTrees[i].Copy(trees[i * 2]);

      delete [] trees;

      treeCount = newTreeCount;
      trees = newTrees;
      treeBits++;
      }

//   printf("Changes in LLK: Families - %6.3f, Unrelateds - %6.3f\n",
//           families_llk1 - families_llk0, unrel_llk1 - unrel_llk0);

   if (best_llk == 0.0 || llk1 > best_llk)
      best_llk = llk1;
   }

double SparseLikelihood::SetLikelihood(HaplotypeSet * set)
   {
   // printf("Entering SetLikelihood ...\n");

   double lk = 0.0;

   factorMatrix.Dimension(set->graphCount, treeCount);

   set->scale = 0;
   for (int i = 0; i < set->graphCount; i++)
      {
      HaplotypeGraph * graph = set->graphs[i];

      // Skip impossible graphs
      if (secondPass && graph->likelihood == 0.0) continue;

      factors = &factorMatrix[i];

      GraphLikelihood(graph);

      // Do some gymnastics to avoid underflows ...
      if (i == 0 || set->scale - 1 > graph->scale)
         {
         lk = graph->likelihood * graph->weight;
         set->scale = graph->scale;
         }
      else if (set->scale == graph->scale)
         lk += graph->likelihood * graph->weight;
      else if (set->scale + 1 == graph->scale)
         lk += graph->likelihood * SCALE_DOWN * graph->weight;
      else if (set->scale - 1 == graph->scale)
         {
         lk = lk * SCALE_DOWN + graph->likelihood * graph->weight;
         set->scale --;
         }
      }

   // Do some gymnastics to ensure all graphs have
   // likelihoods measured on the same scale ...
   for (int i = 0; i < set->graphCount; i++)
      {
      HaplotypeGraph * graph = set->graphs[i];

      if (graph->scale > set->scale)
         if ( graph->scale > set->scale + 1)
            graph->likelihood = 0.0;
         else
            {
            graph->likelihood *= SCALE_DOWN;
            graph->scale --;
            }
      }

   // printf("Leaving SetLikelihood ...\n");
   return lk;
   }

double SparseLikelihood::GraphLikelihood(HaplotypeGraph * graph)
   {
   // printf("  Entering GraphLikelihood ...\n");
   if (!secondPass) factors->Zero();

   for (int tree = 0; tree < treeCount; tree++)
      {
      last.Dimension(graph->count);
      last.Set(0);

      marginals.Dimension(graph->count);

      for (int i = 0; i < graph->count; i++)
         {
         index[i].Dimension(1);
         index[i][0] = 0;

         if (!UpdateIndex(graph, tree, i))
            {
            graph->likelihood = 0.0;
            graph->scale = 0;
            return 0.0;
            }
         }

      double factor = GraphLikelihood(graph, tree, 0);

      if (!secondPass) (*factors)[tree] = factor;
      }

   graph->likelihood = 1.0;
   graph->scale = 0;

   for (int tree = 0; tree < treeCount; tree++)
      {
      graph->likelihood *= (*factors)[tree];

      if (graph->likelihood < SCALE_DOWN)
         {
         graph->likelihood *= SCALE_UP;
         graph->scale++;
         }
      }

   // printf("  Leaving GraphLikelihood ...\n");
   return graph->likelihood * graph->weight;
   }

double SparseLikelihood::GraphLikelihood(HaplotypeGraph * graph, int tree, int i)
   {
   Unknown * u;

   while (true)
      {
      if (i == graph->unknowns.Length())
         return CurrentLikelihood(graph, tree);

      u = graph->unknowns.GetUnknown(i);

      if (u->marker >> treeBits == tree)
         break;

      i++;
      }

   if (u->type == U_MISSING)
      {
      MissingAllele * m = (MissingAllele *) u;

      if (graph->lastGenotype[m->haplotype] != -1)
         if (!ExpandIndex(graph, tree, m->haplotype))
            return 0.0;

      return GraphLikelihood(graph, tree, i + 1);
      }

   /*** u->type == U_GRAPH ***/
   AlleleGraph * g = (AlleleGraph *) u;

   IntArray pos = last;
   IntArray * save = new IntArray [graph->count];

   for (int j = 0; j < graph->count; j++)
      save[j] = index[j];

   bool possible = true;
   for (int j = 0; j < g->haplotypes.Length(); j++)
      {
      possible &= SelectAllele(tree, g->haplotypes[j], g->alleles[j * 2]);
      possible &= UpdateIndex(graph, tree, g->haplotypes[j]);
      if (!possible) break;
      }
   double lk = possible ? GraphLikelihood(graph, tree, i + 1) : 0.0;

   for (int j = 0; j < graph->count; j++)
      index[j] = save[j];
   delete [] save;
   last = pos;

   possible = true;
   for (int j = 0; j < g->haplotypes.Length(); j++)
      {
      possible &= SelectAllele(tree, g->haplotypes[j], g->alleles[j * 2 + 1]);
      possible &= UpdateIndex(graph, tree, g->haplotypes[j]);
      if (!possible) break;
      }
   return possible ? lk + GraphLikelihood(graph, tree, i + 1) : lk;
   }


double SparseLikelihood::CurrentLikelihood(HaplotypeGraph * graph, int tree)
   {
   for (int i = 0; i < graph->count; i++)
      if (last[i] == 0)
         marginals[i] = 1.0;
      else if (last[i] == trees[tree].depth)
         {
         marginals[i] = trees[tree].freqs[index[i][0]];

         for (int j = 1; j < index[i].Length(); j++)
            marginals[i] += trees[tree].freqs[index[i][j]];
         }
      else
         error("Unexpected state while processing haplotype graph");

   double likelihood = 1.0;

   for (int i = 0; i < graph->count; i++)
      likelihood *= marginals[i];

   if (secondPass)
      {
      int   scale = graph->scale;
      likelihood *= graph->weight;

      for (int i = 0; i < treeCount; i++)
         if (i != tree)
            {
            likelihood *= (*factors)[i];
            if (likelihood < SCALE_DOWN && scale)
               scale--, likelihood *= SCALE_UP;
            }

      while (scale--)
         likelihood *= SCALE_UP;

      if ((likelihood / setLikelihood) < 1e-7)
         return 0.0;

      for (int i = 0; i < graph->count; i++)
         {
         double scale = likelihood / (setLikelihood * marginals[i]);

         if (last[i] == 0)
            for (int j = 0; j < trees[tree].new_freqs.Length(); j++)
               trees[tree].new_freqs[j] += trees[tree].freqs[j] * scale;
         else
            for (int j = 0; j < index[i].Length(); j++)
               trees[tree].new_freqs[index[i][j]] += trees[tree].freqs[index[i][j]] * scale;
         }

      return likelihood;
      }

   return likelihood;
   }

bool SparseLikelihood::UpdateIndex(HaplotypeGraph * graph, int tree, int i)
   {
   int marker = last[i] + (tree << treeBits);

   while (last[i] < trees[tree].depth && graph->haplotypes[i][marker] != 0)
      if (!SelectAllele(tree, i, graph->haplotypes[i][marker++]))
         return false;

   return index[i].Length() != 0;
   }

bool SparseLikelihood::ExpandIndex(HaplotypeGraph * graph, int tree, int i)
   {
   int count = 0;

   for (int j = 0; j < index[i].Length(); j++)
      count += trees[tree].entries[index[i][j]];

   new_index.Dimension(count);

   count = 0;

   for (int j = 0; j < index[i].Length(); j++)
      for (int k = 0; k < trees[tree].branches[index[i][j]]->Length(); k++)
         if ((*trees[tree].branches[index[i][j]])[k] != -1)
            new_index[count++] = (*trees[tree].branches[index[i][j]])[k];

   index[i] = new_index;
   last[i]++;

   return UpdateIndex(graph, tree, i);
   }

bool SparseLikelihood::SelectAllele(int tree, int i, int allele)
   {
   int out = 0;

   for (int j = 0; j < index[i].Length(); j++)
      {
      int branch = trees[tree].PeekBranch(index[i][j], allele - 1);

      if (branch != -1) index[i][out++] = branch;
      }

   index[i].Dimension(out);
   last[i]++;

   return index[i].Length() != 0;
   }

void SparseLikelihood::PrintExtras(double minfreq)
   {
   // Variables used for traversing haplotype tree
   IntArray pointer, state;
   int      leaf;

   // Calculate information for this distribution
   double info = 0;

   trees->SetupTraversal(pointer, state);

   while ((leaf = trees->Traverse(pointer, state, minfreq)) != -1)
      info += trees->freqs[leaf] * log(trees->freqs[leaf]);

   // Calculate D' between markers 1..n-1 and marker n
   StringMap haplos, alleles;
   String    haplo, al;

   trees->SetupTraversal(pointer, state);

   while ((leaf = trees->Traverse(pointer, state, minfreq)) != -1)
      {
      haplo.Clear();

      for (int i = 0; i < trees->depth - 1; i++)
         haplo += (char) ('0' + state[i]);

      haplos.Add(haplo);
      alleles.Add(al = state[trees->depth - 1]);
      }

   Vector freqs1(haplos.Length()); freqs1.Zero();
   Vector freqs2(alleles.Length()); freqs2.Zero();

   trees->SetupTraversal(pointer, state);

   while ((leaf = trees->Traverse(pointer, state, minfreq)) != -1)
      {
      haplo.Clear();

      for (int i = 0; i < trees->depth - 1; i++)
         haplo += (char) ('0' + state[i]);

      freqs1[haplos.Find(haplo)] += trees->freqs[leaf];
      freqs2[alleles.Find(al = state[trees->depth - 1])] += trees->freqs[leaf];
      }

   double D = 0;
   double missing_haplotypes = 1.0;

   trees->SetupTraversal(pointer, state);

   while ((leaf = trees->Traverse(pointer, state, minfreq)) != -1)
      {
      haplo.Clear();

      for (int i = 0; i < trees->depth - 1; i++)
         haplo += (char) ('0' + state[i]);

      double f1 = freqs1[haplos.Find(haplo)];
      double f2 = freqs2[alleles.Find(al = state[trees->depth - 1])];

      double _D = trees->freqs[leaf] - f1 * f2;

      if (fabs(_D) < 1e-7) continue;

      double _Dmax = _D < 0 ?
         min(f1 * f2, (1 - f1) * (1 - f2)) :
         min((1 - f1) * f2, f1 * (1 - f2));

      D += fabs(_D/_Dmax) * f1 * f2;
      missing_haplotypes -= f1 * f2;
      }

   // Need to account for expected, but zero frequency haplotypes
   D += missing_haplotypes;

   // Calculate top 5 allele frequencies
   Vector top5;

   top5.Dimension(5);
   top5.Zero();

   trees->SetupTraversal(pointer, state);

   while ((leaf = trees->Traverse(pointer, state, minfreq)) != -1)
      {
      double freq = trees->freqs[leaf];

      int pos = 4;

      while (pos >= 0 && top5[pos] < freq)
         {
         if (pos > 0)
            top5[pos] = top5[pos - 1];
         pos--;
         }

      if (pos < 4)
         top5[pos + 1] = freq;
      }

   for (int i = 1; i < 5; i++)
      top5[i] += top5[i - 1];

   // Print out the statistics
   printf("\n\nAdditional Information\n");
   printf("%10s %10s %10s %10s %10s %10s %10s\n",
          "Info", "D(last,*)", "Top1", "Top2", "Top3", "Top4", "Top5");
   printf("%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
          info, D, top5[0], top5[1], top5[2], top5[3], top5[4]);
   }

void SparseLikelihood::SampleState(HaplotypeGraph * graph)
   {
   // First we calculate the overall probability of current
   // haplotype set ...
   for (int i = 0; i < graph->count; i++)
      if (last[i] == 0)
         marginals[i] = 1.0;
      else if (last[i] == trees->depth)
         {
         marginals[i] = trees->freqs[index[i][0]];

         for (int j = 1; j < index[i].Length(); j++)
            if ( trees->freqs[index[i][j]] > marginals[i] )
               marginals[i] += trees->freqs[index[i][j]];
         }


   double likelihood = marginals.Product();

   // Then we check whether we should sample from the current set
   // or retain the previous configuration
   setLikelihood += likelihood;
   if (globalRandom.Next() > likelihood / setLikelihood)
      return;

   // Finally, we sample a haplotype from the current configuration
   for (int i = 0; i < graph->count; i++)
      if (last[i] == 0)
         selectedState[i] = -1;
      else if (last[i] == trees->depth)
          {
          double alternative = trees->freqs[index[i][0]];
          selectedState[i]   = index[i][0];

          for (int j = 1; j < index[i].Length(); j++)
             {
             double current = trees->freqs[index[i][j]];
             alternative += current;

             if (globalRandom.Next() <  current / alternative )
                selectedState[i] = index[i][j];
             }
          }
   }

void SparseLikelihood::SelectState(HaplotypeGraph * graph)
   {
   // Retrieve haplotype frequencies for the most likely configuration
   for (int i = 0; i < graph->count; i++)
      if (last[i] == 0)
         marginals[i] = 1.0;
      else if (last[i] == trees->depth)
         {
         marginals[i] = trees->freqs[index[i][0]];

         for (int j = 1; j < index[i].Length(); j++)
            if ( trees->freqs[index[i][j]] > marginals[i] )
               marginals[i] = trees->freqs[index[i][j]];
         }

   // Multiply them together
   double likelihood = marginals.Product();

   // If the new state is better than the previous one ...
   if (likelihood < setLikelihood)
      return;
   setLikelihood = likelihood;

   // Then update the current set of haplotypes
   for (int i = 0; i < graph->count; i++)
      if (last[i] == 0)
         selectedState[i] = -1;
      else if (last[i] == trees->depth)
         {
         double current = trees->freqs[index[i][0]];
         selectedState[i] = index[i][0];

         for (int j = 1; j < index[i].Length(); j++)
            if ( trees->freqs[index[i][j]] > current )
               {
               current = trees->freqs[index[i][j]];
               selectedState[i] = index[i][j];
               }
         }
   }

void SparseLikelihood::ExtractConfiguration(HaplotypeGraph * graph, bool sample, int i)
   {
   if (i == graph->unknowns.Length())
      {
      if (sample)
         SampleState(graph);
      else
         SelectState(graph);
      return;
      }

   Unknown * u = graph->unknowns.GetUnknown(i++);

   if (u->type == U_MISSING)
      {
      MissingAllele * m = (MissingAllele *) u;

      if (graph->lastGenotype[m->haplotype] != -1)
         if (!ExpandIndex(graph, 0, m->haplotype))
            return;

      ExtractConfiguration(graph, sample, i);
      return;
      }

   /*** u->type == U_GRAPH ***/
   AlleleGraph * g = (AlleleGraph *) u;

   IntArray pos = last;
   IntArray * save = new IntArray [graph->count];

   for (int j = 0; j < graph->count; j++)
      save[j] = index[j];

   bool possible = true;
   for (int j = 0; j < g->haplotypes.Length(); j++)
      {
      possible &= SelectAllele(0, g->haplotypes[j], g->alleles[j * 2]);
      possible &= UpdateIndex(graph, 0, g->haplotypes[j]);
      if (!possible) break;
      }
   if (possible) ExtractConfiguration(graph, sample, i);

   for (int j = 0; j < graph->count; j++)
      index[j] = save[j];
   delete [] save;
   last = pos;

   possible = true;
   for (int j = 0; j < g->haplotypes.Length(); j++)
      {
      possible &= SelectAllele(0, g->haplotypes[j], g->alleles[j * 2 + 1]);
      possible &= UpdateIndex(graph, 0, g->haplotypes[j]);
      if (!possible) break;
      }
   if (possible) ExtractConfiguration(graph, sample, i);
   }

void SparseLikelihood::ExtractConfiguration
    (HaplotypeGraph * graph, IntArray & configuration, bool sample)
   {
   setLikelihood = 0;

   marginals.Dimension(graph->count);

   selectedState.Dimension(graph->count);
   selectedState.Set(-1);

   last.Dimension(graph->count);
   last.Set(0);

   bool valid = true;
   for (int i = 0; i < graph->count; i++)
      {
      index[i].Dimension(1);
      index[i][0] = 0;

      if (!UpdateIndex(graph, 0, i))
         {
         valid = false;
         break;
         }
      }

   if (valid) ExtractConfiguration(graph, sample);

   configuration = selectedState;
   }

void SparseLikelihood::RetrieveHaplotypes(IntArray & haplos, Vector & freqs, double minFreq)
   {
   haplos.Clear();
   freqs.Clear();

   if (trees->depth != markers)
      error("Attempting to retrieve haplotype frequency information, but\n"
            "none available\n");

   IntArray pointer, state;
   trees->SetupTraversal(pointer, state);

   int    leaf;
   while ((leaf = trees->Traverse(pointer, state, minFreq)) != -1)
      {
      freqs.Push(trees->freqs[leaf]);

      for (int i = 0; i < markers; i++)
         haplos.Push(state[i] + 1);
      }
   }

double SparseLikelihood::GetLogOffset(HaplotypeSet * set)
   {
   return set->scale * log(SCALE_DOWN);
   }


 
