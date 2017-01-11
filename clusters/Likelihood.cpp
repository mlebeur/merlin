////////////////////////////////////////////////////////////////////// 
// clusters/Likelihood.cpp 
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
 
#include "Likelihood.h"
#include "MathConstant.h"
#include "Random.h"
#include "Error.h"

#include <math.h>

Likelihood::Likelihood()
   {
   markers = 0;
   index = NULL;
   tol = 1e-7;
   ztol = 1e-10;
   pseudoMutation = 0.0;
   gradual_smoothing = two_pass_smoothing = false;
   failText = "no failure";
   }

void Likelihood::Initialize(const IntArray & counts)
   {
   markers = counts.Length();
   alleleCounts = counts;

   double number_of_haplotypes = alleleCounts.DoubleProduct();

   if (number_of_haplotypes > (1 << 26))
      error("The number of possible haplotypes (~2^%.1f) is too large",
             log(number_of_haplotypes) / log(2.0));

   haplos = (int) number_of_haplotypes;

//   double memory = 3. * sizeof(double) * number_of_haplotypes * pow(2,-20);
//
//   printf("%d haplotype frequencies will be estimated\n"
//          "  [~%.0f Mb of memory required]\n",
//          haplos, memory);

   frequencies.Dimension(haplos);
   new_frequencies.Dimension(haplos);
   best_frequencies.Clear();

   if (index != NULL) delete [] index;
   index = new IntArray[HaplotypeSet::maxHaploCount];

   best_llk = 0;
   }

void Likelihood::Initialize(const IntArray & counts, const IntArray & haplotypes, const Vector & freqs)
   {
   markers = counts.Length();
   alleleCounts = counts;

   double number_of_haplotypes = 1;

   for (int i = 0; i < markers; i++)
      number_of_haplotypes *= alleleCounts[i];

   if (number_of_haplotypes > (1 << 26))
      error("The number of possible haplotypes (~2^%.1f) is too large",
             log(number_of_haplotypes) / log(2.0));

   haplos = (int) number_of_haplotypes;

   frequencies.Dimension(haplos);
   frequencies.Zero();

   for (int base = 0, f = 0; base < haplotypes.Length(); base += markers)
      {
      int haplotype = haplotypes[base];

      for (int j = 1; j < markers; j++)
         {
         haplotype *= alleleCounts[j];
         haplotype += haplotypes[base + j] - 1;
         }

      frequencies[haplotype] = freqs[f++];
      }
   
   if (index != NULL) delete [] index;
   index = new IntArray[HaplotypeSet::maxHaploCount];

   best_llk = 0;
   }

void Likelihood::SaveSummary()
   {
   if (best_frequencies.Length() == 0) return;

   FILE * f = fopen("first.freq", "at");
   fprintf(f, "%.5f\n", best_frequencies[0]);
   fclose(f);
   }

void Likelihood::RetrieveHaplotypes(IntArray & haplos, Vector & freqs, double minFreq)
   {
   haplos.Clear();
   freqs.Clear();

   for (int i = 0; i < best_frequencies.Length(); i++)
      if (best_frequencies[i] > minFreq)
         {
         int base  = haplos.Length();
         int haplo = i;

         freqs.Push(best_frequencies[i]);
         haplos.Dimension(base + markers);

         for (int j = markers - 1; j >= 0; j--)
            {
            haplos[base + j] = haplo % alleleCounts[j] + 1;
            haplo /= alleleCounts[j];
            }
         }
   }

void Likelihood::Print(double minfreq, HaplotypeSets * sets)
   {
   String label;

   printf("Haplotypes with estimated frequency > %g\n", minfreq);

   int count = 0;
   double total = 0.0;

   for (int i = 0; i < best_frequencies.dim; i++)
      if (best_frequencies[i] > minfreq)
         {
         label = "";

         int haplo = i;
         for (int j = markers - 1; j >= 0; j--)
            {
            int allele = haplo % alleleCounts[j] + 1;
            haplo /= alleleCounts[j];

            if (sets == NULL)
               label = allele + label;
            else
               label = sets->GetAlleleLabel(j, allele) + label;
            }

         printf("%6.2f%% %s\n", best_frequencies[i] * 100.0, (const char *) label);

         total += best_frequencies[i];
         count++;
         }

   printf("\nThese %d haplotypes represent %.2f%% of total probability\n",
          count, total * 100);
   printf("The logLikelihood of the data is %.4f\n\n", best_llk);
   }

void Likelihood::RandomEM(HaplotypeSets * sets)
   {
//   printf("Starting with random allele frequencies...\n");
   for (int i = 0; i < frequencies.dim; i++)
      frequencies[i] = globalRandom.Next();

   frequencies /= frequencies.Sum();

   CoreEM(sets);
//   printf("\n\n");
   }

void Likelihood::CheckSmoothing()
   {
   if (two_pass_smoothing || gradual_smoothing)
      for (int i = 0; i < markers; i++)
         if (alleleCounts[i] != 2)
            error("PseudoMutation model only valid for diallelic markers\n");
   }

void Likelihood::EM(HaplotypeSets * sets)
   {
//   printf("Starting with equal allele frequencies...\n");
   frequencies.Set(1.0 / frequencies.dim);

   CoreEM(sets);

   if (two_pass_smoothing)
      {
      pseudoMutation = 0.05;
      PseudoCoalescent(0, frequencies.Length() - 1);
      CoreEM(sets);
      }

//   printf("\n\n");
   }

void Likelihood::EditedEM(HaplotypeSets * sets)
   {
   int removed = 0;

//   printf("Removing very rare haplotypes...\n");
   for (int i = 0; i < frequencies.dim; i++)
      if (frequencies[i] < 1e-4)
         {
         frequencies[i] = 0;
         removed++;
         }
//   printf("Removed %d haplotypes\n\n", removed);

   frequencies /= frequencies.Sum();

   CoreEM(sets);
//   printf("\n\n");
   }

void Likelihood::CoreEM(HaplotypeSets * sets)
   {
   int count = 0;
   for (int i = 0; i < sets->Length(); i++)
      {
      HaplotypeSet * set = (HaplotypeSet *) sets->Object(i);
      count += set->haploCount * set->copies;
      set->StandardizeWeights();
      }

   if (count == 0) return;

   double llk0 = 0.0, llk1;
   int    pass = 0;

   while (true)
      {
      llk1 = 0.0;

      new_frequencies.Zero();
      missingHaplotypes = 0.0;

      for (int i = 0; i < sets->Length(); i++)
         {
         HaplotypeSet * set = (HaplotypeSet *) sets->Object(i);

         secondPass = false;
         setLikelihood = SetLikelihood(set) * count / set->copies;

         if (setLikelihood == 0.0)
            // We shouldn't get here unless important haplotypes
            // get triaged ... if that happens, this line increments
            // all haplotype frequencies a little so we can try to
            // recover and avoid a crash
            // TODO -- There must be a better solution!
            {
            new_frequencies += 1.0 / new_frequencies.Length() * count;
            continue;
            }

         secondPass = true;
         SetLikelihood(set);

         llk1 += log(setLikelihood * set->copies / count) * set->copies;
         }

      // We have gone through one more round of maximization
      pass++;

      // Update frequencies
      if (missingHaplotypes < 0.99999999)
         {
         frequencies = new_frequencies;
         frequencies.Multiply(1.0 / (1.0 - missingHaplotypes));
         }

//      printf("Pass %d -- lnLikelihood %.3f, Freqs = %.3f, Missing Data = %.3f\n",
//             pass, llk1, frequencies.Sum(), missingHaplotypes);

      if (pass < 8 && gradual_smoothing)
         {
//         best_frequencies = frequencies;
//         Print(0.01);

         pseudoMutation = pow(0.5, pass + 2);
         PseudoCoalescent(0, frequencies.Length() - 1);

//         best_frequencies = frequencies;
//         Print(0.01);
         }

//      if (pass < 3 && pseudoMutation)
//         PseudoCoalescent(0, frequencies.Length() - 1);

//      printf("Pass %2d, log(lk) = %.3f    \r", pass, llk1);
//      fflush(stdout);

      double precision = max(fabs(llk1 + llk0) * tol, ztol);

      // Added a check for futility -- so that we stop after
      // 100 many rounds are carried out ...
      if (fabs(llk1 - llk0) < precision || pass > 100)
         break;

      llk0 = llk1;
      }

   if (best_llk == 0.0 || llk1 > best_llk)
      {
      best_llk = llk1;
      best_frequencies = frequencies;
      }
   }

double Likelihood::SetLikelihood(HaplotypeSet * set)
   {
   double lk = 0.0;

   for (int i = 0; i < set->graphCount; i++)
      lk += GraphLikelihood(set->graphs[i]);

   return lk;
   }

double Likelihood::GraphLikelihood(HaplotypeGraph * graph)
   {
   last.Dimension(graph->count);
   last.Set(0);

   marginals.Dimension(graph->count);

   for (int i = 0; i < graph->count; i++)
      {
      index[i].Dimension(1);
      index[i][0] = 0;

      UpdateIndex(graph, i);
      }

   double likelihood = GraphLikelihood(graph, 0);
   graph->likelihood = likelihood / graph->weight ;
   return likelihood;
   }

double Likelihood::GraphLikelihood(HaplotypeGraph * graph, int i)
   {
   if (i == graph->unknowns.Length())
      return CurrentLikelihood(graph);

   Unknown * u = graph->unknowns.GetUnknown(i);

   if (u->type == U_MISSING)
      {
      MissingAllele * m = (MissingAllele *) u;

      if (graph->lastGenotype[m->haplotype] != -1)
         ExpandIndex(graph, m->haplotype);

      return GraphLikelihood(graph, i + 1);
      }

   /*** u->type == U_GRAPH ***/
   AlleleGraph * g = (AlleleGraph *) u;

   IntArray len(graph->count);
   IntArray pos = last;

   for (int j = 0; j < graph->count; j++)
      len[j] = index[j].Length();

   for (int j = 0; j < g->haplotypes.Length(); j++)
      {
      index[g->haplotypes[j]] *= alleleCounts[g->marker];
      index[g->haplotypes[j]] += g->alleles[j * 2] - 1;
      last[g->haplotypes[j]]++;
      UpdateIndex(graph, g->haplotypes[j]);
      }
   double lk = GraphLikelihood(graph, i + 1);

   for (int j = 0; j < graph->count; j++)
      {
      index[j].Dimension(len[j]);
      ReverseIndex(graph, j, pos[j]);
      }

   for (int j = 0; j < g->haplotypes.Length(); j++)
      {
      index[g->haplotypes[j]] *= alleleCounts[g->marker];
      index[g->haplotypes[j]] += g->alleles[j * 2 + 1] - 1;
      last[g->haplotypes[j]]++;
      UpdateIndex(graph, g->haplotypes[j]);
      }

   return lk + GraphLikelihood(graph, i + 1);
   }

double Likelihood::CurrentLikelihood(HaplotypeGraph * graph)
   {
   for (int i = 0; i < graph->count; i++)
      if (last[i] == 0)
         marginals[i] = 1.0;
      else if (last[i] == markers)
         {
         marginals[i] = frequencies[index[i][0]];

         for (int j = 1; j < index[i].Length(); j++)
            marginals[i] += frequencies[index[i][j]];
         }
      else
         error("Unexpected state while processing haplotype graph");

   double likelihood = graph->weight;

   for (int i = 0; i < graph->count; i++)
      likelihood *= marginals[i];

   if (secondPass && ((likelihood / setLikelihood) > 1e-9))
      {
      for (int i = 0; i < graph->count; i++)
         {
         double scale = likelihood / (setLikelihood * marginals[i]);

         if (last[i] == 0)
            missingHaplotypes += scale;
            // for (int j = 0; j < new_frequencies.dim; j++)
            //   new_frequencies[j] += frequencies[j] * scale;
         else
            for (int j = 0; j < index[i].Length(); j++)
               new_frequencies[index[i][j]] += frequencies[index[i][j]] * scale;
         }

      return likelihood;
      }
   else if (secondPass)
      return 0.0;

   return likelihood;
   }

void Likelihood::UpdateIndex(HaplotypeGraph * graph, int i)
   {
   while (last[i] < markers && graph->haplotypes[i][last[i]] != 0)
      {
      index[i] *= alleleCounts[last[i]];
      index[i] += graph->haplotypes[i][last[i]] - 1;
      last[i]++;
      }
   }

void Likelihood::ReverseIndex(HaplotypeGraph * graph, int i, int pos)
   {
   if (pos == last[i]) return;

   int factor = alleleCounts[--last[i]];

   while (last[i] > pos)
      factor *= alleleCounts[--last[i]];

   index[i] /= factor;
   }

void Likelihood::ExpandIndex(HaplotypeGraph * graph, int i)
   {
   int factor = alleleCounts[last[i]];
   int original = index[i].Length();

   index[i] *= factor;
   index[i].Dimension(index[i].Length() * factor);

   for (int j = 1; j < factor; j++)
      for (int k = original * j; k < original * (j + 1); k++)
         index[i][k] = index[i][k - original] + 1;

   last[i]++;

   UpdateIndex(graph, i);
   }

void Likelihood::PrintExtras()
   {
   // Calculate information for this distribution
   double info = 0;

   for (int i = 0; i < best_frequencies.dim; i++)
      if (frequencies[i] > 1e-7)
         info += best_frequencies[i] * log(best_frequencies[i]);

   // Calculate D' between last marker and haplotype of other markers
   double D = 0.0;

   marginals.Dimension( alleleCounts[markers - 1] );
   marginals.Zero();

   int alleles = alleleCounts[markers - 1];

   for (int i = 0; i < best_frequencies.dim; i++)
      if (frequencies[i] > 1e-7)
         marginals[i % alleles] += frequencies[i];

   for (int i = 0; i < best_frequencies.dim; i += alleles)
      {
      double freq = 0.0;

      for (int j = 0; j < alleles; j++)
         freq += best_frequencies[i + j];

      if (freq < 1e-7) continue;

      for (int j = 0; j < alleles; j++)
         {
         double _D = best_frequencies[i + j] - freq * marginals[j];

         if (fabs(_D) < 1e-7) continue;

         double _Dmax = _D < 0 ?
                min(freq * marginals[j], (1. - freq) * (1. - marginals[j])) :
                min((1. - freq) * marginals[j], freq * (1. - marginals[j]));

         D += fabs(_D/_Dmax) * freq * marginals[j];
         }
      }

   // Calculate cumulative frequencies
   frequencies = best_frequencies;
   frequencies.Sort();

   double top1 = frequencies.dim > 0 ? frequencies[frequencies.dim - 1] : 1.0;
   double top2 = frequencies.dim > 1 ? frequencies[frequencies.dim - 2] + top1 : 1.0;
   double top3 = frequencies.dim > 2 ? frequencies[frequencies.dim - 3] + top2 : 1.0;
   double top4 = frequencies.dim > 3 ? frequencies[frequencies.dim - 4] + top3 : 1.0;
   double top5 = frequencies.dim > 4 ? frequencies[frequencies.dim - 5] + top4 : 1.0;

   printf("\n\nAdditional Information\n");
   printf("%10s %10s %10s %10s %10s %10s %10s\n",
          "Info", "D(last,*)", "Top1", "Top2", "Top3", "Top4", "Top5");
   printf("%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
          info, D, top1, top2, top3, top4, top5);
   }

bool Likelihood::PairwiseLD(double & dprime, double & deltasq, int & coupling)
   {
   coupling = PAIRWISE_EQUILIBRIUM;

   if (markers != 2 || alleleCounts[0] != 2 || alleleCounts[1] != 2)
      {
      // Markers are not biallelic
      dprime = deltasq = 0.0;
      failText = markers == 2 ? "markers not biallelic" : "not a two marker haplotype";
      return false;
      }

   if (best_frequencies.Length() == 0)
      {
      // Haplotype frequency estimation failed
      dprime = deltasq = 0.0;
      failText = "no suitable data [check inheritance and family size?]";
      return false;
      }

  double p = best_frequencies[0] + best_frequencies[2];
  double r = best_frequencies[0] + best_frequencies[1];

  if (min(p, 1 - p) < 1e-8 || min(r, 1 - r) < 1e-8)
     {
     // One of the markers is actually monomorphic
     failText = "one of the markers is monomorphic";
     dprime = deltasq = 0.0;
     return false;
     }

  double D = best_frequencies[0] - p * r;
  dprime = 0;

  if (D > 1e-8)
     {
     dprime = D / min((1-p)*r, (1-r)*p);
     coupling = PAIRWISE_COUPLING;
     }
  else if (D < -1e-8)
     {
     dprime = D / min((1-p)*(1-r), p*r);
     coupling = PAIRWISE_REPULSION;
     }
  else
     coupling = PAIRWISE_EQUILIBRIUM;

  dprime = fabs(dprime);

  deltasq = (best_frequencies[0] * best_frequencies[3] -
             best_frequencies[1] * best_frequencies[2]);
  deltasq *= deltasq;
  deltasq /= p * (1-p) * r * (1-r);

  return true;
  }

void Likelihood::PseudoCoalescent(int start, int stop)
   {
   if (start == stop) return;

   int mid = (start + stop) / 2;

   PseudoCoalescent(start, mid);
   PseudoCoalescent(mid + 1, stop);

   int delta = mid - start;
   for (int i = start; i < mid; i++)
      {
      double a = frequencies[i];
      double b = frequencies[i + delta];

      frequencies[i] = (1.0 -  pseudoMutation) * a + pseudoMutation * b;
      frequencies[i + delta] = (1.0 - pseudoMutation) * b + pseudoMutation * a;
      }
   }

 
