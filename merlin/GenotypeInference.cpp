////////////////////////////////////////////////////////////////////// 
// merlin/GenotypeInference.cpp 
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
 
#include "GenotypeInference.h"
#include "MerlinCluster.h"
#include "MerlinCore.h"
#include "Houdini.h"
#include "Error.h"

#include <math.h>

bool GenotypeInference::inferBest = false;
bool GenotypeInference::inferExpected = false;
bool GenotypeInference::inferProbabilities = false;

bool GenotypeInference::offline = false;

void GenotypeInference::AllocateMemory(Pedigree & ped)
   {
   probability[0].Dimension(ped.count * markers);
   probability[1].Dimension(ped.count * markers);

   // Setting probabilities to -1.0 let us check for individuals whose
   // genotypes were not adequately infered, for whatever reason
   probability[0].Set(-1.0);

#ifdef __CHROMOSOME_X__
   isFemale.Dimension(ped.count);

   for (int i = 0; i < ped.count; i++)
      isFemale = ped[i].sex == 2;
#endif
   }

double GenotypeInference::GetExpectedGenotype(Pedigree & ped, int individual, int marker)
   {
   if (offline)
      return ped[individual].traits[markerKey[marker]] != _NAN_ ?
                ped[individual].traits[markerKey[marker]] :
                2.0 * frequencies[marker];
   else
      return GetExpectedGenotype(individual, marker);
   }

double GenotypeInference::GetExpectedGenotype(int individual, int marker)
   {
   int slot = GetSlot(individual, marker);

   if (probability[0][slot] == -1.0)
      return 2.0 * frequencies[markerKey[marker]];

   return probability[0][slot] * 2 + (1.0 - probability[0][slot] - probability[1][slot]);
   }

double GenotypeInference::GetProbability(int individual, int marker, int genotype)
   {
   int slot = GetSlot(individual, marker);

   if (probability[0][slot] == -1.0)
      {
      double freq = frequencies[markerKey[marker]];

#ifdef __CHROMOSOME_X__
      /* This section applies to X chromosome genotypes for males */
      if (!isFemale[individual])
         switch (genotype)
            {
            case 0  :
               return freq;
            case 1  :
               return 0;
            default : /* case 2 */
               return (1.0 - freq);
            }
#endif

      /* This sections applies to autosomes and the Y chromosome females */
      switch (genotype)
         {
         case 0  :
            return freq * freq;
         case 1  :
            return (1.0 - freq) * freq * 2.0;
         default : /* case 2 */
            return (1.0 - freq) * (1.0 - freq);
         }
      }

   switch (genotype)
      {
      case 0:
         return probability[0][slot];
      case 1:
         return 1.0 - probability[0][slot] - probability[1][slot];
      default: /* case 2 */
         return probability[1][slot];
      }
   }

int GenotypeInference::GetSlot(int individual, int marker)
   {
   marker = markerKey[marker];

   if (marker < 0) return -1;

   return individual * markers + marker;
   }

void GenotypeInference::InferGenotypes(
    Mantra & mantra, TreeInfo & stats,
    Tree & withMarker, Tree & without, Tree & single,
    int marker)
   {
   Pedigree & ped  = *mantra.pedigree;
   Family * family = mantra.family;

   MarkerCluster * cluster = clusters.Enabled() ? clusters.markerToCluster[marker] : NULL;

   FuzzyInheritanceTree alternative;

   double multipoint_likelihood = stats.GetMean(withMarker);

   int markers = cluster == NULL ? 1 : cluster->markerIds.Length();

   for (int m = 0; m < markers; m++)
      {
      double singlepoint_likelihood;

      if (cluster != NULL)
         {
         marker = cluster->markerIds[m];
         mantra.SelectMarker(marker);
         alternative.FuzzyScoreVectors(mantra);

         singlepoint_likelihood = stats.GetMean(alternative);
         }
      else
         singlepoint_likelihood = stats.GetMean(single);

      for (int id = family->first; id <= family->last; id++)
         {
         int slot = GetSlot(id, marker);
         if (slot < 0) continue;

         // If the genotype is known, there is nothing to calculate
         int genotype = ped[id].markers[marker].BinaryCoded();

         if (genotype >= 0 /* and the error rate is zero */ )
            switch (genotype)
               {
               case 1: probability[0][slot] = 1.0; probability[1][slot] = 0.0; continue;
               case 2: probability[0][slot] = 0.0; probability[1][slot] = 1.0; continue;
               case 3: probability[0][slot] = 0.0; probability[1][slot] = 0.0; continue;
               }

         // Otherwise, we first try to set the genotype as if it were homozygous
         // for allele 1
         ped[id].markers[marker][0] = ped[id].markers[marker][1] = 1;

         // Then evaluate the likelihood under that setting
         mantra.SelectMarker(marker);
         alternative.FuzzyScoreVectors(mantra);
         double alternative_likelihood = stats.GetMean(alternative);

         // If the other family members fix the genotype, then we are done ...
         if (alternative_likelihood == singlepoint_likelihood)
            {
            probability[0][slot] = 1.0;
            probability[1][slot] = 0.0;

            // We need to restore the old missing genotype so as
            // not to disrupt missing data patterns for other analyses
            ped[id].markers[marker][0] = ped[id].markers[marker][1] = 0;

            continue;
            }

         // Otherwise, check if we need to update multipoint likelihood
         if (alternative_likelihood == 0.0)
            // Genotype is impossible
            probability[0][slot] = 0.0;
         else
            {
            // If we get here we must recalculate the full multipoint likelihood
            if (cluster != NULL)
               cluster->ScoreLikelihood(mantra, alternative);

            double alternative_multipoint = alternative.MeanProduct(without);

            // This is the conditional probability of the current genotype
            probability[0][slot] = alternative_multipoint / multipoint_likelihood *
                                   exp(alternative.logOffset - single.logOffset);
            }

         // Now, repeat the process for the other homozygous genotype
         ped[id].markers[marker][0] = ped[id].markers[marker][1] = 2;

         // Then evaluate the likelihood under that setting
         mantra.SelectMarker(marker);
         alternative.FuzzyScoreVectors(mantra);
         alternative_likelihood = stats.GetMean(alternative);

         // If the other family members fix the genotype, then we are done ...
         if (alternative_likelihood == singlepoint_likelihood)
            {
            probability[0][slot] = 0.0;
            probability[1][slot] = 1.0;

            // We need to restore the old missing genotype so as
            // not to disrupt missing data patterns for other analyses
            ped[id].markers[marker][0] = ped[id].markers[marker][1] = 0;

            continue;
            }

         // Otherwise, check if we need to update multipoint likelihood
         if (alternative_likelihood == 0.0)
            // Genotype is impossible
            probability[1][slot] = 0.0;
         else
            {
            // If we get here we must recalculate the full multipoint likelihood
            if (cluster != NULL)
               cluster->ScoreLikelihood(mantra, alternative);

            double alternative_multipoint = alternative.MeanProduct(without);

            // This is the conditional probability of the current genotype
            probability[1][slot] = alternative_multipoint / multipoint_likelihood *
                                   exp(alternative.logOffset - single.logOffset);
            }

         ped[id].markers[marker][0] = ped[id].markers[marker][1] = 0;
         }
      }
   }

int GenotypeInference::SelectMarkers(IntArray & markerList, Vector & markerPositions, Vector & analysisPositions)
   {
   markers = 0;

   markerKey.Dimension(Pedigree::markerCount);
   markerKey.Set(-1);

   frequencies.Clear();

   for (int i = 0; i < markerList.Length(); i++)
      {
      MarkerInfo * info = Pedigree::GetMarkerInfo(markerList[i]);

      if (info->freq.Length() != 3) continue;
      if (analysisPositions.FastFind(markerPositions[i]) < 0) continue;

      frequencies.Push(info->freq[1]);
      markerKey[markerList[i]] = markers++;
      }

   return markers;
   }

bool GenotypeInference::resultsAvailable(int marker)
   {
   return markerKey[marker] >= 0;
   }

void GenotypeInference::OutputGenotypes(Pedigree & ped, IntArray & markerList)
   {
   FILE * datfile = fopen(MerlinCore::filePrefix + "-infer.dat", "wt");
   FILE * pedfile = fopen(MerlinCore::filePrefix + "-infer.ped", "wt");

   if (datfile == NULL || pedfile == NULL)
      error("Error opening output files for inferred genotypes\n");

   bool * flips = new bool [markers];

   for (int i = 0, index = 0; i < markerList.Length(); i++)
      if (resultsAvailable(markerList[i]))
         {
         MarkerInfo * info =  ped.GetMarkerInfo(markerList[i]);

         flips[index++] = info->GetAlleleLabel(1) > info->GetAlleleLabel(2);
         }

   String * alleles[2];

   alleles[0] = new String[markers];
   alleles[1] = new String[markers];

   for (int i = 0, index = 0; i < markerList.Length(); i++)
      if (resultsAvailable(markerList[i]))
         {
         const char * label = ped.markerNames[markerList[i]];

         MarkerInfo * info =  ped.GetMarkerInfo(markerList[i]);

         alleles[0][index] = info->GetAlleleLabel(flips[index] ? 2 : 1);
         alleles[1][index] = info->GetAlleleLabel(flips[index] ? 1 : 2);

         const char * label1 = alleles[0][index];
         const char * label2 = alleles[1][index];

         if (inferBest)
            fprintf(datfile, "M %s\n", label);

         if (inferExpected)
            fprintf(datfile, "T COUNT(%s,%s)\n", label1, label);

         if (inferProbabilities)
            {
            fprintf(datfile, "C P(%s=%s/%s)\n", label, label1, label1);
            fprintf(datfile, "C P(%s=%s/%s)\n", label, label1, label2);
            fprintf(datfile, "C P(%s=%s/%s)\n", label, label2, label2);
            }

         index++;
         }

   for (int i = 0; i < ped.count; i++)
      {
      fprintf(pedfile, "%s\t%s\t%s\t%s\t%d",
             (const char *) ped[i].famid, (const char *) ped[i].pid,
             (const char *) ped[i].fatid, (const char *) ped[i].motid, ped[i].sex);

      for (int j = 0, index = 0; j < markerList.Length(); j++)
         if (resultsAvailable(markerList[j]))
            {
            int slot  = GetSlot(i, markerList[j]);
            double z0 = probability[0][slot];
            double z2 = probability[1][slot];
            double z1 = 1.0 - z0 - z2;

            if (flips[index] && z0 >= 0.0)
               { double swap = z0; z0 = z2; z2 = swap; }

            if (z0 < 0.0)
               // These genotypes were not inferred, perhaps because no
               // genotyped relatives were available ... we use their
               // population frequencies instead
               {
               double frequency = frequencies[markerKey[markerList[j]]];

#ifdef __CHROMOSOME_X__
               if (ped[i].sex == SEX_MALE)
                  {
                  z0 = frequency;
                  z1 = 0.0;
                  z2 = 1.0 - frequency;
                  }
               else
#endif
                  {
                  z0 = frequency * frequency;
                  z1 = 2.0 * frequency * (1.0 - frequency);
                  z2 = (1.0 - frequency) * (1.0 - frequency);
                  }
               }

            int best = z0 > z1 ? (z0 > z2 ? 1 : 3) : (z1 > z2 ? 2 : 3);

            if (inferBest)
               if (best == 1 && z0 > .95 || best == 2 && z1 > 0.95 || best == 3 && z2 > 0.95)
                  fprintf(pedfile, " %s/%s",
                    (const char *) (best == 3 ? alleles[1][index] : alleles[0][index]),
                    (const char *) (best == 1 ? alleles[0][index] : alleles[1][index]));
               else
                  fprintf(pedfile, " ./.");

            if (inferExpected)
               fprintf(pedfile, " %.3f", 2 * z0 + z1);

            if (inferProbabilities)
               fprintf(pedfile, " %.3f %.3f %.3f",  z0, z1, z2);

            index++;
            }

      fprintf(pedfile, "\n");
      }

   printf("Infered genotypes stored in [%s-infer.ped]\n",
          (const char *) MerlinCore::filePrefix);

   delete [] flips;
   delete [] alleles[0];
   delete [] alleles[1];

   fclose(datfile);
   fclose(pedfile);
   }

void GenotypeInference::CreateOfflineMarkers(Pedigree & ped)
   {
   for (int i = 0; i < Pedigree::traitCount; i++)
      {
      String & label = Pedigree::traitNames[i];

      if (label.Left(6) != "COUNT(") continue;

      int openp  = label.FindChar('(');
      int closep = label.FindChar(')');
      int comma  = label.FindChar(',');

      if (openp < 0 || comma <= openp || closep <= comma) continue;

      String allele = label.Mid(openp + 1, comma - 1);
      String marker = label.Mid(comma + 1, closep - 1);

      Pedigree::GetMarkerID(marker);

      MarkerInfo * info = Pedigree::GetMarkerInfo(marker);

      if (info->CountAlleles() > 2)
         error("MERLIN can't handle inferred allele counts for marker '%s', since it is not bi-allelic\n",
              (const char *) info->name);

      int allele_number = info->GetAlleleNumber(allele);

      if (allele_number != 1 && allele_number != 2)
         error("Your file with inferred genotypes refers to allele '%s' at marker '%s'\n"
                "However, this allele is not present in the pedigree file or allele frequency file\n\n",
                (const char *) allele, (const char *) marker);

      if (allele_number == 2)
         {
         // Flip reference allele
         for (int j = 0; j < ped.count; j++)
            if (ped[j].traits[i] != _NAN_)
               ped[j].traits[i] = 2.0 - ped[j].traits[i];

         Pedigree::traitLookup.Delete(label);
         label.printf("COUNT(%s,%s)", (const char *) info->alleleLabels[1],
                                      (const char *) marker);
         Pedigree::traitLookup.Add(label, i);
         }
      }
   }

void GenotypeInference::CreateOfflineLookup()
   {
   markerKey.Dimension(Pedigree::markerCount);
   markerKey.Set(-1);

   frequencies.Dimension(Pedigree::markerCount);

   for (int i = 0; i < Pedigree::traitCount; i++)
      {
      String & label = Pedigree::traitNames[i];

      if (label.Left(6) != "COUNT(") continue;

      int openp  = label.FindChar('(');
      int closep = label.FindChar(')');
      int comma  = label.FindChar(',');

      if (openp < 0 || comma <= openp || closep <= comma) continue;

      String allele = label.Mid(openp + 1, comma - 1);
      String marker = label.Mid(comma + 1, closep - 1);

      int markerid = Pedigree::GetMarkerID(marker);
      MarkerInfo * info = Pedigree::GetMarkerInfo(markerid);

      if (info->freq.Length() < 1)
         error("Missing allele frequency information for marker '%s'",
               (const char *) marker);

      markerKey[markerid] = i;
      frequencies[markerid] = info->freq[1];
      }
   }

void GenotypeInference::CheckFrequencies(Pedigree & ped)
   {
   bool checkFailed = false;

   if (!offline) return;

   printf("Checking that allele doses and frequencies are consistent ...\n");

   double avg_diff = 0.0, favg_diff = 0.0;
   double avg_diff_count = 0, favg_diff_count = 0;

   for (int i = 0; i < Pedigree::markerCount; i++)
      {
      MarkerInfo * info = Pedigree::GetMarkerInfo(i);
      double freq = info->freq[1];

      if (frequencies[i] != freq)
         {
         printf("  Internal inconsistency in allele frequencies for marker %s\n",
                (const char *) Pedigree::markerNames[i]);
         checkFailed = true;
         continue;
         }

      int dosage = markerKey[i];

      // Estimate frequency among all individuals ...

      int counts = 0;
      double avg = 0.0;

      for (int j = 0; j < ped.count; j++)
         if (ped[j].traits[i] != _NAN_)
            {
            avg += ped[j].traits[dosage];
            counts++;
            }

      if (counts == 0) continue;
      avg /= 2 * counts;
      avg_diff += fabs(avg - freq);
      avg_diff_count++;

      // Estimate frequency among founders ...

      int fcounts = 0;
      double favg = 0.0;

      for (int j = 0; j < ped.count; j++)
         if (ped[j].traits[dosage] != _NAN_ && ped[j].isFounder())
            {
            favg += ped[j].traits[dosage];
            fcounts++;
            }

      if (fcounts != 0)
         {
         favg /= 2 * fcounts;
         favg_diff += fabs(favg - freq);
         favg_diff_count++;
         }

      // Check for problems and report summary

      if ( (freq - 0.5) * (favg - 0.5) < 0.0 || (freq - 0.5) * (avg - 0.5) < 0.0)
         printf("  Potential flip for %s [freq: %.3f, all: %.3f, founders: %.3f]\n",
                (const char *) Pedigree::markerNames[i], freq, avg, favg);
      else if (fabs(freq - avg) > 0.02 || fabs(freq - favg) > 0.02)
         printf("  Large differency for %s [freq: %.3f, all: %.3f, founders: %.3f]\n",
                (const char *) Pedigree::markerNames[i], freq, avg, favg);
      }

   if (checkFailed) error("Sanity check failed -- aborting due to internal inconsistency");

   if (favg_diff_count) favg_diff /= favg_diff_count;
   if (avg_diff_count) avg_diff /= avg_diff_count;

   printf("  Average difference between doses and frequencies %.3f [all], %.3f [founders]\n\n",
          avg_diff, favg_diff);

   return;
   }

void GenotypeInference::UpdateFrequencies(Pedigree & ped, bool foundersOnly)
   {
   if (!offline) return;

   printf("Updating SNP allele frequencies using estimated genotype counts ...\n");

   double avg_diff = 0.0;
   double avg_diff_count = 0;

   for (int i = 0; i < Pedigree::markerCount; i++)
      {
      MarkerInfo * info = Pedigree::GetMarkerInfo(i);

      if (info->CountAlleles() != 2) continue;

      double freq = info->freq[1];

      int dosage = markerKey[i];

      // Estimate frequency among founders ...

      int counts = 0;
      double avg = 0.0;

      for (int j = 0; j < ped.count; j++)
         if (ped[j].traits[dosage] != _NAN_)
            if (!foundersOnly || ped[j].isFounder())
            {
            avg += ped[j].traits[dosage];
            counts++;
            }

      if (counts == 0) continue;

      avg /= 2 * counts;

      info->freq[1] = frequencies[i] = avg;
      info->freq[2] = 1.0 - avg;

      avg_diff += fabs(avg - freq);
      avg_diff_count++;
      }

   if (avg_diff_count) avg_diff /= avg_diff_count;

   printf("  Adjusted allele frequencies by an average of %.3f\n\n", avg_diff);

   return;
   }

 
