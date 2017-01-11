////////////////////////////////////////////////////////////////////// 
// merlin/MerlinCluster.cpp 
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
 
#include "MerlinCluster.h"
#include "MerlinCore.h"
#include "MathConstant.h"
#include "HaploFamily.h"
#include "SparseLikelihood.h"
#include "Likelihood.h"
#include "Error.h"

#include <math.h>

MarkerClusters clusters;

int MarkerCluster::FindHaplotype(IntArray & haplo)
   {
   if (markerIds.Length() == 0)
      error("MarkerCluster -- Cannot locate haplotypes in empty clusters");

   if (markerIds.Length() != haplo.Length())
      error("MarkerCluster -- Unexpected haplotype length");

   for (int h = 0; h < alleles.Length(); h += markerIds.Length())
      {
      int match = true;

      for (int i = 0; i < markerIds.Length(); i++)
         if (haplo[i] != alleles[h + i])
            {
            match = false;
            break;
            }

      if (match) return h / markerIds.Length();
      }

   return -1;
   }

void MarkerCluster::AddHaplotype(IntArray & haplo, double freq)
   {
   if (markerIds.Length() == 0)
      error("MarkerCluster -- Cannot add haplotypes to empty cluster");

   if (markerIds.Length() != haplo.Length())
      error("MarkerCluster -- Unexpected haplotype length");

   int index = FindHaplotype(haplo);

   if (index >= 0)
      freqs[index] += freq;
   else
      {
      for (int i = 0; i < markerIds.Length(); i++)
         alleles.Push(haplo[i]);

      freqs.Push(freq);
      }
   }

void MarkerCluster::CheckFrequencies()
   {
   double sum = freqs.Sum();

   if ( fabs(sum - 1.0) > 1e-5 )
      {
      String label = "CLUSTER: ";

      for (int i = 0; i < markerIds.Length(); i++)
         label += Pedigree::markerNames[markerIds[i]] + " ";

      if ( fabs(sum) < ZEPS )
         error("%s\n  Haplotype frequencies sum to zero\n",
               (const char *) label);

      printf("%s\n  Haplotype frequencies sum to %f, adjusted to 1.0\n",
             (const char *) label, sum);
      }

   if ( sum != 1.0)
      freqs *= 1.0 / sum;
   }

void MarkerCluster::UpdateAlleleCounts()
   {
   int markers = markerIds.Length();

   alleleCounts.Dimension(markers);

   for (int i = 0; i < markers; i++)
      {
      MarkerInfo * info = Pedigree::GetMarkerInfo(markerIds[i]);

      alleleCounts[i] = info->freq.Length() - 1;

      MerlinCore::ValidateMarker(info);
      }
   }

bool MarkerCluster::EstimateFrequencies(Pedigree & ped, String * messages)
   {
   // Keep user specified frequencies, if these are available
   if (freqs.Length() != 0)
      return false;

   UpdateAlleleCounts();

   FamilyHaplos engine;
   HaplotypeSets  sets;

   engine.markers = markerIds;

   // sets->AllocateAlleleLabels(info.markers);
   // for (int m = 0; m < info.markers; m++)
   //   sets->SetAlleleLabels(m, Pedigree::GetMarkerInfo(info.markerIds[m])->alleleLabels);

   int families = 0;

   for (int f = 0; f < ped.familyCount; f++)
      try
         {
         if (engine.Haplotype(ped, *(ped.families[f])))
            {
            sets.LoadFromMemory(&engine, markerIds.Length());
            families++;
            }
         else if (engine.mendel_errors.Length() || engine.obligate_recombinant)
            families++;
         }
      catch (const TreesTooBig & problem)
         {
         sets.discarded.Push("");
         sets.discarded.Last().printf("%s -  SKIPPED: >%d megabytes needed",
                       (const char *) ped.families[f]->famid, problem.memory_request);
         }

   // Don't try to estimate frequencies if there are no informative families
   if (sets.Length() == 0)
      error("All families have an obligate recombinant or Mendelian inconsistency\n"
            "In %d-marker cluster starting with marker %s\n",
            markerIds.Length(), (const char *) ped.markerNames[markerIds[0]]);

   if (sets.Length() < (int) ceil(families * 0.95))
      {
      if (markerIds.Length() == 1)
         printf("WARNING -- %d of %d families ha%s a Mendelian inconsistency for marker %s\n",
             families - sets.Length(), families,
             families - sets.Length() == 1 ? "s" : "ve",
             (const char *) ped.markerNames[markerIds[0]]);
      else
         printf("WARNING -- %d of %d families ha%s an obligate recombinant or Mendelian inconsistency\n"
            "           In %d-marker cluster starting with marker %s\n",
            families - sets.Length(), families,
            families - sets.Length() == 1 ? "s" : "ve",
            markerIds.Length(), (const char *) ped.markerNames[markerIds[0]]);
      }

   // Check how many haplotype frequencies have to be estimates
   double haploCount = alleleCounts.DoubleProduct();

   if (haploCount < 16384)
      {
      sets.quiet = true;
      sets.Filter(alleleCounts, 20);

      if (sets.discarded.Length())
         {
         // printf("  Exclusions, see [%s] for details.\n", logname);
         messages->catprintf("To speed up calculations, the following %d famil%s not used to estimate %s frequencies %s marker %s:\n",
                            sets.discarded.Length(),
                            sets.discarded.Length() == 1 ? "y was" : "ies were",
                            markerIds.Length() > 1 ? "haplotype" : "allele",
                            markerIds.Length() > 1 ? "in cluster with" : "for",
                            (const char *) Pedigree::markerNames[markerIds[0]]);

         // List skipped families
         for (int i = 0; i < sets.discarded.Length(); i++)
            messages->catprintf("  %s\n", (const char *) sets.discarded[i]);
         }

      Likelihood likelihood;
      likelihood.Initialize(alleleCounts);
      likelihood.EM(&sets);
      likelihood.RetrieveHaplotypes(alleles, freqs, 1e-8);

      #ifdef __DEBUG__
         likelihood.Print(1e-5);
      #endif
      }
   else
      {
      SparseLikelihood likelihood;

      likelihood.Initialize(alleleCounts);
      likelihood.EM(&sets, true /* quiet */);
      likelihood.RetrieveHaplotypes(alleles, freqs, 1e-8);

      #ifdef __DEBUG__
         likelihood.Print(1e-5);
      #endif
      }

   // Haplotype frequency estimation succeeded
   return true;
   }

void MarkerCluster::UpdateAlleleFrequencies()
   {
   // Update marker allele frequencies to reflect the ML
   // estimates of haplotype frequencies
   for (int i = 0; i < markerIds.Length(); i++)
      {
      MarkerInfo * info = Pedigree::GetMarkerInfo(markerIds[i]);

      info->freq.Zero();
      for (int j = 0; j * markerIds.Length() < alleles.Length(); j++)
         info->freq[alleles[j * markerIds.Length() + i]] += freqs[j];

      // info->freq /= info->freq.Sum() + 1e-30;
      }
   }

void MarkerCluster::ScoreLikelihood(Mantra & m, Tree & tree)
   {
   const char * msg = "Cluster at marker %s dropped [%s]";
   const char * mrk = Pedigree::markerNames[markerIds[0]];

   if (freqs.Length() == 0)
      error("MarkerCluster -- Missing haplotype frequency information\n");

   if (alleleCounts.Length() == 0)
      UpdateAlleleCounts();

   FamilyHaplos  engine;
   HaplotypeSet  set;

   engine.markers = markerIds;
   if (!engine.Haplotype(tree, *m.pedigree, *m.family))
      {
      // Bad inheritance or obligate recombinant within cluster
      if (engine.mendel_errors.Length())
         errormsg.printf(msg, mrk, m.family->count > 1 ? "BAD INHERITANCE" : "HETEROZYGOUS MALE");
      else if (engine.obligate_recombinant)
         errormsg.printf(msg, mrk, "OBLIGATE RECOMBINANT");

      tree.MakeMinimalTree(0.0, m.bit_count);
      return;
      }
   else
      errormsg.Clear();

   int markers = markerIds.Length();

   set.LoadFromMemory(&engine, markers);

   SparseLikelihood likelihood;

   likelihood.Initialize(alleleCounts, alleles, freqs);
   likelihood.SetLikelihood(&set);

   engine.RetrieveLikelihoods(tree);
   tree.logOffset = likelihood.GetLogOffset(&set);

   switch (InheritanceTree::mergingStrategy)
      {
      case MERGE_ALL :
         tree.Trim();
         break;
      case MERGE_ZEROS :
         tree.TrimNoMerge();
         break;
      case MERGE_BOOLEAN :
         error("MarkerCluster::ScoreLikelihood encountered an inconsistent state");
         break;
      }

   if (!tree.FindNonZero())
      errormsg.printf(msg, mrk, "UNKNOWN HAPLOTYPE");
   }

void MarkerCluster::GetHaplotypes(Mantra & m, IntArray & vector, IntArray & haplotypes, bool sample)
   {
   FamilyHaplos  engine;
   HaplotypeSet  set;

   engine.markers = markerIds;
   if (engine.Haplotype(*m.pedigree, *m.family, vector))
      {
      set.LoadFromMemory(&engine, markerIds.Length());

      SparseLikelihood likelihood;

      likelihood.Initialize(alleleCounts, alleles, freqs);
      likelihood.ExtractConfiguration(set.graphs[0], haplotypes, sample);
      }
   else
      {
      haplotypes.Dimension(m.two_f);
      haplotypes.Set(-1);
      }

   // Calling functions will usually expect Mantra to be updated to
   // reflect the requested inheritance vector, so we should oblige
   // Setup the inheritance vector
   for (int i = 0; i < m.bit_count; i++)
      {
      int bit = m.bits[i];

      m.vector[bit] = vector[i] == 0 ? m.vector[bit] & ~1 : m.vector[bit] | 1;
      }

   // Identify founder alleles for each individual
   for (int i = m.two_f; i < m.two_n; i++)
      m.state[i] = m.state[m.vector[i]];
   }

void MarkerCluster::RetrieveHaplotype(IntArray & al, int index)
   {
   int markers = markerIds.Length();
   int base = markers * index;

   al.Dimension(markers);
   for (int i = 0; i < markers; i++)
       al[i] = alleles[base++];
   }

MarkerClusters::MarkerClusters()
   {
   count = 0;
   head = NULL;
   markerToCluster = NULL;
   haveOutput = false;
   }

MarkerClusters::~MarkerClusters()
   {
   DeleteList();
   }

void MarkerClusters::DeleteList()
   {
   while (head != NULL)
      {
      MarkerClusterList * old = head;
      head = head->next;
      delete old;
      }

   if (markerToCluster != NULL)
      delete [] markerToCluster;
   }

void MarkerClusters::LoadFromFile(const char * filename, const char * logname)
   {
   if (filename == NULL || filename[0] == 0)
      return;

   IFILE input = ifopen(filename, "rb");

   if (input == NULL)
     {
     printf("WARNING! - Could not open file '%s' with clustering information!\n",
          filename);
     return;
     }

   LoadFromFile(input, logname);

   ifclose(input);
   }

void MarkerClusters::AllocateCrossReference()
   {
   // Initialize lookup table between markers and clusters
   if (markerToCluster == NULL)
      {
      markerToCluster = new MarkerClusterList * [Pedigree::markerCount];

      for (int i = 0; i < Pedigree::markerCount; i++)
         markerToCluster[i] = NULL;
      }
   }

void MarkerClusters::LoadFromFile(IFILE & input, const char * logname)
   {
   AllocateCrossReference();

   // Error message
   const char * errormsg = "Error reading marker cluster information:\n"
      "  *** %s ***\n\n"
      "The following line was encounterd:\n%s\n\n"
      "Each cluster should have a header line listing markers in the cluster\n"
      "Followed by multiple lines declaring haplotypes and their frequencies\n\n"
      "An example declaration follows:\n"
      "CLUSTER MARKER1 MARKER2 MARKER3\n"
      "HAPLO 0.15  1 1 1\n"
      "HAPLO 0.85  2 2 2\n\n";

   // Input buffers
   String      buffer, cluster("CLUSTER"), haplabel("HAPLO");
   StringArray tokens;
   IntArray    haplo, markerIds;
   bool        buffered = false;

   // Log file for additional debugging messages
   FILE * logfile = NULL;
   bool   logOpen = false, logMap = false, logMissing = 0;

   // Loop until input is exhausted
   while (!ifeof(input) || buffered)
      {
      if (!buffered) buffer.ReadLine(input);
      tokens.Clear();
      tokens.AddTokens(buffer);
      buffered = false;

      if (tokens.Length() == 0) continue;

      if (cluster.SlowCompareToStem(tokens[0]) != 0)
         error(errormsg, "Invalid cluster header line", (const char *) buffer);

      // Retrieve marker ids for each marker
      markerIds.Clear();
      for (int i = 1; i < tokens.Length(); i++)
         markerIds.Push(Pedigree::LookupMarker(tokens[i]));

      // Automatically skip haplotype information if none of the markers
      // is in the currently loaded pedigree
      if (markerIds.Length() == 0 || markerIds.Max() < 0)
         {
         while (!ifeof(input))
            {
            buffer.ReadLine(input);
            tokens.Clear();
            tokens.AddTokens(buffer);

            if (haplabel.SlowCompareToStem(tokens[0]) != 0)
               {
               buffered = true;
               break;
               }
            }

         continue;
         }

      // Add new cluster to chain
      head = new MarkerClusterList(head);
      count++;

      // Retrieve marker info for each marker
      MarkerInfo ** info = new MarkerInfo * [markerIds.Length()];
      MarkerInfo * first = NULL;

      bool differentPositions = false;
      for (int i = 0; i < markerIds.Length(); i++)
         if (markerIds[i] >= 0)
            {
            // Retrieve corresponding marker info structures
            info[i] = Pedigree::GetMarkerInfo(markerIds[i]);
            info[i]->UpdateSerial();

            // Check that the marker has been assigned to a chromosome
            if (info[i]->chromosome == -1)
               {
               if (logOpen == false)
                  {
                  logOpen = true;
                  logfile = logname == NULL ? stderr : fopen(logname, "wt");
                  }

               if (logMissing == false)
                  {
                  logMissing = true;
                  printf("MARKER CLUSTERS: Unmapped markers were dropped see [%s]\n",
                         logname == NULL ? "STDERR" : logname);
                  }

               fprintf(logfile, "DROPPED marker %s, map position unknown\n",
                                (const char *) info[i]->name);

               // Discard unmapped markers
               markerIds[i] = -1;
               continue;
               }

            // Check that marker was not previously clustered
            if (markerToCluster[markerIds[i]] != NULL)
               error("Marker %s appears in two different clusters\n",
                     (const char *) info[i]->name);

            // Check that all markers within cluster map to the same chromosome
            if (first == NULL)
               first = info[i];
            else
               if (info[i]->chromosome != first->chromosome)
                   error("Markers %s and %s are in the same cluster\n"
                         "but map to different chromosomes\n",
                         (const char *) info[i]->name, (const char *) first->name);
               else if (info[i]->position != first->position)
                  differentPositions = true;

            // Update cluster information
            head->markerIds.Push(markerIds[i]);
            markerToCluster[markerIds[i]] = head;
            }

      if (differentPositions)
         {
         if (logOpen == false)
            {
            logOpen = true;
            logfile = logname == NULL ? stderr : fopen(logname, "wt");
            }

         if (logMap == false)
            {
            logMap = true;
            printf("MARKER CLUSTERS: Marker map changed, see [%s]\n",
                   logname == NULL ? "STDERR" : logname);
            }

         if (logfile != NULL)
            {
            fprintf(logfile, "CLUSTER:");
            for (int i = 0; i < head->markerIds.Length(); i++)
               fprintf(logfile, " %s",
                       (const char *) Pedigree::markerNames[head->markerIds[i]]);
            AdjustPositions(head->markerIds);
            fprintf(logfile, "\n         Map positions differ, adjusted to %.3f cM\n",
                    Pedigree::GetMarkerInfo(head->markerIds[0])->position * 100);
            }
         }

      while (!ifeof(input))
         {
         buffer.ReadLine(input);
         tokens.Clear();
         tokens.AddTokens(buffer);

         if (tokens.Length() == 0) continue;

         if (haplabel.SlowCompareToStem(tokens[0]) != 0)
            {
            buffered = true;
            break;
            }

         if (tokens.Length() != markerIds.Length() + 2)
            error(errormsg, "Incorrect number of allele states in haplotype", (const char *) buffer);

         double freq = tokens[1].AsDouble();

         if (freq < 0.0)
            error(errormsg, "Haplotype frequencies should be positive", (const char *) buffer);

         haplo.Clear();
         for (int i = 0; i < markerIds.Length(); i++)
            if (markerIds[i] >= 0)
               {
               int allele = Pedigree::LoadAllele(markerIds[i], tokens[i + 2]);

               if (allele <= 0)
                  {
                  error(errormsg, "Invalid allele label encountered", (const char *) buffer);

                  // printf("  Allele label for marker %s is invalid\n",
                  //           (const char *) Pedigree::markerNames[markerIds[i]]);
                  }
               else if (allele >= info[i]->freq.Length())
                  {
                  int previous = info[i]->freq.Length();

                  info[i]->freq.Dimension(allele + 1);
                  for (int j = previous; j <= allele; j++)
                     info[i]->freq[j] = 0.0;

                  // printf("  Allele size range expanded to %d for marker %s\n",
                  //        allele, (const char *) Pedigree::markerNames[markerIds[i]]);
                  }

               haplo.Push(allele);
               }

         head->AddHaplotype(haplo, freq);
         }
      }

   if (logfile != NULL)
      fclose(logfile), haveOutput = true;

   if (count)
      {
      printf("MARKER CLUSTERS: User supplied file defines %d cluster%s\n",
             count, count == 1 ? "" : "s");
      haveOutput = true;
      }
   }

int MarkerClusters::SaveToFile(const char * filename)
   {
   if (count == 0)
      return 0;

   FILE * output = fopen(filename, "wt");

   if (output == NULL)
      {
      haveOutput = true;

      printf("File access error while saving marker clusters...\n"
             "  The file '%s' cannot be opened.\n"
             "Marker clusters have not been saved.\n\n",
             filename);

      return 0;
      }

   SaveToFile(output);

   fclose(output);

   return count;
   }

int MarkerClusters::SaveToFile(FILE * output)
   {
   if (count == 0)
      return 0;

   MarkerClusterList * current = head;

   while (current != NULL)
      {
      int markers = current->CountMarkers();

      // At the beginning of each cluster, list marker names
      fprintf(output, "CLUSTER ");
      for (int i = 0; i < markers; i++)
         fprintf(output, "%s ", (const char *) Pedigree::markerNames[current->markerIds[i]]);
      fprintf(output, "\n");

      // Then, list haplotypes and their frequencies
      for (int h = 0; h < current->CountHaplotypes(); h++)
         {
         fprintf(output, "HAPLO %5.4f ", current->freqs[h]);
         for (int i = 0; i < markers; i++)
            fprintf(output, "%3s ",
               (const char *)
               Pedigree::GetMarkerInfo(current->markerIds[i])->GetAlleleLabel(
                  current->alleles[h * markers + i]));
         fprintf(output, "\n");
         }

      current = current->next;
      }

   // printf("MARKER CLUSTERS: %d cluster%s saved to file\n", count, count == 1 ? "" : "s");
   // haveOutput = true;

   return count;
   }

void MarkerClusters::ClusterByDistance(double minDistance)
   {
   if (minDistance == _NAN_ * 0.01)
      return;

   if (count)
      {
      printf("  Distance clustering criterium ignored\n");
      return;
      }

   IntArray markers;

   AllocateCrossReference();

   int chromosome = -2;

   do {
      markers.Dimension(0);

      chromosome = PedigreeGlobals::SortMarkersInMapOrder(markers, chromosome);
      int markerCount = markers.Length();

      if (markerCount == 0)
         break;

      int clusterStart = 0;
      double lastPosition = Pedigree::GetMarkerInfo(markers[0])->position;

      for (int i = 1; i <= markers.Length(); i++)
         {
         double position = (i == markers.Length()) ? position + minDistance + 1:
                     Pedigree::GetMarkerInfo(markers[i])->position;

         // Whenever two consecutive markers are separated by more
         // than a certain distance, we start a new cluster
         if (position - lastPosition > minDistance)
            {
            // If the previous cluster, had multiple markers
            // then we track it. Trivial clusters can be safely
            // ignored
            if (clusterStart != i - 1)
               CreateCluster(markers, clusterStart, i - 1);

            clusterStart = i;
            }

         lastPosition = position;
         }
   } while (chromosome);

   printf("MARKER CLUSTERS: %d cluster%s where intermarker distance is < %3g cM.\n",
      count, count == 1 ?  "" : "s", minDistance * 100);
   haveOutput = true;
   }

void MarkerClusters::EstimateFrequencies(Pedigree & ped, const char * logname)
   {
   String messages;
   int  done = 1;
   bool nothing_to_do = true;
   MarkerClusterList * current = head;

   while (current != NULL)
      {
      if (current->EstimateFrequencies(ped, &messages))
         {
         if (!MerlinCore::quietOutput)
            printf("MARKER CLUSTERS: Estimated haplotype frequencies for cluster %d of %d\r",
                   done, count);
         fflush(stdout);
         nothing_to_do = false;
         }
      current->UpdateAlleleFrequencies();
      current = current->next;
      done++;
      }

   if (!nothing_to_do)
      {
      printf("MARKER CLUSTERS: Estimated haplotype frequencies "
                              "for %d clusters%20s\n", count, "");
      haveOutput = true;
      }

   if (messages.Length())
      {
      printf("MARKER CLUSTERS: Some families skipped due to computing limitations, see [%s]\n", logname);

      FILE * f = fopen(logname, "wt");

      if (f != NULL)
         {
         fprintf(f, "%s", (const char *) messages);
         fclose(f);
         }
      else
         printf("MARKER CLUSTERS: Error opening log file [%s]!\n", logname);
      }
   }

#ifndef min
#define min(a,b)   ((a) < (b) ? (a) : (b))
#endif

void MarkerClusters::ClusterByRsquared(Pedigree & ped, double threshold)
   {
   if (threshold == _NAN_)
      return;

   if (count)
      {
      printf("  R-squared clustering criterium ignored\n");
      return;
      }

   AllocateCrossReference();

   IntArray markers;
   int chromosome = -2;

   Likelihood     likelihood;
   FamilyHaplos   engine;
   HaplotypeSets  sets;
   IntArray       alleleCounts;

   engine.markers.Dimension(2);
   alleleCounts.Dimension(2);

   bool   firstPass = true;
   double rsq, dprime;

   int    discards = 0;

   do {
      markers.Dimension(0);

      chromosome = PedigreeGlobals::SortMarkersInMapOrder(markers, chromosome);
      int markerCount = markers.Length();

      if (markerCount == 0)
         break;

      if (firstPass)
         printf("MARKER CLUSTERS: Finding clusters of <%d markers where\n"
                "                 r2 between end-point markers is >%.3f\n",
                CLUSTER_SEARCH_MAX, threshold, firstPass = false);

      for (int clusterStart = 0; clusterStart < markerCount; clusterStart++)
         {
         engine.markers[0] = markers[clusterStart];
         alleleCounts[0] = Pedigree::GetMarkerInfo(markers[clusterStart])->freq.Length() - 1;

         if (alleleCounts[0] != 2) continue;

         int clusterEnd = min(clusterStart + 20, markerCount - 1);

         if (!MerlinCore::quietOutput)
            printf("MARKER CLUSTERS: Checking marker %d of %d (chromosome %d)      \r",
               clusterStart, markerCount, Pedigree::GetMarkerInfo(markers[0])->chromosome);
         fflush(stdout);

         while (clusterEnd > clusterStart)
            {
            engine.markers[1] = markers[clusterEnd];
            alleleCounts[1] = Pedigree::GetMarkerInfo(markers[clusterEnd])->freq.Length() - 1;

            if (alleleCounts[1] != 2)
               {
               clusterEnd--;
               continue;
               }

            sets.ClearContents();
            for (int f = 0; f < ped.familyCount; f++)
               if (engine.Haplotype(ped, *(ped.families[f])))
                  sets.LoadFromMemory(&engine, 2);

            sets.quiet = true;
            sets.Filter(alleleCounts, 20);

            if (sets.discarded.Length())
               if (sets.discarded.Length() > discards)
                  discards = sets.discarded.Length();

            likelihood.Initialize(alleleCounts);
            likelihood.EM(&sets);
            likelihood.PairwiseLD(dprime, rsq);

            #ifdef __DEBUG__
               likelihood.Print(1e-5);
            #endif

            if (rsq >= threshold)
               {
               // When we identify a potential cluster, check the number of
               // obligate recombinants within the cluster:
               // if fewer than 5% of families have an obligate recombinant, proceed;
               // if more than 5% of families have an obligate recombinant, do not create this cluster

               engine.markers.Dimension(clusterEnd - clusterStart + 1);
               for (int i = 1; i <= clusterEnd - clusterStart; i++)
                  engine.markers[i] = markers[clusterStart + i];

               int useful_families = 0, recombinant_families = 0;
               for (int f = 0; f < ped.familyCount; f++)
                  try
                     {
                     if (engine.Haplotype(ped, *(ped.families[f])))
                        useful_families++;
                     else if (engine.obligate_recombinant)
                        recombinant_families++;
                     }
                  catch (const TreesTooBig & problem)
                     {
                     }

               engine.markers.Dimension(2);

               // To create a cluster we require that half of families should
               // contribute to allele frequency estimates *and* that no more
               // than 5% of families should have an obligate recombinant.
               if (useful_families > ped.familyCount / 2 &&
                   recombinant_families - 1 < ped.familyCount * 0.05)
                   break;
               }

            clusterEnd--;
            }

      if (clusterEnd > clusterStart)
         CreateCluster(markers, clusterStart, clusterEnd);

      clusterStart = clusterEnd;
      }
   } while (chromosome);

   if (discards)
      printf("  To speed up calculations for some marker pairs, up to\n"
             "  %d families were discarded during haplotype frequency\n"
             "  estimation\n", discards);

   printf("MARKER CLUSTERS: %d cluster%s where pairwise r-squared exceeds %.3f\n",
          count, count == 1 ? "" : "s", threshold);
   haveOutput = true;
   }


void MarkerClusters::AdjustPositions(IntArray & markers)
   {
   if (markers.Length() <= 0) return;

   // Calculate the midpoint for the cluster of markers
   MarkerInfo * first = Pedigree::GetMarkerInfo(markers[0]);
   MarkerInfo * last  = Pedigree::GetMarkerInfo(markers.Last());

   double clusterPosition  = (first->position + last->position) * 0.5;
   double clusterPosMale   = (first->positionMale + last->positionMale) * 0.5;
   double clusterPosFemale = (first->positionFemale + last->positionFemale) * 0.5;

   // Process each marker in the cluster
   for (int j = 0; j < markers.Length(); j++)
      {
      MarkerInfo * info = Pedigree::GetMarkerInfo(markers[j]);
      
      // Assign each marker to the cluster midpoint position
      info->position = clusterPosition;
      info->positionMale = clusterPosMale;
      info->positionFemale = clusterPosFemale;

      info->UpdateSerial();
      }
   }   

void MarkerClusters::CreateCluster(IntArray & markers, int clusterStart, int clusterEnd)
   {
   // Allocate a new cluster
   head = new MarkerClusterList(head);
   count++;
   
   // Calculate the midpoint for the cluster of markers
   MarkerInfo * first = Pedigree::GetMarkerInfo(markers[clusterStart]);
   MarkerInfo * last  = Pedigree::GetMarkerInfo(markers[clusterEnd]);

   double clusterPosition  = (first->position + last->position) * 0.5;
   double clusterPosMale   = (first->positionMale + last->positionMale) * 0.5;
   double clusterPosFemale = (first->positionFemale + last->positionFemale) * 0.5;

   // Process each marker in the cluster
   for (int j = clusterStart; j <= clusterEnd; j++)
      {
      int marker = markers[j];

      MarkerInfo * info = Pedigree::GetMarkerInfo(markers[j]);
      
      // Assign each marker to the cluster midpoint position
      info->position = clusterPosition;
      info->positionMale = clusterPosMale;
      info->positionFemale = clusterPosFemale;

      info->UpdateSerial();

      // Update the cluster
      head->markerIds.Push(marker);
      markerToCluster[marker] = head;
      }
   }

void MarkerClusters::FinishOutput()
   {
   if (haveOutput)
      printf("\n");
   haveOutput = false;
   }

void MarkerCluster::EstimateAlleleFrequencies(Pedigree & ped, const char * logname)
   {
   // Catch error and warning messages during frequency estimation
   String messages;

   // List alleles in the pedigree
   ped.LumpAlleles(0.0);

   // Variables used for pretty output
   bool estimated = false;
   int  line = 3;
   bool condensed = Pedigree::markerCount > 100;
   int  grain = Pedigree::markerCount / 50;

   // Loop through markers in the pedigree
   for (int i = 0; i < Pedigree::markerCount; i++)
      {
      if (!estimated)
         printf("Estimating allele frequencies... [using maximum likelihood]\n   "),
         estimated = true;

      // If there are no alleles, create a dummy allele
      MarkerInfo * info = ped.GetMarkerInfo(i);
      if (info->CountAlleles() == 0)
          info->NewAllele("1");

      // Initialize dummy allele frequencies
      int alleles = info->CountAlleles();
      info->freq.Dimension(alleles + 1);
      info->freq.Set(1.0 / (alleles + 1e-30));
      info->freq[0] = 0.0;

      // Create a dummy cluster
      MarkerCluster marker;
      marker.markerIds.Dimension(1);
      marker.markerIds[0] = i;
      marker.EstimateFrequencies(ped, &messages);
      marker.UpdateAlleleFrequencies();

      if (!condensed)
         {
         if ( line + Pedigree::markerNames[marker.markerIds[0]].Length() + 1 > 79)
            printf("\n   "), line = 3;

         printf("%s ", (const char *) Pedigree::markerNames[marker.markerIds[0]]);
         line += Pedigree::markerNames[marker.markerIds[0]].Length() + 1;
         }
      else
         if (marker.markerIds[0] % grain == 0)
            {
            printf(".");
            fflush(stdout);
            }
      }

   if (messages.Length())
      {
      printf("\n   Some families skipped due to computing limitations, see [%s]", logname);

      FILE * f = fopen(logname, "wt");

      if (f != NULL)
         {
         fprintf(f, "%s", (const char *) messages);
         fclose(f);
         }
      else
         printf("\n   Error opening log file [%s]!", logname);
      }

   if (estimated)
      printf(condensed ? "\nDone estimating frequencies for %d markers\n\n" : "\n\n",
             Pedigree::markerCount);
   }
 
