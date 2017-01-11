////////////////////////////////////////////////////////////////////// 
// merlin/MerlinSimwalk2.cpp 
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
 
#include "MerlinFamily.h"
#include "MerlinSimwalk2.h"

#include <math.h>

void NPLforSimwalk2::OpenFiles()
   {
   String filename(MerlinCore::filePrefix);
   filename += ".npl";

   output = fopen(filename, "wt");

   if (output == NULL)
      error("Error opening file [%s] where intermediate output for "
            "Simwalk2 would be stored", (const char *) filename);

   fprintf(output, "NPLXF version 1.0\n");
   }

void NPLforSimwalk2::CloseFiles()
   {
   fclose(output);
   printf("Intermediate output for Simwalk2 stored in file [%s.npl]\n", (const char *) MerlinCore::filePrefix);
   }

void NPLforSimwalk2::Output()
   {
   KongAndCox & kac = *family.kac;

   if ((kac.scales[0] == 0.0) && (kac.scales[1] == 0.0))
      {
      Skip();
      return;
      }

   // Header line has pedigree name (cols 1-8),
   // and no. of markers (cols. 11 - 15)
   fprintf(output, "%8.8s%7d\n",
          (const char *) family.family->famid, family.markerCount);

   // Headings for the table to follow
   fprintf(output, "%-8s  %9s %9s %9s %9s %9s %9s\n",
          "Marker",
          "Obs_Pair", "Mean_Pair", "SDev_Pair",
          "Obs_All", "Mean_All", "SDev_All");

   // Mean and variance for each statistic
   double pairsMean = kac.means[1];
   double pairsDev = kac.scales[1] ? 1.0/kac.scales[1] : 0.0;

   double allMean = kac.means[0];
   double allDev = kac.scales[0] ? 1.0/kac.scales[0] : 0.0;

   // When there are multiple markers at the same location, MERLIN
   // only calculates one score, but SIMWALK2 expects one score per
   // marker.
   //
   // We iterate over all markers i and output scores for the same
   // analysis position multiple times if required
   //
   int pos = 0;
   int serial = family.serial;
   Vector ** scores = kac.lm.scores;

   for (int i = 0; i < family.markers.Length(); i++)
      {
      // Find position corresponding to current marker
      while (family.analysisPositions[pos] < family.markerPositions[i])
         pos++;

      // Output scores
      fprintf(output, "%-8s ", (const char *) family.ped.markerNames[family.markers[i]]);
      OutputValue(scores[1][serial][pos] * pairsDev + pairsMean);
      OutputValue(pairsMean);
      OutputValue(pairsDev);
      OutputValue(scores[0][serial][pos] * allDev + allMean);
      OutputValue(allMean);
      OutputValue(allDev);
      fprintf(output, "\n");
      }
   }

void NPLforSimwalk2::Skip()
   {
   // Header line has pedigree name (cols 1-8) and zero to denote skipping
   fprintf(output, "%8.8s%7d\n", (const char *) family.family->famid, 0);
   }


void NPLforSimwalk2::OutputValue(double value)
   {
   if (fabs(value) < 1e-10)
      fprintf(output, "         0");
   else
      fprintf(output, "% 10.*G", (value > 1000000 || value < 1.0) ? 4 : 7, value);
   }

 
 
