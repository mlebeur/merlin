////////////////////////////////////////////////////////////////////// 
// regress/Regress.cpp 
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
 
#include "TraitTransformations.h"
#include "RegressParameters.h"
#include "MerlinSimulator.h"
#include "RegressAnalysis.h"
#include "MerlinCluster.h"
#include "MerlinError.h"
#include "MerlinSort.h"
#include "Manners.h"
#include "Error.h"

#include "LongInt.h"

void ShowBanner();

int main(int argc, char * argv[])
   {
   ShowBanner();

   SetupCrashHandlers();

   RegressionParameters pl;

   pl.Read(argc, argv);
   pl.Status();
   pl.Check();

   Pedigree ped;

   ped.Prepare(pl.dataFile);
   ped.Load(pl.pedigreeFile);
   SortFamilies(ped);
   ped.LoadMarkerMap(pl.mapFile, true /* filter out markers not in pedigree */);
   ped.LoadAlleleFrequencies(pl.freqFile);

   // Estimate allele frequencies
   if (pl.alleleFrequencies == FREQ_ML)
      MarkerCluster::EstimateAlleleFrequencies(ped, MerlinCore::filePrefix + "-freqs.log");
   else
      ped.EstimateFrequencies(pl.alleleFrequencies),
      ped.LumpAlleles(0.0);

   if (ped.traitCount < 1)
      error("The data set includes no quantitative trait data\n\n");

   if (pl.trimPedigree)
      ped.Trim();

   if (!MerlinCore::twopoint)
      {
      clusters.LoadFromFile(pl.clusterFile, MerlinCore::filePrefix + "-clusters.log");
      clusters.ClusterByDistance(pl.clusterDistance * 0.01);
      clusters.ClusterByRsquared(ped, pl.clusterRsquared);
      clusters.EstimateFrequencies(ped, MerlinCore::filePrefix + "-clusters-freq.log");
      clusters.FinishOutput();
      }

   if (pl.writeFrequencyFile)
      ped.WriteFreqFile(MerlinCore::filePrefix + ".freq");

   if (pl.inverseNormal)
      InverseNormalTransform(ped);

   if (pl.fitErrorRate)
       {
       ErrorRateEstimator errorEstimator(ped);
       errorEstimator.Estimate();
       }

   int  next_chromosome  = -1;
   bool many_chromosomes = false;

   int run = 0;
   do {
      if (pl.reruns > 1)
         MerlinCore::filePrefix.printf("merlin-%08d-%05d", pl.random_seed, run + 1);

      RegressionAnalysis engine(ped);

      engine.SetupGlobals();

      if (ped.markerCount < 1)
         error("The data set includes no genetic markers\n\n");

      do {
         next_chromosome = engine.SetupMap(next_chromosome);
         many_chromosomes = many_chromosomes || (next_chromosome != 0);

         if (many_chromosomes)
            printf("\nAnalysing Chromosome %d\n\n",
                   ped.GetMarkerInfo(engine.markers[0])->chromosome);

         if (pl.simulateNull)
            Simulator::Simulate(ped, engine.markers);

         for (int i = 0; i < ped.familyCount; i++)
            if (engine.SelectFamily(ped.families[i]))
               engine.Analyse();

         engine.PrintScores();
      } while (next_chromosome);

      engine.CleanupGlobals();

   } while (++run < pl.reruns);
        

   if (pl.writeFrequencyFile)
      {
      ped.WriteFreqFile("merlin.freq");
      printf("Allele frequencies written to file [merlin.freq]\n\n");
      }

   printf("\n");
   }

void ShowBanner()
   {
#ifndef VERSION
   printf("MERLIN - (c) 2000-2003 Goncalo Abecasis");
#else
   printf("MERLIN " VERSION " - (c) 2000-2003 Goncalo Abecasis");
#endif

#ifdef __HUGE_SWAP__
   printf("\nModifications: ");
#else
#ifdef __USE_LONG_INT
   printf("\nModifications: ");
#endif
#endif

#ifdef __HUGE_SWAP__
   printf("HUGE_SWAP ");
#endif

#ifdef __USE_LONG_INT
   printf("64-ALLELES ");
#endif

   printf("\n");
   }


 
