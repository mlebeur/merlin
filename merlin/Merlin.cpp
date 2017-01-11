////////////////////////////////////////////////////////////////////// 
// merlin/Merlin.cpp 
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
 
#include "MerlinParameters.h"
#include "TraitTransformations.h"
#include "MerlinSimulator.h"
#include "MerlinCluster.h"
#include "MerlinFamily.h"
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

   MerlinParameters pl;

   pl.Read(argc, argv);
   pl.Status();

#ifdef __CHROMOSOME_X__
   PedigreeGlobals::chromosomeX = true;
#endif

   Pedigree ped;

   ped.Prepare(pl.dataFile);
   ped.Load(pl.pedigreeFile);
   SortFamilies(ped);
   ped.LoadMarkerMap(pl.mapFile, true /* filter out markers not in pedigree */);
   ped.LoadAlleleFrequencies(pl.freqFile);

   if (pl.inverseNormal)
      InverseNormalTransform(ped);

   // Estimate allele frequencies
   if (pl.alleleFrequencies == FREQ_ML)
      MarkerCluster::EstimateAlleleFrequencies(ped, MerlinCore::filePrefix + "-freqs.log");
   else
      {
      ped.EstimateFrequencies(pl.alleleFrequencies);

      // Recode alleles so more frequent alleles have lower allele numbers internally
      ped.LumpAlleles(0.0);
      }

   if (pl.trimPedigree)
      ped.Trim();

   if (!MerlinCore::twopoint)
      {
      clusters.LoadFromFile(pl.clusterFile, MerlinCore::filePrefix + "-clusters.log");
      clusters.ClusterByDistance(pl.clusterDistance * 0.01);
      clusters.ClusterByRsquared(ped, pl.clusterRsquared);
      clusters.EstimateFrequencies(ped, MerlinCore::filePrefix + "-cluster-freqs.log");
      clusters.FinishOutput();

      if (pl.writeClusters && clusters.Enabled())
         {
         int count = clusters.SaveToFile(MerlinCore::filePrefix + "-clusters.freq");
         printf("Haplotypes for %d clusters written to file [%s-clusters.freq]\n\n", count, (const char *) MerlinCore::filePrefix);
         }
      }

   VarianceComponents::customModels.LoadFromFile(pl.qtlModelFile);

   if (ped.markerCount < 1)
      error("The data set includes no genetic markers\n\n");

   if (ped.sexSpecificMap || KongAndCox::nplExtras)
      Mantra::ignoreCoupleSymmetries = true;

   int  next_chromosome  = -1;
   bool many_chromosomes = false;

   if (pl.fitErrorRate)
      {
      ErrorRateEstimator errorEstimator(ped);
      errorEstimator.Estimate();
      }

   int run = 0;
   do {
      if (pl.reruns > 1)
         MerlinCore::filePrefix.printf("merlin-%08d-%05d", pl.random_seed, run + 1);

      FamilyAnalysis engine(ped);

      engine.SetupGlobals();
      engine.SetupFiles();

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

         printf("\n");
         engine.ShowLikelihood();
         engine.ShowLODs();
      } while (next_chromosome);

      if (pl.saveReplicate)
         {
         printf("Saving simulated pedigree with prefix '%s-replicate' ...\n\n",
                (const char *) MerlinCore::filePrefix );
                
         ped.WriteDataFile((const char *) (MerlinCore::filePrefix + "-replicate.dat"));
         ped.WritePedigreeFile((const char *) (MerlinCore::filePrefix + "-replicate.ped"));
         ped.WriteMapFile((const char *) (MerlinCore::filePrefix + "-replicate.map"));
         ped.WriteFreqFile((const char *) (MerlinCore::filePrefix + "-replicate.freq"));
         }

      printf("\n");

      engine.CloseFiles();
      engine.CleanupGlobals();
   } while (++run < pl.reruns);

   if (pl.writeFrequencyFile)
      {
      ped.WriteFreqFile((const char *) (MerlinCore::filePrefix + ".freq"));
      printf("Allele frequencies written to file [%s.freq]\n", (const char *) MerlinCore::filePrefix);
      }

   printf("\n");
   }

void ShowBanner()
   {
#ifndef VERSION
   printf("MERLIN - (c) 2000-2007 Goncalo Abecasis");
#else
   printf("MERLIN " VERSION " - (c) 2000-2007 Goncalo Abecasis");
#endif

#ifdef __HUGE_SWAP__
   printf("\nModifications: ");
#else
#ifdef __USE_LONG_INT
   printf("\nModifications: ");
#else
#ifdef __CHROMOSOME_X__
   printf("\nModifications: ");
#endif
#endif
#endif

#ifdef __HUGE_SWAP__
   printf("HUGE_SWAP ");
#endif

#ifdef __USE_LONG_INT
   printf("64-ALLELES ");
#endif

#ifdef __CHROMOSOME_X__
   printf("CHROMOSOME-X ");                       
#endif

   printf("\n\n");

#ifdef __VERSION__
   printf("References for this version of Merlin:\n\n"
          "   Abecasis et al (2002) Nat Gen 30:97-101        [original citation]\n"
          "   Fingerlin et al (2004) AJHG 74:432-43          [case selection for association studies]\n"
          "   Abecasis and Wigginton (2005) AJHG 77:754-67   [ld modeling, parametric analyses]\n"
          "   Fingerlin et al (2006) Gen Epidemiol 30:384-96 [sex-specific maps]\n"
          "   Chen and Abecasis (2007) AJHG 81:913-26        [qtl association analysis, qtl simulation]\n\n");
#endif
   }



 
