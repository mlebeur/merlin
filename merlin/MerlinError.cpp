////////////////////////////////////////////////////////////////////// 
// merlin/MerlinError.cpp 
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
 
#include "MerlinError.h"
#include "MerlinCore.h"
#include "MathGold.h"
#include "Houdini.h"

// Estimates error rates using maximum likelihood
//

void ErrorRateEstimator::Estimate()
   {
   // First we try a genotyping a model with genotyping error
   //
   printf("Fitting Maximum Likelihood Estimate for Error Rate...\n");

   alleleError = true;
   FuzzyInheritanceTree::perGenotypeError = 0.0;

   EstimateModel();

   double perAlleleEstimate = min;
   double perAlleleLk = fmin;

   alleleError = false;
   FuzzyInheritanceTree::perAlleleError = 0.0;

   EstimateModel();

   if (perAlleleLk < fmin)
      {
      printf("The per allele error model fits best.\n"
              "    Estimated error rate = %.4f per allele\n\n",
              perAlleleEstimate);
      FuzzyInheritanceTree::perGenotypeError = 0.0;
      FuzzyInheritanceTree::perAlleleError = perAlleleEstimate;
      }
   else
      {
      printf("The per genotype error model fits best.\n"
              "    Estimated error rate = %.4f per genotype\n\n",
              min);
      FuzzyInheritanceTree::perGenotypeError = min;
      FuzzyInheritanceTree::perAlleleError = 0.0;
      }
   }

void ErrorRateEstimator::EstimateModel()
   {
   trace = false;

   printf("   Assuming per %s error model\n", alleleError ? "allele" : "genotype");

   a = 0.00001;
   fa = f(a);

   printf("      lnLikelihood at lower bound (%.3g) is %.3f\n", a, -fa);

   c = 0.10;
   fc = f(c);

   printf("      lnLikelihood at upper bound (%.3g) is %.3f\n", c, -fc);

   b = (c - a) * 0.4;
   fb = f(b);

   trace = !MerlinCore::quietOutput;
   Brent(0.00001);

   printf("      lnLikelihood at MLE (%.3g) is %.3f           \n", min, -fmin);
   }

// Calculates likelihood of data for some arbitrary error rate
//

double ErrorRateEstimator::f(double error_rate)
   {
   if (alleleError)
      FuzzyInheritanceTree::perAlleleError = error_rate;
   else
      FuzzyInheritanceTree::perGenotypeError = error_rate;

   if (trace)
      {
      printf("      lnLikelihood at best guess (%.3g) is %.3f     \r", min, -fmin);
      fflush(stdout);
      }

   return -CalculateLikelihood();
   }

double ErrorRateEstimator::CalculateLikelihood()
   {
   double lnLikelihood = 0.0;

   int  next_chromosome  = 0;
   bool many_chromosomes = false;

   MerlinCore engine(ped);

   engine.SetupGlobals();

   bool quietOutput = engine.quietOutput;
   engine.quietOutput  = true;

   do {
      next_chromosome = engine.SetupMap(next_chromosome);

      many_chromosomes = many_chromosomes || (next_chromosome != 0);

      if (many_chromosomes && !quietOutput)
         {
         printf("C%02d\r", ped.GetMarkerInfo(engine.markers[0])->chromosome);
         fflush(stdout);
         }
         
      for (int i = 0; i < ped.familyCount; i++)
         if (engine.SelectFamily(ped.families[i], false))
            lnLikelihood += engine.CalculateLikelihood();

   } while (next_chromosome);

   if (many_chromosomes && !quietOutput)
      {
      printf("     \r");
      fflush(stdout);
      }

   engine.quietOutput = quietOutput;

   return lnLikelihood;
   } 
