////////////////////////////////////////////////////////////////////// 
// regress/AutoFit.cpp 
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
 
#include "AutoFit.h"

#include "MathStats.h"
#include "MathNormal.h"
#include "Kinship.h"

#ifdef __CHROMOSOME_X__
#include "KinshipX.h"
#endif

#include <math.h>

bool AutoFit::useCovariates = false;
bool AutoFit::useProbands = false;

AutoFit::AutoFit() { }

void AutoFit::FitPolygenicModels(Pedigree & ped)
   {
   // Reset previous estimates
   means.Dimension(ped.traitCount); means.Set(_NAN_);
   variances.Dimension(ped.traitCount); variances.Set(_NAN_);
   heritabilities.Dimension(ped.traitCount); heritabilities.Set(_NAN_);

   // Array of informative individuals for each family
   IntArray * pheno = new IntArray[ped.familyCount];

   // Array of probands for each family
   IntArray probands(ped.familyCount);

   Kinship kin;

#ifdef __CHROMOSOME_X__
   KinshipX kinX;
#endif

   int probandStatus = useProbands ? ped.affectionNames.SlowFind("proband") : -1;

   // Print header line
   printf("FITTED MODELS%s\n",
          useCovariates ? " (for covariate adjusted residuals)" : "");
   printf("======================================================\n"
          "%15s %12s %12s %12s\n",
          "Trait", "Mean", "Variance", "Heritability");

   // Loop through traits in the pedigree
   for (int t = 0; t < ped.traitCount; t++)
      {
      // Reset proband list
      probands.Set(-1);

      // First list phenotyped individuals for each family
      // If we are using covariates, we only consider individuals
      // for which all covariates have been recorded
      for (int f = 0; f < ped.familyCount; f++)
         {
         pheno[f].Dimension(0);
         for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            if (ped[i].isPhenotyped(t) && (!useCovariates || ped[i].isFullyControlled()))
               {
               pheno[f].Push(i);

               if (!useProbands) continue;
               if (ped[i].pid.SlowCompare("proband") == 0 ||
                   probandStatus != -1 && ped[i].affections[probandStatus] == 2)
                   if (probands[f] == -1)
                      probands[f] = i;
                   else
                      error("There are multiple probands in family %s\n\n"
                            "Currently, only a simple correction for single\n"
                            "ascertainment is supported.\n",
                            (const char *) ped.families[f]->famid);
               }
         }

      // Number of families with non-null kinships and phenotypes
      int families = 0;

      // Count useful families
      for (int f = 0; f < ped.familyCount; f++)
         if (pheno[f].Length())
            families++;

      if (families == 0)
         {
         printf("Trait: %s -- No informative families\n\n",
                (const char *) ped.traitNames[t]);
         continue;
         }

      // Count probands
      int probandCount = probands.CountIfGreaterOrEqual(0);

      if (useProbands)
         printf("ASCERTAINMENT for trait %s\n"
                "    %d families ascertained through single proband\n"
                "    %d randomly selected families\n\n",
                (const char *) ped.traitNames[t],
                probandCount, families - probandCount);

      int useful = families + probandCount;

#ifdef __CHROMOSOME_X__
      int vc_count = 3;
#else
      int vc_count = 2;
#endif

      // The normal set class is our workhorse
      NormalSet mvn;
      mvn.Dimension(useful, vc_count);

      // Index to the normal equations within NormalSet
      int index = 0;

      // Build model with polygenes and environment only
      for (int f = 0; f < ped.familyCount; f++)
         if (pheno[f].Length())
            {
            int count = pheno[f].Length();

            // Phenotypes
            mvn[index].scores.Dimension(count);
            for (int i = 0; i < count; i++)
               mvn[index].scores[i] = ped[pheno[f][i]].traits[t];

            // Total number of fixed effects
            int fixedEffects = 1 + (useCovariates ? ped.covariateCount : 0);

            // Setup matrix of fixed effects for each family
            mvn[index].linearModel.Dimension(count, fixedEffects);
            for (int i = 0; i < count; i++)
               {
               // Constant for regressing grand mean
               mvn[index].linearModel[i][0] = 1.0;

               // User specified covariates
               if (useCovariates)
                  for (int j = 1; j <= ped.covariateCount; j++)
                     mvn[index].linearModel[i][j] = ped[pheno[f][i]].covariates[j-1];
               }

            // Setup non-shared variances
            mvn[index].varComponents[0].Dimension(count, count);
            mvn[index].varComponents[0].Identity();

            // Setup polygenic variances
            mvn[index].varComponents[1].Dimension(count, count);
            kin.Setup(*ped.families[f]);

            for (int i = 0; i < count; i++)
               for (int j = i; j < count; j++)
                  mvn[index].varComponents[1][i][j] =
                  mvn[index].varComponents[1][j][i] =
                     2.0 * kin(ped[pheno[f][i]], ped[pheno[f][j]]);

#ifdef __CHROMOSOME_X__
            // Setup polygenic variances for X chromosome
            mvn[index].varComponents[2].Dimension(count, count);
            kinX.Setup(*ped.families[f]);

            for (int i = 0; i < count; i++)
               for (int j = i; j < count; j++)
                  mvn[index].varComponents[2][i][j] =
                  mvn[index].varComponents[2][j][i] =
                     2.0 * kinX(ped[pheno[f][i]], ped[pheno[f][j]]);
#endif

            // Got the basic model set up!
            index++;

            if (probands[f] == -1) continue;

            Person & proband = ped[probands[f]];

            // The overall likelihood should be divided by this individuals'
            // likelihood
            mvn.operators[index] = NORMAL_DIV_LK;

            // Phenotype
            mvn[index].scores.Dimension(1);
            mvn[index].scores[0] = proband.traits[t];

            // Setup matrix of fixed effects, with grand mean and covariates
            mvn[index].linearModel.Dimension(1, fixedEffects);
            mvn[index].linearModel[0][0] = 1.0;
            if (useCovariates)
               for (int j = 1; j <= ped.covariateCount; j++)
                  mvn[index].linearModel[0][j] = proband.covariates[j-1];
            for (int j = ped.covariateCount + 1; j < fixedEffects; j++)
               mvn[index].linearModel[0][j] = 0.0;

            // Setup non-shared variances
            mvn[index].varComponents[0].Dimension(1, 1);
            mvn[index].varComponents[0].Identity();

            // Setup polygenic variances
            mvn[index].varComponents[1].Dimension(1, 1);
            mvn[index].varComponents[1][0][0] =
               2.0 * kin(proband, proband);

#ifdef __CHROMOSOME_X__
            // Setup polygenic variances for X chromosome
            mvn[index].varComponents[2].Dimension(1, 1);
            mvn[index].varComponents[2][0][0] =
               2.0 * kinX(proband, proband);
#endif
            index++;
            }

      // Fit polygenic model
      mvn.Solve();

      double mean = useCovariates ? 0.0 : mvn.means[0];
      double variance = mvn.variances.Sum();
      double heritability = mvn.variances[vc_count - 1] / variance;

      printf("%15s %12.5f %12.5f %12.5f\n",
             (const char *) ped.traitNames[t], mean, variance, heritability);

      means[t] = mean;
      variances[t] = variance;
      heritabilities[t] = heritability;

      if (!useCovariates) continue;

      // Evaluate model to update residuals for each family
      mvn.Evaluate();

      index = 0;

      // First wipe out all the phenotypes, to ensure that we
      // don't get residuals and covariates mixed up ...
      for (int i = 0; i < ped.count; i++)
         ped[i].traits[t] = _NAN_;

      // Replace phenotypes with residuals
      for (int f = 0; f < ped.familyCount; f++)
         if (pheno[f].Length())
            {
            int count = pheno[f].Length();

            // Replace phenotypes with residuals
            for (int i = 0; i < count; i++)
               ped[pheno[f][i]].traits[t] = mvn[index].residuals[i];

            // Proceed to next family
            index++;

            if (probands[f] == -1) continue;

            // Skip over probands
            index++;
            }
      }

   delete [] pheno;
   }


 
