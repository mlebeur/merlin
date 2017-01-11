////////////////////////////////////////////////////////////////////// 
// merlin/VarianceComponents.cpp 
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
 
#include "VarianceComponents.h"
#include "MerlinFamily.h"

#include "MathStats.h"
#include "MathNormal.h"
#include "Kinship.h"

#ifdef __CHROMOSOME_X__
#include "KinshipX.h"
#endif

#include <math.h>
#include <ctype.h>

bool VarianceComponents::useCovariates = false;
bool VarianceComponents::useProbands = false;
double VarianceComponents::unlinkedFraction = 0.0;

QtlModel VarianceComponents::customModels;

VarianceComponents::VarianceComponents()
   {
   perFamily = NULL;
   }

void VarianceComponents::Analyse(FamilyAnalysis & engine)
   {
   Pedigree & ped = engine.ped;
   StringArray & labels = engine.labels;
   MerlinKinship & markerKin = engine.kinship;

   int positions = labels.Length();

   // Array of informative individuals for each family
   IntArray * pheno = new IntArray[ped.familyCount];

   // Array of probands for each family
   IntArray probands(ped.familyCount);

   Kinship kin;

#ifdef __CHROMOSOME_X__
   KinshipX kinX;
#endif

   // Summarize covariates ...
   if (customModels.HaveModels())
      customModels.PrintSummary();
   else
      {
      ListCovariates();
      PrintCovariates();
      }

   int chr = 0;
   String tablename;
   FILE * tablefile = NULL;

   if (MerlinCore::tabulate)
      {
      chr = positions > 0 ? Pedigree::GetMarkerInfo(engine.markers[0])->chromosome : 0;
      tablename.printf("%s-vc-chr%02d.tbl", (const char *) MerlinCore::filePrefix, chr > 0 ? chr : 0);
      tablefile = fopen(tablename, "wt");
      }

   if (tablefile != NULL)
      fprintf(tablefile, "CHR\tPOS\tLABEL\tTRAIT\tH2\tLOD\tPVALUE\n");

   int probandStatus = useProbands ? ped.affectionNames.SlowFind("proband") : -1;

   // Index for repeated measures
   IntArray traitIndex;

   // Label for analysis of this trait
   String   traitLabel;

   // Figure out how many models we have to analyse
   int models = customModels.HaveModels() ? customModels.modelCount : ped.traitCount;

   // Loop through traits in the pedigree
   for (int m = 0; m < models; m++)
      {
      // Retrieve the appropriate trait
      int t = customModels.HaveModels() ? customModels.traitList[m] : m;

      // Retrieve the appropriate set of covariates
      if (customModels.HaveModels()) covariates = customModels.covariateList[m];

      // Check if this trait is repeated
      if (SkipTrait(t)) continue;

      // Default trait label
      traitLabel = ped.traitNames[t];

      // Check for repeat measurements and adjust trait label, if appropriate
      bool repeatedMeasures = IndexMeasurements(traitIndex, traitLabel, t);
      int averagedMeasures = ped.LookupCovariate(ped.traitNames[t] + "_repeats");

      // Reset proband list
      probands.Set(-1);

      // First list phenotyped individuals for each family
      // If we are using covariates, we only consider individuals
      // for which all covariates have been recorded
      for (int f = 0; f < ped.familyCount; f++)
         {
         pheno[f].Dimension(0);
         if (markerKin.CheckFamily(f))
            for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
               if ( IsPhenotyped(ped[i], traitIndex, t) && CheckCovariates(ped[i]) &&
                   (averagedMeasures == -1 || ped[i].isControlled(averagedMeasures)))
                  {
                  if (!repeatedMeasures)
                     pheno[f].Push(i);
                  else for (int j = 0; j < traitIndex.Length(); j++)
                     if (ped[i].isPhenotyped(traitIndex[j]))
                        pheno[f].Push(i);

                  if (!useProbands) continue;
                  if (ped[i].pid.SlowCompare("proband") == 0 ||
                      probandStatus != -1 && ped[i].affections[probandStatus] == 2)
                      if (probands[f] == -1)
                         probands[f] = i;
                      else
                         error("Found multiple probands in family %s\n"
                               "preparing variance component analysis for trait %s\n\n"
                               "Currently, only a simple correction for single\n"
                               "ascertainment is supported.\n",
                               (const char *) ped.families[f]->famid,
                               (const char *) traitLabel);
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
         printf("Trait: %s (No informative families)\n\n",
                (const char *) traitLabel);
         continue;
         }

      // Count probands
      int probandCount = probands.CountIfGreaterOrEqual(0);

      if (useProbands)
         printf("ASCERTAINMENT for trait %s\n"
                "    %d families ascertained through single proband\n"
                "    %d randomly selected families\n\n",
                (const char *) traitLabel,
                probandCount, families - probandCount);

      bool fitMeasurementError = repeatedMeasures ||
                            CheckAveragedMeasures(ped, pheno, averagedMeasures);

      if (repeatedMeasures) SummarizeRepeatedMeasures(ped, traitLabel, pheno);

      bool heterogeneity = unlinkedFraction > 0.0;

      // Setup PDF file
      if (engine.writePDF)
         SetupPDF(engine.pdf, traitLabel, probandCount);

      int useful = families + probandCount;

      if (heterogeneity)
         {
         useful *= 2;

         printf("HETEROGENEITY: Model allows for %.2f%% unlinked families\n\n",
                unlinkedFraction * 100.0);
         }

#ifdef __CHROMOSOME_X__
      int vc_count = 3 + fitMeasurementError;
#else
      int vc_count = 2 + fitMeasurementError;
#endif
      int vc_count2 = heterogeneity ? vc_count * 2 : vc_count;

      // The normal set class is our workhorse
      NormalSet mvn;
      mvn.Dimension(useful, vc_count2);

      // Per family likelihoods under the null and alternative hypothesis
      Vector nullPerFamily;

      // Index to the normal equations within NormalSet
      int index = 0;

      // Build model with polygenes and environment only
      for (int f = 0; f < ped.familyCount; f++)
         if (pheno[f].Length())
            {
            // Track per family contributions to likelihood, if requested
            if (index && (perFamily != NULL))
               mvn.operators[index - 1] |= NORMAL_LAST_OP(NORMAL_RECORD_LLK);

            int count = pheno[f].Length();

            // Use these variable to check whether multiple measurements
            // are extracted from the same individual
            int last_individual = -1;
            int measurement = 0;

            // Phenotypes
            mvn[index].scores.Dimension(count);
            for (int i = 0; i < count; i++)
               {
               if (pheno[f][i] != last_individual) measurement = 0;
               mvn[index].scores[i] = RetrieveMeasurement(ped[pheno[f][i]], traitIndex, measurement);
               last_individual = pheno[f][i];
               }

            // Total number of fixed effects
            int fixedEffects = 1 + covariates.Length();

            if (heterogeneity) fixedEffects *= 2;

            // Setup matrix of fixed effects for each family
            mvn[index].linearModel.Dimension(count, fixedEffects);
            for (int i = 0; i < count; i++)
               {
               // Constant for regressing grand mean
               mvn[index].linearModel[i][0] = 1.0;

               // User specified covariates
               for (int j = 1; j <= covariates.Length(); j++)
                  mvn[index].linearModel[i][j] =
                     ped[pheno[f][i]].covariates[covariates[j-1]];

               // Additional parameters in models with heterogeneity
               for (int j = covariates.Length() + 1; j < fixedEffects; j++)
                  mvn[index].linearModel[i][j] = 0.0;
               }

            // Setup non-shared variances
            mvn[index].varComponents[0].Dimension(count, count);
            if (!repeatedMeasures)
               mvn[index].varComponents[0].Identity();
            else
               for (int i = 0; i < count; i++)
                  for (int j = 0; j < count; j++)
                     mvn[index].varComponents[0][i][j] =
                     mvn[index].varComponents[0][j][i] =
                        (pheno[f][i] == pheno[f][j]) ? 1.0 : 0.0;

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

            // Measurement error component, using averaged measurements
            if (fitMeasurementError && !repeatedMeasures)
               {
               mvn[index].varComponents[vc_count - 1].Dimension(count, count);
               mvn[index].varComponents[vc_count - 1].Zero();

               for (int i = 0; i < count; i++)
                  mvn[index].varComponents[vc_count - 1][i][i] =
                     1.0 / ped[pheno[f][i]].covariates[averagedMeasures];
               }

            // Measurement error component, using actual repeat measurements
            if (fitMeasurementError && repeatedMeasures)
               {
               // The design matrix for the measurement error component is
               // simply the identify matrix ...
               mvn[index].varComponents[vc_count - 1].Dimension(count, count);
               mvn[index].varComponents[vc_count - 1].Identity();
               }

            // Additional variance covariance matrices set to zero
            for (int j = vc_count; j < vc_count2; j++)
               mvn[index].varComponents[j].Dimension(0, 0);

            // Got the basic model set up!
            index++;

            if (heterogeneity && probands[f] == -1)
               {
               // Construct the second likelihood term
               DuplicateModel(mvn[index], mvn[index - 1], vc_count, fixedEffects/2);

               // L(A) = L(A|unlinked)P(unlinked) + L(A|linked)P(linked)

               // To calculate the overall likelihood,
               // scale the first likelihood ...
               mvn.operators[index-1] = NORMAL_SCALE_LLK;
               mvn.weights[index-1] = log(unlinkedFraction);

               // Then scale the second likelihood and sum the two alternative
               // likelihoods before proceeding to the rest of the sample...
               mvn.operators[index] =
                  NORMAL_OP(NORMAL_SCALE_LLK, NORMAL_SUM_LK, NORMAL_MUL_LK, 0);
               mvn.weights[index] = log(1.0 - unlinkedFraction);

               index++;
               }

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
            for (int j = 1; j <= covariates.Length(); j++)
               mvn[index].linearModel[0][j] =
                  proband.covariates[covariates[j-1]];
            for (int j = covariates.Length() + 1; j < fixedEffects; j++)
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

            // Setup measurement error variance
            mvn[index].varComponents[vc_count - 1].Dimension(1, 1);
            mvn[index].varComponents[vc_count - 1][0][0] =
                  1.0 / proband.covariates[averagedMeasures];

            index++;

            if (!heterogeneity) continue;

            // L = L(A|unlinked)/L(P|unlinked)*P(unlinked) +
            //     L(A|linked)/L(P|linked)*P(unlinked)

            mvn.operators[index-2] = NORMAL_NOP;
            mvn.operators[index-1] = NORMAL_OP(NORMAL_DIV_LK, NORMAL_SCALE_LLK, 0, 0);
            mvn.weights[index-1] = log(unlinkedFraction);

            // Construct the second likelihood term
            DuplicateModel(mvn[index], mvn[index - 2], vc_count, fixedEffects/2);
            mvn.operators[index] = NORMAL_NOP;

            index++;

            DuplicateModel(mvn[index], mvn[index - 2], vc_count, fixedEffects/2);
            mvn.operators[index] =
               NORMAL_OP(NORMAL_DIV_LK, NORMAL_SCALE_LLK, NORMAL_SUM_LK, NORMAL_MUL_LK);
            mvn.weights[index] = log(1.0 - unlinkedFraction);

            index++;
            }

      // Record the last log-likelihood of the bunch
      if (index && (perFamily != NULL))
         mvn.operators[index - 1] |= NORMAL_LAST_OP(NORMAL_RECORD_LLK);

      // Fit polygenic model
      mvn.Solve();

      // Track Key Values
      double lkNull = mvn.Evaluate();
      double sampleVar = TotalVariance(mvn.variances, vc_count);
      double sampleH2  = mvn.variances[heterogeneity ? vc_count + 1 : 1];
#ifdef __CHROMOSOME_X__
      double sampleH2X = mvn.variances[heterogeneity ? vc_count + 2 : 2];
#endif

      String errorString;
      if (fitMeasurementError)
         errorString.printf(", error = %.2f%%",
            mvn.variances[vc_count2 - 1] * 100.0 / sampleVar);

      // Record per family likelihoods
      nullPerFamily = mvn.recordedLikelihoods;

      // Print out header for results
#ifdef __CHROMOSOME_X__
      printf("Phenotype: %s [VC] (%d famil%s, h2 = %.2f%% [AUTO] + %.2f%% [X]%s)\n"
             "====================================================================================\n",
            (const char *) traitLabel, families, families == 1 ? "y" : "ies",
            sampleH2 * 100. / sampleVar, sampleH2X * 100. / sampleVar,
            (const char *) errorString);
#else
      printf("Phenotype: %s [VC] (%d famil%s, h2 = %.2f%%%s)\n"
             "==================================================================\n",
            (const char *) traitLabel, families, families == 1 ? "y" : "ies",
            sampleH2 * 100. / sampleVar,
            (const char *) errorString);
#endif
      printf("%20s %8s %7s %7s %7s\n",
             "Position", "H2 ", "ChiSq", "LOD", "pvalue");

      // Add an additional variance component for linked major gene
      mvn.Dimension(useful, vc_count2 + 1);

      // Loop through analysis at individual marker locations
      for (int pos = 0; pos < positions; pos++)
         {
         index = 0;

         // Fill major gene matrix with observed coefficients
         for (int f = 0; f < ped.familyCount; f++)
            if (pheno[f].Length())
               {
               // Likelihood for unlinked portion of sample does not change
               if (heterogeneity)
                  mvn[index++].varComponents[vc_count2].Dimension(0, 0);

               if (heterogeneity && probands[f] != -1)
                  mvn[index++].varComponents[vc_count2].Dimension(0, 0);

               int count = mvn[index].linearModel.rows;

               mvn[index].varComponents[vc_count2].Dimension(count, count);
               for (int i = 0; i < count; i++)
                  for (int j = i; j < count; j++)
                     mvn[index].varComponents[vc_count2][i][j] =
                     mvn[index].varComponents[vc_count2][j][i] =
                                  2.0 * markerKin.Retrieve(f, pos,
                                        ped[pheno[f][i]].traverse,
                                        ped[pheno[f][j]].traverse);

               index++;

               if (probands[f] == -1) continue;

               // Setup kinship coefficient at marker for proband
               mvn[index].varComponents[vc_count2].Dimension(1, 1);
               mvn[index].varComponents[vc_count2][0][0] =
                     2.0 * markerKin.Retrieve(f, pos,
                           ped[probands[f]].traverse,
                           ped[probands[f]].traverse);

               index++;
               }

         // Fit the alternative model
         mvn.Solve();

         // Evaluate the likelihood
         double lk = mvn.Evaluate();

         // Get key values
         double chisq = lkNull > lk ? (lkNull - lk) * 2.0 : 0.0;
         double pvalue = chidist(chisq, 1) * 0.5;
         double h2 = mvn.variances[vc_count2];
         double var = TotalVariance(mvn.variances, vc_count + 1);

         int digits = pvalue < 1.5e-4 ? 5 : int(log(pvalue*0.06666666)*-0.4343);

         double lod = chisq / (2*log(10.0));

         // Print out a one-line summary of results
         printf("%20.20s %7.2f%% %7.2f %7.2f %7.*f\n",
                (const char *) labels[pos],
                h2 * 100 / var, chisq, lod, digits, pvalue);

         if (tablefile != NULL)
            fprintf(tablefile, "%d\t%.3f\t%s\t%s\t%.3f\t%.3f\t%.4g\n",
                     chr, engine.analysisPositions[pos] * 100., (const char *) labels[pos],
                     (const char *) traitLabel, h2 * 100. / var, lod, pvalue);

         if (perFamily != NULL)
            WritePerFamilyLOD(ped, pheno, (const char *) labels[pos],
                              nullPerFamily, mvn.recordedLikelihoods);

         if (engine.writePDF)
            engine.pdf.y[0][pos] = lod;
         }

     if (engine.writePDF)
        engine.pdf.DrawChart();

      printf("\n");
      }

   if (tablefile != 0)
      {
      fclose(tablefile);
      printf("Variance component analysis tabulated to [%s]\n\n", (const char *) tablename);
      }

   delete [] pheno;
   }

void VarianceComponents::ListCovariates()
   {
   const char suffix[] = "_repeats";

   covariates.Clear();

   if (!useCovariates)
      return;

   for (int i = 0; i < Pedigree::covariateCount; i++)
      {
      String & label = Pedigree::covariateNames[i];

      if (label.Length() <= 8)
         covariates.Push(i);
      else
         for (int j = 7; j >= 0; j--)
            if (tolower(label[label.Length() - 8 + j]) != suffix[j])
               {
               covariates.Push(i);
               break;
               }
      }
   }

bool VarianceComponents::CheckCovariates(Person & p)
   {
   for (int i = 0; i < covariates.Length(); i++)
      if (p.covariates[covariates[i]] == _NAN_)
         return false;
   return true;
   }

void VarianceComponents::PrintCovariates()
   {
   if (covariates.Length() == 0)
      return;

   printf("COVARIATES: Variance components model includes the following covariate%s:\n",
          covariates.Length() == 1 ? "" : "s" );

   for (int i = 0, column = 0; i < covariates.Length(); i++)
      {
      String & label = Pedigree::covariateNames[covariates[i]];

      if (column + label.Length() > 78)
         printf("\n", column = 0);

      if (column == 0)
         printf("   ", column = 3);

      printf("%s%s ", (const char *) label,
             (i + 1) != covariates.Length() ? "," : "");

      column += label.Length() + 2;
      }

   printf("\n\n");
   }

void VarianceComponents::DuplicateModel(NormalEquations & dest, NormalEquations & src,
                                        int vc_count, int beta_count)
   {
   // With heterogeneity, the two portions of the sample have
   // similar models, but possibly different parameter values.
   // We achieve this by duplicating variance components and
   // fixed effects ...
   //
   for (int i = 0, j = vc_count; i < vc_count; i++, j++)
      {
      dest.varComponents[i] = src.varComponents[j];
      dest.varComponents[j] = src.varComponents[i];
      }

   int count = src.scores.Length();

   dest.linearModel.Dimension(count, beta_count * 2);
   for (int i = 0; i < count; i++)
      for (int j = 0, k = beta_count; j < beta_count; j++, k++)
         {
         dest.linearModel[i][j] = src.linearModel[i][k];
         dest.linearModel[i][k] = src.linearModel[i][j];
         }

   dest.scores = src.scores;
   }

double VarianceComponents::TotalVariance(Vector & components, int vc_count)
   {
   // When modelling heterogeneity, the first set of variance components
   // refers to unlinked sample. So we total the last vc_count components
   double var = 0.0;

   for (int i = 0, j = components.Length() - 1; i < vc_count; i++, j--)
      var += components[j];

   return var;
   }

void VarianceComponents::OpenPerFamilyFile()
   {
   String filename(MerlinCore::filePrefix);
   filename += ".vc";

   perFamily = fopen(filename, "wt");

   if (perFamily == NULL)
      error("Opening file %s for storing per family contributions to VC lod score\n",
            (const char *) filename);

   fprintf(perFamily, "%20s %10s %10s %10s %10s\n",
           "FAMILY", "POSITION", "LLK_NULL", "LLK_ALT", "LOD");
   }

void VarianceComponents::ClosePerFamilyFile()
   {
   printf("LOD score contributions for individual families stored in file [%s.vc].\n",
          (const char *) (MerlinCore::filePrefix));
   fclose(perFamily);
   }

void VarianceComponents::WritePerFamilyLOD(Pedigree & ped, IntArray * pheno,
     const char * position, Vector & nullLLK, Vector & fullLLK)
   {
   double lastNull = 0.0;
   double lastFull = 0.0;
   int    index = 0;

   for (int f = 0; f < ped.familyCount; f++)
      if (pheno[f].Length())
          {
          double null = nullLLK[index] - lastNull;
          double full = fullLLK[index] - lastFull;
          double lod = 2 * (full - null) / (2*log(10.0));

          fprintf(perFamily, "%20s %10s %10.5f %10.5f %10.5f\n",
                  (const char *) ped.families[f]->famid, position,
                  null, full, lod);

          lastNull = nullLLK[index];
          lastFull = fullLLK[index++];
          }
   }

void VarianceComponents::SummarizeRepeatedMeasures(Pedigree & ped, String & label, IntArray * pheno)
   {
   int measurements = 0;
   int individuals = 0;

   for (int i = 0; i < ped.familyCount; i++)
      if (pheno[i].Length())
         {
         measurements += pheno[i].Length();

         individuals++;

         for (int j = 1; j < pheno[i].Length(); j++)
            if (pheno[i][j] != pheno[i][j - 1])
               individuals++;
         }

   printf("Repeated measurements found for trait %s\n"
          "   An average of %.1f measurements per subject were taken\n"
          "   A measurement error component will be fitted\n\n",
          (const char *) label, measurements / (individuals + 1e-30));
   }

bool VarianceComponents::CheckAveragedMeasures(Pedigree & ped, IntArray * pheno, int repeatedMeasures)
   {
   bool fitMeasurementError = false;

   if (repeatedMeasures != -1)
      {
      printf("Information on repeated measurements found\n");

      int first = 0;
      double sum = 0.0, count = 1e-30;

      for (int f = 0; f < ped.familyCount; f++)
         for (int i = 0; i < pheno[f].Count(); i++)
            {
            Person & p   = ped[pheno[f][i]];
            int measures = (int) ped[pheno[f][i]].covariates[repeatedMeasures];

            if (measures <= 0)
               error("Invalid number of repeated measurements for\n"
                     "individual %s in family %s\n",
                     (const char *) p.pid, (const char *) p.famid);

            if (first == 0)
               first = measures;
            else if (first != measures)
               fitMeasurementError = true;

            sum += measures;
            count += 1.0;
            }

      printf("   An average of %.1f measurements per subject were taken\n", sum / count);

      printf(fitMeasurementError ?
             "   A measurement error component will be fitted\n\n" :
             "   A measurement error component will not be fitted\n"
             "   because everyone was measured the same number of times\n\n");
      }

   return fitMeasurementError;
   }

void VarianceComponents::SetupPDF(MerlinPDF & pdf, String & trait, int probandCount)
   {
   bool heterogeneity = unlinkedFraction > 0.0;

   pdf.PrepareChart();
   pdf.chart.yAxis.SetMin(0.0);
   pdf.chart.yAxis.minMax = 1.0;
   pdf.chart.yAxis.label = "LOD score";
   pdf.chart.title = trait;
   pdf.chart.title += " [VC";

   if (useCovariates)
      {
      pdf.chart.title += ", ";
      pdf.chart.title += covariates.Length();
      pdf.chart.title += " covars";
      }

   if (useProbands)
      {
      pdf.chart.title += ", ";
      pdf.chart.title += probandCount;
      pdf.chart.title += " probands";
      }

   if (heterogeneity)
      {
      pdf.chart.title += ", ";
      pdf.chart.title += (int) (unlinkedFraction * 100.0);
      pdf.chart.title += "% unlinked";
      }

   pdf.chart.title += "]";
   }

bool VarianceComponents::ExtractMeasurementNumber(const String & label, int & prefix_len, int & number)
   {
   char suffix[] = "_measurement";

   prefix_len = 0;
   number = -1;

   if (label.Length() < 13)
      return false;

   // First extract the measurement number
   int number_start = label.Length() - 1;

   if (!isdigit(label[number_start]))
      return false;

   while (number_start > 12 && isdigit(label[number_start - 1]))
      number_start--;

   number = atoi(((const char *) label) + number_start);

   // Next check whether suffix matches
   for (int i = 0; i < 12; i++)
      if (suffix[i] != tolower(label[number_start - 12 + i]))
         return false;

   prefix_len = number_start - 12;

   return true;
   }

int  VarianceComponents::FindMeasurement(const String & label, int measurement)
   {
   String lookup = label + "_measurement" + measurement;

   return Pedigree::LookupTrait(lookup);
   }

bool VarianceComponents::SkipTrait(int trait)
   {
   int prefix = 0, measurement_number = 0;

   if (ExtractMeasurementNumber(Pedigree::traitNames[trait], prefix, measurement_number) &&
       measurement_number > 1)
       return true;

   return false;
   }

double VarianceComponents::RetrieveMeasurement(Person & person, IntArray & measurements, int & measurement)
   {
   while (person.traits[measurements[measurement]] == _NAN_)
      measurement++;

   return person.traits[measurements[measurement++]];
   }

bool VarianceComponents::IndexMeasurements(IntArray & index, String & label, int trait)
   {
   // Clear index of repeated measurements
   index.Clear();

   // Then check whether we should build one for the current trait
   int measurement = 1, prefix = 0;

   if (!ExtractMeasurementNumber(Pedigree::traitNames[trait], prefix, measurement) ||
       measurement != 1)
       {
       index.Push(trait);
       return false;
       }

   label = Pedigree::traitNames[trait].Left(prefix);

   do
      index.Push(FindMeasurement(label, measurement++));
   while (index.Last() != -1);

   index.Pop();

   return index.Length() > 1;
   }

bool VarianceComponents::IsPhenotyped(Person & person, IntArray & measurements, int trait)
   {
   if (person.isPhenotyped(trait))
      return true;

   for (int i = 1; i < measurements.Length(); i++)
      if (person.isPhenotyped(measurements[i]))
         return true;

   return false;
   }

void VarianceComponents::PrintPvalue(double pvalue)
   {
   if (pvalue > 1.5e-4)
      printf("%7.*f", int(log(pvalue*0.06666666)*-0.4343)+1, pvalue);
   else
      printf("%7.1e", pvalue);
   }

 
