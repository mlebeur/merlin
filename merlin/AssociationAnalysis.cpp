////////////////////////////////////////////////////////////////////// 
// merlin/AssociationAnalysis.cpp 
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
 
#include "AssociationAnalysis.h"
#include "MerlinFamily.h"
#include "MathStats.h"
#include "Kinship.h"

#ifdef  __CHROMOSOME_X__
#include "KinshipX.h"
#endif

RefinedQtlModel AssociationAnalysis::refinedAssociationModels;

void AssociationAnalysis::AnalyseAssociation(FamilyAnalysis & engine)
   {
   Pedigree & ped = engine.ped;
   MerlinKinship & markerKin = engine.kinship;

   int positions = engine.analysisPositions.Length();
   int markers = engine.markers.Length();

   // Array of phenotyped individuals for each family
   IntArray * pheno = new IntArray[ped.familyCount];

   // Summarize covariates ...
   if (customModels.HaveModels())
      customModels.PrintSummary();
   else
      {
      ListCovariates();
      PrintCovariates();
      }

   // Dummy vector used for some calls to variance components code
   IntArray emptyList;

   // Label for analysis of this trait
   String   traitLabel;

   // Figure out how many models we have to analyse
   int models = customModels.HaveModels() ? customModels.modelCount : ped.traitCount;

   int chr = 0;
   String tablename;
   FILE * tablefile = NULL;

   if (MerlinCore::tabulate)
      {
      chr = positions > 0 ? Pedigree::GetMarkerInfo(engine.markers[0])->chromosome : 0;
      tablename.printf("%s-assoc-chr%02d.tbl", (const char *) MerlinCore::filePrefix, chr > 0 ? chr : 0);
      tablefile = fopen(tablename, "wt");
      }

   if (tablefile != NULL)
      fprintf(tablefile, "CHR\tSNP\tALLELE\tFREQ\tTRAIT\tEFFECT\tH2\tLOD\tPVALUE\n");


   // Loop through traits in the pedigree
   for (int m = 0; m < models; m++)
      {
      // Retrieve the appropriate trait
      int t = customModels.HaveModels() ? customModels.traitList[m] : m;

      // Retrieve the appropriate set of covariates
      if (customModels.HaveModels()) covariates = customModels.covariateList[m];

      // Default trait label
      traitLabel = ped.traitNames[t];

      // First list phenotyped individuals for each family
      // If we are using covariates, we only consider individuals
      // for which all covariates have been recorded
      for (int f = 0; f < ped.familyCount; f++)
         {
         pheno[f].Dimension(0);
         if (markerKin.CheckFamily(f))
            for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
               if ( ped[i].isPhenotyped(t) && CheckCovariates(ped[i]))
                  pheno[f].Push(i);
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

      // When analysing custom models, we can skip SNPs and positions to save time
      Vector importantPositions;

      if (customModels.HaveModels() && customModels.candidateList[m].Length())
         {
         importantPositions.Dimension(customModels.candidateList[m].Length());
         importantPositions.Clear();

         for (int i = 0; i < markers; i++)
            if (!customModels.SkipMarker(m, engine.markers[i]))
               importantPositions.Push(engine.markerPositions[i]);

         if (importantPositions.Length() == 0)
            continue;
         }

      // Setup PDF file
      if (engine.writePDF)
         SetupPDF(engine.pdf, traitLabel);

#ifdef __CHROMOSOME_X__
      int vc_count = 3;
#else
      int vc_count = 2;
#endif

      // The normal set class is our workhorse
      NormalSet mvn;
      mvn.Dimension(families, vc_count);

      // Per family likelihoods under the null and alternative hypothesis
      Vector nullPerFamily;

      PolygenicModel(ped, t, mvn, pheno);

      // Track Key Values
      double lkNull = mvn.Evaluate();
      double sampleVar = TotalVariance(mvn.variances, vc_count);
      double sampleH2  = mvn.variances[1];
#ifdef __CHROMOSOME_X__
      double sampleH2X = mvn.variances[2];
#endif

      // Print out header for results
#ifdef __CHROMOSOME_X__
      printf("Phenotype: %s [ASSOC] (%d famil%s, h2 = %.2f%% [AUTO] + %.2f%% [X])\n"
             "===============================================================================\n",
            (const char *) traitLabel, families, families == 1 ? "y" : "ies",
            sampleH2 * 100. / sampleVar, sampleH2X * 100. / sampleVar);
#else
      printf("Phenotype: %s [ASSOC] (%d famil%s, h2 = %.2f%%)\n"
             "===============================================================================\n",
            (const char *) traitLabel, families, families == 1 ? "y" : "ies",
            sampleH2 * 100. / sampleVar);
#endif

      printf("---- LINKAGE TEST RESULTS ----  ----------- ASSOCIATION TEST RESULTS ----------\n");

      printf("%9s %6s %5s %7s %10s %6s %7s %7s %6s %7s\n",
             "Position", "H2 ", "LOD", "pvalue", "Marker", "Allele", "Effect", "H2", "LOD", "pvalue");

      // Add an additional variance component for linked major gene
      mvn.Dimension(families, vc_count + 1);

      // These are use to track the position of the most interesting result
      int    peak_marker = -1, peak_digits = 0, peak_pos = 0;
      double peak_lod = -1., peak_effect = 0., peak_pvalue = 1., peak_h2 = 0.;

      // Loop through analysis at individual marker locations
      for (int pos = 0, marker = 0, pdfOffset = 0; pos < positions; pos++)
         {
         int index = 0;

         if (importantPositions.Length() &&
             importantPositions.FastFind(engine.analysisPositions[pos]) < 0)
             continue;

         // Fill major gene matrix with observed coefficients
         for (int f = 0; f < ped.familyCount; f++)
            if (pheno[f].Length())
               {
               int count = mvn[index].linearModel.rows;

               mvn[index].varComponents[vc_count].Dimension(count, count);
               for (int i = 0; i < count; i++)
                  for (int j = i; j < count; j++)
                     mvn[index].varComponents[vc_count][i][j] =
                     mvn[index].varComponents[vc_count][j][i] =
                                  2.0 * markerKin.Retrieve(f, pos,
                                        ped[pheno[f][i]].traverse,
                                        ped[pheno[f][j]].traverse);

               index++;
               }

         // Fit the alternative model
         mvn.Solve();

         // Evaluate the likelihood assuming a linked major gene
         double lk = mvn.Evaluate();

         // Record per family likelihoods
         nullPerFamily = mvn.recordedLikelihoods;

         // Get key values
         double chisq = lkNull > lk ? (lkNull - lk) * 2.0 : 0.0;
         double pvalue = chidist(chisq, 1) * 0.5;
         double h2 = mvn.variances[vc_count];
         double var = TotalVariance(mvn.variances, vc_count + 1);

         double lod = chisq / (2*log(10.0));

         // Print out a short summary of linkage results
         printf("%9.9s %5.1f%% %5.2f ",
                (const char *) engine.labels[pos], h2 * 100 / var, lod);

         // Print out p-value for linkage tests
         PrintPvalue(pvalue);

         if (engine.writePDF)
            engine.pdf.y[0][pos + pdfOffset] = lod;

         // Check if current analysis position corresponds to an exact
         // marker location ...
         while (marker < markers &&
                engine.markerPositions[marker] < engine.analysisPositions[pos])
            marker++;

         // If it doesn't continue
         if (marker >= markers || engine.markerPositions[marker] > engine.analysisPositions[pos])
            {
            printf("           ---\n");
            continue;
            }

         // Loop to analyze all genotyped markers at current position
         bool first_marker = true;
         while (marker < markers && engine.markerPositions[marker] == engine.analysisPositions[pos])
            {
            // Absolute marker id
            int markerId = engine.markers[marker];

            // Only analyzed markers for which genotypes were inferred
            if (engine.imputationEngine.resultsAvailable(markerId) == false)
               {
               marker++;
               continue;
               }

            // With custom models, only a subset of markers is analysed
            if (customModels.HaveModels() && customModels.SkipMarker(m, markerId))
               {
               marker++;
               continue;
               }

            // Total number of fixed effects in the model
            int parameters = mvn[0].linearModel.cols;
            index = 0;

            for (int f = 0; f < ped.familyCount; f++)
               if (pheno[f].Length())
                  {
                  int count = mvn[index].linearModel.rows;

                  mvn[index].linearModel.Dimension(count, parameters + 1);
                  for (int i = 0; i < count; i++)
                     mvn[index].linearModel[i][parameters] =
                        engine.imputationEngine.GetExpectedGenotype(pheno[f][i], markerId);

                  index++;
                  }

            // Fit allelic association model
            mvn.Solve();

            // Evaluate the likelihood assuming the current marker is associated
            double lkAssoc = mvn.Evaluate();

            // Get key values
            double assoc_chisq  = lk > lkAssoc ? (lk - lkAssoc) * 2.0 : 0.0;
            double assoc_lod = assoc_chisq / (2*log(10.0));
            double assoc_pvalue = chidist(assoc_chisq, 1);
            double assoc_effect = mvn[0].linearModel.cols == parameters ? 0.0 : mvn.means.Last();

            // Calculate variance explained by SNP
            double freq = ped.GetMarkerInfo(markerId)->freq[1];
            double genotypeVariance = 2.0 * freq * (1.0 - freq);
            double assoc_h2 = assoc_effect * assoc_effect * genotypeVariance / var;

            // Update PDF file -- must improve this to allow multiple scores at
            // the same position ...
            if (engine.writePDF)
               {
               // Gymnastics to allow multiple markers at the same position
               if (!first_marker)
                  {
                  pdfOffset++;

                  engine.pdf.x.Insert(pos + pdfOffset, engine.pdf.x[pos + pdfOffset - 1]);
                  engine.pdf.y.Dimension(engine.pdf.y.rows, engine.pdf.y.cols + 1);
                  engine.pdf.y[0][pos + pdfOffset] = lod;

                  for (int i = 0; i < engine.pdf.y.rows; i++)
                     engine.pdf.y[i].Last() = _NAN_;
                  }

               engine.pdf.y[1][pos + pdfOffset] = assoc_lod;
               }

            // Remove marker genotypes from the model before moving on to the
            // next marker
            index = 0;
            for (int f = 0; f < ped.familyCount; f++)
               if (pheno[f].Length())
                  {
                  int count = mvn[index].linearModel.rows;

                  mvn[index].linearModel.Dimension(count, parameters);

                  index++;
                  }

            // Print out the effect size in a pretty fashion ...
            int digits = (assoc_effect < 999.9995 && assoc_effect > -99.9995) ? 3 :
                         (assoc_effect < 9999.995 && assoc_effect > -999.995) ? 2 :
                         (assoc_effect < 99999.95 && assoc_effect > -9999.95) ? 1 : 0;

            // Check for weird results that can occur with uninformative genotypes
            // or colinear variables
            if (assoc_h2 > 5.0)
               printf("%*s %10.10s %6s %7s %7s %6s %7s",
                   first_marker ? 0 : 30, "",
                   (const char *) ped.markerNames[engine.markers[marker]],
                   (const char *) ped.GetMarkerInfo(engine.markers[marker])->GetAlleleLabel(1),
                   "-", "-", "-", "-");
            else
               {
               // Print out results of association analysis
               printf("%*s %10.10s %6s %7.*f %6.2f%% %6.2f ",
                   first_marker ? 0 : 30, "",
                   (const char *) ped.markerNames[engine.markers[marker]],
                   (const char *) ped.GetMarkerInfo(engine.markers[marker])->GetAlleleLabel(1),
                   digits, assoc_effect, assoc_h2 * 100.0, assoc_lod);

               PrintPvalue(assoc_pvalue);
               }

            if (tablefile != NULL)
               {
               fprintf(tablefile, "%d\t%s\t%s\t%.3f\t%s\t",
                   chr, (const char *) ped.markerNames[engine.markers[marker]],
                   (const char *) ped.GetMarkerInfo(engine.markers[marker])->GetAlleleLabel(1),
                   freq,
                   (const char *) traitLabel);

               if (assoc_h2 > 5.0)
                  fprintf(tablefile, "-\t-\t-\t-\n");
               else
                  fprintf(tablefile, "%.3f\t%.3f\t%.3f\t%.4g\n",
                          assoc_effect, assoc_h2 * 100., assoc_lod, assoc_pvalue);
               }

            if (perFamily != NULL)
               WritePerFamilyLOD(ped, pheno, (const char *) ped.markerNames[engine.markers[marker]],
                                 nullPerFamily, mvn.recordedLikelihoods);

            // Check if we found a new peak of allelic association
            if (assoc_lod > peak_lod)
               {
               peak_pos = pos + pdfOffset;
               peak_marker = marker;
               peak_digits = digits;
               peak_lod = assoc_lod;
               peak_effect = assoc_effect;
               peak_h2 = assoc_h2;
               peak_pvalue = assoc_pvalue;
               }

            printf("\n");

            marker++;

            first_marker = false;
            }

         if (first_marker)
            printf("           ---\n");
         }

      if (peak_marker >= 0 && engine.imputationEngine.availableMarkers() > 1)
         {
         printf("%30s %10.10s %6s %7.*f %6.2f%% %6.2f ",
                "Peak -->",
                (const char *) ped.markerNames[engine.markers[peak_marker]],
                (const char *) ped.GetMarkerInfo(engine.markers[peak_marker])->GetAlleleLabel(1),
                peak_digits, peak_effect, peak_h2 * 100.0, peak_lod);

         PrintPvalue(peak_pvalue);

         printf("\n");

         if (engine.writePDF)
            {
            engine.pdf.chart.SetSeriesLabel(2, ped.markerNames[engine.markers[peak_marker]]);
            engine.pdf.y[2][peak_pos] = peak_lod;
            }
         }

      if (peak_marker >= 0 && customModels.HaveModels())
         {
         if (customModels.modelCount != refinedAssociationModels.modelCount)
            refinedAssociationModels.Dimension(customModels.modelCount);

         if (refinedAssociationModels.CheckForImprovement(m, peak_lod))
            {
            refinedAssociationModels.scores[m] = peak_lod;
            refinedAssociationModels.labels[m] = ped.markerNames[engine.markers[peak_marker]];
            refinedAssociationModels.comments[m].printf(
               " Adding %s to the model produced a LOD of %.3f, p-value %#.3g",
               (const char *) ped.markerNames[engine.markers[peak_marker]],
               peak_lod, peak_pvalue);

            refinedAssociationModels.values[m].Dimension(ped.count);
            for (int i = 0; i < ped.count; i++)
               refinedAssociationModels.values[m][i] =
                  engine.imputationEngine.GetExpectedGenotype(i, engine.markers[peak_marker]);
            }
         }

      if (engine.writePDF)
         engine.pdf.DrawChart();

      printf("\n");
      }


   if (tablefile != 0)
      {
      fclose(tablefile);
      printf("Association results tabulated in [%s]\n\n", (const char *) tablename);
      }

   delete [] pheno;
   }

void AssociationAnalysis::SetupPDF(MerlinPDF & pdf, String & trait)
   {
   pdf.PrepareChart(3);
   pdf.chart.yAxis.SetMin(0.0);
   pdf.chart.yAxis.minMax = 1.0;
   pdf.chart.yAxis.label = "LOD score";
   pdf.chart.title = trait;
   pdf.chart.title += " [Association Analysis]";
   pdf.chart.SetSeriesLabel(0, "Linkage");
   pdf.chart.SetSeriesLabel(1, "Association");
   pdf.chart.SetSeriesLabel(2, "Peak Marker");
   pdf.chart.SetMarkerShape(1, lmCircle);
   pdf.chart.ShowLine(1, false);
   pdf.chart.ShowMarker(1, true);
   pdf.chart.SetMarkerShape(2, lmCross);
   pdf.chart.SetSeriesColor(2, 1.0, 0.0, 0.0);
//   pdf.chart.SetLineColor(2, 1.0, 0.0, 0.0);
   pdf.chart.ShowLine(2, false);
   pdf.chart.ShowMarker(2, true);
   pdf.chart.useLegend = true;
   pdf.y.Set(_NAN_);
   }

void AssociationAnalysis::PolygenicModel(
   Pedigree & ped, int t, NormalSet & mvn, IntArray * pheno)
   {
   // Kinship matrices used for scoring polygenic component
   Kinship kin;
#ifdef __CHROMOSOME_X__
   KinshipX kinX;
#endif

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

         // Phenotypes
         mvn[index].scores.Dimension(count);
         for (int i = 0; i < count; i++)
            mvn[index].scores[i] = ped[pheno[f][i]].traits[t];

         // Total number of fixed effects
         int fixedEffects = 1 + covariates.Length();

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

         index++;
         }

   // Record the last log-likelihood of the bunch
   if (index && (perFamily != NULL))
      mvn.operators[index - 1] |= NORMAL_LAST_OP(NORMAL_RECORD_LLK);

   // Fit polygenic model
   mvn.Solve();
   }

void AssociationAnalysis::OutputRefinedModels(const String & prefix, Pedigree & ped)
   {
   if (customModels.HaveModels())
      {
      printf("Refined association models stored in [%s-assoc-covars.*]\n",
             (const char *) prefix);

      refinedAssociationModels.WriteFiles(prefix + "-assoc-covars", ped, customModels);
      }
   }

void AssociationAnalysis::OpenPerFamilyFile(int chromosome)
   {
   String filename(MerlinCore::filePrefix);
   filename += "-chr";
   filename += chromosome;
   filename += ".assoc";

   perFamily = fopen(filename, "wt");

   if (perFamily == NULL)
      error("Opening file %s for storing per family contributions to association score\n",
            (const char *) filename);

   fprintf(perFamily, "%20s %10s %10s %10s %10s\n",
           "FAMILY", "POSITION", "LLK_NULL", "LLK_ALT", "LOD");
   }

void AssociationAnalysis::ClosePerFamilyFile(int chromosome)
   {
   printf("Association score contributions for individual families stored in file [%s-chr%d.assoc].\n",
          (const char *) (MerlinCore::filePrefix), chromosome);
   fclose(perFamily);
   }

 
