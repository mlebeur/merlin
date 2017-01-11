////////////////////////////////////////////////////////////////////// 
// merlin/FastAssociation.cpp 
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
 
#include "FastAssociation.h"
#include "MerlinFamily.h"
#include "MathStats.h"

double FastAssociationAnalysis::fastFilter = _NAN_;

RefinedQtlModel FastAssociationAnalysis::refinedFastModels;

void FastAssociationAnalysis::AnalyseAssociation(FamilyAnalysis & engine)
   {
   Pedigree & ped = engine.ped;

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
      tablename.printf("%s-fastassoc-chr%02d.tbl", (const char *) MerlinCore::filePrefix, chr > 0 ? chr : 0);
      tablefile = fopen(tablename, "wt");
      }

   if (tablefile != NULL)
      fprintf(tablefile, "CHR\tSNP\tAL1\tAL2\tFREQ1\tTRAIT\tEFFECT\tSE\tH2\tLOD\tPVALUE\n");

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
         for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            if (ped[i].isPhenotyped(t) && CheckCovariates(ped[i]))
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

      // Fit a base polygenic model
      PolygenicModel(ped, t, mvn, pheno);

      // Evaluate the likelihood at the fitted point to ensure that
      // polygenic model residuals are available
      mvn.Evaluate();

      // Track Key Values
      double sampleVar = TotalVariance(mvn.variances, vc_count);
      double sampleH2  = mvn.variances[1];
#ifdef __CHROMOSOME_X__
      double sampleH2X = mvn.variances[2];
#endif

      // Record per family likelihoods
      nullPerFamily = mvn.recordedLikelihoods;

      if (sampleVar == 0.0)
         {
         printf("WARNING: Phenotype %s is constant and unsuitable for association testing\n\n",
                (const char *) traitLabel);
         continue;
         }

      // Print out header for results
#ifdef __CHROMOSOME_X__
      printf("Phenotype: %s [FAST-ASSOC] (%d famil%s, h2 = %.2f%% [AUTO] + %.2f%% [X])\n"
             "==============================================================================\n",
            (const char *) traitLabel, families, families == 1 ? "y" : "ies",
            sampleH2 * 100. / sampleVar, sampleH2X * 100. / sampleVar);
#else
      printf("Phenotype: %s [FAST-ASSOC] (%d famil%s, h2 = %.2f%%)\n"
             "==============================================================================\n",
            (const char *) traitLabel, families, families == 1 ? "y" : "ies",
            sampleH2 * 100. / sampleVar);
#endif

      printf("%10s %13s %7s %7s %7s %7s %7s %7s\n",
             "Position", "Marker", "Allele", "Effect", "StdErr", "H2", "LOD", "pvalue");

      // These are use to track the position of the most interesting result
      int    peak_marker = -1, peak_digits = -1, peak_sdigits = -1, peak_pos = -1, tests = 0;
      double peak_lod = -1., peak_effect = 0., peak_stderr = 0., peak_h2 = 0., peak_pvalue = 1.;

      // Loop through analysis at individual marker locations
      for (int pos = 0, marker = 0, pdfOffset = 0; pos < positions; pos++)
         {
         // Check if current analysis position corresponds to an exact
         // marker location ...
         while (marker < markers &&
                engine.markerPositions[marker] < engine.analysisPositions[pos])
            marker++;

         // If it doesn't continue
         if (marker >= markers || engine.markerPositions[marker] > engine.analysisPositions[pos])
            continue;

         // We may need special handling of PDF output when there are
         // multiple markers in the same position
         bool first_marker = true;

         // Loop to analyze all genotyped markers at current position
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

            // Only analyze markers that are not filter
            if (customModels.HaveModels() && customModels.SkipMarker(m, markerId))
               {
               marker++;
               continue;
               }

            // Vector of genotypes for each family
            Vector genotypes;

            // Allele frequency
            double freq = ped.GetMarkerInfo(markerId)->freq[1];

            // Expected genotype for each individual
            double expectedGenotype = 2.0 * freq;

            // Variance genotype scores
            double genotypeVariance = 2.0 * freq * (1.0 - freq);

            // Z statistic evaluating evidence for association
            double numerator = 0.0, denominator = 1e-10; // SSxy = 0.0, SSxx = 0.0;

            // Loop through families and calculate score statistic to evaluate
            // evidence for association
            for (int f = 0, index = 0, count; f < ped.familyCount; f++)
               if ((count = pheno[f].Length()) > 0)
                  {
                  // First calculate a vector of expected genotypes
                  genotypes.Dimension(count);
                  for (int i = 0; i < count; i++)
                     genotypes[i] =
                        engine.imputationEngine.GetExpectedGenotype(ped, pheno[f][i], markerId) -
                        expectedGenotype;

                  // Next, calculate Genotype * SIGMA^-1
                  mvn[index].cholesky.BackSubst(genotypes);

                  // Numerator of test statistic is Genotype * SIGMA^-1 * Phenotype
                  numerator += mvn[index].cholesky.x.InnerProduct(mvn[index].residuals);

                  // Denominator of test statistic is Genotype * SIGMA^-1 * Genotype
                  denominator += mvn[index].cholesky.x.InnerProduct(genotypes);

                  // printf("n = %.5f, d = %.5f\n", numerator * mvn.variances.Sum(), denominator * mvn.variances.Sum());

                  // Update residual sum of squares for estimating effect size
                  // SSxx += genotypes.SumSquares();
                  // SSxy += mvn[index].residuals.InnerProduct(genotypes);

                  // Update family index
                  index++;
                  }

            // Get key values
            double assoc_chisq  = numerator * numerator / denominator;
            double assoc_lod    = assoc_chisq / (2*log(10.0));
            double assoc_pvalue = chidist(assoc_chisq, 1);
            double assoc_effect = -numerator / denominator; // -SSxy / SSxx;
            double assoc_stderr = 1.0 / sqrt(denominator);
            double assoc_h2     = assoc_effect * assoc_effect * genotypeVariance / sampleVar;
            tests++;

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

                  for (int i = 0; i < engine.pdf.y.rows; i++)
                     engine.pdf.y[i].Last() = _NAN_;
                  }

               engine.pdf.y[1][pos + pdfOffset] = assoc_lod;
               first_marker = false;
               }

            // Print out the effect size in a pretty fashion ...
            int digits = (assoc_effect < 999.9995 && assoc_effect > -99.9995) ? 3 :
                         (assoc_effect < 9999.995 && assoc_effect > -999.995) ? 2 :
                         (assoc_effect < 99999.95 && assoc_effect > -9999.95) ? 1 : 0;

            int sdigits = (assoc_stderr < 999.9995 && assoc_stderr > -99.9995) ? 3 :
                          (assoc_stderr < 9999.995 && assoc_stderr > -999.995) ? 2 :
                          (assoc_stderr < 99999.95 && assoc_stderr > -9999.95) ? 1 : 0;

            // Check if we found a new peak of allelic association
            if (assoc_lod > peak_lod)
               {
               peak_pos = pos + pdfOffset;
               peak_marker = marker;
               peak_digits = digits;
               peak_lod = assoc_lod;
               peak_effect = assoc_effect;
               peak_sdigits = sdigits;
               peak_stderr = assoc_stderr;
               peak_h2 = assoc_h2;
               peak_pvalue = assoc_pvalue;
               }

            if (fastFilter != _NAN_ && assoc_pvalue > fastFilter)
               {
               marker++;
               continue;
               }

            // Check for weird results that can occur with uninformative genotypes
            // or colinear variables
            if (assoc_h2 > 5.0)
               printf("%10.10s %13.13s %7s %7s %7s %7s %7s %7s",
                   (const char *) engine.labels[pos],
                   (const char *) ped.markerNames[engine.markers[marker]],
                   (const char *) ped.GetMarkerInfo(engine.markers[marker])->GetAlleleLabel(1),
                   "-", "-", "-", "-", "-");
            else
               {
               // Print out results of association analysis
               printf("%10.10s %13.13s %7s %7.*f %7.*f %6.2f%% %7.3f ",
                   (const char *) engine.labels[pos],
                   (const char *) ped.markerNames[engine.markers[marker]],
                   (const char *) ped.GetMarkerInfo(engine.markers[marker])->GetAlleleLabel(1),
                   digits, assoc_effect, sdigits, assoc_stderr,
                   assoc_h2 * 100.0, assoc_lod);

               PrintPvalue(assoc_pvalue);
               }


            if (tablefile != NULL)
               {
               fprintf(tablefile, "%d\t%s\t%s\t%s\t%.3f\t%s\t",
                   chr, (const char *) ped.markerNames[engine.markers[marker]],
                   (const char *) ped.GetMarkerInfo(engine.markers[marker])->GetAlleleLabel(1),
                   (const char *) ped.GetMarkerInfo(engine.markers[marker])->GetAlleleLabel(2),
                   freq,
                   (const char *) traitLabel);

               if (assoc_h2 > 5.0)
                  fprintf(tablefile, "-\t-\t-\t-\t-\n");
               else
                  fprintf(tablefile, "%.3f\t%.3f\t%.3f\t%.3f\t%.4g\n",
                          assoc_effect, assoc_stderr, assoc_h2 * 100., assoc_lod, assoc_pvalue);
               }

            printf("\n");

            marker++;
            }
         }

      if (peak_marker >= 0 && tests > 1)
         {
         printf("%10s %13.13s %7s %7.*f %7.*f %6.2f%% %7.3f ",
                "Peak -->",
                (const char *) ped.markerNames[engine.markers[peak_marker]],
                (const char *) ped.GetMarkerInfo(engine.markers[peak_marker])->GetAlleleLabel(1),
                peak_digits, peak_effect, peak_sdigits, peak_stderr,
                peak_h2 * 100.0, peak_lod);

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
         if (customModels.modelCount != refinedFastModels.modelCount)
            refinedFastModels.Dimension(customModels.modelCount);

         if (refinedFastModels.CheckForImprovement(m, peak_lod))
            {
            refinedFastModels.scores[m] = peak_lod;
            refinedFastModels.labels[m] = ped.markerNames[engine.markers[peak_marker]];
            refinedFastModels.comments[m].printf(
               " Adding %s to the model produced a LOD of %.3f, p-value %#.3g",
               (const char *) ped.markerNames[engine.markers[peak_marker]],
               peak_lod, peak_pvalue);

            refinedFastModels.values[m].Dimension(ped.count);
            for (int i = 0; i < ped.count; i++)
               refinedFastModels.values[m][i] =
                  engine.imputationEngine.GetExpectedGenotype(ped, i, engine.markers[peak_marker]);
            }
         }

      if (engine.writePDF)
         engine.pdf.DrawChart();

      printf("\n");
      }

   if (tablefile != 0)
      {
      fclose(tablefile);
      printf("Score test association results tabulated in [%s]\n\n", (const char *) tablename);
      }

   delete [] pheno;
   }

void FastAssociationAnalysis::OutputFastModels(const String & prefix, Pedigree & ped)
   {
   if (customModels.HaveModels())
      {
      printf("Refined fast association models stored in [%s-fast-covars.*]\n",
             (const char *) prefix);

      refinedFastModels.WriteFiles(prefix + "-fast-covars", ped, customModels);
      }
   }


 
