////////////////////////////////////////////////////////////////////// 
// merlin/KongAndCox.cpp 
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
 
#include "KongAndCox.h"

#include "MerlinFamily.h"
#include "MathStats.h"
#include "NPL-ASP.h"
#include "NPL-QTL.h"

#include <math.h>

// User defined variables
//

bool KongAndCox::nplAll = false;
bool KongAndCox::nplPairs = false;
bool KongAndCox::nplQtl = false;
bool KongAndCox::nplDeviates = false;
bool KongAndCox::nplExponential = false;
bool KongAndCox::rawOutput = false;

// Additional non-parametric linkage tests

bool KongAndCox::nplExtras = false;

// Routines for evaluating LOD scores calculated with the linear model
//

LinearModel::LinearModel()
   {
   statistics = families = positions = 0;
   min = max = NULL;
   scores = NULL;
   }

LinearModel::~LinearModel()
   {
   FreeMemory();
   }

void LinearModel::FreeMemory()
   {
   if (min != NULL) delete [] min;
   if (max != NULL) delete [] max;

   if (scores != NULL)
      {
      for (int s = 0; s < statistics; s++)
         delete [] scores[s];
      delete [] scores;
      }

   min = max = NULL;
   scores = NULL;
   }

void LinearModel::AllocateMemory()
   {
   FreeMemory();

   min = new Vector[statistics];
   max = new Vector[statistics];
   scores = new Vector * [statistics];

   for (int s = 0; s < statistics; s++)
      scores[s] = new Vector [families];
   }

void LinearModel::Dimension(int Statistics, int Families, int Positions)
   {
   if (statistics != Statistics || families != Families)
      {
      statistics = Statistics;
      families = Families;
      AllocateMemory();
      }

   positions = Positions;
   for (int i = 0; i < statistics; i++)
      {
      min[i].Dimension(families);
      min[i].Zero();
      max[i].Dimension(families);
      max[i].Zero();

      for (int j = 0; j < families; j++)
         {
         scores[i][j].Dimension(positions);
         scores[i][j].Zero();
         }
      }
   }

double LinearModel::f(double delta)
   {
   double sum = 0.0;

   for (int i = 0; i < families; i++)
      sum += 2.0 * log(1.0 + delta * scores[stat][i][pos]);

   return -sum;
   }

void LinearModel::UpdateDeltaBounds()
   {
   // Calculate the range of Z-scores and d-hat
   double lo = 0.0, hi = 0.0;

   // Count the number of informative families
   informativeFamilies = 0;

   for (int i = 0; i < families; i++)
      if (min[stat][i] != max[stat][i])
         {
         if (min[stat][i] < lo) lo = min[stat][i];
         if (max[stat][i] > hi) hi = max[stat][i];
         informativeFamilies++;
         }

   if (lo == hi)
      {
      minDelta = maxDelta = scale = 0.0;
      return;
      }

   minDelta = -1.0 / hi;
   maxDelta = -1.0 / lo;
   scale = 1.0 / sqrt((double) (informativeFamilies));
   }

void LinearModel::CalculateMinScore(double & Z, double & delta, double & chisq)
   {
   Z = chisq = 0.0;
   delta = minDelta + TINY;

   // Calculate the minimum & maximum possible scores
   for (int i = 0; i < families; i++)
      {
      Z     += min[stat][i];
      chisq += 2.0 * log(1 + delta * min[stat][i]);
      }
   Z *= scale;
   }

void LinearModel::CalculateMaxScore(double & Z, double & delta, double & chisq)
   {
   Z = chisq = 0.0;
   delta = maxDelta - TINY;

   // Calculate the minimum & maximum possible scores
   for (int i = 0; i < families; i++)
      {
      Z     += max[stat][i];
      chisq += 2.0 * log(1 + delta * max[stat][i]);
      }
   Z *= scale;
   }

void LinearModel::CalculateScore(double & Z, double & delta, double & chisq)
   {
   // First, calculate the Z score
   Z = 0.0;
   for (int i = 0; i < families; i++)
      Z += scores[stat][i][pos];
   Z *= scale;

   // Setup boundaries for minimization
   a = minDelta + TINY;
   c = maxDelta - TINY;
   b = a + (c - a) * GOLD;
   fb = f(b);

   // Optimize function -- i.e. find value of delta that maximizes
   // likelihood and gives the single parameter Kong and Cox LOD score
   Brent();

   delta = ScalarMinimizer::min;
   chisq = -ScalarMinimizer::fmin;
   }

// Routines for evaluating LOD scores calculated with the exponential model
//

ExponentialModel::ExponentialModel()
   {
   statistics = families = positions = 0;
   scores = probs = NULL;
   }

ExponentialModel::~ExponentialModel()
   {
   FreeMemory();
   }

void ExponentialModel::AllocateMemory()
   {
   FreeMemory();

   scores = new Vector * [statistics];
   probs = new Vector * [statistics];

   for (int s = 0; s < statistics; s++)
      {
      scores[s] = new Vector [families];
      probs[s] = new Vector [families];
      }
   }

void ExponentialModel::FreeMemory()
   {
   if (scores != NULL)
      {
      for (int s = 0; s < statistics; s++)
         {
         delete [] scores[s];
         delete [] probs[s];
         }

      delete [] scores;
      delete [] probs;
      }
   }

void ExponentialModel::Dimension(int Statistics, int Families, int Positions)
   {
   if (statistics != Statistics || families != Families)
      {
      statistics = Statistics;
      families = Families;
      AllocateMemory();
      }

   positions = Positions;
   for (int i = 0; i < statistics; i++)
      for (int j = 0; j < families; j++)
         {
         scores[i][j].Dimension(0);
         probs[i][j].Dimension(0);
         }
   }

void ExponentialModel::NullDistribution(int statistic, int family, Vector & z, Vector & p)
   {
   // Count distinct NPL scores
   int distinct = z.Length();

   // Save list of possible NPL scores
   scores[statistic][family] = z;

   // Ignore families where NPL score is invariant
   if (distinct <= 1)
      return;

   // Save null distribution of NPL scores
   probs[statistic][family].Dimension(z.Length() * (positions + 1));
   for (int i = 0; i < distinct; i++)
      probs[statistic][family][i] = p[i];
   }

void ExponentialModel::UpdateDistribution(int statistic, int family, int position, Vector & p)
   {
   int distinct = scores[statistic][family].Length();
   int current  = distinct * (position + 1);

   // Ignore families where NPL score is invariant
   if (distinct <= 1)
      return;

   for (int i = 0; i < distinct; i++, current++)
      probs[statistic][family][current] = p[i];
   }

double ExponentialModel::f(double x)
   {
   double score = 0.0;

   for (int f = 0; f < families; f++)
      if (scores[stat][f].Length() > 1)
         {
         double scale = 0.0;
         double sum = 0.0;

         // Number of possible values for NPL statistic
         int distinct = scores[stat][f].Length();

         // Calculate the location of the first NPL statistic
         int current  = (pos + 1) * distinct;

         for (int i = 0; i < distinct; i++, current++)
            {
            // Calculate updated probability for each NPL score
            double factor = exp(x * scores[stat][f][i]);

            // This summation updates standardization constant
            scale += probs[stat][f][i] * factor;

            // This summation updates the observed statistic
            sum   += probs[stat][f][current] * factor;
            }

         score += log(sum / scale);
         }

   return -2.0 * score;
   }

void ExponentialModel::UpdateDeltaBounds()
   {
   // Calculate the range of Z-scores and d-hat
   double lo = 0.0, hi = 0.0;

   for (int i = 0; i < families; i++)
      if (scores[stat][i].Length() > 1)
         {
         if (scores[stat][i][0] < lo)
            lo = scores[stat][i][0];

         if (scores[stat][i].Peek() > hi)
            hi = scores[stat][i].Peek();
         }

   minDelta = lo < -10 ? 100 / fabs(lo) : -9.999;
   maxDelta = hi > 10  ? 100 / fabs(hi) : +9.999;
   }

void ExponentialModel::CalculateMinScore(double & delta, double & chisq)
   {
   delta = minDelta;
   chisq = 0.0;

   for (int f = 0; f < families; f++)
      if (scores[stat][f].Length() > 1)
         {
         double scale = 0.0, sum = 0.0, norm = 0.0;

         // Number of possible values for NPL statistic
         int distinct = scores[stat][f].Length();

         for (int i = 0; i < distinct; i++)
            {
            // Calculate updated probability for each NPL score
            double factor = exp(delta * scores[stat][f][i]);

            // This summation updates standardization constant
            scale += probs[stat][f][i] * factor;

            // This summation updates the observed statistic
            sum   += factor * factor;
            norm  += factor;
            }

         chisq += log((sum / norm) / scale);
         }

   chisq *= 2.0;
   }

void ExponentialModel::CalculateMaxScore(double & delta, double & chisq)
   {
   delta = maxDelta;
   chisq = 0.0;

   for (int f = 0; f < families; f++)
      if (scores[stat][f].Length() > 1)
         {
         double scale = 0.0, sum = 0.0, norm = 0.0;

         // Number of possible values for NPL statistic
         int distinct = scores[stat][f].Length();

         for (int i = 0; i < distinct; i++)
            {
            // Calculate updated probability for each NPL score
            double factor = exp(delta * scores[stat][f][i]);

            // This summation updates standardization constant
            scale += probs[stat][f][i] * factor;

            // This summation updates the observed statistic
            sum   += factor * factor;
            norm  += factor;
            }

         chisq += log((sum / norm) / scale);
         }

   chisq *= 2.0;
   }

void ExponentialModel::CalculateScore(double & delta, double & chisq)
   {
   // Setup boundaries for minimization
   a = minDelta;
   c = maxDelta;
   b = a + (c - a) * GOLD;
   fb = f(b);

   // Optimize function -- i.e. find value of delta that maximizes
   // likelihood and gives the single parameter Kong and Cox LOD score
   Brent();

   delta = ScalarMinimizer::min;
   chisq = -ScalarMinimizer::fmin;
   }

// Routines for interface with the user and internal merlin routines
//

KongAndCox::KongAndCox(Pedigree & ped)
   {
   // Initialize key variables
   file = NULL;
   tablefile = NULL;
   rawOutputFile = NULL;

   // Initialize pointer variables
   index = NULL;
   npl = NULL;

   // How many NPL statistics to calculate
   statistics = (nplAll + nplPairs + nplExtras * 2) * ped.affectionCount +
                (nplQtl + nplDeviates) * ped.traitCount;

   AllocateMemory();
   LabelAnalyses(ped);
   }

void KongAndCox::AllocateMemory()
   {
   if (statistics > 0)
      {
      means.Dimension(statistics);
      scales.Dimension(statistics);

      npl = new Tree[statistics];
      if (nplExponential) index = new IndexTree[statistics];
      }
   }

void KongAndCox::LabelAnalyses(Pedigree & ped)
   {
   pheno.Clear();
   phenoTag.Clear();

   if (nplAll)
      for (int a = 0; a < ped.affectionCount; a++)
         pheno.Add(ped.affectionNames[a] + " [ALL]"),
         phenoTag.Add(ped.affectionNames[a] + ":all");

   if (nplPairs)
      for (int a = 0; a < ped.affectionCount; a++)
         pheno.Add(ped.affectionNames[a] + " [Pairs]"),
         phenoTag.Add(ped.affectionNames[a] + ":pairs");

   if (nplQtl)
      for (int t = 0; t < ped.traitCount; t++)
         pheno.Add(ped.traitNames[t] + " [QTL]"),
         phenoTag.Add(ped.traitNames[t] + ":qtl");

   if (nplDeviates)
      for (int t = 0; t < ped.traitCount; t++)
         pheno.Add(ped.traitNames[t] + " [QTL, Mean = 0.0]"),
         phenoTag.Add(ped.traitNames[t] + ":deviates");

   if (nplExtras)
      for (int a = 0; a < ped.affectionCount; a++)
         {
         pheno.Add(ped.affectionNames[a] + " [Maternal Pairs]");
         phenoTag.Add(ped.affectionNames[a] + ":maternal_pairs");

         pheno.Add(ped.affectionNames[a] + " [Paternal Pairs]");
         phenoTag.Add(ped.affectionNames[a] + ":paternal_pairs");
         }
   }

KongAndCox::~KongAndCox()
   {
   if (npl != NULL) delete [] npl;
   if (index != NULL) delete [] index;
   }

void KongAndCox::Setup(AnalysisInfo & info)
   {
   families = info.families;
   positions = info.positions;

   lm.Dimension(statistics, families, positions);

   if (nplExponential)
      em.Dimension(statistics, families, positions);
   }

void KongAndCox::SetupFamily(AnalysisInfo & info)
   {
   int discrete = Pedigree::affectionCount;
   int quantitative = Pedigree::traitCount;
   int s = 0;

   if (nplAll)
      {
      NPL_ALL engine;

      for (int a = 0; a < discrete; a++)
         {
         engine.ScoreNPL(*info.m, a);
         npl[s].Copy(engine);
         ScoreNPLBounds(s++, info.famno);
         }
      }

   if (nplPairs)
      {
      NPL_ASP_Tree * engine = FamilyAnalysis::simwalk2 ?
         (NPL_ASP_Tree *) new NPL_Pairs_For_Simwalk2 :
         (NPL_ASP_Tree *) new NPL_Pairs;

      for (int a = 0; a < discrete; a++)
         {
         engine->ScoreNPL(*info.m, a);
         npl[s].Copy(*engine);
         ScoreNPLBounds(s++, info.famno);
         }

      delete engine;
      }

   if (nplQtl)
      {
      NPL_QTL engine;

      for (int t = 0; t < quantitative; t++)
         {
         engine.ScoreWithSampleMean(*info.m, t);
         npl[s].Copy(engine);
         ScoreNPLBounds(s++, info.famno);
         }
      }

   if (nplDeviates)
      {
      NPL_QTL engine;

      for (int t = 0; t < quantitative; t++)
         {
         engine.ScoreWithMeanZero(*info.m, t);
         npl[s].Copy(engine);
         ScoreNPLBounds(s++, info.famno);
         }
      }

   if (nplExtras)
      {
      NPL_Pairs_Maternal maternal;
      NPL_Pairs_Paternal paternal;

      for (int a = 0; a < discrete; a++)
         {
         maternal.ScoreNPL(*info.m, a);
         npl[s].Copy(maternal);
         ScoreNPLBounds(s++, info.famno);

         paternal.ScoreNPL(*info.m, a);
         npl[s].Copy(paternal);
         ScoreNPLBounds(s++, info.famno);
         }
      }
   }

bool KongAndCox::CalculateNPL(Pedigree & ped)
   {
   return (
           (nplAll || nplPairs || nplExtras) && ped.affectionCount ||
           (nplQtl || nplDeviates) && ped.traitCount
          );
   }

void KongAndCox::ScoreNPLBounds(int s, int f)
   {
   stats.GetBounds(npl[s]);

   means[s] = stats.mean;
   if (stats.var > TINY)
      {
      scales[s] = 1.0 / sqrt(stats.var);
      lm.min[s][f] = (stats.min - stats.mean) * scales[s];
      lm.max[s][f] = (stats.max - stats.mean) * scales[s];
      }
   else
      scales[s] = 0.0;

   if (nplExponential)
     {
      index[s].MakeIndex(npl[s]);
      index[s].uniqValues -= stats.mean;
      index[s].uniqValues *= scales[s];

     index[s].EvenFrequencies();
      em.NullDistribution(s, f, index[s].uniqValues, index[s].freqs);
      }

   npl[s].Pack();
   }

void KongAndCox::SkipFamily(AnalysisInfo & info)
   {
   for (int s = 0; s < statistics; s++)
      {
      scales[s] = 0.0;
      means[s] = 0.0;

      lm.min[s][info.famno] = 0.0;
      lm.max[s][info.famno] = 0.0;
      lm.scores[s][info.famno].Zero();

      if (!nplExponential) continue;

      em.scores[s][info.famno].Dimension(0);
      em.probs[s][info.famno].Dimension(0);
      }
   }

void KongAndCox::Analyse(AnalysisInfo & info, Tree & tree, int position)
   {
   if (info.lk == -1.0)
      info.lk = tree.WeightedSum();

   for (int s = 0; s < statistics; s++)
      // Only analyse informative families (eg. multiple related affecteds)
      if (scales[s] != 0.0)
         {
         // Calculate Z-score for linear model
         npl[s].UnPack();
         lm.scores[s][info.famno][position] =
            (tree.MeanProduct(npl[s])/info.lk - means[s]) * scales[s];
         npl[s].RePack();

         if (!nplExponential) continue;

         // Calculate conditional distribution of Z-scores for exponential model
         index[s].UpdateFrequencies(tree, 1.0 / info.lk);
         em.UpdateDistribution(s, info.famno, position, index[s].freqs);
         }
   }

void KongAndCox::Report(AnalysisInfo & info)
   {
   // Create an alias for convenience
   MerlinPDF & pdf = *info.pdf;

   // Cycle through each trait and position
   for (int s = 0; s < statistics; s++)
      {
      // Select the appropriate statistic
      em.stat = lm.stat = s;

      // Calculate linear model bounds for delta parameter
      lm.UpdateDeltaBounds();

      if (lm.minDelta == lm.maxDelta)
         {
         printf("No informative families for %s\n\n",
            (const char *) pheno[s]);
         continue;
         }

      // Update bounds for delta parameter in exponential model
      if (nplExponential)
         em.UpdateDeltaBounds();

      printf("Phenotype: %s (%d famil%s)\n"
             "======================================================%s\n"
             "            Pos   Zmean  pvalue %selta    LOD  pvalue%s\n",
            (const char *) pheno[s], lm.informativeFamilies,
            lm.informativeFamilies == 1 ? "y" : "ies",
            nplExponential ? "=========================" : "",
            nplExponential ? "linD" : "   d",
            nplExponential ? " expDelta    LOD  pvalue" : "");

      if (pdf.isOpen())
         {
         pdf.PrepareChart(nplExponential ? 2 : 1);
         pdf.chart.yAxis.SetMin(0.0);
         pdf.chart.yAxis.minMax = 1.0;
         pdf.chart.yAxis.label = "LOD score";
         pdf.chart.title = pheno[s];

         if (nplExponential)
            {
            pdf.chart.SetSeriesLabel(0, "K&C lin");
            pdf.chart.SetSeriesLabel(1, "K&C exp");
            pdf.chart.useLegend = true;
            }
         }

      double Z, delta, chisq;

      // Calculate the minimum and maximum possible LOD scores for each model
      lm.CalculateMinScore(Z, delta, chisq);

      OutputLabel(0, _NAN_, pheno[s], "min", Z);
      OutputLOD(delta, chisq);

      if (nplExponential)
         {
         em.CalculateMinScore(delta, chisq);
         OutputLOD(delta, chisq);
         }

      if (tablefile != NULL) fprintf(tablefile, "\n");
      printf("\n");

      lm.CalculateMaxScore(Z, delta, chisq);

      OutputLabel(0, _NAN_, pheno[s], "max", Z);
      OutputLOD(delta, chisq);

      if (nplExponential)
         {
         em.CalculateMaxScore(delta, chisq);
         OutputLOD(delta, chisq);
         }

      if (tablefile != NULL) fprintf(tablefile, "\n");
      printf("\n");

      // Calculate LOD scores at each analysis position
      for (int pos = 0; pos < positions; pos++)
         {
         // Update position indicators in each subclass
         lm.pos = em.pos = pos;

         // Calculation using linear model
         lm.CalculateScore(Z, delta, chisq);

         OutputLabel(info.chromosome, (*info.analysisPositions)[pos],
                    pheno[s], (*info.labels)[pos], Z);
         OutputLOD(delta, chisq);

         // Output scores for individual families
         if (file)
            OutputPerFamily(*info.m->pedigree, (*info.labels)[pos], delta);

         // Update PDF file
         if (info.pdf->isOpen())
            info.pdf->y[0][pos] =
               chisq * (delta > 0 ? 0.2171472409516 : -0.2171472409516);

         if (!nplExponential)
            {
            if (tablefile != NULL) fprintf(tablefile, "\n");
            printf("\n");
            continue;
            }

         // Calculation using exponential model
         em.CalculateScore(delta, chisq);

         OutputLOD(delta, chisq);

         // Update PDF file
         if (info.pdf->isOpen())
            info.pdf->y[1][pos] =
               chisq * (delta > 0 ? 0.2171472409516 : -0.2171472409516);

         if (tablefile != NULL) fprintf(tablefile, "\n");
         printf("\n");
         }

      if (pdf.isOpen())
         pdf.DrawChart();

      printf("\n");
      }

   // Output raw information for downstream processors
   if (rawOutputFile)
      OutputRawScores(*info.m->pedigree, info.chromosome, (*info.labels));
   }

void KongAndCox::OpenFiles(AnalysisInfo & info, String & prefix)
   {
   if (MerlinCore::tabulate)
      {
      tablename = prefix + "-nonparametric.tbl";
      tablefile = fopen(tablename, "wt");

      if (tablefile == NULL)
         error("Opening file [%s] for tabulating NPL scores\n", (const char *) tablename);

      fprintf(tablefile,"CHR\tPOS\tLABEL\tANALYSIS\tZSCORE\tDELTA\tLOD\tPVALUE%s\n",
                         nplExponential ? "\tExDELTA\tExLOD\tPVALUE" : "");
      }

   if (info.perFamily)
      {
      filename = prefix + ".lod";

      file = fopen((const char *) filename, "wt");

      if (file == NULL)
         error("Opening file [%s] for storing NPL scores for each family\n",
               (const char *) filename);

      fprintf(file, "%10s %15s %10s %10s %10s %10s %10s\n",
              "FAMILY", "TRAIT", "LOCATION", "Z-SCORE", "pLOD", "DELTA", "LOD");
      }

   if (rawOutput)
      {
      rawOutputFilename = prefix + ".zscore";

      rawOutputFile = fopen((const char *) rawOutputFilename, "wt");

      if (rawOutputFile == NULL)
         error("Opening file [%s] for storing nonparametric Z-scores\n",
               (const char *) rawOutputFilename);
      }
   }

void KongAndCox::OutputLabel(int chr, double position,
                             const char * method, const char * poslabel, double Z)
   {
   double pvalue = ndist(Z);
   int digits = pvalue < 1.5e-4 ? 5 : int(log(pvalue * 0.06666666)*-0.4343);

   printf("%15s %7.2f %7.*f ", poslabel, Z, digits, pvalue);

   if (tablefile != NULL)
      if (chr == 0 && position == _NAN_)
         fprintf(tablefile, "na\tna\t%s\t%s\t%.3f", poslabel, method, Z);
      else
         fprintf(tablefile, "%d\t%.3f\t%s\t%s\t%.3f", chr, position * 100., poslabel, method, Z);
   }

void KongAndCox::OutputLOD(double delta, double chisq)
   {
   double LOD    = chisq * 0.2171472409516;
   double pvalue = chisq > 0 ? chidist(chisq, 1.0) * 0.5 : 0.5;

   if (delta < 0)
      {
      pvalue = 1.0 - pvalue;
      LOD = -LOD;
      }

   int ldigits = fabs(LOD) < 100 ? 2 : fabs(LOD) < 1000 ? 1 : 0;
   int digits = pvalue < 1.5e-4 ? 5 : int(log(pvalue*0.06666666)*-0.4343);

   printf("%8.3f %6.*f %7.*f ", delta, ldigits, LOD, digits, pvalue);

   if (tablefile != NULL)
      fprintf(tablefile, "\t%.3f\t%.3f\t%.4g", delta, LOD, pvalue);
   }

void KongAndCox::CloseFiles()
   {
   if (tablefile != NULL)
      {
      fclose(tablefile);
      printf("NPL scores tabulated in [%s]\n", (const char *) tablename);
      }

   if (file != NULL)
      {
      fclose(file);
      printf("NPL scores for individual families stored in [%s]\n",
             (const char *) filename);
      }

   if (rawOutputFile != NULL)
      {
      fclose(rawOutputFile);
      printf("Nonparametric Z-scores for individual families stored in [%s]\n",
             (const char *) rawOutputFilename);
      }
   }

const char * KongAndCox::TaskDescription()
   {
   return "Scoring Nonparametric Statistics";
   }

void KongAndCox::OutputPerFamily(Pedigree & ped, const char * label, double delta)
   {
   int stat = lm.stat;
   int pos = lm.pos;

   for (int i = 0; i < families; i++)
      if (lm.min[stat][i] != lm.max[stat][i])
         {
         // Calculate delta maximized within this family
         double delta_i = lm.scores[stat][i][pos] < 0 ?
            -1.0 / lm.max[stat][i] : -1.0 / lm.min[stat][i];

         // k = log10(e)
         const double k = 0.4342944819032;

         // Calculate LOD score maximized within this family
         double lod_i = k * log(1.0 + delta_i * lm.scores[stat][i][pos]);

         // Calculate this families contribution to overall LOD
         double partial = k * log(1.0 + delta * lm.scores[stat][i][pos]);

         fprintf(file, "%10s %15s %10s %10f %10f %10f %10f\n",
            (const char *) ped.families[i]->famid,
            (const char *) pheno[stat], label,
            lm.scores[stat][i][pos], (delta >= 0) ? partial : -partial,
            delta_i, (delta_i >= 0) ? lod_i : -lod_i);
         }
   }

void KongAndCox::OutputRawScores(Pedigree & ped, int chromosome, StringArray & labels)
   {
   fprintf(rawOutputFile, "CHROMOSOME %d\n", chromosome);

   fprintf(rawOutputFile, "POSITIONS ");
   for (int i = 0; i < positions; i++)
      fprintf(rawOutputFile, "%s ", (const char *) labels[i]);
   fprintf(rawOutputFile, "\n");

   for (int i = 0; i < statistics; i++)
      {
      fprintf(rawOutputFile, "ANALYSIS %s\n", (const char *) phenoTag[i]);

      for (int j = 0; j < families; j++)
         if (lm.min[i][j] != lm.max[i][j])
            {
            fprintf(rawOutputFile, "FAMILY %s ZMIN %.6f ZMAX %.6f SCORES ",
                    (const char *) ped.families[j]->famid,
                    lm.min[i][j], lm.max[i][j]);

            for (int k = 0; k < positions; k++)
               fprintf(rawOutputFile, "%.6f ", lm.scores[i][j][k]);

            fprintf(rawOutputFile, "\n");
            }
      }
   }
 
