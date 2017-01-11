////////////////////////////////////////////////////////////////////// 
// merlin/MerlinModel.cpp 
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
 
#include "MerlinModel.h"
#include "MerlinCore.h"
#include "MathStats.h"
#include "MerlinPDF.h"

#include <math.h>

String MerlinModel::modelsFile("param.tbl");

MerlinModel::MerlinModel()
   {
   modelCount = 0;
   models = NULL;
   }

MerlinModel::~MerlinModel()
   {
   if (models != NULL)
      delete [] models;
   }

void MerlinModel::LoadModels()
   {
   // Delete previously allocated models
   if (models != NULL)
      {
      delete [] models;
      labels.Clear();
      modelCount = 0;
      models = NULL;
      }

   labels.Clear();

   // Load models file into memory
   StringArray input;
   input.Read(modelsFile);

   if (input.Length() == 0)
      return;

   models = new DisModel[input.Length()];

   printf("Analysis models will be retrieved from file [%s]...\n",
         (const char *) modelsFile);

   IFILE f = ifopen(modelsFile, "rb");
   while (!ifeof(f))
      if (models[modelCount].ReadFromFile(f))
         {
         labels.Push(models[modelCount].label);
         modelCount++;
         }

   ifclose(f);

   printf("Table processed. %d models recognized\n\n", modelCount);
   }

int MerlinModel::CountModels()
   {
   return modelCount;
   }

AnalysisTask * MerlinModel::CreateTask()
   {
   return new ParametricAnalysis(models, labels);
   }

double ParametricAnalysis::scale = 1.0 / (log(10.0));

ParametricAnalysis::ParametricAnalysis(DisModel * m, StringArray & l)
   {
   modelCount = l.Length();
   labels = l;

   scores = new Vector [modelCount];
   rawScores = new Matrix [modelCount];
   impossible = new IntArray [modelCount];
   models = new DisModel [modelCount];
   likelihoods = new Tree [modelCount];

   for (int i = 0; i < modelCount; i++)
      models[i] = m[i];

   baseline.Dimension(modelCount);

   scoreFile = NULL;
   scoreTable = NULL;
   }

ParametricAnalysis::~ParametricAnalysis()
   {
   delete [] impossible;
   delete [] likelihoods;
   delete [] models;
   delete [] scores;
   delete [] rawScores;
   }

void ParametricAnalysis::Setup(AnalysisInfo & info)
   {
   for (int i = 0; i < modelCount; i++)
      {
      scores[i].Dimension(info.positions);
      scores[i].Zero();

      impossible[i].Dimension(info.positions);
      impossible[i].Zero();

      rawScores[i].Dimension(info.positions, info.families);
      rawScores[i].Set(1.0);
      }
   }

void ParametricAnalysis::SetupFamily(AnalysisInfo & info)
   {
   ParametricLikelihood tree;

   for (int i = 0; i < modelCount; i++)
      {
      tree.model = &(models[i]);
      tree.FuzzyScoreVectors(*info.m);
      likelihoods[i].Copy(tree);

      // Check that likelihood is not zero ...
      baseline[i] = tree.WeightedSum();

      if (baseline[i] == 0.0)
         error("Disease model %s is too unlikely, family %s has zero likelihood\n",
                (const char *) labels[i], (const char *) (*info.famid));

      baseline[i] = log(baseline[i]);
      }
   }

void ParametricAnalysis::AnalyseUninformative(AnalysisInfo & info, Tree & tree)
   {
   if (scoreFile != NULL)
      for (int i = 0; i < modelCount; i++)
         for (pos = 0; pos < info.positions; pos++)
            fprintf(scoreFile, "%15s %15s %15s %10.4f\n",
                    (const char *) labels[i],
                    (const char *) info.m->family->famid,
                    (const char *) (*info.labels)[pos],
                    0.00);
   }

void ParametricAnalysis::Analyse(AnalysisInfo & info, Tree & tree, int position)
   {
   if (info.lk == -1.0)
      info.lk = tree.WeightedSum();

   double offset = -log(info.lk);

   for (int i = 0; i < modelCount; i++)
      {
      double lkPos = tree.MeanProduct(likelihoods[i]);
      double lkRatio = 0.0; // Initialization avoids compiler warning

      if (lkPos <= 0.0)
         {
         impossible[i][position] = true;
         rawScores[i][position][info.famno] = 1E-50;
         }
      else
         {
         lkRatio = log(lkPos) + offset - baseline[i];
         scores[i][position] += lkRatio;
         rawScores[i][position][info.famno] = exp(lkRatio);
         }

      if (scoreFile != NULL)
         {
         fprintf(scoreFile, "%15s %15s %15s ",
                 (const char *) labels[i],
                 (const char *) (*info.famid),
                 (const char *) (*info.labels)[position]);

         if (lkPos <= 0.0)
            fprintf(scoreFile, "%10s\n", "-INF");
         else
            fprintf(scoreFile, "%10.4f\n", lkRatio * scale);
         }
      }
   }

void ParametricAnalysis::ReportFamily()
   {
   for (int i = 0; i < modelCount; i++)
      likelihoods[i].Discard();
   }

void ParametricAnalysis::Report(AnalysisInfo & info)
   {
   for (model = 0; model < modelCount; model++)
      {
      if (info.drawPDF)
         {
         info.pdf->PrepareChart(4);
         info.pdf->chart.yAxis.minMin = -10.0;
         info.pdf->chart.yAxis.minMax = 4.0;
         info.pdf->chart.SetSeriesColor(0, 0.70, 0.50, 0.50);
         info.pdf->y[0].Set(0.0);
         info.pdf->chart.SetSeriesColor(1, 0.75, 0.75, 0.75);
         info.pdf->y[1].Set(3.0);
         info.pdf->chart.SetSeriesColor(2, 0.75, 0.75, 0.75);
         info.pdf->y[2].Set(-2.0);
         info.pdf->chart.SetSeriesColor(3, 0.0, 0.0, 0.0);
         info.pdf->chart.title = "Parametric Analysis for " + labels[model];
         }

      printf("Parametric Analysis, Model %s\n", (const char *) labels[model]);
      printf("=======================================================\n");
      printf("%15s %10s %10s %10s\n", "POSITION", "LOD", "ALPHA", "HLOD");

      for (pos = 0; pos < scores[model].Length(); pos++)
         {
         EstimateAlpha();

         if (impossible[model][pos])
            {
            printf("%15s %10s %10.3f %10.3f\n",
                   (const char *) (*info.labels)[pos], "-INFINITY", alpha, hlod);

            if (scoreTable != NULL)
               fprintf(scoreTable, "%d\t%.4f\t%s\t%s\t%s\t%.4f\t%.4f\n",
                       info.chromosome, (*info.analysisPositions)[pos],
                       (const char *) (*info.labels)[pos],
                       (const char *) labels[model],
                       "-INF", alpha, hlod);

            if (info.drawPDF)
               info.pdf->y[3][pos] = -999.0;
            }
         else
            {
            double lod = scores[model][pos] * scale;

            printf("%15s %10.3f %10.3f %10.3f\n",
                   (const char *) (*info.labels)[pos], lod, alpha, hlod);

            if (scoreTable != NULL)
               fprintf(scoreTable, "%d\t%.4f\t%s\t%s\t%.4f\t%.4f\t%.4f\n",
                       info.chromosome, (*info.analysisPositions)[pos],
                       (const char *) (*info.labels)[pos],
                       (const char *) labels[model],
                       lod, alpha, hlod);

            if (info.drawPDF)
               info.pdf->y[3][pos] = lod;
            }
         }

      if (info.drawPDF)
         info.pdf->DrawChart();

      printf("\n\n");
      }
   }

void ParametricAnalysis::OpenFiles(AnalysisInfo & info, String & prefix)
   {
   if (MerlinCore::tabulate)
      {
      tablename = prefix + "-parametric.tbl";
      scoreTable = fopen(tablename, "wt");

      if (scoreTable == NULL)
         error("Could not open table [%s] for parametric analysis results\n",
              (const char *) tablename);

      fprintf(scoreTable, "CHR\tPOS\tLABEL\tMODEL\tLOD\tALPHA\tHLOD\n");
      }

   if (info.perFamily == false)
      return;

   filename = prefix + ".par";
   scoreFile = fopen(filename, "wt");

   if (scoreFile == NULL)
      error("Could not open file [%s] for storing results of parametric analysis\n",
            (const char *) scoreFile);

   fprintf(scoreFile, "%15s %15s %15s %10s\n", "MODEL", "FAMILY", "POSITION", "LOD");
   }

void ParametricAnalysis::CloseFiles()
   {
   if (scoreTable != NULL)
      {
      printf("Parametric LOD scores tabulated in [%s]\n", (const char *) tablename);

      fclose(scoreTable);
      scoreTable = NULL;
      }

   if (scoreFile == NULL)
      return;

   printf("Parametric LOD scores for individual families stored in file [%s].\n",
          (const char *) filename);

   fclose(scoreFile);
   scoreFile = NULL;
   }

void ParametricAnalysis::EstimateAlpha()
   {
   // for (alpha = 0.00; alpha < 1.00; alpha += 0.05)
   //   printf("%7.3f %7.3f\n", alpha, - f(alpha) * scale);

   // Setup range for minimization
   a = 0.0;
   c = 1.0;

   // As a starting point, assume that most families are unlinked
   b = 1.0;
   fb = f(b);

   // Call optimization routine
   Brent(0.0001);

   if (min < 1E-5 && fmin > 0.00 && fmin < 0.0001)
      alpha = hlod = 0.0;
   else
      {
      alpha = min;
      hlod = -fmin * scale;
      }
   }

double ParametricAnalysis::f(double alpha)
   {
   double complement = 1.0 - alpha;
   double score = 0.0;

   for (int j = 0; j < rawScores[model].cols; j++)
      score += log(complement + alpha * rawScores[model][pos][j]);

   return -score;
   }

const char * ParametricAnalysis::TaskDescription()
   {
   return ("Scoring Parametric Models");
   }
 
