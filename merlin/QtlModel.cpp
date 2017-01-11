////////////////////////////////////////////////////////////////////// 
// merlin/QtlModel.cpp 
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
 
#include "QtlModel.h"
#include "StringArray.h"
#include "Pedigree.h"

QtlModel::QtlModel()
   {
   modelCount = 0;
   covariateList = NULL;
   }

QtlModel::~QtlModel()
   {
   if (covariateList != NULL)
      {
      delete [] covariateList;
      delete [] candidateList;
      }
   }

void QtlModel::LoadFromFile(const char * filename)
   {
   IFILE f = ifopen(filename, "rb");

   if (f == NULL)
      return;

   LoadFromFile(f);

   ifclose(f);
   }

void QtlModel::LoadFromFile(IFILE & f)
   {
   // Free previously allocated model
   if (covariateList != NULL)
      delete [] covariateList;

   // Read input in one fell swoop
   StringArray input, tokens;
   input.Read(f);

   // Count the maximum number of models we are dealing with
   int maxModels = 0;
   for (int i = 0; i < input.Length(); i++)
      {
      input[i].Trim();

      if (input[i].Left(5).SlowCompare("TRAIT") == 0)
         maxModels++;
      }
   if (maxModels == 0) maxModels = 1;

   // Initialize
   modelCount = 0;
   covariateList = new IntArray[maxModels];
   candidateList = new IntArray[maxModels];

   int trait = -1;
   bool warned = false;

   // Print a message
   printf("Parsing file with QTL models ...\n");
   printf("   There can be up to three lines per model, each starting with\n"
          "   the keyword TRAIT, COVARIATES, or CANDIDATES.\n");

   // Loop through models and update stuff
   for (int i = 0; i < input.Length(); i++)
      {
      bool valid = false;

      // The first few lines ensure we skip comments and empty lines
      input[i].Trim();

      if (input[i][0] == '#' || input[i][0] == '/' && input[i][1] == '/' || input[i].IsEmpty())
         continue;

      tokens.ReplaceTokens(input[i]);

      // First check whether we have new trait information
      if (tokens[0].SlowCompare("TRAIT") == 0 && tokens.Length() > 1)
         {
         trait = Pedigree::LookupTrait(tokens[1]);

         if (trait != -1)
            {
            tokens.Delete(0);
            tokens.Delete(0);

            traitList.Push(trait);
            modelCount++;

            valid = true;
            }
         }

      if (tokens.Length() == 0)
         continue;

      // Next, check whether we have covariate information
      if (tokens[0].SlowCompare("COVARIATES") == 0 && tokens.Length() > 1 && trait != -1)
         {
         tokens.Delete(0);

         while (tokens.Length() && tokens[0].SlowCompare("CANDIDATES") != 0)
            {
            int covariate = Pedigree::LookupCovariate(tokens[0]);

            if (covariate >= 0)
               covariateList[modelCount - 1].Push(covariate);

            tokens.Delete(0);
            }

         valid = true;
         }

      if (tokens.Length() == 0)
         continue;

      // Finally, check whether we have candidate information
      if (tokens[0].SlowCompare("CANDIDATES") == 0 && tokens.Length() > 1 && trait != -1)
         {
         while (tokens.Length() > 1)
            {
            int candidate = Pedigree::LookupMarker(tokens[1]);

            if (candidate >= 0)
               candidateList[modelCount - 1].Push(candidate);

            tokens.Delete(1);
            }

         valid = true;
         }

      if (valid || warned) continue;

      warned = true;

      printf("  WARNING: I am having trouble parsing your model file.\n"
             "           Problem entries will be ignored!\n");
      }

   for (int i = 0; i < modelCount; i++)
      candidateList[i].Sort();

   printf("%d models for QTL analysis parsed\n\n", modelCount);
   }

void QtlModel::PrintSummary()
   {
   if (modelCount == 0)
      {
      printf("No custom models loaded\n");
      return;
      }

   int covariates = 0;
   int candidates = 0;

   for (int i = 0; i < modelCount; i++)
      {
      covariates += covariateList[i].Length();
      candidates += candidateList[i].Length();
      }

   printf("%d custom quantitative trait models loaded, the average model lists...\n", modelCount);
   printf("   %.2f covariates\n", covariates / (double) modelCount);
   printf("   %.2f candidate SNPs for association analysis\n", candidates / (double) modelCount);

   if (modelCount > 3)
      {
      printf("\n");
      return;
      }

   printf("\n");

   for (int i = 0; i < modelCount; i++)
      {
      printf("CUSTOM QUANTITATIVE TRAIT MODEL #%d\n", i+1);
      printf("===================================\n\n");

      printf("TRAIT: %s\n", (const char *) Pedigree::traitNames[traitList[i]]);

      if (covariateList[i].Length())
         {
         printf("  COVARIATES: ");

         int togo = 60;

         for (int j = 0; j < covariateList[i].Length(); j++)
            {
            togo -= Pedigree::covariateNames[covariateList[i][j]].Length() + 2;

            if (togo < 0)
               {
               printf("...");
               break;
               }

            printf("%s%s", (const char *) Pedigree::covariateNames[covariateList[i][j]],
                           j == covariateList[i].Length() - 1 ? "" : ", ");
            }

         printf("\n");
         }
      else
         printf("  No covariates\n");

      if (candidateList[i].Length())
         {
         printf("  CANDIDATE SNPS: ");

         int togo = 58;

         for (int j = 0; j < candidateList[i].Length(); j++)
            {
            togo -= Pedigree::markerNames[candidateList[i][j]].Length() + 2;

            if (togo < 0)
               {
               printf("...");
               break;
               }

            printf("%s%s", (const char *) Pedigree::markerNames[candidateList[i][j]],
                           j == candidateList[i].Length() - 1 ? "" : ", ");
            }

         printf("\n");
         }

      printf("\n");
      }
   }

bool QtlModel::SkipMarker(int model, int marker)
   {
   if (candidateList[model].Length() == 0)
      return false;

   return (candidateList[model].FastFind(marker) < 0);
   }

RefinedQtlModel::RefinedQtlModel()
   {
   values = NULL;
   }

RefinedQtlModel::~RefinedQtlModel()
   {
   if (values != NULL)
      delete [] values;
   }

void RefinedQtlModel::Dimension(int count)
   {
   modelCount = count;

   labels.Clear();
   labels.Dimension(count);

   comments.Clear();
   comments.Dimension(count);

   scores.Dimension(count);
   scores.Zero();

   if (values != NULL)
      delete [] values;
   values = new Vector[count];
   }

void RefinedQtlModel::PrintDataFile(const char * filename)
   {
   FILE * f = fopen(filename, "wt");

   if (f == NULL)
      {
      printf("  WARNING: Error opening file [%s].\n", filename);
      return;
      }

   PrintDataFile(f);

   fclose(f);
   }

void RefinedQtlModel::PrintPedigreeFile(const char * filename, Pedigree & ped)
   {
   FILE * f = fopen(filename, "wt");

   if (f == NULL)
      {
      printf("  WARNING: Error opening file [%s].\n", filename);
      return;
      }

   PrintPedigreeFile(f, ped);

   fclose(f);
   }

void RefinedQtlModel::PrintModelFile(const char * filename, QtlModel & base)
   {
   FILE * f = fopen(filename, "wt");

   if (f == NULL)
      {
      printf("  WARNING: Error opening file [%s].\n", filename);
      return;
      }

   PrintModelFile(f, base);

   fclose(f);
   }

void RefinedQtlModel::PrintDataFile(FILE * f)
   {
   for (int i = 0; i < modelCount; i++)
      if (scores[i] > 0.0)
         fprintf(f, "C %s\n", (const char *) labels[i]);
   }

void RefinedQtlModel::PrintPedigreeFile(FILE * f, Pedigree & ped)
   {
   for (int i = 0; i < ped.count; i++)
         {
         fprintf(f, "%s %s %s %s %d",
                  (const char *) ped[i].famid, (const char *) ped[i].pid,
                  (const char *) ped[i].fatid, (const char *) ped[i].motid,
                  ped[i].sex);

         for (int j = 0; j < modelCount; j++)
            if (scores[j] > 0.0)
               fprintf(f, " %.3f", values[j][i]);

         fprintf(f, "\n");
         }
   }

void RefinedQtlModel::PrintModelFile(FILE * f, QtlModel & base)
   {
   for (int i = 0; i < modelCount; i++)
      if (scores[i] > 0.0)
         {
         fprintf(f, "////////////////////////////////////////////\n"
                    "// %s\n"
                    "//\n", (const char *) comments[i]);

         fprintf(f, "TRAIT %s\n", (const char *) Pedigree::traitNames[base.traitList[i]]);

         fprintf(f, "COVARIATES ");
         for (int j = 0; j < base.covariateList[i].Length(); j++)
            fprintf(f, "%s ", (const char *) Pedigree::covariateNames[base.covariateList[i][j]]);
         fprintf(f, "%s\n", (const char *) labels[i]);

         if (base.candidateList[i].Length())
            {
            fprintf(f, "CANDIDATES ");
            for (int j = 0; j < base.candidateList[i].Length(); j++)
               if (Pedigree::markerNames[base.candidateList[i][j]] != labels[i])
                  fprintf(f, "%s ", (const char *) Pedigree::markerNames[base.candidateList[i][j]]);
            fprintf(f, "\n");
            }
         }
      else
         {
         fprintf(f, "////////////////////////////////////////////\n"
                    "// No refined model available for %s ...\n"
                    "//\n", (const char *) Pedigree::traitNames[base.traitList[i]]);
         }
   }

void RefinedQtlModel::WriteFiles(const String & prefix, Pedigree & ped, QtlModel & base)
   {
   PrintDataFile(prefix + ".dat");
   PrintPedigreeFile(prefix + ".ped", ped);
   PrintModelFile(prefix + ".tbl", base);
   }

bool RefinedQtlModel::CheckForImprovement(int model, double score)
   {
   return (score > scores[model]);
   }



 
