////////////////////////////////////////////////////////////////////// 
// regress/RegressAnalysis.cpp 
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
 
#include "RegressAnalysis.h"
#include "AutoFit.h"

// Trait model parameters
//

double RegressionAnalysis::mean = 0.0;
double RegressionAnalysis::variance = 1.0;
double RegressionAnalysis::heritability = 0.5;
double RegressionAnalysis::testRetestCorrel = 1.0;
String RegressionAnalysis::modelsFile("models.tbl");
bool   RegressionAnalysis::autoFit = false;

// Parameters for model fitting engine
//

bool   RegressionAnalysis::unrestricted = false;

// Output flags
//

bool RegressionAnalysis::rankFamilies = false;
bool RegressionAnalysis::writePDF = false;

// Constructor and destructor
//

RegressionAnalysis::RegressionAnalysis(Pedigree & pedigree)
   : MerlinCore(pedigree), pdf(analysisPositions), kinship(mantra)
   {
   regress = NULL;
   }

RegressionAnalysis::~RegressionAnalysis()
   {
   if (regress != NULL)
      delete [] regress;
   }

// Setup global storage
//

bool RegressionAnalysis::ReadModelsFromFile()
   {
   StringArray models;
   models.Read(modelsFile);

   if (models.Length() == 0)
      return false;

   regress = new FancyRegression[models.Length()];

   printf("Retrieving analysis models from file [%s]...\n",
          (const char *) modelsFile);

   modelCount = 0;

   StringArray tokens;
   for (int i = 0, line = 0; i < models.Length(); i++)
      {
      models[i].Trim();

      // Skip comments
      if (models[i][0] == '#') continue;

      // Divide each line into tokens
      tokens.Clear();
      tokens.AddTokens(models[i]);

      // Skip blank lines
      if (tokens.Length() == 0) continue;

      // Print message for tracing...
      printf("   Input: %s\n", (const char *) models[i], line++);

      // Need a minimum of four tokens per line
      if (tokens.Length() < 4)
         {
         printf(" Skipped: Trait name, mean, variance and heritability required.\n");
         continue;
         }

      regress[modelCount].trait = ped.LookupTrait(tokens[0]);

      if (regress[modelCount].trait < 0)
         {
         printf(line == 1 ? " Skipped: Appears to be a header line\n" :
                            " Skipped: Trait %s not listed in the data file\n",
                            (const char *) tokens[0]);
         continue;
         }

      // First check that mean, variance and heritability are valid numbers
      bool fail = false;
      for (int j = 1; j <= 3; j++)
         {
         char * ptr = NULL;
         strtod(tokens[j], &ptr);
         fail |= ptr[0] != 0;
         }

      // If one of the values is not a valid number, skip
      if (fail)
         {
         printf(line == 1 ? " Skipped: Appears to be a header line\n" :
                            " Skipped: Invalid numeric format\n");
         continue;
         }

      regress[modelCount].mean = tokens[1];
      regress[modelCount].variance = tokens[2];
      regress[modelCount].heritability = tokens[3];

      if (tokens.Length() > 4)
         {
         regress[modelCount].label = tokens[4];

         for (int j = 5; j < tokens.Length(); j++)
            {
            regress[modelCount].label += " ";
            regress[modelCount].label += tokens[j];
            }
         }
      else
         regress[modelCount].label.printf("Model %d", modelCount + 1);

      regress[modelCount].shortLabel = regress[modelCount].label;
      regress[modelCount].testRetestCorrel = testRetestCorrel;
      regress[modelCount].bounded = !unrestricted;

      printf("        Model loaded and labelled %s\n", (const char *) regress[modelCount].label);

      modelCount++;
      }

   if (modelCount == 0)
      {
      printf("No valid models, default model will be used\n\n");
      return false;
      }

   printf("Table processed. %d models recognized\n\n", modelCount);

   return true;
   }

void RegressionAnalysis::SetupDefaultModels()
   {
   modelCount = ped.traitCount;

   regress = new FancyRegression[modelCount];

   for (int i = 0; i < modelCount; i++)
      {
      regress[i].shortLabel = ped.traitNames[i];
      regress[i].label.printf("Trait: %s", (const char *) ped.traitNames[i]);
      regress[i].trait = i;
      regress[i].mean = mean;
      regress[i].variance = variance;
      regress[i].heritability = heritability;
      regress[i].testRetestCorrel = testRetestCorrel;
      regress[i].bounded = !unrestricted;
      }
   }

void RegressionAnalysis::AutoFitModels()
   {
   // Fit trait models
   AutoFit engine;
   engine.FitPolygenicModels(ped);

   // Allocate models, up to one per trait
   modelCount = 0;
   regress = new FancyRegression[ped.traitCount];

   for (int i = 0; i < ped.traitCount; i++)
      if (engine.means[0] != _NAN_)
         {
         regress[modelCount].shortLabel = ped.traitNames[i];
         regress[modelCount].label.printf("Trait: %s", (const char *) ped.traitNames[i]);
         regress[modelCount].trait = i;
         regress[modelCount].mean = engine.means[i];
         regress[modelCount].variance = engine.variances[i];
         regress[modelCount].heritability = engine.heritabilities[i];
         regress[modelCount].testRetestCorrel = 1.0;
         modelCount++;
         }

   printf("\n");
   }

void RegressionAnalysis::SetupGlobals()
   {
   MerlinCore::SetupGlobals();

   if (regress != NULL)
      delete [] regress;

   if (autoFit)
      AutoFitModels();
   else if (!ReadModelsFromFile())
      SetupDefaultModels();

   if (rankFamilies)
      RankFamilies();

   betweenMarkers = scanChromosome = multipoint = true;

   if (writePDF)
      pdf.OpenFile("merlin-regress.pdf");
   }

void RegressionAnalysis::RankFamilies()
   {
   printf("Family Informativeness\n"
       "======================\n"
       "%15s %15s %7s %7s %7s %7s %7s\n",
       "Family", "Trait", "People", "Phenos", "Pairs", "Info", "ELOD20");

   IntArray totalCount(modelCount); totalCount.Zero();
   IntArray totalPheno(modelCount); totalPheno.Zero();
   IntArray totalPairs(modelCount); totalPairs.Zero();

   for (int trait = 0; trait < modelCount; trait++)
      regress[trait].ResetMaxInfo();

   for (int i = 0; i < ped.familyCount; i++)
      {
      Family * family = ped.families[i];

      mantra.Prepare(ped, *family);

      if (mantra.bit_count > maxBits)
         {
         printf("%15s *** %d-bit family skipped ***\n",
            (const char *) family->famid, mantra.bit_count);
         continue;
         }

      mantra.PrepareIBD();
      tree.MakeMinimalTree(1.0, mantra.bit_count);
      kinship.SelectFamily(&ped, family);
      kinship.Calculate(tree);

      for (int trait = 0; trait < modelCount; trait++)
         {
         regress[trait].SetupFamily(kinship);

         int phenotypes = 0;
         for (int j = family->first; j <= family->last; j++)
            if (ped[j].isPhenotyped(regress[trait].trait))
               phenotypes++;

         totalCount[trait] += family->count;
         totalPheno[trait] += phenotypes;
         totalPairs[trait] += regress[trait].CountPairs();

         printf("%15s %15s %7d %7d %7d %7.3f %7.3f\n",
            (const char *) family->famid,
            (const char *) regress[trait].shortLabel,
            family->count, phenotypes, regress[trait].CountPairs(),
            regress[trait].maxInfo,
            regress[trait].maxInfo * 0.04 * 0.2171472409516);
         }
      }

   for (int trait = 0; trait < modelCount; trait++)
      printf("%15s %15s %7d %7d %7d %7.3f %7.3f\n",
         "TOTAL", (const char *) regress[trait].shortLabel,
         totalCount[trait], totalPheno[trait], totalPairs[trait],
         regress[trait].sumMaxInfo,
         regress[trait].sumMaxInfo * 0.04 * 0.2171472409516);

   printf("\n");
   }

// Setup storage that depends on map-length
//

int RegressionAnalysis::SetupMap(int chromosome)
   {
   chromosome = MerlinCore::SetupMap(chromosome);

   // Setup chromosome labels for PDF output
   if (writePDF)
      pdf.SelectChromosome(ped.GetMarkerInfo(markers[0])->chromosome);

   for (int i = 0; i < modelCount; i++)
      {
      regress[i].SetPositionCount(analysisPositions.Length());
      regress[i].ResetMaxInfo();
      }

   return chromosome;
   }

// Setup analysis for each family
//

bool RegressionAnalysis::SelectFamily(Family * f, bool warnOnSkip)
   {
   if (MerlinCore::SelectFamily(f, warnOnSkip))
      {
      kinship.SelectFamily(&ped, f);

      ProgressReport("Preparing Matrices", 0, modelCount + 1);
      tree.MakeMinimalTree(1.0, mantra.bit_count);
      kinship.Calculate(tree);

      for (int i = 0; i < modelCount; i++)
         {
         ProgressReport("Preparing Matrices", i+1, modelCount + 1);
         regress[i].SetupFamily(kinship);
         }

      return true;
      }

   return false;
   }

// Carry out analysis for a family
//

void RegressionAnalysis::AnalyseLocation(int pos, Tree & inheritance)
   {
   kinship.Calculate(inheritance);

   for (int i = 0; i < modelCount; i++)
      regress[i].Analyse(kinship, pos);
   }

// Put it all together
//

void RegressionAnalysis::Analyse()
   {
   if (analysisPositions.Length() == 0)
      error("The list of positions to analyze is empty.");

   singlepoint = NULL;
   right = NULL;

   try
      {
      BasicTree::FreeSwap();

      ScoreSinglepoint();

      if (informativeCount != 0)
         {
         if (ScoreConditionals())
            ScanChromosome();
         else
            PrintMessage("  SKIPPED: Requires impossible recombination pattern");

         if (!twopoint) delete [] right;
         }

      CleanMessages();
      delete [] singlepoint;
      }
   catch (const TreesTooBig & problem)
      {
      PrintMessage("  SKIPPED: At least %d megabytes needed", problem.memory_request);
      CleanMessages();

      if (singlepoint != NULL) delete [] singlepoint;
      if (right != NULL) delete [] right;
      }
   catch (const OutOfTime & problem)
      {
      PrintMessage("  SKIPPED: Analysis would require more than %d minutes\n", maxMinutes);
      CleanMessages();

      if (singlepoint != NULL) delete [] singlepoint;
      if (right != NULL) delete [] right;
      };
   };

// Output LOD scores at each location
//

void RegressionAnalysis::PrintScores()
   {
   int chr = 0;
   String tablename;
   FILE * tablefile = NULL;

   if (MerlinCore::tabulate)
      {
      chr = markers.Length() > 0 ? Pedigree::GetMarkerInfo(markers[0])->chromosome : 0;
      tablename.printf("%s-regress-chr%02d.tbl", (const char *) MerlinCore::filePrefix, chr > 0 ? chr : 0);
      tablefile = fopen(tablename, "wt");
      }

   if (tablefile != NULL)
      fprintf(tablefile, "CHR\tPOS\tPHENOTYPE\tH2\tSD\tINFO\tLOD\tPVALUE\n");

   for (int i = 0; i < modelCount; i++)
      regress[i].PrintScores(tablefile, chr, pdf, labels);

   if (tablefile != NULL)
      {
      printf("Linkage results stored in file [%s]\n\n", (const char *) tablename);

      fclose(tablefile);
      tablefile = NULL;
      }
   }

 
