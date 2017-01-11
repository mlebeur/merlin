////////////////////////////////////////////////////////////////////// 
// merlin/MerlinFamily.cpp 
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
 
#include "MerlinHaplotype.h"
#include "MerlinCluster.h"
#include "MerlinFamily.h"
#include "MerlinModel.h"
#include "AssociationAnalysis.h"
#include "FastAssociation.h"
#include "InformationContent.h"
#include "MapFunction.h"
#include "MathConstant.h"
#include "KongAndCox.h"
#include "Houdini.h"

// General analyses
bool FamilyAnalysis::estimateIBD = false;
bool FamilyAnalysis::estimateKinship = false;
bool FamilyAnalysis::estimateKinship15 = false;
bool FamilyAnalysis::estimateInformation = false;
bool FamilyAnalysis::estimateMatrices = false;
bool FamilyAnalysis::selectCases = false;
bool FamilyAnalysis::storeKinshipForVc = false;
bool FamilyAnalysis::storeKinshipForAssoc = false;
bool FamilyAnalysis::fastAssociationAnalysis = false;
bool FamilyAnalysis::calcLikelihood = false;

// Haplotyping
bool FamilyAnalysis::bestHaplotype = false;
int  FamilyAnalysis::sampledHaplotypes = 0;
bool FamilyAnalysis::allHaplotypes = false;

// Output options
bool FamilyAnalysis::simwalk2 = false;
bool FamilyAnalysis::perFamily = false;
bool FamilyAnalysis::writePDF = false;

// Error detection and genotype inference
bool FamilyAnalysis::inferGenotypes = false;
bool FamilyAnalysis::outputInferredGenotypes = false;
bool FamilyAnalysis::findErrors = false;

FamilyAnalysis::FamilyAnalysis(Pedigree & pedigree) :
   MerlinCore(pedigree), matrix(mantra),
   kinship(mantra), kinship15(mantra),
   ibd(mantra), hybrid(*this),
   pdf(analysisPositions)
   {
   errorfile = NULL;
   taskList = NULL;
   }

FamilyAnalysis::~FamilyAnalysis()
   {
   // Delete generic analysis tasks
   AnalysisTask * current;

   while (taskList != NULL)
      {
      current = taskList;
      taskList = taskList->next;
      if (current != &kinship) delete current;
      }
   }

bool FamilyAnalysis::SelectFamily(Family * f, bool warnOnSkip)
   {
   if (MerlinCore::SelectFamily(f, warnOnSkip))
      // If analysis for this family is proceeding ...
      {
      if (estimateMatrices) matrix.SelectFamily(&ped, family);
      if (estimateIBD) ibd.SelectFamily(&ped, family);
      if (estimateKinship15) kinship15.SelectFamily(&ped, family);

      return true;
      }

   // Tell simwalk2 we have skipped this family
   if (simwalk2) hybrid.Skip();

   return false;
   }

void FamilyAnalysis::SetupGlobals()
   {
   MerlinCore::SetupGlobals();

   // Information content calculations requested?
   if (KongAndCox::CalculateNPL(ped))
      NewTask(kac = new KongAndCox(ped));

   // Parametric analyses requested?
   MerlinModel parametricModels;

   parametricModels.LoadModels();

   if (parametricModels.CountModels())
      NewTask(parametricModels.CreateTask());

   // Information Content Calculations requested?
   if (estimateInformation)
      NewTask(new InformationContent);

   // Record how kinship information will be used
   if (estimateKinship || selectCases || storeKinshipForVc || storeKinshipForAssoc)
      {
      kinship.write = estimateKinship;
      kinship.store = storeKinshipForVc || storeKinshipForAssoc;
      kinship.selectCases = selectCases;

      NewTask(&kinship);
      }

   // Carry out between marker analyses?
   betweenMarkers = estimateIBD || estimateKinship15 || estimateMatrices ||
      taskList != NULL;

   // Carry out analyses at each marker?
   scanChromosome = betweenMarkers || findErrors || inferGenotypes;

   // Perform right conditional calculations?
   multipoint     = scanChromosome || bestHaplotype && zeroRecombination ||
                    (sampledHaplotypes != 0) || allHaplotypes || calcLikelihood;

   // Stop right conditional calculations early?
   lazyRightConditional = !(bestHaplotype && zeroRecombination || allHaplotypes ||
                            sampledHaplotypes || calcLikelihood);

   // Merge gene flow outcomes with identical likelihoods, but possibly
   // implying different allelic states
   InheritanceTree::mergingStrategy = !allHaplotypes ? MERGE_ALL : MERGE_ZEROS;

   // Check that we have the correct type of data for simwalk2 analyses
   if (simwalk2 && ped.affectionCount != 1)
      error("Joint analysis with SimWalk2 requires pedigree with exactly one discrete trait\n"
            "The current pedigree includes %d discrete traits\n", ped.traitCount);

   // Information for generic tasks
   taskInfo.m = &mantra;
   taskInfo.drawPDF = writePDF;
   taskInfo.pdf = &pdf;
   taskInfo.perFamily = perFamily;
   }

int FamilyAnalysis::SetupMap(int chromosome)
   {
   int next_chromosome = MerlinCore::SetupMap(chromosome);

   // Setup chromosome labels for PDF output
   if (writePDF)
#ifdef __CHROMOSOME_X__
      pdf.SelectChromosomeX();
#else
      pdf.SelectChromosome(ped.GetMarkerInfo(markers[0])->chromosome);
#endif

   // Setup genotype inference engine
   if (inferGenotypes)
      {
      imputationEngine.SelectMarkers(markers, markerPositions, analysisPositions);
      imputationEngine.AllocateMemory(ped);
      }

   // Initialize overall likelihood for this chromosome
   lkSum = 0.0;
   lkCount = 0;

   // Chromosome specific setup information for generic tasks
   taskInfo.chromosome = ped.GetMarkerInfo(markers[0])->chromosome;
   taskInfo.positions = analysisPositions.Length();
   taskInfo.analysisPositions = &analysisPositions;
   taskInfo.families = ped.familyCount;
   taskInfo.labels = &labels;

   for (AnalysisTask * task = taskList; task != NULL; task = task->next)
      task->Setup(taskInfo);

   return next_chromosome;
   }

void FamilyAnalysis::Analyse()
   {
   if (analysisPositions.Length() == 0)
      error("List of positions to analyze is empty\n");

   singlepoint = NULL;
   right = NULL;

   // Update task information
   taskInfo.famno = mantra.family->serial;
   taskInfo.famid = &(mantra.family->famid);

   try
      {
      FreeSwap();

      ScoreSinglepoint();

      if (informativeCount == 0)
         SkipFamily();
      else
         {
         // Setup for generic tasks
         for (AnalysisTask * task = taskList; task != NULL; task = task->next)
            {
            ProgressReport(task->TaskDescription());
            task->SetupFamily(taskInfo);
            }

         if (multipoint)
            {
            if (ScoreConditionals())
               {
               if (scanChromosome) ScanChromosome();
               if (calcLikelihood && perFamily) PrintMessage("  lnLikelihood = %.3f    ", likelihood);
               if (simwalk2) hybrid.Output();
               if (allHaplotypes) haplo.All(*this);
               for (int i = 0; i < sampledHaplotypes; i++) haplo.Sample(*this);
               if (bestHaplotype && zeroRecombination) haplo.MostLikely(*this);
               lkSum += likelihood;
               lkCount ++;
               }
            else
               {
               PrintMessage("  SKIPPED: Requires impossible recombination pattern");
               AbortAnalysis();
               }
            if (!twopoint)
               {
               delete [] right;
               right = NULL;
               }
            }
         if (bestHaplotype && !zeroRecombination) haplo.MostLikely(*this);
         }

      CleanMessages();
      delete [] singlepoint;
      }
   catch (const TreesTooBig & problem)
      {
      PrintMessage("  SKIPPED: At least %d megabytes needed", problem.memory_request);
      CleanMessages();

      FreeMemory();
      AbortAnalysis();
      }
   catch (const OutOfTime & problem)
      {
      PrintMessage("  SKIPPED: Analysis would require more than %d minutes\n", maxMinutes);
      CleanMessages();

      FreeMemory();
      AbortAnalysis();
      };
   };

void FamilyAnalysis::FreeMemory()
   {
   if (singlepoint != NULL) delete [] singlepoint;
   if (right != NULL) delete [] right;

   right = NULL;
   singlepoint = NULL;
   }

void FamilyAnalysis::AbortAnalysis()
   {
   for (AnalysisTask * task = taskList; task != NULL; task = task->next)
      task->SkipFamily(taskInfo);

   if (simwalk2) nplInv[0] = nplInv[1] = 0.0;
   };

void FamilyAnalysis::SkipFamily()
   {
   // Print out overall gene-flow likelihood
   if (calcLikelihood && perFamily)
      PrintMessage("  lnLikelihood = %.3f", likelihood);
      
   // And update overall likelihood for the sample
   lkSum += likelihood;
   lkCount++;

   // If there is no information on this family
   // IBD and kinship coefficients / matrices are the same at all markers
   Tree dummy;
   dummy.MakeMinimalTree(1.0, mantra.bit_count);

   for (AnalysisTask * task = taskList; task != NULL; task = task->next)
      task->AnalyseUninformative(taskInfo, dummy);

   if (estimateIBD || estimateKinship15 || estimateMatrices)
      {
      double sum = 0.0; // Initialization avoids compiler warning

      double scale = pow(2.0, mantra.bit_count);

      if (estimateIBD) sum = ibd.Calculate(dummy, labels[0]) * scale;
      if (estimateMatrices ) sum = matrix.Calculate(dummy, labels[0]) * scale;
      if (estimateKinship15) sum = kinship15.Calculate(dummy, labels[0]) * scale;

      for (int i = 1; i < labels.Length(); i++)
         {
         if (estimateIBD) ibd.Output(labels[i], sum);
         if (estimateKinship15) kinship15.Output(labels[i], sum);
         if (estimateMatrices) matrix.Output(labels[i], sum);
         }
      }

   // Output dummy haplotype -- something to look at -- keeps people happy
   if (allHaplotypes || bestHaplotype)
      haplo.UninformativeFamily(*this, false);

   for (int i = 0; i < sampledHaplotypes; i++)
      haplo.UninformativeFamily(*this, true);

   // And tell Simwalk2 this family isn't very interesting (to us!)
   if (simwalk2) hybrid.Skip();
   }

void FamilyAnalysis::AnalyseLocation(int pos, Tree & inheritance)
   {
   // The probability of the data summed all possible inheritance vectors
   double & lk = taskInfo.lk = -1.0;

   if (estimateKinship15) kinship15.Calculate(inheritance, labels[pos]);
   if (estimateMatrices) lk = matrix.Calculate(inheritance, labels[pos]);
   if (estimateIBD) lk = ibd.Calculate(inheritance, labels[pos]);

   // Execute generic analysis hooks
   for (AnalysisTask * task = taskList; task != NULL; task = task->next)
      task->Analyse(taskInfo, inheritance, pos);
   }

void FamilyAnalysis::GenotypeAnalysisHook
   (Tree & withMarker, Tree & withoutMarker, int informativeMarker)
   {
   if (findErrors)
      ErrorCheck(withMarker, withoutMarker, informativeMarker);

   if (inferGenotypes)
      imputationEngine.InferGenotypes(
         mantra, stats, withMarker, withoutMarker, singlepoint[informativeMarker],
         markers[informativeMarkers[informativeMarker]]);
   }

void FamilyAnalysis::ErrorCheck
     (Tree & withMarker, Tree & without, int informativeMarker)
   {
   int marker = markers[informativeMarkers[informativeMarker]];
   MarkerCluster * cluster = clusters.Enabled() ? clusters.markerToCluster[marker] : NULL;

   if (cluster == NULL && informativeCount <= 1) return;

   FuzzyInheritanceTree alternative;

   double multipoint_likelihood = stats.GetMean(withMarker);

   int markers = cluster == NULL ? 1 : cluster->markerIds.Length();
   for (int m = 0; m < markers; m++)
      {
      double singlepoint_likelihood;
      if (cluster != NULL)
         {
         marker = cluster->markerIds[m];
         mantra.SelectMarker(marker);
         alternative.FuzzyScoreVectors(mantra);

         singlepoint_likelihood = stats.GetMean(alternative);
         }
      else
         singlepoint_likelihood = stats.GetMean(singlepoint[informativeMarker]);

      for (int error = family->first; error <= family->last; error++)
         if (ped[error].markers[marker].isKnown())
            {
            Alleles saved_genotype = ped[error].markers[marker];
            ped[error].markers[marker][0] = ped[error].markers[marker][1] = 0;

            mantra.SelectMarker(marker);
            alternative.FuzzyScoreVectors(mantra);

            double alternative_likelihood = stats.GetMean(alternative);

            if (alternative_likelihood != singlepoint_likelihood)
               {
               if (cluster != NULL)
                  cluster->ScoreLikelihood(mantra, alternative);

               double alternative_multipoint = alternative.MeanProduct(without);

               double score =
                  (multipoint_likelihood / alternative_multipoint) /
                  (singlepoint_likelihood / alternative_likelihood);

               if (cluster != NULL)
                  score *= exp(singlepoint[informativeMarker].logOffset -
                               alternative.logOffset);

               if (score < 0.025)
                  {
                  PrintMessage("  %s genotype for individual %s is unlikely [%.3e]",
                              (const char *) ped.markerNames[marker],
                              (const char *) ped[error].pid, score);
                  fprintf(errorfile, "%10s %10s %10s %#10.3g\n",
                        (const char *) family->famid,
                        (const char *) ped[error].pid,
                        (const char *) ped.markerNames[marker], score);
                  }
               }

            ped[error].markers[marker] = saved_genotype;
            }
      }
   }

void FamilyAnalysis::ShowLODs()
   {
   if (storeKinshipForVc)
      vc.Analyse(*this);

   if (storeKinshipForAssoc)
      {
      AssociationAnalysis engine;

      if (perFamily) engine.OpenPerFamilyFile(taskInfo.chromosome);
      engine.AnalyseAssociation(*this);
      if (perFamily) engine.ClosePerFamilyFile(taskInfo.chromosome);
      }

   if (fastAssociationAnalysis)
      {
      FastAssociationAnalysis engine;

      engine.AnalyseAssociation(*this);
      }

   // Generic analysis hooks
   for (AnalysisTask * task = taskList; task != NULL; task = task->next)
      task->Report(taskInfo);

   // Genotype inference
   if (outputInferredGenotypes)
      imputationEngine.OutputGenotypes(ped, markers);
   }

void FamilyAnalysis::OpenErrorFile()
   {
   errorfile = fopen((const char *) (MerlinCore::filePrefix + ".err"), "wt");
   if (errorfile == NULL)
      error("Can't open error file [%s.err]\n", (const char *) MerlinCore::filePrefix);
   fprintf(errorfile, "%10s %10s %10s %10s\n",
      "FAMILY", "PERSON", "MARKER", "RATIO");
   }

void FamilyAnalysis::SetupFiles()
   {
   if (findErrors) OpenErrorFile();
   if (perFamily && storeKinshipForVc) vc.OpenPerFamilyFile();
   if (estimateMatrices) matrix.OpenFile();
   if (estimateKinship15) kinship15.OpenFile();
   if (estimateIBD) ibd.OpenFile();
   if (bestHaplotype || (sampledHaplotypes != 0) || allHaplotypes) haplo.OpenFile();
   if (simwalk2) hybrid.OpenFiles();
   if (writePDF) pdf.OpenFile((const char *) (MerlinCore::filePrefix + ".pdf"));

   for (AnalysisTask * task = taskList; task != NULL; task = task->next)
      task->OpenFiles(taskInfo, MerlinCore::filePrefix);
   }

void FamilyAnalysis::CloseFiles()
   {
   if (findErrors)
      {
      printf("Unlikely genotypes listed in file [%s.err]\n", (const char *) MerlinCore::filePrefix);
      fclose(errorfile);
      }
   if (perFamily && storeKinshipForVc) vc.ClosePerFamilyFile();
   if (estimateMatrices) matrix.CloseFile();
   if (estimateKinship15) kinship15.CloseFile();
   if (estimateIBD) ibd.CloseFile();
   if (bestHaplotype || (sampledHaplotypes != 0) || allHaplotypes) haplo.CloseFile();
   if (simwalk2) hybrid.CloseFiles();
   if (writePDF) pdf.CloseFile();
   if (storeKinshipForAssoc) AssociationAnalysis::OutputRefinedModels(MerlinCore::filePrefix, ped);
   if (fastAssociationAnalysis) FastAssociationAnalysis::OutputFastModels(MerlinCore::filePrefix, ped);

   for (AnalysisTask * task = taskList; task != NULL; task = task->next)
      task->CloseFiles();
   }

void FamilyAnalysis::ShowLikelihood()
   {
   if (!calcLikelihood) return;

   printf("lnLikelihood for %d families (%d skipped) = %.3f\n\n",
          lkCount, ped.familyCount - lkCount, lkSum);
   }

void FamilyAnalysis::NewTask(AnalysisTask * task)
   {
   task->next = taskList;
   taskList = task;
   }

bool FamilyAnalysis::GenotypeAnalysisEnabled(int informativeMarker, bool global)
   {
   if (findErrors || inferGenotypes)
      {
      // Remap the current marker
      int marker = global ? informativeMarker : informativeMarkers[informativeMarker];

      // Check if we are in a region that is actually being analysed
      // if (start != _NAN_ && markerPositions[marker] < start * 0.01)
      //   return false;
      if (analysisPositions.Length() == 0 ||
          markerPositions[marker] < analysisPositions[0] ||
          markerPositions[marker] > analysisPositions.Last() )
         return false;

      // Find out if the marker belongs to a cluster of markers in LD
      MarkerCluster * cluster = clusters.Enabled() ? clusters.markerToCluster[markers[marker]] : NULL;

      if (cluster != NULL)
         {
         // Check whether we are trying to impute genotypes for at least one marker in cluster
         if (!findErrors)
            {
            bool imputationEnabled = false;

            for (int m = 0; m < cluster->markerIds.Length(); m++)
               if (imputationEngine.resultsAvailable(cluster->markerIds[m]))
                  {
                  imputationEnabled = true;
                  break;
                  }

            if (!imputationEnabled) return false;
            }

         // If clusters are enabled we need to check whether any of several markers is typed
         for (int m = 0; m < cluster->markerIds.Length(); m++)
            for (int i = family->first, marker = cluster->markerIds[m]; i <= family->last; i++)
               if (ped[i].markers[marker].isKnown())
                  return true;
         }
      else
         {
         // Map marker globally
         marker = markers[marker];

         // Check whether we are trying to impute genotypes for this marker
         if (!findErrors && !imputationEngine.resultsAvailable(marker))
            return false;

         // Check if there is at least one genotyped individual
         for (int i = family->first; i <= family->last; i++)
            if (ped[i].markers[marker].isKnown())
               return true;
         }
      }

   return false;
   }

 
