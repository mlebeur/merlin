////////////////////////////////////////////////////////////////////// 
// offline/Main.cpp 
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
 
#include "Parameters.h"
#include "MerlinFamily.h"
#include "FastAssociation.h"
#include "TraitTransformations.h"

int main(int argc, char ** argv)
   {
   printf("MERLIN -- Offline Association Analysis\n"
          "          (c) 2006-2007 Goncalo Abecasis\n\n");

   String datfile("merlin.dat"), pedfile("merlin.ped");
   String freqfile("merlin.freq"), mapfile("merlin.map");
   String datinfer("merlin-infer.dat"), pedinfer("merlin-infer.ped");
   String covfile("covars.tbl");
   bool   inverseNormal = false;
   bool   sanityCheck = false;
   bool   updateFrequencies = false;

   BEGIN_LONG_PARAMETERS(additional)
      LONG_PARAMETER_GROUP("Inferred Genotypes")
         LONG_STRINGPARAMETER("datinfer", &datinfer)
         LONG_STRINGPARAMETER("pedinfer", &pedinfer)
      LONG_PARAMETER_GROUP("Analysis Options")
         LONG_PARAMETER("inverseNormal", &inverseNormal)
         LONG_PARAMETER("useCovariates", &VarianceComponents::useCovariates)
         LONG_DOUBLEPARAMETER("filter", &FastAssociationAnalysis::fastFilter)
         LONG_STRINGPARAMETER("custom", &covfile)
      LONG_PARAMETER_GROUP("Output Files")
         LONG_STRINGPARAMETER("prefix", &MerlinCore::filePrefix)
         LONG_PARAMETER("pdf", &FamilyAnalysis::writePDF)
         LONG_PARAMETER("tabulate", &MerlinCore::tabulate)
   BEGIN_LEGACY_PARAMETERS()
         LONG_PARAMETER("sanityCheck", &sanityCheck)
         LONG_PARAMETER("updateFrequencies", &updateFrequencies)
   END_LONG_PARAMETERS();

   // Process parameters
   ParameterList pl;

   pl.Add(new StringParameter('d', "Data File", datfile));
   pl.Add(new StringParameter('p', "Pedigree File", pedfile));
   pl.Add(new StringParameter('m', "Map File", mapfile));
   pl.Add(new StringParameter('f', "Frequency File", freqfile));
   pl.Add(new LongParameters("Additional Options", additional));
   pl.Read(argc, argv);
   pl.Status();

   // Load pedigree files, one with inferred genotypes, the other with
   // phenotype data

   Pedigree            ped;
   PedigreeDescription genotypes;

   ped.Prepare(datfile);
   int realTraits = Pedigree::traitCount;
   int realCovariates = Pedigree::covariateCount;

   genotypes.Load(datinfer);

   ped.LoadAlleleFrequencies(freqfile, true);
   ped.Load(pedfile);

   ped.multiFileCount = 1;
   ped.pd = genotypes;
   ped.Load(pedinfer);

   VarianceComponents::customModels.LoadFromFile(covfile);

   FamilyAnalysis engine(ped);

   engine.imputationEngine.offline = true;
   engine.imputationEngine.CreateOfflineMarkers(ped);
   engine.imputationEngine.CreateOfflineLookup();

   ped.LoadMarkerMap(mapfile, true /* filter out markers not in pedigree */);

   Pedigree::traitCount = realTraits;
   Pedigree::covariateCount = realCovariates;

   if (inverseNormal)
      InverseNormalTransform(ped);

   int  next_chromosome  = -1;
   bool many_chromosomes = false;

   if (updateFrequencies)
      engine.imputationEngine.UpdateFrequencies(ped);

   if (sanityCheck)
      engine.imputationEngine.CheckFrequencies(ped);

   engine.SetupFiles();
   do {
      next_chromosome = engine.SetupMap(next_chromosome);

      many_chromosomes = many_chromosomes || (next_chromosome != 0);

      if (many_chromosomes)
         printf("\nAnalysing Chromosome %d\n\n",
                ped.GetMarkerInfo(engine.markers[0])->chromosome);

      FastAssociationAnalysis modeler;

      modeler.AnalyseAssociation(engine);
   } while (next_chromosome);
   engine.CloseFiles();

   FastAssociationAnalysis::OutputFastModels(MerlinCore::filePrefix, ped);
   }
 
