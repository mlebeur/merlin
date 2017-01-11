////////////////////////////////////////////////////////////////////// 
// regress/RegressParameters.cpp 
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
 
#include "RegressParameters.h"
#include "RegressAnalysis.h"
#include "MerlinCache.h"
#include "Pedigree.h"
#include "Houdini.h"
#include "AutoFit.h"
#include "Random.h"
#include "Error.h"

// Memory limit for gene flow trees
int  RegressionParameters::maxMegabytes = 0;
bool RegressionParameters::trimPedigree = false;

// Carry out gene dropping simulations
bool RegressionParameters::simulateNull = 0;
int  RegressionParameters::reruns = 0;

// Output allele frequencies
bool RegressionParameters::writeFrequencyFile = 0;

// Estimate an error rate from the data
bool RegressionParameters::fitErrorRate = false;

// Apply inverse normal transform
bool RegressionParameters::inverseNormal = false;

// Marker Clustering
String RegressionParameters::clusterFile;
double RegressionParameters::clusterDistance = _NAN_;
double RegressionParameters::clusterRsquared = _NAN_;

BEGIN_LONG_PARAMETERS(longParameters)
   LONG_PARAMETER_GROUP("User Model")
      LONG_DOUBLEPARAMETER("mean", &RegressionAnalysis::mean)
      LONG_DOUBLEPARAMETER("variance", &RegressionAnalysis::variance)
      LONG_DOUBLEPARAMETER("heritability", &RegressionAnalysis::heritability)
      LONG_DOUBLEPARAMETER("testRetest", &RegressionAnalysis::testRetestCorrel)
   LONG_PARAMETER_GROUP("Automatic Models")
      LONG_PARAMETER("randomSample", &RegressionAnalysis::autoFit)
      LONG_PARAMETER("useCovariates", &AutoFit::useCovariates)
      LONG_PARAMETER("sexAsCovariate", &Pedigree::sexAsCovariate)
      LONG_PARAMETER("inverseNormal", &RegressionParameters::inverseNormal)
   LONG_PARAMETER_GROUP("Errors")
      LONG_DOUBLEPARAMETER("perAllele", &FuzzyInheritanceTree::perAlleleError)
      LONG_DOUBLEPARAMETER("perGenotype", &FuzzyInheritanceTree::perGenotypeError)
      LONG_PARAMETER("fit", &RegressionParameters::fitErrorRate)
//      LONG_PARAMETER("autofit", &RegressionParameters::autofit)
   LONG_PARAMETER_GROUP("Recombination")
      EXCLUSIVE_PARAMETER("zero", &MerlinCore::zeroRecombination)
      EXCLUSIVE_PARAMETER("one", &MerlinCore::oneRecombination)
      EXCLUSIVE_PARAMETER("two", &MerlinCore::twoRecombinations)
      EXCLUSIVE_PARAMETER("three", &MerlinCore::threeRecombinations)
      EXCLUSIVE_PARAMETER("singlepoint", &MerlinCore::twopoint)
   LONG_PARAMETER_GROUP("Positions")
      LONG_INTPARAMETER("steps", &MerlinCore::steps_per_interval)
      LONG_DOUBLEPARAMETER("maxStep", &MerlinCore::maxDistance)
      LONG_DOUBLEPARAMETER("minStep", &MerlinCore::minDistance)
      LONG_DOUBLEPARAMETER("grid", &MerlinCore::grid)
      LONG_DOUBLEPARAMETER("start", &MerlinCore::start)
      LONG_DOUBLEPARAMETER("stop", &MerlinCore::stop)
   LONG_PARAMETER_GROUP("Marker Clusters")
      LONG_STRINGPARAMETER("clusters", &RegressionParameters::clusterFile)
      LONG_DOUBLEPARAMETER("distance", &RegressionParameters::clusterDistance)
      LONG_DOUBLEPARAMETER("rsq", &RegressionParameters::clusterRsquared)
   LONG_PARAMETER_GROUP("Basic Performance")
      LONG_INTPARAMETER("bits", &MerlinCore::maxBits)
      LONG_INTPARAMETER("megabytes", &RegressionParameters::maxMegabytes)
      LONG_INTPARAMETER("minutes", &MerlinCore::maxMinutes)
      LONG_PARAMETER("trim", &RegressionParameters::trimPedigree)
   LONG_PARAMETER_GROUP("Performance")
      LONG_PARAMETER("noCoupleBits", &Mantra::ignoreCoupleSymmetries)
      LONG_PARAMETER("swap", &MerlinCore::useSwap)
      LONG_STRINGPARAMETER("cache", &MerlinCache::directory)
   LONG_PARAMETER_GROUP("Output")
      LONG_STRINGPARAMETER("prefix", &MerlinCore::filePrefix)
      LONG_PARAMETER("pdf", &RegressionAnalysis::writePDF)
      LONG_PARAMETER("tabulate", &MerlinCore::tabulate)
      LONG_PARAMETER("quiet", &MerlinCore::quietOutput)
      LONG_PARAMETER("markerNames", &MerlinCore::useMarkerNames)
   LONG_PARAMETER_GROUP("Others")
      LONG_PARAMETER("simulate", &RegressionParameters::simulateNull)
      LONG_INTPARAMETER("reruns", &RegressionParameters::reruns)
      LONG_PARAMETER("rankFamilies", &RegressionAnalysis::rankFamilies)
      LONG_PARAMETER("unrestriced", &RegressionAnalysis::unrestricted)
END_LONG_PARAMETERS();

RegressionParameters::RegressionParameters() :
   pedigreeFile("merlin.ped"),
   dataFile("merlin.dat"),
   mapFile("merlin.map"),
   freqFile(""),
   alleleFrequencies(FREQ_ALL),
   random_seed(123456)
   {
   Add(new StringParameter('d', "Data File", dataFile));
   Add(new StringParameter('p', "Pedigree File", pedigreeFile));
   Add(new StringParameter('x', "Missing Value Code", Pedigree::missing, false));
   Add(new StringParameter('m', "Map File", mapFile));
   Add(new StringParameter('t', "Trait Models File", RegressionAnalysis::modelsFile));
   Add(new AlleleFrequencyParameter('f', "Allele Frequencies", alleleFrequencies, freqFile));
   Add(new IntParameter('r', "Random Seed", random_seed));
   Add(new LongParameters("Regression Analysis Options", longParameters));

   // The -b switch is obsolete and is thus hidden (use --bits instead)
   Add(new HiddenInteger('b', "Maximum Bits", MerlinCore::maxBits));
   // The -i switch is obsolete and is thus hidden (use --steps instead)
   Add(new HiddenInteger('i', "Steps Per Interval", MerlinCore::steps_per_interval));
   }

void RegressionParameters::Read(int argc, char ** argv, int start)
   {
   ParameterList::Read(argc, argv, start);

   globalRandom.Reset(random_seed);

   // The number of analysis steps between markers must be positive
   if (MerlinCore::steps_per_interval < 0)
      MerlinCore::steps_per_interval = 0;

   // Memory usage limit should be positive
   if (maxMegabytes < 0)
      maxMegabytes = 0;

   // The bit limit for analysis should be positive
   if (MerlinCore::maxBits < 0)
      MerlinCore::maxBits = 0;

   // Limited range of options in two-point mode
   if (MerlinCore::twopoint)
      {
      // No between marker analyses
      MerlinCore::steps_per_interval = 0;

      // Positions along chromosome not meaningful
      MerlinCore::useMarkerNames = true;
      }

   if (AutoFit::useCovariates)
      Enforce(RegressionAnalysis::autoFit, true,
         "The --randomSamples option is required for modelling covariate effects\n");

   if (RegressionParameters::inverseNormal)
      Enforce(RegressionAnalysis::autoFit, true,
         "The --randomSamples option is required for the inverse normal transform\n");

   if (RegressionAnalysis::autoFit)
      Enforce(RegressionAnalysis::testRetestCorrel, 1.0,
         "The --testRetest option is disabled for automatically fitted models\n");

   if (reruns > 1)
      Enforce(simulateNull, true, "The --reruns option automatically enables the --simulate option\n");
      
   BasicTree::maxNodes = (1024 * 1024) / sizeof(TreeNode) * maxMegabytes;
   }

void RegressionParameters::Check()
   {
   if (RegressionAnalysis::variance <= 0.0)
      error("Trait variance must be positive");

   if (RegressionAnalysis::heritability <= 0.0 || RegressionAnalysis::heritability >= 1.0)
      error("Trait heritability must be between 0.0 and 1.0");

   if (RegressionAnalysis::testRetestCorrel <= 0.0 || RegressionAnalysis::testRetestCorrel > 1.0)
      error("Test-retest correlation must be between 0.0 and 1.0");

   if (RegressionAnalysis::testRetestCorrel <= RegressionAnalysis::heritability)
      error("Heritability must be smaller than the test-retest correlation");
   }

AlleleFrequencyParameter::AlleleFrequencyParameter
   (char c, const char * desc, int & how, String & file) :
   Parameter(c, desc, &how), status("ALL INDIVIDUALS")
   {
   how = FREQ_ALL;
   filename = &file;
   }

void AlleleFrequencyParameter::Translate(const char * str)
   {
   String value(str);
   int * estimator = (int *) var;

   if (value == "a" || value == "A")
      {
      * estimator = FREQ_ALL;
      * filename = "";
      status = "ALL INDIVIDUALS";
      }
   else if (value == "e" || value == "E")
      {
      * estimator = FREQ_EQUAL;
      * filename = "";
      status = "ASSUMED EQUAL";
      }
   else if (value == "f" || value == "F")
      {
      * estimator = FREQ_FOUNDERS;
      * filename = "";
      status = "FOUNDERS ONLY";
      }
   else if  (value == "m" || value == "M")
      {
      * estimator = FREQ_ML;
      * filename = "";
      status = "MAX LIKELIHOOD";
      }
   else
      {
      * estimator = FREQ_ALL;
      * filename = value;
      status = value;
      }
   }

void AlleleFrequencyParameter::Status()
   {
   printf("%*s : %*s (-%c[a|e|f|m|file])\n", nameCol, description,
          statusCol, (const char *) status, ch);
   }

   
 
