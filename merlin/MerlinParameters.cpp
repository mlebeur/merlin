////////////////////////////////////////////////////////////////////// 
// merlin/MerlinParameters.cpp 
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
// Wednesday October 17, 2007
// 
 
#include "MerlinParameters.h"
#include "MerlinSimulator.h"
#include "FastAssociation.h"
#include "MerlinFamily.h"
#include "MerlinModel.h"
#include "MerlinCache.h"
#include "HaploFamily.h"
#include "KongAndCox.h"
#include "Pedigree.h"
#include "Houdini.h"
#include "Random.h"

int  MerlinParameters::maxMegabytes = 0;
bool MerlinParameters::trimPedigree = false;

// Specialized analyses
bool MerlinParameters::fitErrorRate = false;
bool MerlinParameters::varianceComponents = false;
bool MerlinParameters::associationAnalysis = false;
String MerlinParameters::qtlModelFile("cov.tbl");

// Simulation parameters
bool MerlinParameters::simulateNull = false;
int  MerlinParameters::reruns = 0;
bool MerlinParameters::saveReplicate = false;

// Output files
bool MerlinParameters::writeFrequencyFile = false;

// Marker Clustering
String MerlinParameters::clusterFile;
double MerlinParameters::clusterDistance = _NAN_;
double MerlinParameters::clusterRsquared = _NAN_;
bool   MerlinParameters::writeClusters = false;

// Trait transformations
bool   MerlinParameters::inverseNormal = false;

BEGIN_LONG_PARAMETERS(longParameters)
   LONG_PARAMETER_GROUP("General")
      LONG_PARAMETER("error", &FamilyAnalysis::findErrors)
      LONG_PARAMETER("information", &FamilyAnalysis::estimateInformation)
      LONG_PARAMETER("likelihood", &FamilyAnalysis::calcLikelihood)
      LONG_STRINGPARAMETER("model", &MerlinModel::modelsFile)
//   LONG_PARAMETER_GROUP("Errors")
//      LONG_PARAMETER("flag", &FamilyAnalysis::findErrors)
//      LONG_DOUBLEPARAMETER("perAllele", &FuzzyInheritanceTree::perAlleleError)
//      LONG_DOUBLEPARAMETER("perGenotype", &FuzzyInheritanceTree::perGenotypeError)
//      LONG_PARAMETER("fit", &MerlinParameters::fitErrorRate)
   LONG_PARAMETER_GROUP("IBD States")
      LONG_PARAMETER("ibd", &FamilyAnalysis::estimateIBD)
      LONG_PARAMETER("kinship", &FamilyAnalysis::estimateKinship)
      LONG_PARAMETER("matrices", &FamilyAnalysis::estimateMatrices)
#ifndef __CHROMOSOME_X__
      LONG_PARAMETER("extended", &FamilyAnalysis::estimateKinship15)
#endif
      LONG_PARAMETER("select", &FamilyAnalysis::selectCases)
   LONG_PARAMETER_GROUP("NPL Linkage")
      LONG_PARAMETER("npl", &KongAndCox::nplAll)
      LONG_PARAMETER("pairs", &KongAndCox::nplPairs)
      LONG_PARAMETER("qtl", &KongAndCox::nplQtl)
      LONG_PARAMETER("deviates", &KongAndCox::nplDeviates)
//      LONG_PARAMETER("extras", &KongAndCox::nplExtras)
      LONG_PARAMETER("exp", &KongAndCox::nplExponential)
   LONG_PARAMETER_GROUP("VC Linkage")
      LONG_PARAMETER("vc", &MerlinParameters::varianceComponents)
      LONG_PARAMETER("useCovariates", &VarianceComponents::useCovariates)
      LONG_PARAMETER("ascertainment", &VarianceComponents::useProbands)
      LONG_DOUBLEPARAMETER("unlinked", &VarianceComponents::unlinkedFraction)
   LONG_PARAMETER_GROUP("Association")
      LONG_PARAMETER("infer", &FamilyAnalysis::inferGenotypes)
      LONG_PARAMETER("assoc", &MerlinParameters::associationAnalysis)
      LONG_PARAMETER("fastAssoc", &FamilyAnalysis::fastAssociationAnalysis)
      LONG_DOUBLEPARAMETER("filter", &FastAssociationAnalysis::fastFilter)
      LONG_STRINGPARAMETER("custom", &MerlinParameters::qtlModelFile)
   LONG_PARAMETER_GROUP("Haplotyping")
      LONG_PARAMETER("best", &FamilyAnalysis::bestHaplotype)
      LONG_SMARTINTPARAMETER("sample", &FamilyAnalysis::sampledHaplotypes)
      LONG_PARAMETER("all", &FamilyAnalysis::allHaplotypes)
      LONG_PARAMETER("founders", &MerlinHaplotype::outputFounders)
      LONG_PARAMETER("horizontal", &MerlinHaplotype::outputHorizontal)
//      LONG_PARAMETER("descent", &MerlinHaplotype::outputGraph)
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
   LONG_PARAMETER_GROUP("LD Clusters")
      LONG_STRINGPARAMETER("clusters", &MerlinParameters::clusterFile)
      LONG_DOUBLEPARAMETER("distance", &MerlinParameters::clusterDistance)
      LONG_DOUBLEPARAMETER("rsq", &MerlinParameters::clusterRsquared)
      LONG_PARAMETER("cfreq", &MerlinParameters::writeClusters)
   LONG_PARAMETER_GROUP("Limits")
      LONG_INTPARAMETER("bits", &MerlinCore::maxBits)
      LONG_INTPARAMETER("megabytes", &MerlinParameters::maxMegabytes)
      LONG_INTPARAMETER("minutes", &MerlinCore::maxMinutes)
   LONG_PARAMETER_GROUP("Performance")
      LONG_PARAMETER("trim", &MerlinParameters::trimPedigree)
#ifndef __CHROMOSOME_X__
      LONG_PARAMETER("noCoupleBits", &Mantra::ignoreCoupleSymmetries)
#endif
      LONG_PARAMETER("swap", &FamilyAnalysis::useSwap)
      LONG_PARAMETER("smallSwap", &MerlinCore::smallSwap)
//      LONG_STRINGPARAMETER("cache", &MerlinCache::directory)
   LONG_PARAMETER_GROUP("Output")
      LONG_PARAMETER("quiet", &MerlinCore::quietOutput)
      LONG_PARAMETER("markerNames", &MerlinCore::useMarkerNames)
      LONG_PARAMETER("frequencies", &MerlinParameters::writeFrequencyFile)
      LONG_PARAMETER("perFamily", &FamilyAnalysis::perFamily)
      LONG_PARAMETER("pdf", &FamilyAnalysis::writePDF)
      LONG_PARAMETER("tabulate", &MerlinCore::tabulate)
      LONG_STRINGPARAMETER("prefix", &MerlinCore::filePrefix)
   LONG_PARAMETER_GROUP("Simulation")
      LONG_PARAMETER("simulate", &MerlinParameters::simulateNull)
      LONG_INTPARAMETER("reruns", &MerlinParameters::reruns)
      LONG_PARAMETER("save", &MerlinParameters::saveReplicate)
      LONG_STRINGPARAMETER("trait", &Simulator::model)
   BEGIN_LEGACY_PARAMETERS()
      LONG_PARAMETER("inferBest", &GenotypeInference::inferBest)
      LONG_PARAMETER("inferExpected", &GenotypeInference::inferExpected)
      LONG_PARAMETER("inferProbabilities", &GenotypeInference::inferProbabilities)
      LONG_PARAMETER("sexAsCovariate", &Pedigree::sexAsCovariate)
      LONG_PARAMETER("inverseNormal", &MerlinParameters::inverseNormal)
      LONG_STRINGPARAMETER("positions", &MerlinCore::positionList)
      LONG_DOUBLEPARAMETER("iTheta", &MerlinCore::imaginaryTheta)
      LONG_PARAMETER("zscores", &KongAndCox::rawOutput)
      LONG_PARAMETER("simwalk2", &FamilyAnalysis::simwalk2)
      LONG_PARAMETER("error", &FamilyAnalysis::findErrors)
      LONG_STRINGPARAMETER("reuse", &Simulator::reuseLabel)
END_LONG_PARAMETERS();

MerlinParameters::MerlinParameters() :
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
   Add(new AlleleFrequencyParameter('f', "Allele Frequencies", alleleFrequencies, freqFile));
   Add(new IntParameter('r', "Random Seed", random_seed));
   Add(new LongParameters("Data Analysis Options", longParameters));

   // The -b switch is obsolete and is thus hidden (use --bits instead)
   Add(new HiddenInteger('b', "Maximum Bits", FamilyAnalysis::maxBits));
   // The -i switch is obsolete and is thus hidden (use --steps instead)
   Add(new HiddenInteger('i', "Steps Per Interval", FamilyAnalysis::steps_per_interval));
   }

void MerlinParameters::Read(int argc, char ** argv, int start)
   {
   ParameterList::Read(argc, argv, start);

   globalRandom.Reset(random_seed);

   // The number of analysis steps between markers must be positive
   if (FamilyAnalysis::steps_per_interval < 0)
      FamilyAnalysis::steps_per_interval = 0;

   // Memory usage limit should be positive
   if (maxMegabytes < 0)
      maxMegabytes = 0;

   // The bit limit for analysis should be positive
   if (FamilyAnalysis::maxBits < 0)
      FamilyAnalysis::maxBits = 0;

   // The number of sampled haplotypes should be positive
   if (FamilyAnalysis::sampledHaplotypes < 0)
      FamilyAnalysis::sampledHaplotypes = 0;

   // Clustering distance must be positive
   if (clusterDistance < 0.0)
      clusterDistance = _NAN_;

   // Clustering r-squared must be between 0 and 1
   if (clusterRsquared < 0.0)
      clusterRsquared = 0.0;

   if (clusterRsquared > 1.0)
      clusterRsquared = 1.0;

   // Limited range of options available with error modelling
   if (FuzzyInheritanceTree::perAlleleError != 0.0 || 
       FuzzyInheritanceTree::perGenotypeError != 0.0 ||
       fitErrorRate)
      {
      char * reason = "The --%s option is disabled when error modeling is enabled";

      Enforce(clusterDistance, _NAN_, reason, "distance");
      Enforce(clusterRsquared, _NAN_, reason, "rsq");
      Enforce(clusterFile, "", reason, "clusters");
      }

   // Limited range of options in two-point mode
   if (FamilyAnalysis::twopoint)
      {
      char * reason = "The %s option is disabled during "
                      "singlepoint analysis.\n";

      // No between marker analyses
      Enforce(FamilyAnalysis::steps_per_interval, 0, reason, "--steps");
      Enforce(FamilyAnalysis::minDistance, _NAN_, reason, "--minStep");
      Enforce(FamilyAnalysis::maxDistance, _NAN_, reason, "--maxStep");
      Enforce(FamilyAnalysis::grid, _NAN_, reason, "--grid");

      // No error checking
      Enforce(FamilyAnalysis::findErrors, false, reason, "--flag");

      // No haplotyping
      Enforce(FamilyAnalysis::bestHaplotype, false, reason, "--best");
      Enforce(FamilyAnalysis::sampledHaplotypes, 0, reason, "--sample");
      Enforce(FamilyAnalysis::allHaplotypes, false, reason, "--all");

      // Positions along chromosome not meaningful
      Enforce(FamilyAnalysis::useMarkerNames, true,
              "The --markerNames option is automatically enabled during singlepoint analysis\n");
      }

   // Limited range of options when using imaginary components for recombination fraction
   if (MerlinCore::imaginaryTheta < 0.0)
      MerlinCore::imaginaryTheta = 0.0;

   if (MerlinCore::imaginaryTheta > 0.0)
      {
      const char * reason = "The %s option is disabled when complex "
                            "recombination fractions are in use\n";

      // No error checking
      Enforce(FamilyAnalysis::findErrors, false, reason, "--flag");

      // No haplotyping
      Enforce(FamilyAnalysis::bestHaplotype, false, reason, "--best");
      Enforce(FamilyAnalysis::sampledHaplotypes, 0, reason, "--sample");
      Enforce(FamilyAnalysis::allHaplotypes, false, reason, "--all");
      }

   if (FamilyAnalysis::allHaplotypes)
      {
      char * reason = "The --%s option is disabled when listing all haplotypes\n";

      Enforce(clusterDistance, _NAN_, reason, "distance");
      Enforce(clusterRsquared, _NAN_, reason, "rsq");
      Enforce(clusterFile, "", reason, "clusters");
      }

   // Limited range of options when using simwalk2
   if (FamilyAnalysis::simwalk2)
      {
      const char * reason = "The %s option is disabled when creating "
                            "intermediate results for SimWalk2\n";

      // No between marker analyses
      Enforce(FamilyAnalysis::steps_per_interval, 0, reason, "--steps");
      Enforce(FamilyAnalysis::minDistance, _NAN_, reason, "--minStep");
      Enforce(FamilyAnalysis::maxDistance, _NAN_, reason, "--maxStep");
      Enforce(FamilyAnalysis::grid, _NAN_, reason, "--grid");
      Enforce(FamilyAnalysis::start, _NAN_, reason, "--start");
      Enforce(FamilyAnalysis::stop, _NAN_, reason, "--stop");

      const char * reason2 = " %s option enabled when creating "
                             "intermediate results for SimWalk2\n";

      // Simwalk2 might want the marker names
      Enforce(FamilyAnalysis::useMarkerNames, true, reason2, "--markerNames");

      if (KongAndCox::nplPairs)
         messages += "The --pairs option will ignore inbreeding "
                     "to allow joint analysis with Simwalk2\n";

      // Non-parametric statistics are required
      Enforce(KongAndCox::nplAll, true, reason2, "--npl");
      Enforce(KongAndCox::nplPairs, true, reason2, "--pairs");
      }

   // Checked that proportion of unlinked families is within range
   if (VarianceComponents::unlinkedFraction < 0.0 ||
       VarianceComponents::unlinkedFraction >= 1.0)
      VarianceComponents::unlinkedFraction = 0.0;

   // Limited range of options zero recombination is assumed
   if (MerlinCore::zeroRecombination)
      {
      const char * reason = "The --%s option is incompatible with the --zero flag\n"
                            "and will be disabled\n";

      Enforce(associationAnalysis, false, reason, "--assoc options");
      Enforce(FamilyAnalysis::findErrors, false, reason, "--flag / --error option");
      }

   if (reruns > 1)
      Enforce(simulateNull, true, "The --reruns option automatically enables the --simulate option\n");

   FamilyAnalysis::storeKinshipForVc = varianceComponents;
   FamilyAnalysis::storeKinshipForAssoc = associationAnalysis;

   // With the --infer parameter, we actively output all information
   if (FamilyAnalysis::inferGenotypes)
      GenotypeInference::inferBest = GenotypeInference::inferExpected =
      GenotypeInference::inferProbabilities = true;

   // With the specific options, we enable genotype inference but don't
   // touch the individual flags
   if (GenotypeInference::inferBest || GenotypeInference::inferExpected ||
       GenotypeInference::inferProbabilities)
       FamilyAnalysis::inferGenotypes = true;

   // The infer genotypes flag controls inference, whereas the output
   // flag determines whether results get saved to disk
   if (FamilyAnalysis::inferGenotypes)
      FamilyAnalysis::outputInferredGenotypes = true;

   if (associationAnalysis || FamilyAnalysis::fastAssociationAnalysis)
      FamilyAnalysis::inferGenotypes = true;

   // Some basic calculations to convert megabyte memory usage
   // limit into maximum TreeNode count
   BasicTree::maxNodes = (1024 * 1024) / sizeof(TreeNode) * maxMegabytes;

   // If we are saving replicated data, then we should also
   // be simulating it!
   simulateNull |= saveReplicate;

   // The max bit parameters for processing haplotype clusters
   // should be the same as for other analyses
   FamilyHaplos::maxBits = MerlinCore::maxBits;
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


 
