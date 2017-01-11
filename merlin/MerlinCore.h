////////////////////////////////////////////////////////////////////// 
// merlin/MerlinCore.h 
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
 
#ifndef __MERLINCORE_H__
#define __MERLINCORE_H__

///////////////////////////////////////////////////////////////
// This class handles some fundamental pedigree analysis tasks
// that allow likelihoods to be calculated at arbitrary genomic
// locations.  Derived classes can use these likelihoods to
// calculate linkage scores, information content, ...
//

#include "Conquer.h"
#include "TreeInfo.h"

class MerlinCore
   {
   public:
      MerlinCore(Pedigree & ped);
      virtual ~MerlinCore();

      virtual void SetupGlobals();
      virtual void CleanupGlobals();

      virtual int  SetupMap(int chromosome = -1);
      virtual void SexSpecificMap();

      virtual bool SelectFamily(Family * f, bool warnOnSkip = true);

      // Iterators for likelihood calculations
      void ScoreSinglepoint();
      void CalculateThetas();
      bool ScoreConditionals();
      void ScanChromosome();

      // Routines for managing intermediate singlepoint results
      void StoreSinglepoint(int position);
      void RetrieveSinglepoint(int position);
      void ReStoreSinglepoint(int position);

      // Routines for managing intermediate multipoint results
      void StoreMultipoint(int position);
      void RetrieveMultipoint(int position);
      void ReStoreMultipoint(int position);

      // Calculates likelihood of marker data for currently selected family
      double CalculateLikelihood();

      // Functions for handling maps with zero recombination
      void SortInformation(int start, int end);
      bool ZeroRecombination();

      // Limits to stop Merlin from attempting problems that are too tough
      static int  maxBits;
      static int  maxMinutes;

      // Analysis positions
      static int    steps_per_interval;
      static double maxDistance;
      static double minDistance;
      static double grid;
      static double stop;
      static double start;
      static String positionList;

      // General flags for controlling output
      static bool quietOutput;
      static bool useMarkerNames;
      static bool tabulate;
      static String filePrefix;

      // Flags for controlling multipoint calculations
      static bool zeroRecombination;
      static bool oneRecombination;
      static bool twoRecombinations;
      static bool threeRecombinations;
      static bool twopoint;

      // Flags for controlling error model
      static double imaginaryTheta;

      // Flags for controlling performance
      static bool useSwap;
      static bool smallSwap;

      // Internal flags to minimize useless calculations
      static bool multipoint;
      static bool scanChromosome;
      static bool betweenMarkers;
      static bool lazyRightConditional;

      // Bounds for rescaling inheritance vectors
      double lowBound, rescale, lnScale;

      // These functions provide for buffered output and reduce clutter
      void PrintHeader();
      void PrintMessage(const char * msg, ...);
      void ProgressReport(const char * msg, int done = 0, int total = 1);
      void CleanMessages();

      // Index for locations where analysis will be carried out and
      // corresponding descriptive labels (e.g. marker names)
      StringArray labels;
      Vector      analysisPositions;

      // Index between markers in pedigree and markers in the chromosome
      // being analysed
      IntArray markers;
      int      markerCount;
      Vector   markerPositions, markerMalePositions, markerFemalePositions;

      // Information on subset of selected markers that is informative
      // in the currently selected family
      IntArray informativeMarkers;
      Vector   informativePositions;
      Vector   informationScores;
      Vector   informationBounds;
      IntArray informationIndex;
      int      informativeCount;

      // Between marker recombination fractions along sex-specific map
      Vector   femaleMarkerTheta;
      Vector   maleMarkerTheta;

      // These array store between marker recombination fractions
      // for informative markers in each family
      bool     keyThetas;
      Vector   keyMaleTheta;
      Vector   keyFemaleTheta;

      // Recombination fractions between analysis positions
      Vector   femalePositions;
      Vector   malePositions;

      // Trees with likelihoods at each marker
      InheritanceTree * singlepoint;
      Tree            * right;
      Vector          rightScale, leftScale;

      // These variables track the currently selected family
      Mantra     mantra;
      Pedigree & ped;
      Family   * family;
      int        serial;

      // Used for testing marker informativeness and other summary statistics
      // on inheritance trees
      TreeInfo   stats;

      // Flag indicating whether to try fast calculations for very sparse trees
      bool useSparse;

      // Check if marker fits allele count limitations
      static void ValidateMarker(MarkerInfo * info);

   protected:
      // Used for managing output buffering
      int  printHeader, cleanOutput;

      // Likelihood for the current family
      double likelihood;

      // Allows derived classes to perform custom analyses at each location
      virtual void AnalyseLocation(int pos, Tree & inheritance);

      // Allows derived classes to perform genotype inference
      virtual bool GenotypeAnalysisEnabled(int marker, bool globalMapping);
      virtual void GenotypeAnalysisHook(Tree & withMarker, Tree & without, int marker);
      virtual void GenotypeAnalysisHandler(Multipoint & left, Multipoint & full, int marker);

      // Swap file space
      TreeManager  singlepointSwap;
      TreeManager  multipointSwap;
      TreeManager  scaffoldSwap;

      // Swap file management information
      int          lastSinglepointSwap;
      int          multipointGrid;
      int          multipointSwapStart, multipointSwapEnd;

      // Swap file helper functions
      void         FreeSwap();

      bool         memoryManagement;

   private:
      // Update left conditional probabilities at the edge of map
      // (used when skipping through positions where we don't want to evaluate
      //  lod scores, etc.)
      void RunLeft(Multipoint & left, Multipoint & leftMarker, int marker);

      // Calculate left conditional & multipoint probabilities
      void WalkLeft(Multipoint & left, Multipoint & leftMarker,
                    Multipoint & inheritance, int marker);
   };

// This empty class is used to signal OutOfTime exceptions
class OutOfTime
   {
   public:
      OutOfTime() { };
      OutOfTime(const OutOfTime &) { };
   };

#endif



 
