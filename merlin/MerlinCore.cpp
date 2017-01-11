////////////////////////////////////////////////////////////////////// 
// merlin/MerlinCore.cpp 
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
 
#include "MerlinCore.h"
#include "MapFunction.h"
#include "MerlinCache.h"
#include "MerlinCluster.h"
#include "MathStats.h"
#include "Houdini.h"
#include "Error.h"

#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>
#include <ctype.h>

// Limits to stop Merlin from attempting problems that are too tough
int MerlinCore::maxBits = 24;
int MerlinCore::maxMinutes = 0;

// Positions to analyse
int    MerlinCore::steps_per_interval = 0;
double MerlinCore::maxDistance = _NAN_;
double MerlinCore::minDistance = _NAN_;
double MerlinCore::grid  = _NAN_;
double MerlinCore::start = _NAN_;
double MerlinCore::stop  = _NAN_;
String MerlinCore::positionList;

// Output options
bool MerlinCore::quietOutput = false;
bool MerlinCore::useMarkerNames = false;
bool MerlinCore::tabulate = false;
String MerlinCore::filePrefix = "merlin";

// Swap options
bool MerlinCore::useSwap = false;
bool MerlinCore::smallSwap = false;

// Internal flags to minimize useless calculations
bool MerlinCore::multipoint = false;
bool MerlinCore::scanChromosome = false;
bool MerlinCore::betweenMarkers = false;
bool MerlinCore::lazyRightConditional = false;

// Flags for controlling recombination
bool MerlinCore::zeroRecombination = false;
bool MerlinCore::oneRecombination = false;
bool MerlinCore::twoRecombinations = false;
bool MerlinCore::threeRecombinations = false;
bool MerlinCore::twopoint = false;

// Error model
double MerlinCore::imaginaryTheta = 0.00;

// Constructor and Destructor
//

MerlinCore::MerlinCore(Pedigree & pedigree) : mantra(), ped(pedigree)
   {
   memoryManagement = false;
   printHeader = false;
   lowBound = pow(2.0, -129);
   rescale  = pow(2.0, 258);
   lnScale  = 258 * log(0.5);
   }

MerlinCore::~MerlinCore()
   {
   MerlinCache::FreeBuffers();
   }

// Global setup that must be carried out before calculating likelihods
//

void MerlinCore::SetupGlobals()
   {
   // Set the maximum number of recombinants in multipoint calculations
   Multipoint::maximum_recombinants =
      threeRecombinations * 3 + twoRecombinations * 2 + oneRecombination;

   if (useSwap || smallSwap)
      {
      BasicTree::SetupSwap();
      singlepointSwap.OpenFiles();
      multipointSwap.OpenFiles();
      scaffoldSwap.OpenFiles();
      }
   }

void MerlinCore::CleanupGlobals()
   {
   if (useSwap || smallSwap)
      {
      bool   validSize = BasicTree::SwapFileSize() >= 0 &&
                         singlepointSwap.GetFileSize() >= 0 &&
                         multipointSwap.GetFileSize() >= 0 &&
                         scaffoldSwap.GetFileSize() >= 0;

      double fileSize = BasicTree::SwapFileSize() +
                        singlepointSwap.GetFileSize() +
                        multipointSwap.GetFileSize() +
                        scaffoldSwap.GetFileSize();

      BasicTree::CloseSwap(true);
      singlepointSwap.CloseFiles();
      multipointSwap.CloseFiles();
      scaffoldSwap.CloseFiles();

      if (validSize)
         printf("Swap file usage peaked at < %.0fM\n", fileSize * 1e-6 + 1);
      else
         printf("Swap file usage peaked at > 2048M\n");
      }
   }

// Setup map for analysing the next chromosome at currently selected locations
//

int MerlinCore::SetupMap(int chromosome)
   {
   markers.Dimension(0);
   int next = PedigreeGlobals::SortMarkersInMapOrder(markers, chromosome);

   markerCount = markers.Length();

   if (markerCount == 0) error("No markers could be placed");

   // Check that stop position for analysis is greater than start position
   if (MerlinCore::start != _NAN_ && MerlinCore::stop != _NAN_ &&
       MerlinCore::start > MerlinCore::stop)
      warning("Requested position for starting analyses is greater than final position\n");

   // Update swap file information
   BasicTree::UpdateSwapThreshold(markerCount * 2 + 10);
   multipointGrid = introot(markerCount);

   // Clear previous information
   labels.Clear();
   analysisPositions.Clear();

   // Allocate arrays of known size
   markerPositions.Dimension(markerCount);
   markerFemalePositions.Dimension(markerCount);
   markerMalePositions.Dimension(markerCount);

   String label;

   // Retrieve marker positions
   for (int i = 0; i < markerCount; i++)
      {
      MarkerInfo * info = PedigreeGlobals::GetMarkerInfo(markers[i]);
      ValidateMarker(info);
      markerPositions[i] = info->position;
      markerMalePositions[i] = info->positionMale;
      markerFemalePositions[i] = info->positionFemale;
      }

   // For maps with no recombination, only analyse one location
   if (zeroRecombination)
      {
      analysisPositions.Push(markerPositions[0]);
      labels.Push(useMarkerNames ?
                  PedigreeGlobals::GetMarkerInfo(markers[0])->name :
                  label = (markerPositions[0] * 100.0));
      }
   // Singlepoint analysis at marker locations only
   else if (twopoint)
      {
      analysisPositions = markerPositions;

      for (int i = 0; i < markerCount; i++)
         labels.Push(PedigreeGlobals::GetMarkerInfo(markers[i])->name);
      }
   // Analysis at positions specified by the user
   else if (positionList.Length())
      {
      // Split list into tokens
      StringArray tokens;
      tokens.AddTokens(positionList, ',');

      // Find out what chromosome we are analyzing
      int chr = PedigreeGlobals::GetMarkerInfo(markers[0])->chromosome;

      // Each position should either be a marker name or a numeric position
      for (int i = 0; i < tokens.Length(); i++)
         {
         int lookup = ped.LookupMarker(tokens[i]);

         // If it is a marker name, we check the chromosome
         if (lookup >= 0)
            {
            MarkerInfo * info = ped.GetMarkerInfo(lookup);

            if (info->chromosome == chr &&
                (start == _NAN_ || info->position >= start * 0.01) &&
                (stop  == _NAN_ || info->position <= stop  * 0.01))
               {
               analysisPositions.Push(info->position);
               labels.Push(tokens[i]);
               }
            }
         // Otherwise we simply add it to the list
         else if (isdigit(tokens[i][0]) || tokens[i][0] == 0)
            {
            double position = tokens[i].AsDouble();

            if ((start == _NAN_ || position >= start) &&
                (stop  == _NAN_ || position <= stop))
               {
               analysisPositions.Push(position * 0.01);
               labels.Push(tokens[i]);
               }
            }
         }

      // TODO -- More efficient sort would be better here ...
      // Insertion sort will only work well for a small number of positions
      for (int i = 1; i < analysisPositions.Length(); i++)
         for (int j = i - 1, k = i; j >= 0; j--, k--)
            if (analysisPositions[k] < analysisPositions[j])
               {
               labels[j].Swap(labels[k]);
               analysisPositions.Swap(j,k);
               }
            else
               break;
      }
   // Analysis along grid of equally spaced locations
   else if (grid != _NAN_)
      {
      double position = start != _NAN_ ? start * 0.01 :
         PedigreeGlobals::GetMarkerInfo(markers[0])->position;
      double mapEnd = stop != _NAN_ ? stop * 0.01 :
         PedigreeGlobals::GetMarkerInfo(markers.Last())->position;
      double step = grid * 0.01;

      while (true)
         {
         analysisPositions.Push(position);
         labels.Push(label = position * 100.0);
         if (mapEnd - position < 1e-5) break;
         position += step;
         }
      }
   // Analysis at marker positions and equally spaced locations in-between
   else
      {
      double min = minDistance == _NAN_ ? _NAN_ : minDistance * 0.01;
      double max = maxDistance == _NAN_ ? _NAN_ : maxDistance * 0.01;

      double last = markerPositions[0];

      for (int i = 0; i < markerCount; i++)
         {
         if (markerPositions[i] > last)
            {
            double delta = markerPositions[i] - last;
            int    steps = steps_per_interval;
            double step  = delta / (steps + 1);

            if ((max != _NAN_) && (max < step))
               {
               steps = (int) (delta / (max + 1e-6));
               step  = delta / (steps + 1);
               }

            if ((min != _NAN_) && min > step)
               steps = (int) (delta / (minDistance  + 1e-6) - 0.999);

            if (steps < 0) steps = 0;

            for (int j = 1; j <= steps; j++)
               {
               double position = last + delta * j / (steps+1);

               analysisPositions.Push(position);
               labels.Push(label = position * 100.0);
               }
            }

         if (markerPositions[i] > last || i == 0)
            {
            analysisPositions.Push(last = markerPositions[i]);
            labels.Push(useMarkerNames ?
               PedigreeGlobals::GetMarkerInfo(markers[i])->name :
               label = (last * 100.0));
            }
         }

      if (start != _NAN_)
         {
         int firstPosition = analysisPositions.Length() -
                             analysisPositions.CountIfGreaterOrEqual(start * 0.01);

         if (firstPosition > 0)
            {
            for (int i = firstPosition; i < analysisPositions.Length(); i++)
               {
               analysisPositions[i - firstPosition] = analysisPositions[i];
               labels[i - firstPosition] = labels[i];
               }

            analysisPositions.Dimension(analysisPositions.Length() - firstPosition);
            labels.Dimension(analysisPositions.Length());
            }
         }

      if (stop != _NAN_)
         {
         int lastPosition = analysisPositions.Length() -
                            analysisPositions.CountIfGreater(stop * 0.01);

         analysisPositions.Dimension(lastPosition);
         labels.Dimension(lastPosition);
         }
      }

   SexSpecificMap();

   if (analysisPositions.Length() &&
       analysisPositions.Last() - analysisPositions.First() > 300)
       warning("According to the current map, this chromosome is %.0f Morgans long.\n"
               "Check units in your map file.\n",
                (analysisPositions.Last() - analysisPositions.First()));

   return next;
   }

// Construct sex specific map using current sex averaged map as a template
//

void MerlinCore::SexSpecificMap()
   {
   // First calculate recombination fractions between consecutive markers
   femaleMarkerTheta.Dimension(markerCount - 1);
   maleMarkerTheta.Dimension(markerCount - 1);
   keyFemaleTheta.Dimension(markerCount - 1);
   keyMaleTheta.Dimension(markerCount - 1);

   for (int i = 1; i < markerCount; i++)
      {
      MarkerInfo * info1 = ped.GetMarkerInfo(markers[i - 1]);
      MarkerInfo * info2 = ped.GetMarkerInfo(markers[i]);

      femaleMarkerTheta[i - 1] =
         DistanceToRecombination(info2->positionFemale - info1->positionFemale);
      maleMarkerTheta[i - 1] =
         DistanceToRecombination(info2->positionMale - info1->positionMale);
      }

   // Check if there are any analysis positions ...
   if (analysisPositions.Length() == 0) return;

   // Pre-allocate arrays to speed things up
   malePositions.Dimension(analysisPositions.Length());
   femalePositions.Dimension(analysisPositions.Length());

   // Variables to keep track of current position
   int pos = 0;
   int marker = 0;

   // Next calculate recombination fractions between each analysis position an
   // the first marker (using the sex averaged map)
   while (analysisPositions[pos] < markerPositions[marker])
      {
      double delta = analysisPositions[pos] - markerPositions[marker];
      malePositions[pos] = markerMalePositions[marker] + delta;
      femalePositions[pos] = markerFemalePositions[marker] + delta;
      pos++;
      }

   // Calculate recombination fractions for all intervening positions
   while (marker + 1 < markerCount)
      {
      marker++;

      double delta = markerPositions[marker] - markerPositions[marker - 1];
      double maleDelta = markerMalePositions[marker] - markerMalePositions[marker - 1];
      double femaleDelta = markerFemalePositions[marker] - markerFemalePositions[marker - 1];

      while (pos < analysisPositions.Length() &&
             analysisPositions[pos] < markerPositions[marker])
         {
         double f = (analysisPositions[pos] - markerPositions[marker - 1]) / delta;

         malePositions[pos] = markerMalePositions[marker - 1] + maleDelta * f;
         femalePositions[pos] = markerFemalePositions[marker - 1] + femaleDelta * f;
         pos++;
         }
      }

   // Calculate recombination fractions between final analysis positions and
   // the last marker (using the sex averaged map)
   while (analysisPositions.Length() > pos)
      {
      double delta = analysisPositions[pos] - markerPositions[marker];
      malePositions[pos] = markerMalePositions[marker] + delta;
      femalePositions[pos] = markerFemalePositions[marker] + delta;
      pos++;
      }
   }

// Setup analyses for a new family
//

bool MerlinCore::SelectFamily(Family * f, bool warnOnSkip)
   {
   mantra.Prepare(ped, *f);
   mantra.PrepareIBD();

   printHeader = true;
   cleanOutput = false;

   family = f;
   serial = f->serial;

   memoryManagement = useSwap || smallSwap;

   if (!memoryManagement && BasicTree::maxNodes)
      {
      memoryManagement = mantra.bit_count > (int) (CHAR_BIT * sizeof(mantra.bit_count) - 2) ||
                         (2 << mantra.bit_count) > BasicTree::LargeTreeSize() ||
                         BasicTree::LargeTreeSize() <= 256;

      if (memoryManagement)
         PrintMessage("  Advanced Memory Management enabled for large pedigree");
      }

   // if (mantra.couples)
   //   PrintMessage("  %d symmetric founder couple%s detected",
   //               mantra.couples, mantra.couples == 1 ? "" : "s");

   if (mantra.bit_count <= maxBits)
      return true;

   if (warnOnSkip)
      PrintMessage("  SKIPPED: Maximum complexity set at %d bits", maxBits);
   CleanMessages();
   return false;
   }

// Calculate Likelihood of Observed Marker Data
//

double MerlinCore::CalculateLikelihood()
   {
   singlepoint = NULL;
   right = NULL;

   try {
      FreeSwap();

      ScoreSinglepoint();

      if (informativeCount > 0)
         ScoreConditionals();

      CleanMessages();

      if (right != NULL)
         delete [] right;
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

   return likelihood;
   }

// Calculation of single marker likelihoods
//

void MerlinCore::ScoreSinglepoint()
   {
   MerlinCache cache;

   keyThetas = false;

   singlepoint = new InheritanceTree[ped.markerCount];
   FuzzyInheritanceTree tempVectors;

   informativeCount = 0;
   informativeMarkers.Clear();
   informativePositions.Clear();
   informationScores.Clear();
   informationBounds.Clear();

#ifdef __PRINTOUT_STATS__
   Vector nodeCounts, leafCounts, zeroCounts, infoScores;
#endif

   // Initialize a file cache for this family
   cache.OpenCache(mantra);

   // Reset the likelihood
   likelihood = 0;

   for (int m = 0; m < markerCount; m++)
      {
      bool error_message_printed = false;
      int  markerid = markers[m];

      ProgressReport("Singlepoint", m, markerCount);

      if (clusters.Enabled() && clusters.markerToCluster[markerid] != NULL)
         {
         MarkerCluster * cluster = clusters.markerToCluster[markerid];

         // Score using clustering algorithms
         cluster->ScoreLikelihood(mantra, tempVectors);

         // Clear error messages
         if (!cluster->errormsg.IsEmpty())
            PrintMessage("  %s", (const char *) cluster->errormsg,
                                  error_message_printed = true);

         // TO DO -- Code for implementing cache when there are clustered markers

         // Skip subsequent markers that are in the same cluster
         m += cluster->markerIds.Length() - 1;
         }
      else
         {
         mantra.SelectMarker(markerid);

         if (!cache.RetrieveFromCache(mantra, tempVectors, Pedigree::markerNames[markerid]))
            {
            tempVectors.FuzzyScoreVectors(mantra);
            cache.SaveToCache(mantra, tempVectors, Pedigree::markerNames[markerid]);
            }
         }

#ifndef __PRINTOUT_STATS__
      stats.GetTreeInfo(tempVectors);
#else
      stats.GetDetailedInfo(tempVectors);
      nodeCounts.Push(stats.nodeCount);
      leafCounts.Push(stats.leafCount);
      zeroCounts.Push(stats.zeroCount);
      infoScores.Push(stats.information);
#endif

      // if (!tempVectors.FindNonZero())
      //   PrintMessage("Something is wrong ...\n");

      if (stats.mean == 0.0)
         {
         // TODO -- be helpful and try to pin down problem genotypes
         if (!error_message_printed)
            PrintMessage("  Skipping Marker %-10s [BAD INHERITANCE]",
                 (const char *) ped.markerNames[markers[m]]);

         // Replace with uninformative inheritance tree
         tempVectors.MakeMinimalTree(1.0, mantra.bit_count);
         stats.information = 0.0;

         // Erase genotypes for this marker
         for (int i = family->first; i <= family->last; i++)
            ped[i].markers[markers[m]][2] = ped[i].markers[markers[m]][1] = 0;
         }

      if (stats.information > 1e-5 || twopoint ||
           InheritanceTree::mergingStrategy != MERGE_ALL &&
           tempVectors.nodes[0].type != TREE_NODE_ZERO &&
           tempVectors.nodes[0].type != TREE_NODE_LEAF ||
           GenotypeAnalysisEnabled(m, true))
         {
         // The next few lines implement the complex recombination
         // fraction of Terwilliger and Goring
         if (imaginaryTheta)
            {
            Multipoint projection;

            projection.Copy(tempVectors);
            projection.MoveAlong(mantra, imaginaryTheta, 1);
            tempVectors.Copy(projection);
            }

         singlepoint[informativeCount].Copy(tempVectors);
         StoreSinglepoint(informativeCount++);

         informativeMarkers.Push(m);
         informativePositions.Push(markerPositions[m]);
         informationScores.Push(stats.information);
         informationBounds.Push(stats.information > 0.3 ? stats.GetMinInfo(tempVectors) : 0.0);

         likelihood += tempVectors.logOffset;
         }
      else
         {
         // Likelihood for all gene flow patterns is effectively equal
         // if (stats.information > 1e-10 && !quietOutput)
         //   PrintMessage("  Skipping Marker %-10s [Info = %.6f]",
         //       (const char *) ped.markerNames[markers[m]], stats.information);

         // stats.mean as calculated by GetTreeInfo must be scaled
         if (stats.mean)
            likelihood += log(stats.mean) + mantra.bit_count * log(0.5) + tempVectors.logOffset;
         }
      }

   if (Multipoint::maximum_recombinants || twopoint || zeroRecombination ||
       mantra.bit_count < 3)
      useSparse = false;
   else
      {
      // By default we don't assume sparseness
      useSparse = false;

      // Lower information bounds if founder couple symmetries are being used
      if (mantra.couples)
         // TODO -- work out what a better penalty for couples is ...
         informationBounds.Add(-1.0 * mantra.couples / mantra.bit_count);

      // Take advantage of sparseness when there are consecutive
      // very informative markers
      for (int i = 0; i < informativeCount - 1; i++)
         if (informationBounds[i] + informationBounds[i + 1] > 1.0)
            useSparse = true;
      }

#ifdef __PRINTOUT_STATS__
   int median = nodeCounts.dim / 2;
   int ciLo = nodeCounts.dim * 25 / 1000;
   int ciHi = nodeCounts.dim * 975 / 1000;

   char * labels [] = { "Nodes", "Leafs", "Zeros", "Info" };
   Vector * vect [] = { &nodeCounts, &leafCounts, &zeroCounts, &infoScores};

   for (int i = 0; i < 4; i++)
      {
      Vector & vector = *vect[i];

      vector.Sort();
      PrintMessage("  % 8s: %.2f (%.2f to %.2f) Median = %.2f, 95%% CI = %.2f - %.2f",
         labels[i], vector.Average(), vector.Min(), vector.Max(),
         vector[median], vector[ciLo], vector[ciHi]);
      }
#endif

   cache.CloseCache();
   }

// Calculation of between marker recombination fractions
//

void MerlinCore::CalculateThetas()
   {
   if (keyThetas) return;

   // keyFemaleTheta.Dimension(informativeCount - 1);
   // keyMaleTheta.Dimension(informativeCount - 1);

   for (int m = informativeCount - 2; m >= 0; m--)
      {
      double theta[2];
      int    mrk = informativeMarkers[m];

      theta[0] = 1.0 - 2.0 * femaleMarkerTheta[mrk];
      theta[1] = 1.0 - 2.0 * maleMarkerTheta[mrk];

      while (++mrk < informativeMarkers[m + 1])
         {
         theta[0] *= 1.0 - 2.0 * femaleMarkerTheta[mrk];
         theta[1] *= 1.0 - 2.0 * maleMarkerTheta[mrk];
         }

      theta[0] = (1.0 - theta[0]) / 2;
      theta[1] = (1.0 - theta[1]) / 2;

      keyFemaleTheta[m] = theta[0];
      keyMaleTheta[m] = theta[1];
      }

   keyThetas = true;
   }

// Calculation of right conditional likelihoods
//

bool MerlinCore::ScoreConditionals()
   {
   if (zeroRecombination)
      return ZeroRecombination();

   if (twopoint) return true;

   CalculateThetas();

   right = new Tree[informativeCount];
   rightScale.Dimension(informativeCount);
   leftScale.Dimension(informativeCount);

   Multipoint engine;

   if (useSparse) engine.SetupConditioning(mantra);

   int last = informativeCount - 1;

   rightScale[last] = 1.0;

   RetrieveSinglepoint(last);
   right[last].Copy(singlepoint[last]);
   ReStoreSinglepoint(last);

   bool   check_time = (mantra.bit_count >= 16) && (maxMinutes > 0);
   double seconds = 0;  // Initialization avoids compiler warning
   time_t start_time;

   if (check_time)
      {
      seconds = maxMinutes * 60.;
      // It would be good to check for problems in time call, but
      // errno is not actually part of ANSI specification
      // errno = 0;
      time(&start_time);
      // if (errno) check_time = false;
      }

   int rescalings = 0;

   for (int m = last - 1; m >= 0; m--)
      {
      ProgressReport("Right Conditional", last - m, informativeCount);

      // Retrieve recombination fraction
      double theta[2];

      theta[0] = keyFemaleTheta[m];
      theta[1] = keyMaleTheta[m];

      // Retrieve inheritance at m+1 given m + 1, m + 2, ...
      engine.Clear();
      engine.Copy(right[m+1]);
      StoreMultipoint(m+1);

      if (useSparse && (informationBounds[m+1] + informationBounds[m] > 1.0))
         {
         // Try fast conditioning for informative locations?
         RetrieveSinglepoint(m);
         right[m].Copy(singlepoint[m]);
         ReStoreSinglepoint(m);

         // Calculate conditional probability without intermediates
         engine.Condition(right[m], theta, rightScale[m+1]);
         }
      else
         {
         // Calculate inheritance at m given m+1, m+2, ...
         engine.MoveAlong(mantra, theta, rightScale[m+1]);

         // Calculate inheritance at m given *** m ***, m + 1, ...
         RetrieveSinglepoint(m);
         engine.Multiply(singlepoint[m]);
         ReStoreSinglepoint(m);

         // Store it for future use
         right[m].Copy(engine);
         }

      // Check that the likelihood hasn't hit zero
      if (stats.GetMean(right[m]) == 0.0) return false;

      // If the tree is very unlikely, multiply by a constant to avoid underflow
      rightScale[m] = stats.mean < lowBound ? (rescalings+=m!=0, rescale) : 1.0;

      // Check if we have moved past the first analysis position
      if (lazyRightConditional && analysisPositions.First() > informativePositions[m] && m > 0)
         {
         StoreMultipoint(m);
         return true;
         }

      if (check_time)
         if (difftime(time(NULL), start_time) > seconds)
            throw OutOfTime();
      }

   if (informativeCount == 1) stats.GetMean(right[0]);
   likelihood += log(stats.mean) + rescalings * lnScale;

   StoreMultipoint(0);
   return true;
   }

// Likelihood calculation assuming no recombination
//

void MerlinCore::SortInformation(int start_index, int end_index)
   {
   if (start_index >= end_index) return;

   double average = (informationScores[informationIndex[start_index]] +
                     informationScores[informationIndex[end_index]]) * 0.5;

   int a = start_index, b = end_index;

   while (true)
      {
      while ((informationScores[informationIndex[a]] >= average) && a < b)
         a++;

      while ((informationScores[informationIndex[b]] <= average) && a < b)
         b--;

      if (a == b) break;

      informationIndex.Swap(a, b);
      }

   SortInformation(start_index, a == end_index ? a - 1 : a);
   SortInformation(a == start_index ? a + 1 : a, end_index);
   }

bool MerlinCore::ZeroRecombination()
   {
   // Identify the most informative markers
   informationIndex.Dimension(informativeCount);
   informationIndex.SetSequence(0, 1);
   SortInformation(0, informativeCount - 1);

   // Combined tree includes all nodes that are non-zero in all trees
   Tree combined;

   // Allocate storage for final result
   right = new Tree[1];
   RetrieveSinglepoint(informationIndex[0]);
   right[0].Copy(singlepoint[informationIndex[0]]);
   ReStoreSinglepoint(informationIndex[0]);

   // Keep track of how many likelihood rescalings were carried out
   int rescalings = 0;

   // Locate inheritance vectors that are compatible with all markers
   for (int i = 1; i < informativeCount; i++)
      {
      // Incorporate next most informative marker
      combined.Copy(right[0]);
      RetrieveSinglepoint(informationIndex[i]);
      combined.Multiply(singlepoint[informationIndex[i]]);
      ReStoreSinglepoint(informationIndex[i]);
      combined.TrimNoMerge();
      right[0].Copy(combined);

      // Check that is possible given zero recombination
      if (stats.GetMean(right[0]) == 0.0)
         return false;

      // Avoid underflow by multiplying by large rescale constant
      if (stats.mean < lowBound)
         {
         rescalings++;
         right[0].Multiply(0, rescale);
         }
      }

   likelihood += log(stats.mean) + rescalings * lnScale;

   StoreMultipoint(0);
   return true;
   }

// Calculation of left conditional likelihoods
//

void MerlinCore::RunLeft(Multipoint & left, Multipoint & leftMarker, int marker)
   {
   RetrieveSinglepoint(marker + 1);

   if (marker == -1)
      {
      leftMarker.Copy(singlepoint[0]);
      leftScale[0] = 1.0;
      }
   else
      {
      // First, retrieve the appropriate recombination fraction
      double theta[2];

      theta[0] = keyFemaleTheta[marker];
      theta[1] = keyMaleTheta[marker];

      // Load previous left conditioned likelihoods
      left.Copy(leftMarker);

      // Load the singlepoint likelihoods at the next location
      leftMarker.Copy(singlepoint[marker + 1]);

      if (useSparse && (informationBounds[marker] + informationBounds[marker + 1] > 1.0) &&
          !GenotypeAnalysisEnabled(marker + 1, false))
         {
         // Build left conditioned probability up to marker + 1 (inclusive)
         left.Condition(leftMarker, theta, leftScale[marker]);
         leftScale[marker + 1] = stats.GetMean(leftMarker) < lowBound ? rescale : 1.0;
         }
      else
         {
         // Do the multipoint calculation
         left.MoveAlong(mantra, theta, leftScale[marker]);

         // Get the full left conditioned probability
         leftMarker.Multiply(left);
         leftScale[marker + 1] = stats.GetMean(leftMarker) < lowBound ? rescale : 1.0;
         }
       }

   // Free up memory
   ReStoreSinglepoint(marker + 1);
   }

void MerlinCore::WalkLeft
   (Multipoint & left, Multipoint & leftMarker, Multipoint & full, int marker)
   {
   // Load required information to memory
   RetrieveMultipoint(marker + 1);
   RetrieveSinglepoint(marker + 1);

   // Calculate the new left conditioned probabilities
   if (marker != -1)
      {
      // First, retrieve the appropriate recombination fraction
      double theta[2];

      theta[0] = keyFemaleTheta[marker];
      theta[1] = keyMaleTheta[marker];

      // Load previous left conditioned likelihoods
      left.Copy(leftMarker);

      // Load the singlepoint likelihoods at the next location
      leftMarker.Copy(singlepoint[marker + 1]);

      if (useSparse && (informationBounds[marker] + informationBounds[marker + 1] > 1.0) &&
          !GenotypeAnalysisEnabled(marker + 1, false))
         {
         // Build left conditioned probability up to marker + 1 (inclusive)
         left.Condition(leftMarker, theta, leftScale[marker]);
         leftScale[marker + 1] = stats.GetMean(leftMarker) < lowBound ? rescale : 1.0;

         // TODO -- We could choose to condition on left / right ?
         // Get right conditioned probability
         full.Copy(right[marker + 1]);
         left.Condition(full, theta, leftScale[marker]);
         }
      else
         {
         // Do the multipoint calculation
         left.MoveAlong(mantra, theta, leftScale[marker]);

         // Get the full left conditioned probability
         leftMarker.Multiply(left);
         leftScale[marker + 1] = stats.GetMean(leftMarker) < lowBound ? rescale : 1.0;

         // Get the full left and right conditioned probability
         full.Copy(right[marker + 1]);
         full.Multiply(left);
         }
       }
   else
      {
      // This is for the first marker
      leftMarker.Copy(singlepoint[0]);
      full.Copy(right[0]);
      leftScale[0] = 1.0;
      }

   // Check the marker for errors
   if (GenotypeAnalysisEnabled(marker + 1, false))
      GenotypeAnalysisHandler(left, full, marker + 1);

   // Free up memory
   ReStoreSinglepoint(marker + 1);
   ReStoreMultipoint(marker + 1);
   }

// Scan chromosome and calculate likelihoods at each location for analyses
//

void MerlinCore::ScanChromosome()
   {
   // Total number of analysis locations
   int positionCount = analysisPositions.Length();

   if (positionCount == 0)
      return;

   // If no recombination ...
   if (zeroRecombination)
      {
      RetrieveMultipoint(0);
      AnalyseLocation(0, right[0]);
      ReStoreMultipoint(0);

      return;
      }

   // Special case of two point analysis
   if (twopoint)
      {
      Multipoint prior;
      prior.MakeMinimalTree(1.0, mantra.bit_count);

      for (int pos = 0; pos < positionCount; pos++)
         {
         RetrieveSinglepoint(pos);
         AnalyseLocation(pos, singlepoint[pos]);

         if (GenotypeAnalysisEnabled(pos, false))
            GenotypeAnalysisHook(singlepoint[pos], prior, pos);

         ReStoreSinglepoint(pos);
         }
      return;
      }

   // Left-conditioned probabilities up to left-marker
   Multipoint leftMarker;

   // Left and right conditionals up to current position
   Multipoint leftInheritance, inheritance;

   // Setup auxiliary variables for processing very sparse trees
   if (useSparse) leftInheritance.SetupConditioning(mantra);

   // Intermediate variables used for calculating recombination fractions
   double theta[2];

   // First we analyse any locations preceding the first informative marker
   int pos = 0;
   double firstMarkerPos = informativePositions[0];

   while (pos < positionCount && analysisPositions[pos] < firstMarkerPos)
      if (betweenMarkers)
         {
         ProgressReport("Scanning Chromosome", pos, positionCount);

         int marker = informativeMarkers[0];

         theta[0] = DistanceToRecombination(markerFemalePositions[marker] - femalePositions[pos]);
         theta[1] = DistanceToRecombination(markerMalePositions[marker] - malePositions[pos]);

         RetrieveMultipoint(0);
         inheritance.Copy(right[0]);
         inheritance.MoveAlong(mantra, theta);
         ReStoreMultipoint(0);

         AnalyseLocation(pos++, inheritance);
         }
      else
         pos++;

   // Next analyse positions, including marker locations and locations
   // between markers
   double lastMarkerPos = informativePositions.Last();

   // Anchors for starting position
   int leftAnchor  = -1;
   int rightAnchor =  0;

   // Quickly run through markers that are before the first analysis position
   if (pos == 0)
      {
      double position = analysisPositions[pos];

      while (rightAnchor + 1 < informativeCount && informativePositions[rightAnchor + 1] <= position)
         {
         RunLeft(leftInheritance, leftMarker, leftAnchor++);
         rightAnchor++;
         }
      }

   while (pos < positionCount && analysisPositions[pos] <= lastMarkerPos)
      {
      ProgressReport("Scanning Chromosome", pos, positionCount);

      double position = analysisPositions[pos];

      // Move anchors so they flank the current position
      while (rightAnchor < informativeCount &&
             informativePositions[rightAnchor] <= position)
         {
         WalkLeft(leftInheritance, leftMarker, inheritance, leftAnchor++);
         rightAnchor++;
         }

      // If we are analysing at marker, the hard work is done
      if (informativePositions[leftAnchor] == position)
          {
          // nothing to do!
          }
      // For some types of analysis, we don't need between marker information
      else if (!betweenMarkers)
         {
         pos++;
         continue;
         }
      // Otherwise we are really between markers
      else
         {
         int marker = informativeMarkers[leftAnchor];

         theta[0] = DistanceToRecombination(femalePositions[pos] - markerFemalePositions[marker]);
         theta[1] = DistanceToRecombination(malePositions[pos] - markerMalePositions[marker]);

         leftInheritance.Copy(leftMarker);
         leftInheritance.MoveAlong(mantra, theta, leftScale[leftAnchor]);

         marker = informativeMarkers[rightAnchor];

         theta[0] = DistanceToRecombination(markerFemalePositions[marker] - femalePositions[pos]);
         theta[1] = DistanceToRecombination(markerMalePositions[marker] - malePositions[pos]);

         RetrieveMultipoint(rightAnchor);
         inheritance.Copy(right[rightAnchor]);
         inheritance.MoveAlong(mantra, theta, rightScale[rightAnchor]);
         inheritance.Multiply(leftInheritance);
         ReStoreMultipoint(rightAnchor);
         }

      AnalyseLocation(pos++, inheritance);
      }

   if (pos < positionCount)
      while (rightAnchor++ < informativeCount)
         WalkLeft(leftInheritance, leftMarker, inheritance, leftAnchor++);

   // Finalize analyse any locations beyond the last informative marker
   if (betweenMarkers)
      while (pos < positionCount)
         {
         ProgressReport("Scanning Chromosome", pos, positionCount);

         int marker = informativeMarkers.Last();

         theta[0] = DistanceToRecombination(femalePositions[pos] - markerFemalePositions[marker]);
         theta[1] = DistanceToRecombination(malePositions[pos] - markerMalePositions[marker]);

         inheritance.Copy(leftMarker);
         inheritance.MoveAlong(mantra, theta);

         AnalyseLocation(pos++, inheritance);
         }
   }

// Output buffering functions
//

void MerlinCore::PrintHeader()
   {
   printHeader = false;

   printf("Family: %5s - Founders: %-2d - Descendants: %-2d - Bits: %-2d\n",
          (const char *) family->famid, mantra.f, mantra.n - mantra.f,
          mantra.bit_count);
   }

void MerlinCore::PrintMessage(const char * msg, ...)
   {
   if (printHeader) PrintHeader();

   va_list ap;

   va_start (ap, msg);
   vprintf(msg, ap);
   printf("\n");
   va_end(ap);
   }

void MerlinCore::ProgressReport(const char * msg, int done, int total)
   {
   if (quietOutput || msg == NULL) return;

   long check = (0x10000 >> mantra.bit_count) - 1;
   if (mantra.bit_count > 16 || ((done & check ) == check))
      {
      if (printHeader)
         printf("Family %s (%d bits) -",
               (const char *) family->famid, mantra.bit_count);
      printf("  %s: %3d%%%20s\r", msg, done * 100 / total, "");
      fflush(stdout);
      cleanOutput = true;
      }
   }

void MerlinCore::CleanMessages()
   {
   if (cleanOutput)
      {
      printf("%70s\r", "");
      fflush(stdout);
      }

   if (!printHeader)
      printf("\n");
   }

// Enforce maximum number of alleles on markers
//

void MerlinCore::ValidateMarker(MarkerInfo * info)
   {
#ifdef __USE_LONG_INT
   if (info->freq.dim > 64)
       error("Please recode marker %s (and possibly others), which has %d alleles\n"
             "This version of merlin supports up to 63 alleles\n",
             (const char *) info->name, info->freq.dim - 1);
#else
   if (info->freq.dim > 32)
       error("Please recode marker %s (and possibly others), which has %d alleles\n"
             "This version of merlin handles 31 alleles or fewer per marker\n\n"
             "Alternatively, you could recompile Merlin with the __USE_LONG_INT\n"
             "option, which enables support for up to 63 alleles\n\n"
             "See the Merlin Makefile for details.\n",
             (const char *) info->name, info->freq.dim - 1);
#endif
   }

// Stub functions that should be replaced by derived classes

void MerlinCore::AnalyseLocation(int, Tree &) { }
void MerlinCore::GenotypeAnalysisHook(Tree &, Tree &, int ) { }

// This function prepares structures to be handled by the genotype inference
// and error detection handlers, which require of evaluation of different
// likelihoods, before and after small modifications to available genotype
// data
//

void MerlinCore::GenotypeAnalysisHandler
         (Multipoint & left, Multipoint & withMarker, int informativeMarker)
   {
   Multipoint withoutMarker;

   if (++informativeMarker != informativeCount)
      {
      double theta[2];

      theta[0] = keyFemaleTheta[informativeMarker - 1];
      theta[1] = keyMaleTheta[informativeMarker - 1];

      RetrieveMultipoint(informativeMarker);
      withoutMarker.Copy(right[informativeMarker]);
      ReStoreMultipoint(informativeMarker);
      withoutMarker.MoveAlong(mantra, theta, rightScale[informativeMarker]);

      if (informativeMarker != 1)
         withoutMarker.Multiply(left);
      }

   if (informativeCount == 1)
      left.MakeMinimalTree(1.0, left.bit_count);

   GenotypeAnalysisHook(withMarker,
                        informativeMarker == informativeCount ? left : withoutMarker,
                        informativeMarker - 1);
   }

// This function checks whether a genotype based analysis is enabled for the
// current marker -- if it is, we don't skip over uninformative markers

bool MerlinCore::GenotypeAnalysisEnabled(int, bool )
   {
   return false;
   }

// Routines for swapping singlepoint results in and out of memory
// (when the swap file exceeds 1 GB, these routines discard
// and recalculate singlepoint results as needed).

void MerlinCore::StoreSinglepoint(int position)
   {
   // Check if swap file is open
   if (!memoryManagement)
      return;

   // Update list of positions stored in file
   if (position == 0)
      lastSinglepointSwap = (useSwap && !smallSwap) ? markerCount : -1;

   // Check if tree is complex enough to merit swapping
   if (!singlepoint[position].IsSwappable())
      return;

   // Check if swapping is enabled for the present location
   if (position > lastSinglepointSwap)
      {
      singlepoint[position].Discard();
      return;
      }

   // If everything checks out, then swap the current location
   singlepoint[position].Pack(singlepointSwap);

   // Check if swap file is filled up (i.e. exceeds 1 GB)
   double usage = singlepointSwap.GetFileUsage();

   if (usage < 0.0 || usage > 1024 * 1024 * 1024)
      lastSinglepointSwap = position;
   }

void MerlinCore::RetrieveSinglepoint(int position)
   {
   // Check if swap file is open
   if (!memoryManagement)
      return;

   // Check if tree is complex enough to merit swapping
   if (singlepoint[position].IsInMemory())
      return;

   // Check if swapping is enabled for the present location
   if (position > lastSinglepointSwap)
      {
      // If the swap file filled up, we may need to recalculate some numbers
      MerlinCache cache;
      FuzzyInheritanceTree engine;

      // Check if we can retrieve data from a cache ...
      cache.OpenCache(mantra);

      // Lookup the id for the marker of interest
      int  markerid = markers[informativeMarkers[position]];

      // Check if it is part of a cluster of markers in LD
      if (clusters.Enabled() && clusters.markerToCluster[markerid] != NULL)
         {
         MarkerCluster * cluster = clusters.markerToCluster[markerid];

         // Score using clustering algorithms
         cluster->ScoreLikelihood(mantra, engine);
         }
      else
         {
         mantra.SelectMarker(markerid);

         if (!cache.RetrieveFromCache(mantra, engine, Pedigree::markerNames[markerid]))
            engine.FuzzyScoreVectors(mantra);
         }

      // Apply Terwilliger and Goring's complex recombination fraction
      if (imaginaryTheta)
         {
         Multipoint projection;

         projection.Copy(engine);
         projection.MoveAlong(mantra, imaginaryTheta, 1);
         singlepoint[position].Copy(projection);
         }
      else
         singlepoint[position].Copy(engine);
      }
   else
      singlepoint[position].UnPack(singlepointSwap);
   }

void MerlinCore::ReStoreSinglepoint(int position)
   {
   // Check if swap file is open
   if (!memoryManagement)
      return;

   // Check if tree is complex enough to merit swapping
   if (!singlepoint[position].IsSwappable())
      return;

   // Check if swapping is enabled for the present location
   if (position > lastSinglepointSwap)
      singlepoint[position].Discard();
   else
      singlepoint[position].RePack();
   }

void MerlinCore::StoreMultipoint(int position)
   {
   // Check if swap file is open
   if (!memoryManagement)
      return;

   // Update list of positions stored in file
   if (position == informativeCount - 1)
      {
      multipointSwapStart = 1;
      multipointSwapEnd = (useSwap && !smallSwap) ? informativeCount - 1 : multipointGrid - 1;
      }

   // Check if tree is complex enough to merit swapping
   if (!right[position].IsSwappable())
      return;

   // We always swap grid points, so as to keep recalculation
   // times manageable
   if (position % multipointGrid == 0 || position == informativeCount - 1 || zeroRecombination)
      {
      right[position].Pack(scaffoldSwap);
      return;
      }

   // Otherwise, we dynamically enable and disable swapping to preserve memory
   if (position > multipointSwapEnd)
      {
      right[position].Discard();
      return;
      }

   // If everything checks out, then swap the current location
   right[position].Pack(multipointSwap);

   // If we are already in economy mode, nothing else we can do
   if (multipointSwapEnd == multipointGrid - 1)
      return;

   // Check if swap file is filled up (i.e. exceeds 1 GB)
   double usage = multipointSwap.GetFileUsage();

   // If so, we discard contents unless we are "nearly done" ...
   if ((usage < 0.0 || usage > 1024 * 1024 * 1024) && position > multipointGrid * 2)
      {
      // Discard previously swapped results
      multipointSwap.Free();
      for (int i = position; i < informativeCount - 1; i++)
         if (right[i].IsInSwapFile() && !(i % multipointGrid == 0))
            right[i].Discard();

      // Only enable swapping for a small chunk of data
      multipointSwapStart = 0;
      multipointSwapEnd = multipointGrid - 1;
      }
   }

void MerlinCore::RetrieveMultipoint(int position)
   {
   // Check if swap file is open
   if (!memoryManagement)
      return;

   // Check if tree is complex enough to merit swapping
   if (right[position].IsInMemory())
      return;

   // Grid points are always available from a separate swap file
   if (position % multipointGrid == 0 || position == informativeCount - 1 || zeroRecombination)
      {
      right[position].UnPack(scaffoldSwap);
      return;
      }

   // Check if swapping is enabled for the present location
   if (position >= multipointSwapStart && position <= multipointSwapEnd)
      {
      right[position].UnPack(multipointSwap);
      return;
      }

   // Otherwise, we first reset the swap file state
   multipointSwap.Free();

   for (int i = multipointSwapStart; i <= multipointSwapEnd; i++)
      if (right[i].IsSwappable())
         right[i].Discard();

   // And then proceed to recalculate right conditional probabilities
   // for a small window of markers
   bool repack = false;
   multipointSwapStart = position;
   multipointSwapEnd = position + 1;

   while (!right[multipointSwapEnd].IsInMemory() &&
          !right[multipointSwapEnd].IsInSwapFile())
      multipointSwapEnd++;

   if (right[multipointSwapEnd].IsInSwapFile())
      {
      right[multipointSwapEnd].UnPack(scaffoldSwap);
      repack = true;
      }

   // The following section repeats the multipoint calculation as needed
   Multipoint engine;

   if (useSparse) engine.SetupConditioning(mantra);

   for (int m = --multipointSwapEnd; m >= multipointSwapStart; m--)
      {
      // Retrieve recombination fraction
      double theta[2];

      theta[0] = keyFemaleTheta[m];
      theta[1] = keyMaleTheta[m];

      // Retrieve inheritance at m+1 given m + 1, m + 2, ...
      engine.Clear();
      engine.Copy(right[m+1]);
      if (m < multipointSwapEnd)
         right[m+1].Pack(multipointSwap);

      if (useSparse && (informationBounds[m+1] + informationBounds[m] > 1.0))
         {
         // Try fast conditioning for informative locations?
         RetrieveSinglepoint(m);
         right[m].Copy(singlepoint[m]);
         ReStoreSinglepoint(m);

         // Calculate conditional probability without intermediates
         engine.Condition(right[m], theta, rightScale[m+1]);
         }
      else
         {
         // Calculate inheritance at m given m+1, m+2, ...
         engine.MoveAlong(mantra, theta, rightScale[m+1]);

         // Calculate inheritance at m given *** m ***, m + 1, ...
         RetrieveSinglepoint(m);
         engine.Multiply(singlepoint[m]);
         ReStoreSinglepoint(m);

         // Store it for future use
         right[m].Copy(engine);
         }
      }

   // Free memory, if appropriate
   if (repack)
      right[multipointSwapEnd + 1].RePack();

   // The PackOnly function keeps the tree in memory after writing swap file
   right[position].PackOnly(multipointSwap);
   }

void MerlinCore::ReStoreMultipoint(int position)
   {
   // Check if swap file is open
   if (multipointSwap.skeleton == NULL || scaffoldSwap.skeleton == NULL)
      return;

   // Check if tree is complex enough to merit swapping
   if (!right[position].IsSwappable())
      return;

   // Discard information at the current position
   right[position].RePack();
   }

void MerlinCore::FreeSwap()
   {
   BasicTree::FreeSwap();

   singlepointSwap.Free();
   multipointSwap.Free();
   scaffoldSwap.Free();
   }
 
