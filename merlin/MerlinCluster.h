////////////////////////////////////////////////////////////////////// 
// merlin/MerlinCluster.h 
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
 
#ifndef __MERLINCLUSTER_H__
#define __MERLINCLUSTER_H__

#include "Pedigree.h"
#include "Mantra.h"
#include "Tree.h"

#define  CLUSTER_SEARCH_MAX   20

class MarkerCluster
   {
   public:
      // Ids for each marker in the cluster
      IntArray markerIds;
      IntArray alleleCounts;

      // Allele states for each haplotype
      //   states 1 .. m for haplotype 1
      //   states m + 1 .. 2m for haplotype 2
      IntArray alleles;
      Vector   freqs;

      // Error message for haplotyping
      String   errormsg;

      void AddHaplotype(IntArray & al, double freq);
      int  FindHaplotype(IntArray & al);
      void RetrieveHaplotype(IntArray & al, int index);
      int  RetrieveAllele(int index, int marker) { return alleles[index * markerIds.Length() + marker]; }

      int  CountHaplotypes()   { return freqs.Length();     }
      int  CountMarkers()      { return markerIds.Length(); }

      // Check that allele frequencies sum to 1.0
      void CheckFrequencies();

      // Update allele counts for each marker
      void UpdateAlleleCounts();

      // Estimate allele frequencies this cluster
      bool EstimateFrequencies(Pedigree & ped, String * messages);
      void UpdateAlleleFrequencies();
      static void EstimateAlleleFrequencies(Pedigree & ped, const char * logname);

      // Score likelihoods of observed genotypes for a particular cluster
      // using the current set of haplotype frequency estimates
      void ScoreLikelihood(Mantra & m, Tree & tree);

      // Sample one haplotype configuration, or retrieve the most likely
      // configuration
      void GetHaplotypes(Mantra & m, IntArray & vector, IntArray & haplotypes, bool sample);
   };

class MarkerClusterList : public MarkerCluster
   {
   public:
      MarkerClusterList * next;

      MarkerClusterList()
         { next = NULL; }

      MarkerClusterList(MarkerClusterList * nextNode)
         { next = nextNode; }
   };

class MarkerClusters
   {
   public:
      MarkerClusters();
      ~MarkerClusters();

      // Marker clusters are maintained as a linked list
      MarkerClusterList * head;

      // An auxiliary array maps marker ids to clusters
      MarkerClusterList ** markerToCluster;

      // Functions to retrieve cluster information for individual markers
      bool IsClustered(int markerId)
         { return (markerToCluster != NULL) && (markerToCluster[markerId] != NULL); }
      MarkerCluster * GetCluster(int markerId)
         { return (MarkerCluster *) markerToCluster[markerId]; }

      // Delete linked list of clusters
      void DeleteList();

      // Load cluster information from a file
      void LoadFromFile(const char * filename, const char * logname = NULL);
      void LoadFromFile(IFILE & input, const char * logname = NULL);

      // Save cluster information to a file
      int SaveToFile(const char * filename);
      int SaveToFile(FILE * output);

      // Automatically cluster markers according to the distance between them
      void ClusterByDistance(double minDistance);

      // Automatically cluster markers according to pairwise rsq
      void ClusterByRsquared(Pedigree & ped, double minRsquared);

      // Allocate cross reference tables for relating markers to clusters, if necessary
      void AllocateCrossReference();

      // Check if clustering has been activated
      bool Enabled() { return (count > 0); }

      // Estimate allele frequencies for all clusters (if required)
      void EstimateFrequencies(Pedigree & ped, const char * logfile);

      // Prints a blank line after all other messages
      void FinishOutput();

   protected:
      void CreateCluster(IntArray & markers, int clusterStart, int clusterEnd);
      void AdjustPositions(IntArray & markers);

      int  count;
      bool haveOutput;
   };

extern MarkerClusters clusters;

#endif

 
