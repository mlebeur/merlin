////////////////////////////////////////////////////////////////////// 
// clusters/HaploFamily.h 
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
 
#ifndef __HAPLOFAMILY_H__
#define __HAPLOFAMILY_H__

#include "BasicHash.h"
#include "Mantra.h"
#include "Magic.h"

class FounderGraph;

class FamilyHaplos
   {
   public:
      String         label;
      IntArray       markers;
      FounderGraph * graphs;

      // Information on errors
      IntArray       mendel_errors;
      bool           obligate_recombinant;

      // Tuning parameters
      static int     maxBits;

      FamilyHaplos();
      ~FamilyHaplos();

      bool Haplotype(Pedigree & ped, Family & f);
      bool Haplotype(Pedigree & ped, Person & p);

      bool Haplotype(Tree & joint, Pedigree & ped, Family & f);
      bool Haplotype(Tree & joint, Pedigree & ped, Person & p);
      bool HaplotypeMaleX(Tree & joint, Pedigree & ped, Person & p);

      String & RetrieveGenotypes(String & ped, Family & f);
      String & RetrieveGenotypes(String & ped, Person & p);

      bool Haplotype(Pedigree & ped, Family & f, IntArray & vector);

      // Retrieves likelihoods for each node in a tree indexed with ListAllRecursively
      void RetrieveLikelihoods(Tree & tree);
      void RetrieveLikelihoods(Tree & tree, int node, int bit);

   private:
      Mantra m;
      InheritanceTree engine;

      // Memory management
      void DeleteGraphs();

      // Functions for flattening inheritance trees
      void Flatten();

      // Functions used for iterating through inheritance space
      void ListAllRecursively(Tree & tree, int node, int bit, double prob);
      bool RetrieveGraph(bool canFail = false);

      // Temporary arrays used for iterating through inheritance space
      IntArray  * graph, * fixed, vector;
      LongArray * alleles, * alleles2;

      // Memory management functions
      void FreeMemory();
      void AllocateMemory(int haplos, int bits);

      int allocatedMarkers;

      // Lookup table for facilitating processing of large families
      // with many possible founder haplotype sets
      BasicHash table;
   };

class FounderGraph
   {
   public:
      int            markers;
      IntArray     * graph;
      LongArray    * alleles, * alleles2;
      double         weight;
      double       * likelihood;
      FounderGraph * next;

      FounderGraph(int mrk);
      ~FounderGraph();

      FounderGraph * Append(IntArray * graph, IntArray * fixed, LongArray * alleles,
                            LongArray * alleles2, double weight, BasicHash * table = NULL);

      bool Compare(IntArray * graph, LongArray * alleles, LongArray * alleles2);
      int  Hash(IntArray * graph, LongArray * alleles, LongArray * alleles2);

      int  Length();

      int  GetAllele(int mrk, int founder);
      int  GetAllele2(int mrk, int founder);
   };


#endif

 
