////////////////////////////////////////////////////////////////////// 
// merlin/Houdini.h 
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
 
#ifndef __HOUDINI_H__
#define __HOUDINI_H__

#include "MathMatrix.h"
#include "Magic.h"

class FuzzyInheritanceTree : public InheritanceTree
   {
   public:
      // Construct and destructor look after memory allocation
      FuzzyInheritanceTree();
      virtual ~FuzzyInheritanceTree();

      // Scores the likelihoods of all inheritance vectors
      void FuzzyScoreVectors(Mantra & m);

      // Calculates the likelihood of a single inheritance vector
      // and stores the most likely state for each allele in the
      // alleles vector
      double  FuzzyEvaluateInheritanceVector(Mantra & m, int * vector,
              int * graph, longint * alleles, longint * alleles2, int * fixed);

      // The following two variables define the two available genotype
      // error models
      static double perAlleleError, perGenotypeError;

   protected:
      // Short cuts for situations where fast alternatives exist ...
      virtual bool ShortScoreVectors(Mantra & m);
      virtual double ShortEvaluateVector(
              Mantra & m, int * vector, int * graph,
              longint * alleles, longint * alleles2, int * fixed);

      virtual void FillPenetranceMatrix(Mantra & m);

      // This stores the probability of the observed phenotypes for
      // each individual conditional on each possible genotype
      Matrix penetrances;

      // This stores the base likelihood for each family. It is useful
      // to try and guess the modal likelihood for all inheritance vectors
      // so as to avoid recalculating it.
      double priorLikelihood;

      // Recurse throught the inheritance vector space and gradually update
      // founder allele graphs so as to calculate all likelihoods in a timely
      // manner
      int ScoreRecursive(Mantra & m, int current_bit, int last_pivot);
      int FindRedundancy(Mantra & m, int pivot);

      // This function allows us to deal with cases where the founder
      // couple symmetry has been applied to individuals whose phenotypes
      // differ
      void UnravelAssymmetries(Mantra & m);

      // These functions help us to unravel founder couple symmetries when necessary
      virtual bool IdenticalPenetrances(int founder1, int founder2);
      virtual void SwapPenetrances(int founder1, int founder2);

      // These arrays store information about founder allele graphs
      IntArray graph, weight;       // Information about connected components
      IntArray graph_descendants;   // Information required for early evaluation of inheritance graphs
      IntArray graph_journal;       // Information required to undo changes in graph components
      Vector   graph_likelihood;    // Likelihoods for each inheritance graph component

      // A list of phenotyped descendants for each founder allele
      IntArray * carriers;

      // Connections in the founder allele graph are also stored in
      // a matrix representation and as an array of links
      IntArray * matrix, * links;

      // These functions track founder allele graphs
      int  ProcessConnection(int founder1, int founder2);
      void RemoveConnection(int founder1, int founder2);
      int  FindComponent(int allele);
      int  JoinAlleles(int founder1, int founder2);
      void UndoGraph();

      // These arrays track alleles for each founder
      IntArray * possibleAlleles, * alleles;

      // This variable determines whether a meta-allele, grouping all unobserved
      // alleles is created
      bool metaAllele;

      // This array tracks frequencies of undistinguishable alleles for each founder
      Vector metaAlleles;

      // These functions track alleles for each founder
      void ShowAllele(int founder, int allele);
      void HideAllele(int founder, int allele);
      virtual void UpdateAlleles(int founder, int carrier);
      virtual void UndoAlleles(int founder, int carrier);

      // Return the number of possible alleles for the problem at hand
      virtual int CountAlleles(Mantra & m);

      // This function allocates sufficient memory to process the current pedigree
      virtual void AllocateMemory(Mantra & m);

      // Variable tracking the amount of allocated memory
      int allocated_founders;

      // Information about phenotyped descendants for each individual
      IntArray isTyped;          // True for phenotyped individuals, false otherwise
      IntArray typedDescendants; // Number of phenotyped descendants for each allele

      // These functions setup the previous two arrays
      void PreparePhenotypes(Mantra & m);
      virtual bool isPhenotyped(Mantra & m, int who);

      // Likelihood information
      double CalculateLikelihood(Mantra & m);

      // Peeling routines
      bool   PeelComponent(Mantra & m, int component);
      int    PeelFounder(Mantra & m, int founder, int order);
      void   PeelBranch(Mantra & m, int branchStart, bool partial = true);

      // Routine for finding the most likely set of allele states
      void   FindModalState(Mantra & m, int component, longint * alleles);

      // Routines for looping through possible allele sets
      void   InitializeAlleleSet(int branch);
      bool   IncrementAlleleSet(int branch);
      void   ResetAlleleStates(int branch);

      // Variables used during peeling
      IntArray peelSet;       // Current set of alleles to be peeled
      IntArray visited;       // Identifies founder alleles processed so far
      IntArray visitIndex;    // Order in which founder alleles visited
      IntArray phenoVisits;   // Identifies genotypes processed so far
      IntArray phenoList;     // Genotypes currently being processed
      int      peelKey;       // Unique ID for each round of peeling

      // Allele frequency and related information
      Vector * condition;
      Vector * freq;

      // Product of allele frequencies for the first k founder alleles
      Vector product;

      // Allele state information during peeling
      IntArray alleleState;
      IntArray alleleStateIndex;

      // These function retrieve key probabilities
      virtual double GetAlleleFrequency(int founder, int allele);
      virtual double JointProbability(Mantra & m);

      // These functions initialize and destroy the initial allele list
      // for each founder
      virtual void SetupAlleleList(Mantra & );
      virtual void CleanUpAlleleList(Mantra & );

   private:
      // Pointer to allele frequency information
      Vector * allele_frequencies;

      // These functions fill in the penetrance matrix using a
      // specific error model
      void AlleleErrorProbs(Mantra & m);
      void GenotypeErrorProbs(Mantra & m);

      // This function extracts two integer alleles from a bit-coded genotype
      void ExtractAlleles(longint geno, int & allele1, int & allele2);

      // These arrays store information on the observed alleles for
      // each individual
      IntArray observedHi, observedLo;

      // Functions of convenience for helping memory allocator
      void ReallocateArray(IntArray * & ptr, int old_size, int new_size);
      void ReallocateVector(Vector * & ptr, int old_size, int new_size);
   };

#endif
 
