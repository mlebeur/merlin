////////////////////////////////////////////////////////////////////// 
// merlin/Parametric.h 
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
 
#ifndef __PARAMETRIC_H__
#define __PARAMETRIC_H__

#include "Magic.h"

class DiseaseModel
   {
   public:
      // Affection status to which model applies
      int affection;

      // Allele frequency for the disease allele (-)
      double p;

      // The three penetrances (+/+, +/-, -/-)
      double pen[3];

      DiseaseModel() { }
      DiseaseModel(DiseaseModel & init) { Copy(init); }

      DiseaseModel & operator = (DiseaseModel & rhs)
         { Copy(rhs); return *this; }

      void Copy(DiseaseModel & source);
      bool CheckModel();
   };

class LikelihoodTree : public Tree
   {
   public:
      LikelihoodTree();
      ~LikelihoodTree();

      // Scores the likelihoods of all inheritance vectors
      void ScoreModel(Mantra & m);

      // The frequency of the disease allele
      DiseaseModel model;

   private:
      // Stacks for the most commonly used arrays
      int stack_depth;

      // Routines for managing temporary storage
      void AllocateStack(Mantra & m);
      void FreeStack();

      // Recurse throught the pedigree and store the likelihood of all
      // possible inheritance vectors
      int ScoreRecursive(Mantra & m, int current_bit, int last_pivot,
               IntArray * graph, IntArray * graph_phenos,
               int freec);

      // Check if branching is required (while recursing through the pedigree)
      int FindRedundancy(Mantra & m, int pivot);

      // Calculate the likelihood for a single inheritance vector
      void   UpdateLikelihood(Mantra & m, IntArray & graph, IntArray & graph_phenos, int pivot);
      double CalculateLikelihood(Mantra & m, IntArray & graph);

      // Calculate likelihood for a component of the inheritance graph
      void   PrepareComponent(Mantra & m, IntArray & graph, int component, int pivot);
      double ComponentLikelihood(Mantra & m, IntArray & graph, int component);
      double PartialLikelihood(Mantra & m, IntArray & graph, int component);

      // Allelic state information for each component
      IntArray allele;
      IntArray activeComponent, mappedComponent;

      // Connected components in founder allele graph
      IntArray * components;
      IntArray * graphs;

      // Structures to allow early calculation of graph components
      IntArray * graph_phenotypes;
      Vector     graph_likelihood;

      // Temporary storage for ComponentLikelihood and PartialLikelihood
      IntArray cPheno, cState;
   };


#endif
 
