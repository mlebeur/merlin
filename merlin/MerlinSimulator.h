////////////////////////////////////////////////////////////////////// 
// merlin/MerlinSimulator.h 
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
 
#ifndef __MERLIN_SIMULATOR_H__
#define __MERLIN_SIMULATOR_H__

#include "Pedigree.h"
#include "ParametricLikelihood.h"

class Simulator
   {
   public:
      // Simulate
      static void Simulate(Pedigree & ped, IntArray & markers);

      // User interface functions
      static void NullData(Pedigree & ped, IntArray & markers);
      static void Alternative(Pedigree & ped, IntArray & markers, DisModel & model,
                         double positionFemale, double positionMale);
      static void SimulateQTL(Pedigree & ped, Family & family, IntArray & alleles);

      // Location for simulated locus
      static String model;

      // Pedigree containing founder haplotypes
      // to be used in gene-drops
      static String   reuseLabel;
      static Pedigree reuseSource;
      static int      reuseCounter;

   private:
      // Core simulation modules
      static void SampleInheritanceVector(Tree & tree, int node, int bit,
                                   IntArray & current, IntArray & best,
                                   double & sum, double weight = 1.0);
      static void SimulateFamily(Pedigree & ped, Family & family, IntArray & markers, IntArray & vector);

      static void LoadReuseableHaplotypes();

      // Variables used for simulating QTL loci
      static int    qtlMarker;
      static int    qtlTrait;
      static double sQtl, sPolygenes, sEnvironment;
   };


#endif

 
