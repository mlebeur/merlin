////////////////////////////////////////////////////////////////////// 
// regress/RegressParameters.h 
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
 
#ifndef __MERLINPARAMETERS_H__
#define __MERLINPARAMETERS_H__

#include "StringBasics.h"
#include "Parameters.h"

class RegressionParameters : public ParameterList
   {
   public:
      String pedigreeFile;
      String dataFile;
      String mapFile;
      String freqFile;
      int    alleleFrequencies;
      int    random_seed;

      static int  maxMegabytes;
      static bool trimPedigree;

      static bool fitErrorRate;

      static bool simulateNull;
      static int  reruns;

      static bool writeFrequencyFile;
      static bool inverseNormal;

      static double heritability;
      static double variance;
      static double mean;

      // Marker Clustering
      static String clusterFile;
      static double clusterDistance;
      static double clusterRsquared;

      RegressionParameters();

      virtual void Read(int argc, char ** argv, int start = 1);
      void Check();
   };

#define  FREQ_ML        99   // For maximum likelihood allele frequencies

class AlleleFrequencyParameter : public Parameter
   {
   public:
   AlleleFrequencyParameter(char c, const char * desc, int & how, String & file);

   virtual void Status();

   protected:
      String * filename;
      String status;

      virtual void Translate(const char * value);
   };


#endif

 
