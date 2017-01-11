////////////////////////////////////////////////////////////////////// 
// extras/pedmerge.cpp 
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
 
#include "Pedigree.h"
#include "Error.h"

#include "string.h"
#include "stdlib.h"

int main(int argc, char * argv[])
   {
   printf("PedMerge - Pedigree Merge (c) 1999 Goncalo Abecasis\n\n");

   if (argc < 2)
      {
      printf("Usage: pedmerge input1 input2 ... output\n\n"
         "This program will try to merge a set of paired pedigree (.ped)\n"
         "and data (.dat) files into a single composite pedigree.\n\n"
         "For example:\n\n"
         "    > pedmerge a b c\n\n"
         "Will create the files c.dat and c.ped including all the phenotype\n"
         "data and individuals in a.dat, a.ped, b.dat and b.ped.\n\n"
         "WARNING: pedmerge will overwrite output files without checking\n\n");
      exit(0);
      }

   PedigreeDescription * pd = new PedigreeDescription[argc];
   Pedigree ped;
   String   filename;

   for (int i = 1; i < argc - 1; i++)
      {
      filename = argv[i];
      filename += ".dat";
      printf("Reading data file %s ...\n", (const char *) filename);
      ped.Prepare(filename);
      pd[i] = ped.pd;
      }

   for (int i = 1; i < argc - 1; i++)
      {
      filename = argv[i];
      filename += ".ped";
      printf("Reading pedigree file %s ...\n", (const char *) filename);
      ped.pd = pd[i];
      ped.Load(filename);
      }

   if (ped.MarkerPositionsAvailable())
      {
      filename = argv[argc - 1];
      filename += ".map";
      printf("Writing map file %s ...\n"
             "   * CHECK CAREFULLY, map positions may not merge correctly\n",
             (const char *) filename);
      ped.WriteMapFile(filename);
      }

   if (ped.AlleleFrequenciesAvailable())
      {
      filename = argv[argc - 1];
      filename += ".freq";
      printf("Writing frequency file %s ...\n", (const char *) filename);
      ped.WriteFreqFile(filename);
      }

   filename = argv[argc - 1];
   filename += ".dat";
   printf("\nWriting data file %s ...\n", (const char *) filename);
   ped.WriteDataFile(filename);

   filename = argv[argc - 1];
   filename += ".ped";
   printf("Writing pedigree file %s ...\n", (const char *) filename);
   ped.WritePedigreeFile(filename);

   delete [] pd;

   return 0;
   }

 
