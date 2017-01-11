////////////////////////////////////////////////////////////////////// 
// extras/pedwipe.cpp 
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
#include "Parameters.h"
#include "QuickIndex.h"

int main(int argc, char * argv[])
   {
   printf("PedWipe - (c) 2000 Goncalo Abecasis\n"
          "Automatically wipe out genotypes from a pedigree file\n\n");

   String pedfile("merlin.ped");
   String datafile("merlin.dat");
   String errorfile("merlin.err");

   bool showTallies = false;

   ParameterList pl;

   pl.Add(new StringParameter('d', "Data File", datafile));
   pl.Add(new StringParameter('p', "Pedigree File", pedfile));
   pl.Add(new StringParameter('e', "Errors File", errorfile));
   pl.Add(new SwitchParameter('t', "Show Tallies", showTallies));

   pl.Read(argc, argv);
   pl.Status();

   Pedigree ped;

   ped.Prepare(datafile);
   ped.Load(pedfile);

   StringArray errors, tokens;
   errors.Read(errorfile);

   int count = 0;
   StringIntMap perMarker, perFamily, perPerson;

   for (int i = 1; i < errors.Length(); i++)
      {
      tokens.Clear();
      tokens.AddTokens(errors[i]);

      if (tokens.Length() < 3) continue;

      Person * person = ped.FindPerson(tokens[0], tokens[1]);

      int markerid = ped.LookupMarker(tokens[2]);

      if (person == NULL)
         {
         printf("Person %s.%s not found ... \n",
                (const char *) tokens[0], (const char *) tokens[1]);
         continue;
         }

      if (markerid == -1)
         {
         printf("Marker %s not found ... \n",
             (const char *) tokens[2]);
         continue;
         }

      printf("Person %s.%s, marker %s wiped.\n",
         (const char *) tokens[0], (const char *) tokens[1],
         (const char *) tokens[2]);

      person->markers[markerid].one = 0;
      person->markers[markerid].two = 0;

      perPerson.IncrementCount(tokens[0] + "." + tokens[1]);
      perFamily.IncrementCount(tokens[0]);
      perMarker.IncrementCount(tokens[2]);

      count++;
      }

   if (perMarker.Length() == 0)
      {
      printf("No errors found in merlin.err\n");
      }
   else if (showTallies && count)
      {
      printf("\nSummary of Errors\n");
      printf("=================\n\n");

      QuickIndex index;

      printf("Per Marker:  (average = %.2f)\n"
             "-----------------------------\n",
            (double) count / (double) ped.markerCount);
      index.IndexCounts(perMarker);
      index.Reverse();
      for (int i = 0; i < perMarker.Length(); i++)
         printf(" %3d errors for marker %s\n",
                perMarker.GetCount(index[i]),
                (const char *) perMarker[index[i]]);

      printf("\nPer Family: (average = %.2f)\n"
             "----------------------------\n",
            (double) count / (double) ped.familyCount);
      index.IndexCounts(perFamily);
      index.Reverse();
      for (int i = 0; i < perFamily.Length(); i++)
         printf(" %3d errors for family %s\n",
                perFamily.GetCount(index[i]), (const char *) perFamily[index[i]]);

      printf("\nPer Person: (average = %.2f)\n"
             "----------------------------\n",
            (double) count / (double) ped.count);
      index.IndexCounts(perPerson);
      index.Reverse();
      for (int i = 0; i < perPerson.Length(); i++)
         printf(" %3d errors for person %s\n",
               perPerson.GetCount(index[i]), (const char *) perPerson[index[i]]);
      }

   printf("\nWriting out edited files [wiped.*] ...\n\n");

   if (ped.markerInfoCount)
      {
      ped.WriteMapFile("wiped.map");
      ped.WriteFreqFile("wiped.freq");
      }

   ped.WriteDataFile("wiped.dat");
   ped.WritePedigreeFile("wiped.ped");
   }
 
