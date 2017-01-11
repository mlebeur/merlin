////////////////////////////////////////////////////////////////////// 
// extras/hapmapConverter.cpp 
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
 
#include "StringMap.h"
#include "StringArray.h"
#include "IntArray.h"
#include "Error.h"
#include "Parameters.h"

int main(int argc, char ** argv)
   {
   String templateFile, genotypeFile;
   String mapFile("mapfile"), pedFile("pedfile"), datFile("datfile");
   bool coriellIds = false;

   printf("hapmapConverter -- (c) 2004-2007 Goncalo Abecasis\n\n");
   printf("This program converts genotype files downloaded from the HapMap website\n"
          "into MERLIN format. Sample pedigree template and genotype files are \n"
          "included in the MERLIN examples subdirectory\n\n");

   ParameterList pl;

   pl.Add(new StringParameter('t', "Pedigree Template", templateFile));
   pl.Add(new StringParameter('g', "Genotype File", genotypeFile));
   pl.Add(new StringParameter('m', "Output Map File", mapFile));
   pl.Add(new StringParameter('d', "Output Data File", datFile));
   pl.Add(new StringParameter('p', "Output Pedigree File", pedFile));
   pl.Add(new SwitchParameter('c', "Use Coriell Ids", coriellIds));

   pl.Read(argc, argv);
   pl.Status();

   // Input buffers
   String line;
   StringArray tokens;

   // Read in the pedigree template
   StringArray  rowId;
   StringArray  rows;

   FILE * input = fopen(templateFile, "rt");

   if (input == NULL)
      error("Opening template file\n");

   while (!feof(input))
      {
      line.ReadLine(input);

      tokens.Clear();
      tokens.AddTokens(line);

      if (tokens.Length() == 0) continue;

      if (tokens.Length() == 6)
         {
         String tag = tokens[4].Left(7);

         tokens.Clear();
         tokens.Push("");
         tokens[0].printf("%s", (const char *) tag);
         tokens.AddTokens(line);
         }

      if (tokens.Length() != 7)
         error("Unexpected format for pedigree template\n\n"
               "Problem line transcribed below:\n%s\n",
               (const char *) line);

      rowId.Add(tokens[6]);
      rows.Add("");

      if (coriellIds)
         rows[rows.Length() - 1] = tokens[6] + " 1 0 0 1 ";
      else
         for (int i = 0; i < 5; i++)
            rows[rows.Length() - 1] += tokens[i] + "\t";
      }

   fclose(input);

   input = fopen(genotypeFile, "rt");

   if (input == NULL)
      error("Opening genotype file\n");

   tokens.Clear();

   while (tokens.Length() == 0)
      {
      line.ReadLine(input);
      tokens.AddTokens(line);
      }

   // Locate map file columns
   IntArray mapKey(3);

   mapKey[0] = tokens.Find("chrom");
   mapKey[1] = tokens.Find("rs#");
   mapKey[2] = tokens.Find("pos");

   if (mapKey.Min() < 0)
      error("Columns labelled 'chrom', 'rs#', 'pos'\n");

   // Link labels to allele names
   IntArray rowKey(rows.Length());

   for (int i = 0; i < rows.Length(); i++)
      {
      rowKey[i] = tokens.Find(rowId[i]);

      if (rowKey[i] == -1) continue;

      rows[i].LockBuffer(100000);
      rows[i].UnlockBuffer();
      }

   int cols = tokens.Length();

   FILE * map = fopen(mapFile, "wt");
   FILE * dat = fopen(datFile, "wt");
   FILE * ped = fopen(pedFile, "wt");

   if (map == NULL)
      error("Opening map file\n");

   if (dat == NULL)
      error("Opening data file\n");

   if (ped == NULL)
      error("Opening pedigree file\n");

   String previous;

   while (!feof(input))
      {
      line.ReadLine(input);

      tokens.Clear();
      tokens.AddTokens(line);

      if (tokens.Length() == 0) continue;

      if (tokens.Length() < cols)
         error("Marker %s -- Too few columns in genotype file\n",
               (const char *) tokens[0]);

      int ch = 0;
      while (!isdigit(tokens[mapKey[0]][ch])
             && tokens[mapKey[0]][ch]
             && tokens[mapKey[0]][ch] != 'x' && tokens[mapKey[0]][ch] != 'X'
             && tokens[mapKey[0]][ch] != 'y' && tokens[mapKey[0]][ch] != 'Y')
             ch++;

      fprintf(map, "%s\t%s\t%s\n",
              ((const char *) tokens[mapKey[0]]) + ch,
              (const char *) tokens[mapKey[1]],
              (const char *) tokens[mapKey[2]]);

      fprintf(dat, "%s %s\n",
         previous == tokens[mapKey[1]] ? "S2" : (const char *) "M",
         (const char *) tokens[mapKey[1]]);

      previous = tokens[mapKey[1]];

      for (int i = 0; i < rowKey.Length(); i++)
         if (rowKey[i] > -1)
            {
            switch (tokens[rowKey[i]][0])
               {
               case 'A': case 'a' : rows[i] += '1'; break;
               case 'C': case 'c' : rows[i] += '2'; break;
               case 'G': case 'g' : rows[i] += '3'; break;
               case 'T': case 't' : rows[i] += '4'; break;
               default : rows[i] += '.'; break;
               }

            rows[i] += "/";

            switch (tokens[rowKey[i]][1])
               {
               case 'A': case 'a' : rows[i] += '1'; break;
               case 'C': case 'c' : rows[i] += '2'; break;
               case 'G': case 'g' : rows[i] += '3'; break;
               case 'T': case 't' : rows[i] += '4'; break;
               default : rows[i] += '.'; break;
               }

            rows[i] += "\t";
            }
      }

   for (int i = 0; i < rowKey.Length(); i++)
      if (rowKey[i] > -1)
         fprintf(ped, "%s\n", (const char *) rows[i]);

   fclose(dat);
   fclose(map);
   fclose(ped);
   fclose(input);
   }
 
