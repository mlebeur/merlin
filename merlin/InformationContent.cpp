////////////////////////////////////////////////////////////////////// 
// merlin/InformationContent.cpp 
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
 
#include "InformationContent.h"
#include "MerlinCore.h"
#include "MerlinPDF.h"

InformationContent::InformationContent()
   {
   infoFile = NULL;
   infoTable = NULL;
   }

InformationContent::~InformationContent()
   {
   }

void InformationContent::Setup(AnalysisInfo & info)
   {
   // Allocate memory for storing information content
   infoBits.Dimension(info.positions);
   infoBits.Zero();

   // And total bit count at each analysed position
   totalBits.Dimension(info.positions);
   totalBits.Zero();
   }

void InformationContent::Analyse(AnalysisInfo & info, Tree & tree, int pos)
   {
   // Calculate information
   stats.GetTreeInfo(tree);

   infoBits[pos] += stats.information * info.m->bit_count;
   totalBits[pos] += info.m->bit_count;

   if (infoFile != NULL)
      fprintf(infoFile, "%10s %10s %12.5f %5d\n",
         (const char *) (*info.famid),
         (const char *) (*info.labels)[pos],
         stats.information, info.m->bit_count);
   }

void InformationContent::AnalyseUninformative(AnalysisInfo & info, Tree & )
   {
   totalBits += info.m->bit_count;

   if (infoFile != NULL)
      for (int pos = 0; pos < info.positions; pos++)
         fprintf(infoFile, "%10s %10s %12.5f %5d\n",
            (const char *) (*info.famid),
            (const char *) (*info.labels)[pos],
            0.0, info.m->bit_count);
   }


void InformationContent::Report(AnalysisInfo & info)
   {
   printf("%15s %6s\n", "Position", "Info");

   for (int i = 0; i < infoBits.Length(); i++)
      printf("%15s %6.4f\n", (const char *) (*info.labels)[i],
             totalBits[i] == 0.0 ? 0.0 :
             infoBits[i] / totalBits[i]);

   if (infoTable != NULL)
      for (int i = 0; i < infoBits.Length(); i++)
         fprintf(infoTable, "%d\t%.4f\t%s\t%.5f\n",
                info.chromosome, (*info.analysisPositions)[i] * 100.,
                (const char *) (*info.labels)[i],
                totalBits[i] == 0.0 ? 0.0 :
                infoBits[i] / totalBits[i]);

   if (info.drawPDF)
      {
      info.pdf->PrepareChart();

      info.pdf->chart.title = "Information Content";
      info.pdf->chart.yAxis.SetMin(0.0);
      info.pdf->chart.yAxis.SetMax(1.0);
      info.pdf->chart.yAxis.label = "Information Content";

      for (int i = 0; i < infoBits.Length(); i++)
         info.pdf->y[0][i] = totalBits[i] == 0.0 ? 0.0 :
                             infoBits[i] / totalBits[i];

      info.pdf->DrawChart();
      }

   printf("\n");
   }

void InformationContent::OpenFiles(AnalysisInfo & info, String & prefix)
   {
   if (MerlinCore::tabulate)
      {
      tablename = prefix + "-info.tbl";

      infoTable = fopen(tablename, "wt");

      if (infoTable == NULL)
         error("Can't open information table [%s]", (const char *) tablename);

      fprintf(infoTable, "CHR\tPOS\tLABEL\tINFO\n");
      }

   if (info.perFamily == false)
      return;

   filename = prefix + ".inf";

   infoFile = fopen((const char *) filename, "wt");

   if (infoFile == NULL)
      error("Can't open information file [%s]\n", (const char *) filename);

   fprintf(infoFile, "%10s %10s %12s %5s\n", "FAMILY", "POSITION", "INFORMATION", "BITS");
   }

void InformationContent::CloseFiles()
   {
   if (infoTable != NULL)
      {
      printf("Information content tabulated in file [%s]\n", (const char *) tablename);

      fclose(infoTable);
      infoTable = NULL;
      }

   if (infoFile == NULL)
      return;

   printf("Information scores for individual families in file [%s]\n",
          (const char *) filename);

   fclose(infoFile);
   infoFile = NULL;
   }
 
