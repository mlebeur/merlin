////////////////////////////////////////////////////////////////////// 
// merlin/MerlinPDF.cpp 
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
 
#include "MerlinPDF.h"

MerlinPDF::MerlinPDF(Vector & analysisPositions) : positions(analysisPositions)
   {
   }

MerlinPDF::~MerlinPDF()
   {
   if (isOpen()) CloseFile();
   }

void MerlinPDF::OpenFile(const char * pdfFilename)
   {
   PDF::OpenFile(pdfFilename);

   info.author = "Merlin (c) 2000-2003 Goncalo Abecasis";
   info.subject = "Automatically generated summaries";

   page.SetSize(psLetterR);

   filename = pdfFilename;
   }

void MerlinPDF::CloseFile()
   {
   printf("Graphical output in file [%s]\n", (const char *) filename);
   PDF::CloseFile();
   }

void MerlinPDF::PrepareChart(int series)
   {
   chart.Reset();

   chart.xAxis.label = chromosome + "Position (cM)";

   double nudge = positions.First() == positions.Last() ? 1.0 : 0.0;
   chart.xAxis.SetMin(positions.First() * 100.0 - nudge);
   chart.xAxis.SetMax(positions.Last() * 100.0 + nudge);
   
   chart.useLegend = false;

   chart.Dimension(series, positions.Length());
   x.Dimension(positions.Length());
   y.Dimension(series, positions.Length());

   for (int i = 0; i < series; i++)
      chart.ShowMarker(i, positions.Length() == 1);

   for (int i = 0; i < positions.Length(); i++)
      x[i] = positions[i] * 100.0;
   }

void MerlinPDF::SelectChromosome(int chr)
   {
   if (chr == 0)
      chromosome = "";
   else
      {
      chromosome = "Chromosome ";
      chromosome += chr;
      chromosome += " ";
      }
   }

void MerlinPDF::SelectChromosomeX()
   {
   chromosome = "Chromosome X ";
   }

void MerlinPDF::DrawChart()
   {
   chart.SetDataValues(x, y);
   chart.Draw(*this);
   }
 
