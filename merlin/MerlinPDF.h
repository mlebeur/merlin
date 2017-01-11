////////////////////////////////////////////////////////////////////// 
// merlin/MerlinPDF.h 
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
 
#ifndef __MERLINPDF_H__
#define __MERLINPDF_H__

#include "PDF.h"
#include "PDFlinechart.h"

class MerlinPDF : public PDF
   {
   public:
      PDFLineChart chart;
      Vector & positions;

      MerlinPDF(Vector & plotLocations);
      ~MerlinPDF();

      void OpenFile(const char * filename);
      void CloseFile();

      void PrepareChart(int series = 1);
      void DrawChart();

      void SelectChromosome(int chr);
      void SelectChromosomeX();

      Vector x;
      Matrix y;

   private:
      String filename;
      String chromosome;

   };

#endif

 
