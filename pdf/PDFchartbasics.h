////////////////////////////////////////////////////////////////////// 
// pdf/PDFchartbasics.h 
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
 
/* Written by Jan Wigginton */

#ifndef __PDFCHARTBASICS_H__
#define __PDFCHARTBASICS_H__

#include "IntArray.h"
#include "MathMatrix.h"
#include "PDF.h"
#include "PDFchartaxis.h"
#include "PDFchartlegend.h"
#include "PDFchartmarker.h"
#include "StringBasics.h"
#include "Constant.h"

#define PDFCHART_EDGE 0.1

class PDFChartBasics
   {
   public:

      PDFChartAxis   xAxis, yAxis;
      PDFChartLegend legend;

      String      title, subTitle;
      bool        useLegend, useColor, useHeader;
      bool        drawHGrid, drawXY, drawQuadrants;
      bool        drawVConnector;
      bool        symmetricAxes;

      PDFFonts    titleFont;
      double      titleFontSize;

      //Alias for legend.lines pointer
      PDFChartLine * & lines;

      PDFChartBasics();
      virtual ~PDFChartBasics();

      void SetSeriesColor(int series, double red, double green, double blue);
      void SetSeriesGray(int series, double gray);
      void SetSeriesLabel(int series, const char * label);

      void Draw(PDF & pdf);
      void DrawInBox(PDF & pdf, double x_0, double y_0, double x_1, double y_1);
      void DrawInGrid(PDF & pdf, int row_to_use, int col_to_use, int rows, int cols, double spacer = 0.05);
      void DrawInUpperLeft(PDF & pdf, double bump = 0.0);
      void DrawInUpperRight(PDF & pdf, double bump = 0.0);
      void DrawInLowerLeft(PDF & pdf, double bump = 0.0);
      void DrawInLowerRight(PDF & pdf, double bump = 0.0);

   protected:

      // Variables related to the chart's page size/type
      double      height, width;
      double      titleHeight;
      double      x0, y0, x1, y1;
      double      x0Page, y0Page;
      double      xScale, yScale;
      bool        noData;
      bool        keepLegend;

      // Drawing subroutines...
      void DrawAxis(PDF & pdf,  PDFChartAxis & axis);
      void DrawTicks(PDF & pdf,  PDFChartAxis & axis);
      void DrawTickLabels(PDF & pdf, PDFChartAxis & axis);
      void DrawAxisLabels(PDF & pdf, PDFChartAxis & axis);
      void DrawTitle(PDF & pdf);
      void DrawOutline(PDF & pdf);
      virtual void DrawBody(PDF & pdf)= 0;
      virtual void DrawLegend(PDF & pdf) = 0;

      void InitializePage(PDF & pdf);
      void InitializeBlankGraph(PDF & pdf);

      virtual void CloseChart(PDF & pdf);
      virtual bool OpenChart(PDF & pdf) = 0 ;
      void Reset();
      void Init();

      // Scaling, sizing subroutines
      int  ReadData(const char *filename, Matrix &temp_values);
      void SetChartDimensions(PDF & pdf);
      virtual void SetOtherDimensions(PDF & pdf);

      void SetTickSize(PDFChartAxis & axis);
      void SetAxisDimensions(PDF & pdf, PDFChartAxis & axis);
      void SetTitleDimensions(PDF & pdf);
      void SetLegendDimensions(PDF & pdf);
      void SetAxisScaling();

      // Utility functions
      double MapX(double pos);
      double MapY(double pos);
      double Space();

   };

#endif




 
