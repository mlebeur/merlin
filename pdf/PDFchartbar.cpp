////////////////////////////////////////////////////////////////////// 
// pdf/PDFchartbar.cpp 
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
 
#include "PDFchartbar.h"
#include "PDF.h"
#include "Constant.h"

PDFChartBar::PDFChartBar()
  {
  count = 0;
  hasData = false;
  lowerBound = upperBound = _NAN_;
  }

PDFChartBar::~PDFChartBar()
  {
  }

void PDFChartBar::operator=(const PDFChartBar & rhs)
   {
   red = rhs.red;
   green = rhs.green;
   blue = rhs.blue;
   gray = rhs.gray;
   isInitialized = rhs.isInitialized;

   hasData = rhs.hasData;
   lowerBound = rhs.lowerBound;
   upperBound = rhs.upperBound;
   count = rhs.count;
   }

 
 
 
 
