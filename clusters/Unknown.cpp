////////////////////////////////////////////////////////////////////// 
// clusters/Unknown.cpp 
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
 
#include "Unknown.h"

MissingAllele * SetOfUnknowns::GetMissing(const ::String & name)
   {
   int index = Find(name, create_missing);
   return (MissingAllele *) objects[index];
   }

AlleleGraph * SetOfUnknowns::GetGraph(const ::String & name)
   {
   int index = Find(name, create_graph);
   return (AlleleGraph *) objects[index];
   }

void * SetOfUnknowns::create_missing()
   {
   return new MissingAllele;
   }

void * SetOfUnknowns::create_graph()
   {
   return new AlleleGraph;
   }

SetOfUnknowns::~SetOfUnknowns()
   {
   for (int i = 0; i < Length(); i++)
      delete (Unknown *) objects[i];
   }



 
