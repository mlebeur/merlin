////////////////////////////////////////////////////////////////////// 
// clusters/Unknown.h 
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
 
#ifndef __UNKNOWN_H__
#define __UNKNOWN_H__

#include "IntArray.h"
#include "StringMap.h"

#define U_GRAPH      0
#define U_MISSING    1
#define U_INVALID    2

class Unknown
   {
   public:
      int type;
      int marker;

      Unknown()
         {
         marker = -1;
         type = U_INVALID;
         }

      virtual ~Unknown() {}
   };

class MissingAllele : public Unknown
   {
   public:
      int haplotype;

      MissingAllele()
         { type = U_MISSING; };

      virtual ~MissingAllele() {}
   };

class AlleleGraph : public Unknown
   {
   public:
      IntArray haplotypes;
      IntArray alleles;

      AlleleGraph()
         { type = U_GRAPH; }

      virtual ~AlleleGraph() {}

      void Append(int haplo, int allele1, int allele2)
         {
         haplotypes.Push(haplo);
         alleles.Push(allele1);
         alleles.Push(allele2);
         }
   };

class SetOfUnknowns : public StringMap
   {
   public:
      Unknown * GetUnknown(int index)
         { return (Unknown *) Object(index); };

      MissingAllele * GetMissing(const ::String & name);
      AlleleGraph * GetGraph(const ::String & name);

      ~SetOfUnknowns();
      
   private:
      static void * create_missing();
      static void * create_graph();
   };

#endif

 
