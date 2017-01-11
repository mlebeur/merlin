////////////////////////////////////////////////////////////////////// 
// merlin/DiseaseModel.h 
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
 
#ifndef __DISEASEMODEL_H__
#define __DISEASEMODEL_H__

#include "MathMatrix.h"
#include "IntArray.h"

class Person;

class DisModel
   {
   public:
      String   label;

      int      affection;
      double   frequency;

      Matrix   penetrances;

      IntArray variableTypes;
      IntArray variables;
      IntArray conditions;
      Vector   values;

      bool ReadFromFile(IFILE & input);

      bool     ValidatePerson(Person & person);
      Vector & GetPenetrances(Person & person);

      void     Copy(const DisModel & source);
      DisModel & operator = (DisModel & source)  { Copy(source); return *this; }

      bool     CheckModel();

   private:
      int  TranslateCondition(String & token);
      bool TranslateVariable(String & token, int & variableType, int & variableId );
      bool TranslatePenetrances(String & token, Vector & penetrances);

      bool CheckCondition(double value1, int comparison, double value2);
      bool ValidateNumber(const char * string);

      void Clear();
   };

#endif

 
