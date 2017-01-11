////////////////////////////////////////////////////////////////////// 
// merlin/DiseaseModel.cpp 
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
 
#include "DiseaseModel.h"
#include "Pedigree.h"

#define _COMPARISON_INVALID_    -1
#define _COMPARISON_EQUAL_       0
#define _COMPARISON_GREATER_     1
#define _COMPARISON_LESSER_      2
#define _COMPARISON_GE_          3
#define _COMPARISON_LE_          4
#define _COMPARISON_NOT_EQUAL_   5

#define _CONDITION_ON_SEX_           0
#define _CONDITION_ON_COVARIATE_     1
#define _CONDITION_ON_TRAIT_         2
#define _CONDITION_NONE_             3

bool DisModel::ReadFromFile(IFILE & f)
   {
   // First clear away any previous settings
   Clear();

   String      input, buffer;
   StringArray tokens, extra;

   // If we get past the end-of-file, then stop
   if (ifeof(f))
      return false;

   input = buffer.ReadLine(f);
   input.Trim();

   // Skip comments
   if (input[0] == '#')
      return false;

   // Divide each line into tokens
   tokens.Clear();
   tokens.AddTokens(input);

   // Skip blank lines
   if (tokens.Length() == 0)
      return false;

   // Print message for tracing...
   printf("      Input: %s\n", (const char *) buffer);

   // Need a minimum of four tokens per line
   if (tokens.Length() < 4)
      {
      printf("      Error: Trait name, allele frequency, penetrances, and label for analysis required.\n");
      return false;
      }

   // First lookup trait label
   affection = Pedigree::LookupAffection(tokens[0]);
   if (affection < 0)
      {
      printf("      Error: Token '%s' does not correspond to affection status in the pedigree\n",
             (const char *) tokens[0]);
      return false;
      }

   // First check that allele frequency is valid
   bool fail = !ValidateNumber(tokens[1]);

   frequency = tokens[1];

   if (fail || frequency <= 0.0 || frequency >= 1.0)
      {
      printf("      Error: Token '%s' is not a valid allele frequency\n",
             (const char *) tokens[1]);
      return false;
      }

   // If penetrances are missing, we assume that liability class information
   // will follow, with one liability class per line
   if (tokens[2][0] == '*')
      {
      while (true)
         {
         if (ifeof(f))
            {
            printf("      Error: Unexpected reached end-of-file, incomplete liability class matrix\n");
            return false;
            }

         // Read additional lines
         input = buffer.ReadLine(f);
         input.Trim();

         // skip blank lines and comments
         if (input.Length() == 0 || input[0] == '#')
            continue;

         extra.ReplaceTokens(input);

         // In principle, non-blank lines should not give zero tokens, but
         // this serves as safety net
         if (extra.Length() == 0)
            continue;

         // Produce log ...
         printf("      Input: %s\n", (const char *) buffer);

         int variableType, variableId;

         if (!TranslateVariable(extra[0], variableType, variableId))
            {
            printf("      Error: Can't condition on '%s'\n", (const char *) extra[0]);
            return false;
            }

         variableTypes.Push(variableType);
         variables.Push(variableId);

         // If there is no condition, we are done
         if (variableType == _CONDITION_NONE_)
            {
            // Check that line format is correct
            if (extra.Length() != 2)
               {
               printf("      Error: '%s' should be followed by three comma separated penetrances\n", (const char *) extra[0]);
               return false;
               }

            penetrances.Dimension(penetrances.rows + 1, 3);

            // Read final set of penetrances
            if (!TranslatePenetrances(extra[1], penetrances.Last()))
               return false;

            break;
            }

         if (extra.Length() != 4)
            {
            printf("      Error: Line should include conditioning variable, operator, value and penetrances\n"
                   "             For example, AGE > 47 0.01,0.01,1.00\n\n");
            return false;
            }

         int condition = TranslateCondition(extra[1]);

         if (condition < 0)
            return false;

         conditions.Push(condition);

         if (variableType == _CONDITION_ON_SEX_ &&
            (extra[2].SlowCompare("MALE") == 0 || extra[2].SlowCompare("M") == 0))
            extra[2] = "1";

         if (variableType == _CONDITION_ON_SEX_ &&
            (extra[2].SlowCompare("FEMALE") == 0 || extra[2].SlowCompare("F") == 0))
            extra[2] = "2";

         if (!ValidateNumber(extra[2]))
            {
            printf("      Error: '%s' is not a valid number\n", (const char *) extra[2]);
            return false;
            }

         values.Push(extra[2]);

         penetrances.Dimension(penetrances.rows + 1, 3);

         // Read final set of penetrances
         if (!TranslatePenetrances(extra[3], penetrances.Last()))
            return false;
         }
      }
   else
      {
      penetrances.Dimension(1, 3);

      if (!TranslatePenetrances(tokens[2], penetrances[0]))
         return false;
      }

   // Finally load the text description for the model
   label = tokens[3];
   for (int j = 4; j < tokens.Length(); j++)
      {
      label += " ";
      label += tokens[j];
      }
   printf("  Validated: Loaded model '%s'\n", (const char *) label);

   return true;
   }

int DisModel::TranslateCondition(String & token)
   {
   if (token.Length() == 1)
      {
      if (token[0] == '=') return _COMPARISON_EQUAL_;
      if (token[0] == '<') return _COMPARISON_LESSER_;
      if (token[0] == '>') return _COMPARISON_GREATER_;
      }
   else
      {
      if (token == "==" || token == "is" ) return _COMPARISON_EQUAL_;
      if (token == "<=") return _COMPARISON_LE_;
      if (token == ">=") return _COMPARISON_GE_;
      if (token == "<>" || token == "!=") return _COMPARISON_NOT_EQUAL_;
      }

   printf("      Error: '%s' is not a recognize comparison operator\n", (const char *) token);
   return _COMPARISON_INVALID_;
   }

bool DisModel::TranslateVariable(String & token, int & variableType, int & variableId )
   {
   variableId = Pedigree::LookupCovariate(token);
   if (variableId >= 0)
      {
      variableType = _CONDITION_ON_COVARIATE_;
      return true;
      }

   variableId = Pedigree::LookupTrait(token);
   if (variableId >= 0)
      {
      variableType = _CONDITION_ON_TRAIT_;
      return true;
      }

   if (token.SlowCompare("sex") == 0)
      {
      variableType = _CONDITION_ON_SEX_;
      variableId = 0;
      return true;
      }

   if (token.SlowCompare("Otherwise") == 0 || token.SlowCompare("Else") == 0)
      {
      variableType = _CONDITION_NONE_;
      variableId = 0;
      return true;
      }

   return false;
   }

bool DisModel::TranslatePenetrances(String & token, Vector & pen)
   {
   StringArray ptokens;
   bool fail = false;

   // The third token should be a set of penetrances separated by commas
   ptokens.AddColumns(token, ',');

   // Next check that penetrances are valid
   for (int j = 0; j < ptokens.Length(); j++)
      fail |= !ValidateNumber(ptokens[j]);

   if (fail || ptokens.Length() != 3)
      {
      printf("      Error: '%s' is not a comma separated list of three penetrances\n",
             (const char *) token);
      return false;
      }

   for (int j = 0; j < 3; j++)
      pen[j] = ptokens[j];

   if (pen.Min() < 0.0 || pen.Max() > 1.0)
      {
      printf("      Error: Penetrances should be between 0.0 and 1.0\n",
             (const char *) token);
      return false;
      }

   return true;
   }

void DisModel::Clear()
   {
   penetrances.Dimension(0,0);

   variableTypes.Clear();
   variables.Clear();
   conditions.Clear();

   values.Clear();
   }

bool DisModel::ValidatePerson(Person & person)
   {
   if (person.affections[affection] == 0)
      return false;

   for (int i = 0; i < variableTypes.Length(); i++)
      switch (variableTypes[i])
         {
         case _CONDITION_ON_COVARIATE_ :
            if (person.covariates[variables[i]] == _NAN_)
               return false;
            break;
         case _CONDITION_ON_TRAIT_ :
            if (person.traits[variables[i]] == _NAN_)
               return false;
            break;
         case _CONDITION_ON_SEX_ :
            if (person.sex == 0)
               return false;
            break;
         }
   return true;
   }

Vector & DisModel::GetPenetrances(Person & person)
   {
   for (int i = 0; i < variableTypes.Length(); i++)
      switch (variableTypes[i])
         {
         case _CONDITION_ON_COVARIATE_ :
            if (CheckCondition(person.covariates[variables[i]], conditions[i], values[i]))
               return penetrances[i];
            break;
         case _CONDITION_ON_TRAIT_ :
            if (CheckCondition(person.traits[variables[i]], conditions[i], values[i]))
               return penetrances[i];
            break;
         case _CONDITION_ON_SEX_ :
            if (CheckCondition(person.sex, conditions[i], values[i]))
               return penetrances[i];
            break;
         }
   return penetrances.Last();
   }

bool DisModel::CheckCondition(double value1, int comparison, double value2)
   {
   switch (comparison)
      {
      case _COMPARISON_GREATER_ :
         return value1 > value2;
      case _COMPARISON_EQUAL_ :
         return value1 == value2;
      case _COMPARISON_LESSER_ :
         return value1 < value2;
      case _COMPARISON_LE_ :
         return value1 <= value2;
      case _COMPARISON_GE_ :
         return value1 >= value2;
      case _COMPARISON_NOT_EQUAL_ :
         return value1 != value2;
      }
   return false;
   }

void DisModel::Copy(const DisModel & source)
   {
   label = source.label;

   affection = source.affection;
   frequency = source.frequency;

   penetrances = source.penetrances;

   variableTypes = source.variableTypes;
   variables = source.variables;
   conditions = source.conditions;
   values = source.values;
   }

bool DisModel::ValidateNumber(const char * string)
   {
   char * ptr = NULL;

   strtod(string, &ptr);

   return ptr[0] == 0;
   }

bool DisModel::CheckModel()
   {
   if (frequency <= 0.0 || frequency >= 1.0)
      return false;

   if (penetrances.Min() < 0.0 || penetrances.Max() > 1.0)
      return false;

   return true;
   }
 
