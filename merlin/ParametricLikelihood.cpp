////////////////////////////////////////////////////////////////////// 
// merlin/ParametricLikelihood.cpp 
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
 
#include "ParametricLikelihood.h"
#include "Error.h"

ParametricLikelihood::ParametricLikelihood()
   {
   model = NULL;
   }

ParametricLikelihood::~ParametricLikelihood()
   {
   }

int ParametricLikelihood::CountAlleles(Mantra & /* m */)
   {
   return 2;
   }

void ParametricLikelihood::FillPenetranceMatrix(Mantra & m)
   {
   // Check that a disease model is available
   if (model == NULL)
      error("Attempting to score parametric likelihood with disease model\n");

   // Two allele disease model is our default
   // int alleles = 2;
   int genotypes = 3;

   // Allocate matrix for storing probability of observed genotypes
   // conditional on the number of matching alleles in the underlying
   // true genotype.
   penetrances.Dimension(m.n, genotypes);

   // Allocate an array to keep track of which alleles are compatible
   // with which individuals
   compatibleAlleles.Dimension(m.n * 2);

   // Constant factor for likelihood
   priorLikelihood = 1.0;

   // Calculate individual specific penetrances for each individual
   for (int i = 0; i < m.n; i++)
      {
      Person & p = m.GetPerson(i);

      int disease_status = p.affections[model->affection];

      if (!model->ValidatePerson(p))
         {
         // No phenotype data for this invidual ...
         penetrances[i].Set(1.0);
         compatibleAlleles[i * 2] = compatibleAlleles[i * 2 + 1] = 1;
         continue;
         }

      Vector & mp = model->GetPenetrances(p);
      for (int j = 0; j < 3; j++)
         penetrances[i][j] = disease_status == 2 ? mp[j] : 1.0 - mp[j];

      compatibleAlleles[i * 2] = penetrances[i][1] > 0.0 || penetrances[i][0] > 0.0;
      compatibleAlleles[i * 2 + 1] = penetrances[i][1] > 0.0 || penetrances[i][2] > 0.0;
      }
   }

bool ParametricLikelihood::isPhenotyped(Mantra & /* m */, int i)
   {
   i = i / 2;

   return (penetrances[i][0] != 1.0) ||
          (penetrances[i][1] != 1.0) ||
          (penetrances[i][2] != 1.0);

   // Person & p = m.GetPerson(i / 2);
   //
   // return model->ValidatePerson(p);
   }

double ParametricLikelihood::GetAlleleFrequency(int /* founder */, int allele)
   {
   return allele == 0 ? 1.0 - model->frequency : model->frequency;
   }

void ParametricLikelihood::UpdateAlleles(int founder, int carrier)
   {
   carrier = carrier / 2;

   if (compatibleAlleles[carrier * 2] == false)
      HideAllele(founder, 0);

   if (compatibleAlleles[carrier * 2 + 1] == false)
      HideAllele(founder, 1);
   }

void ParametricLikelihood::UndoAlleles(int founder, int carrier)
   {
   carrier = carrier / 2;

   if (compatibleAlleles[carrier * 2] == false)
      ShowAllele(founder, 0);

   if (compatibleAlleles[carrier * 2 + 1] == false)
      ShowAllele(founder, 1);
   }

double ParametricLikelihood::JointProbability(Mantra & m)
   {
   // Retrieve allele frequency portion of likelihood for current state
   double lk = product.Peek();

   // Adjust likelihood according to the probability of observing each phenotype
   for (int i = 0; i < phenoList.Length(); i++)
      {
      int who = phenoList[i];

      int al1 = alleleState[m.state[who]];
      int al2 = alleleState[m.state[who + 1]];

      int genotype = al1 > al2 ? al1 * (al1 + 1) / 2 + al2:
                                 al2 * (al2 + 1) / 2 + al1;

      lk *= penetrances[who / 2][genotype];
      }

   return lk;
   }

bool ParametricLikelihood::ShortScoreVectors(Mantra & /* m */)
   {
   // There is no short cut for scoring likelihoods in this case
   return false;
   }

void ParametricLikelihood::SetupAlleleList(Mantra & m)
   {
   for (int j = 0; j < m.two_f; j++)
       ShowAllele(j, 0),
       ShowAllele(j, 1);
   }

void ParametricLikelihood::CleanUpAlleleList(Mantra & m)
   {
   for (int j = 0; j < m.two_f; j++)
       HideAllele(j, 0),
       HideAllele(j, 1);
   }

bool ParametricLikelihood::IdenticalPenetrances(int founder1, int founder2)
   {
   if (penetrances[founder1] != penetrances[founder2])
      return false;

   if (compatibleAlleles[founder1 * 2] != compatibleAlleles[founder2 * 2] ||
       compatibleAlleles[founder1 * 2 + 1] != compatibleAlleles[founder2 * 2 + 1])
       return false;

   return true;
   }

#define SWAPINT(a,b)    {int tmp = (a); (a) = (b); (b) = tmp; }

void ParametricLikelihood::SwapPenetrances(int founder1, int founder2)
   {
   penetrances.SwapRows(founder1, founder2);

   SWAPINT(compatibleAlleles[founder1 * 2], compatibleAlleles[founder2 * 2]);
   SWAPINT(compatibleAlleles[founder1 * 2 + 1], compatibleAlleles[founder2 * 2 + 1]);
   }


 
