////////////////////////////////////////////////////////////////////// 
// clusters/HaploGraph.cpp 
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
 
#include "HaploGraph.h"
#include "StringArray.h"
#include "MathVector.h"
#include "Error.h"

#include <ctype.h>

void HaplotypeGraph::LoadFromMemory(FounderGraph * founder, int haplos, int markers)
   {
   Dimension(haplos, markers);

   int haplo = 0;

   String index;
   lastGenotype.Set(-1);

   while (haplo < count)
      {
      for (int m = 0; m < markers; m++)
         if (founder->graph[m][haplo] >= 0)
            {
            haplotypes[haplo][m] = 0;

            index.printf("%dG%d", m, founder->graph[m][haplo]);

            AlleleGraph * graph = unknowns.GetGraph(index);

            graph->marker = m;
            graph->Append(haplo, founder->GetAllele(m, haplo), founder->GetAllele2(m, haplo));

            lastGenotype[haplo] = m;
            }
         else if (founder->alleles[m][haplo] == NOTZERO)
            {
            haplotypes[haplo][m] = 0;

            index.printf("%dM%d", m, haplo);

            MissingAllele * missing = unknowns.GetMissing(index);

            missing->haplotype = haplo;
            missing->marker = m;
            }
         else
            haplotypes[haplo][m] = founder->GetAllele(m, haplo),
            lastGenotype[haplo] = m;

      haplo++;
      }

   founder->likelihood = &likelihood;
   }

double HaplotypeGraph::Complexity(IntArray & alleleCounts)
   {
   Vector factors;

   factors.Dimension(count);
   factors.Set(1);

   double c = 1;

   for (int i = 0; i < unknowns.Length(); i++)
      {
      Unknown * u = (Unknown *) unknowns.Object(i);

      if (u->type == U_GRAPH)
         c *= 2;
      else
         {
         MissingAllele * m = (MissingAllele *) u;

         if (lastGenotype[m->haplotype] == -1)
            continue;

         factors[m->haplotype] *= alleleCounts[m->marker];
         }
      }

   return c * factors.Sum();
   }

void HaplotypeGraph::Dimension(int haplos, int markers)
   {
   count = haplos;

   lastGenotype.Dimension(count);
   lastGenotype.Set(-1);

   if (haplotypes != NULL)
      delete [] haplotypes;

   haplotypes = new IntArray[count];

   for (int i = 0; i < count; i++)
      haplotypes[i].Dimension(markers);
   }


 
