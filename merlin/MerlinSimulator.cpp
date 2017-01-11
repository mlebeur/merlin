////////////////////////////////////////////////////////////////////// 
// merlin/MerlinSimulator.cpp 
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
 
#include "MerlinSimulator.h"
#include "MerlinCluster.h"
#include "MapFunction.h"
#include "ParametricLikelihood.h"
#include "Random.h"
#include "Error.h"

#include <math.h>

String Simulator::model;

String   Simulator::reuseLabel;
Pedigree Simulator::reuseSource;
int      Simulator::reuseCounter = -1;

// Qtl parameters
int    Simulator::qtlMarker = -1;
int    Simulator::qtlTrait  = -1;
double Simulator::sQtl = 0.0;
double Simulator::sPolygenes = 0.0;
double Simulator::sEnvironment = 0.0;

void Simulator::Simulate(Pedigree & ped, IntArray & markers)
   {
   int badModel = false;

   LoadReuseableHaplotypes();

   if (!model.IsEmpty())
      {
      badModel = true;

      StringArray  tokens;
      DisModel     parameters;

      tokens.AddColumns(model, ',');

      if (tokens.Length() >= 6 && tokens.Length() <= 7)
         {
         parameters.affection = ped.LookupAffection(tokens[0]);
         parameters.frequency = tokens[1].AsDouble();

         parameters.penetrances.Dimension(1,3);

         for (int i = 0; i < 3; i++)
            parameters.penetrances[0][i] = tokens[i + 2].AsDouble();

         double femalePosition = tokens[5].AsDouble() * 0.01;
         double malePosition = tokens.Last().AsDouble() * 0.01;

         if (parameters.CheckModel())
            {
            Alternative(ped, markers, parameters, femalePosition, malePosition);
            return;
            }
         }

      if (tokens.Length() == 5)
         {
         qtlTrait = ped.LookupTrait(tokens[0]);
         qtlMarker = ped.LookupMarker(tokens[1]);

         double varQtl = tokens[2].AsDouble();
         double varPoly = tokens[3].AsDouble();
         double varEnv = tokens[4].AsDouble();

         double qtlFreq = qtlMarker >= 0 ? ped.GetMarkerInfo(qtlMarker)->freq[1] : 0.0;

         if (qtlTrait < 0 || qtlMarker < 0 || qtlFreq <= 0.0 ||
             varQtl < 0.0 || varPoly < 0.0 || varEnv < 0.0)
            qtlMarker = -1;
         else
            {
            sQtl = sqrt(varQtl / (2.0 * qtlFreq * (1.0  - qtlFreq)));
            sPolygenes = sqrt(varPoly);
            sEnvironment = sqrt(varEnv);

            badModel = false;
            }
         }
      }

   if (badModel)
      {
      printf("SIMULATOR: Model specified incorrectly, simulating data under the null\n\n"
             "  Correct format for simulating conditional on disease status is:\n"
             "    --trait:AFFECTION,FREQ(-),PEN(+/+),PEN(+/-),PEN(-/-),POSITION\n\n"
             "  Correct format for simulating a QTL is:\n"
             "    --trait:TRAIT,SNP,Var(QTL),Var(Polygenes),Var(Environment)\n\n");

      model.Clear();
      }

   NullData(ped, markers);
   }

void Simulator::Alternative
   (Pedigree & ped, IntArray & markers, DisModel & model,
    double femalePosition, double malePosition)
   {
   LoadReuseableHaplotypes();

   Mantra mantra;
   ParametricLikelihood tree;

   // These arrays track inheritance vectors
   IntArray current, best;

   // Separate markers to the left and right of disease locus position
   IntArray left, right;
   for (int i = 0; i < markers.Length(); i++)
      if (ped.GetMarkerInfo(markers[i])->positionFemale < femalePosition ||
          ped.GetMarkerInfo(markers[i])->positionMale < malePosition)
         left.Push(markers[i]);
      else break;

   for (int i = left.Length(); i < markers.Length(); i++)
      right.Push(markers[i]);
   left.Reverse();

   // Check that position is consistent with male and female recombination
   // maps ...
   MarkerInfo * leftMarker = left.Length() ? ped.GetMarkerInfo(left[0]) : NULL;
   MarkerInfo * rightMarker = right.Length() ? ped.GetMarkerInfo(right[0]) : NULL;

   if (right.Length() && (rightMarker->positionFemale < femalePosition ||
                          rightMarker->positionMale < malePosition) ||
       left.Length() && (leftMarker->positionFemale > femalePosition ||
                         leftMarker->positionMale > malePosition))
          error("MerlinSimulator -- Disease locus position not consistent with\n"
                "                   male and female genetic maps.\n");

   for (int f = 0; f < ped.familyCount; f++)
      {
      mantra.Prepare(ped, *ped.families[f]);

      tree.model = &model;
      tree.FuzzyScoreVectors(mantra);

      current.Dimension(tree.bit_count);
      current.Zero();

      double likelihood = 0.0;
      SampleInheritanceVector(tree, 0, tree.bit_count - 1, current, best, likelihood);

      // Then extract an inheritance vector
      for (int i = 0; i < mantra.bit_count; i++)
         {
         int bit = mantra.bits[i];

         mantra.vector[bit] = best[i] == 0 ? mantra.vector[bit] & ~1 : mantra.vector[bit] | 1;
         }

      best.Dimension(mantra.two_n - mantra.two_f);
      for (int i = mantra.two_f; i < mantra.two_n; i++)
         best[i - mantra.two_f] = mantra.vector[i] & 1;

      int savedCounter = reuseCounter;

      if (left.Length())
         {
         current = best;

         double thetaF = DistanceToRecombination(femalePosition - leftMarker->positionFemale);
         double thetaM = DistanceToRecombination(malePosition - leftMarker->positionMale);

         for (int j = 0; j < best.Length(); j += 2)
            {
            current[j] ^= (globalRandom.Next() < thetaF ? 1 : 0);
            current[j + 1] ^= (globalRandom.Next() < thetaM ? 1 : 0);
            }

         SimulateFamily(ped, *ped.families[f], left, current);
         }

      if (right.Length())
         {
         reuseCounter = savedCounter;

         double thetaF = DistanceToRecombination(rightMarker->positionFemale - femalePosition);
         double thetaM = DistanceToRecombination(rightMarker->positionMale - malePosition);

         for (int j = 0; j < best.Length(); j += 2)
            {
            best[j] ^= (globalRandom.Next() < thetaF ? 1 : 0);
            best[j + 1] ^= (globalRandom.Next() < thetaM ? 1 : 0);
            }

         SimulateFamily(ped, *ped.families[f], right, best);
         }
      }
   }

void Simulator::NullData(Pedigree & ped, IntArray & markers)
   {
   IntArray vector;

   for (int f = 0; f < ped.familyCount; f++)
      {
      Family * family = ped.families[f];

      vector.Dimension(family->nonFounders * 2);

      // Create a random inheritance vector for the first marker
      for (int i = 0; i < vector.Length(); i++)
         vector[i] = globalRandom.Binary();

      SimulateFamily(ped, *family, markers, vector);
      }
   }

void Simulator::SimulateFamily(Pedigree & ped, Family & family,
                               IntArray & markers, IntArray & vector)
   {
   IntArray alleles, mzTwins;

   alleles.Dimension(family.count * 2);

   double femalePosition = -1, malePosition = -1;

   for (int i = 0; i < markers.Length(); i++)
      {
      // Get information on the current marker (position, allele frequencies)
      MarkerInfo * info = ped.GetMarkerInfo(markers[i]);
      bool isClustered = clusters.Enabled() && clusters.markerToCluster[markers[i]] != NULL;
      MarkerCluster * cluster = isClustered ? clusters.markerToCluster[markers[i]] : NULL;

      // Drop founder alleles, by sampling each chromosome
      for (int j = 0; j < family.founders; j++)
         for (int chr = 0; chr < 2; chr++)
            if (reuseCounter < 0)
               {
               int  allele = 0;
               double rand = globalRandom.Next(), sum = 0.0;

               // By default, we sample one allele according to marker
               // allele frequencies
               Vector * freq = &info->freq;

               // For clustered markers, we sample a haplotype which
               // will determine alleles at multiple consecutive
               // markers
               if (isClustered)
                  {
                  freq = &cluster->freqs;
                  allele = -1;
                  }

               while (sum < rand && (allele + 1 < freq->Length()))
                  sum += (*freq)[++allele];

               alleles[j * 2 + chr] = allele;
               }
            else
               alleles[j * 2 + chr] =
                  reuseSource[reuseCounter + j].markers[markers[i]][chr];

#ifdef __CHROMOSOME_X__
      // All males are homozygous for X
      for (int j = 0; j < family.founders; j++)
         if (ped[family.path[j]].sex == SEX_MALE)
            alleles[j * 2 + 1] = alleles[j * 2];
#endif

      // Update the inheritance vector using the recombination fraction
      if (femalePosition != -1)
         {
         double thetaF = DistanceToRecombination(fabs(info->positionFemale - femalePosition));
         double thetaM = DistanceToRecombination(fabs(info->positionMale - malePosition));

         for (int j = 0; j < vector.Length(); j += 2)
            {
            vector[j] ^= (globalRandom.Next() < thetaF ? 1 : 0);
            vector[j + 1] ^= (globalRandom.Next() < thetaM ? 1 : 0);
            }
         }
      femalePosition = info->positionFemale;
      malePosition = info->positionMale;

      // Assign offspring genotypes using founder alleles + vector
      for (int j = 0; j < family.nonFounders; j++)
         {
         int index  = j + family.founders;
         int person = family.path[index];
         Person & p = ped[person];

         // Check if we are dealing with a twin...
         if (p.zygosity & 1)
            {
            bool newTwin = true;

            // If so, check if individual is preceded by cotwin
            for (int k = 0; k < mzTwins.Length(); k++)
               // and if we find a co-twin...
               if (p.isTwin(ped[mzTwins[k]]))
                  {
                  // Retrieve the cotwin
                  Person & cotwin = ped[mzTwins[k]];

                  // Make sure the twins have identical inheritance vectors
                  vector[j * 2] =
                     vector[(cotwin.traverse - family.founders) * 2];
                  vector[j * 2 + 1] =
                     vector[(cotwin.traverse - family.founders) * 2 + 1];

                  // This is not a new twinship!
                  newTwin = false;

                  break;
                  }

            // We finish by adding a new twin to the list...
            if (newTwin) mzTwins.Push(p.serial);
            }

         int mother_allele = p.mother->traverse*2 + vector[j*2];
         int father_allele = p.father->traverse*2 + vector[j*2+1];

         index *= 2;
         alleles[index] = alleles[mother_allele];
#ifdef __CHROMOSOME_X__
         // All males are homozygous for X
         if (p.sex == SEX_MALE)
            alleles[index + 1] = alleles[index];
         else
#endif
         alleles[index + 1] = alleles[father_allele];
         }

      // Replicate the pattern of missing data in the pedigree
      if (!isClustered)
         {
         for (int j = 0; j < family.count; j++)
            {
            Person & person    = ped[family.path[j]];
            Alleles & genotype = person.markers[markers[i]];

            if (genotype.isKnown())
               {
               genotype[0] = alleles[j * 2];
               genotype[1] = alleles[j * 2 + 1];
               }
            else
               genotype[0] = genotype[1] = 0;
            }

         if (qtlMarker == markers[i])
            SimulateQTL(ped, family, alleles);
         }
      else
         {
         for (int m = 0; m < cluster->markerIds.Length(); m++)
            {
            int markerid = cluster->markerIds[m];

            for (int j = 0; j < family.count; j++)
               {
               Person & person    = ped[family.path[j]];
               Alleles & genotype = person.markers[markerid];

               if (genotype.isKnown())
                  {
                  genotype[0] = cluster->RetrieveAllele(alleles[j * 2], m);
                  genotype[1] = cluster->RetrieveAllele(alleles[j * 2 + 1], m);
                  }
               else
                  genotype[0] = genotype[1] = 0;
               }
            }

         if (qtlMarker >= 0 && clusters.markerToCluster[qtlMarker] == cluster)
            {
            int m = cluster->markerIds.Find(qtlMarker);

            for (int j = 0; j < family.count; j++)
               {
               alleles[j * 2] = cluster->RetrieveAllele(alleles[j * 2], m);
               alleles[j * 2 + 1] = cluster->RetrieveAllele(alleles[j * 2 + 1], m);
               }

            SimulateQTL(ped, family, alleles);
            }

         i += cluster->markerIds.Length() - 1;
         }
      }

   if (reuseCounter >= 0)
      reuseCounter += family.founders;
   }

void Simulator::SampleInheritanceVector(Tree & tree, int node, int bit,
                                        IntArray & current, IntArray & best,
                                        double & sum, double weight)
   {
   switch (tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         {
         weight *= tree.nodes[node].value;
         sum += weight;
         if (weight > 0 && (globalRandom < weight / sum))
            {
            while (bit >= 0)
               {
               current[bit] = globalRandom.Binary();
               bit--;
               }
            best = current;
            }
         break;
         }
      case TREE_NODE_ONE :
         current[bit] = globalRandom.Binary();
         SampleInheritanceVector(tree, tree.nodes[node].child[0], bit - 1,
                                 current, best, sum, weight);
         break;
      case TREE_NODE_TWO :
         current[bit] = 0;
         SampleInheritanceVector(tree, tree.nodes[node].child[0], bit - 1,
                                 current, best, sum, weight * 0.5);
         current[bit] = 1;
         SampleInheritanceVector(tree, tree.nodes[node].child[1], bit - 1,
                                 current, best, sum, weight * 0.5);
         break;
      }
   }

void Simulator::LoadReuseableHaplotypes()
   {
   if (reuseCounter >= 0)
      reuseCounter = 0;
   else if (!reuseLabel.IsEmpty())
      {
      reuseSource.Prepare(reuseLabel + ".dat");
      reuseSource.Load(reuseLabel + ".ped");

      reuseCounter = 0;

      if (clusters.Enabled())
         error("Cannot gene drop user specified haplotypes when clusters are enabled\n");
      }
   }

void Simulator::SimulateQTL(Pedigree & ped, Family & family, IntArray & alleles)
   {
   Vector  poly(family.count);

   // Loop through all individuals in the family
   for (int i = 0; i < family.count; i++)
      {
      Person & person = ped[family.path[i]];

      // First sample polygenic effects for each individual
      if (person.isFounder())
         poly[i] = globalRandom.Normal() * sPolygenes;
      else
         poly[i] = globalRandom.Normal() * sPolygenes * 0.7071068 +
            (poly[person.father->traverse] + poly[person.mother->traverse]) * 0.5;

      // Check if individual is phenotyped ...
      if (!person.isPhenotyped(qtlTrait))
         continue;

      // Incorporate QTL effects and random error
      person.traits[qtlTrait] = poly[i] +
         sQtl * (int(alleles[i * 2] == 1) + int(alleles[i * 2 + 1] == 1)) +
         globalRandom.Normal() * sEnvironment;
      }
   }
 
