////////////////////////////////////////////////////////////////////// 
// merlin/MerlinBitSet.cpp 
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
 
#include "MerlinBitSet.h"
#include "Random.h"
#include "Error.h"

#include <math.h>

// MerlinBitSet
//

double MerlinBitSet::MaxLikelihood(IntArray & current, IntArray & prev, double theta)
   {
   Update(current, prev, theta);

   double adjust = pow(1.0 - theta, ambiguous);

   best = 0;

   for (int i = 1; i < type; i++)
      if (lk[i] > lk[best])
         best = i;

   return lk[best] * adjust;
   }

double MerlinBitSet::Likelihood(IntArray & current, IntArray & prev, double theta)
   {
   Update(current, prev, theta);

   double likelihood = lk[0];

   for (int i = 1; i < type; i++)
      likelihood += lk[i];

   return likelihood;
   }

void MerlinBitSet::SelectMostLikely(IntArray & current, IntArray & prev)
   {
   if (ambiguous)
      switch (type)
         {
         case BS_DEFAULT :
         case BS_FOUNDER :
            {
            for (int i = 0; i < bits.Length(); i++)
               {
               int j = bits[i];

               if (current[j] == 2)
                  current[j] = prev[j] ^ best;
               }
            }
            break;
         case BS_COUPLE :
            {
            int b1 = (best & 1) ? 1 : 0;
            int b2 = (best & 2) ? 1 : 0;
            int b3 = (best & 4) ? 1 : 0;

            for (int i = 0; i < bits.Length(); i++)
               {
               int j = bits[i];

               if (current[j] == 2)
                  {
                  int k =  b3 ? couple_flips[i] : bits[i];
                  current[j] = prev[k] ^ (p1[k] & b1) ^ (p2[k] & b2) ^ (p3[k] & b3);
                  }
               }
            }
            break;
         }
   }

void MerlinBitSet::SelectBySampling(IntArray & current, IntArray & prev, double theta)
   {
   double omt = 1.0 - theta;

   // Determine outcome of ambiguous meiosis by sampling
   if (ambiguous)
      switch (type)
         {
         case BS_DEFAULT :
            // No special symmetries involved ...
            // Outcomes is either same as previous (P = 1.0 - theta),
            // or a recombination (P = theta)
            {
            for (int i = 0; i < bits.Length(); i++)
               {
               int j = bits[i];

               if (current[j] == 2)
                  {
                  bool recombinant = (globalRandom > theta) ? false : true;
                  current[j] = recombinant ? !prev[j] : prev[j];
                  lk[0] *= recombinant ? theta : omt;
                  }
               }
            }
            break;
         case BS_FOUNDER :
            // Founder meiosis ...
            // Decide whether to set bit to zero or one considering underlying
            // possibilities (that is flipped or unflipped founder).
            {
            for (int i = 0; i < bits.Length(); i++)
               {
               int j = bits[i];

               if (current[j] == 2)
                  {
                  double lk0[2] = { lk[0] * (prev[j] == 0 ? omt : theta),
                                    lk[1] * (prev[j] == 1 ? omt : theta) };
                  double p0 = (lk0[0] + lk0[1]) / (lk[0] + lk[1]);
                  int   bit = (globalRandom < p0) ? 0 : 1;

                  current[j] = bit;
                  lk[0] *= prev[j] == bit ? omt : theta;
                  lk[1] *= prev[j] != bit ? omt : theta;
                  }
               }
            }
            break;
         case BS_COUPLE :
            // Symmetric founder couple
            // This is the hard work ... to work out probability of zero bit
            // need to consider all eight underlying possibilities !
            {
            for (int i = 0; i < bits.Length(); i++)
               {
               int j = bits[i];

               if (current[j] == 2)
                  {
                  int k = couple_flips[i];

                  double lk0[8] = {
                     lk[0] * ((prev[j] == 0) ? omt : theta),
                     lk[1] * ((prev[j] ^ p1[j]) == 0 ? omt : theta),
                     lk[2] * ((prev[j] ^ p2[j]) == 0 ? omt : theta),
                     lk[3] * ((prev[j] ^ p1[j] ^ p2[j]) == 0 ? omt : theta),
                     lk[4] * ((prev[k] ^ p3[k]) == 0 ? omt : theta),
                     lk[5] * ((prev[k] ^ p1[k] ^ p3[k]) == 0 ? omt : theta),
                     lk[6] * ((prev[k] ^ p2[k] ^ p3[k]) == 0 ? omt : theta),
                     lk[7] * ((prev[k] ^ p1[k] ^ p2[k] ^ p3[k]) == 0 ? omt : theta) };
                  double p0 =
                     (lk0[0]+lk0[1]+lk0[2]+lk0[3]+lk0[4]+lk0[5]+lk0[6]+lk0[7])/
                     ( lk[0]+ lk[1]+ lk[2]+ lk[3]+ lk[4]+ lk[5]+ lk[6]+ lk[7]);
                  int   bit = (globalRandom < p0) ? 0 : 1;

                  current[j] = bit;
                  lk[0] = bit == 0 ? lk0[0] : lk[0] - lk0[0];
                  lk[1] = bit == 0 ? lk0[1] : lk[1] - lk0[1];
                  lk[2] = bit == 0 ? lk0[2] : lk[2] - lk0[2];
                  lk[3] = bit == 0 ? lk0[3] : lk[3] - lk0[3];
                  lk[4] = bit == 0 ? lk0[4] : lk[4] - lk0[4];
                  lk[5] = bit == 0 ? lk0[5] : lk[5] - lk0[5];
                  lk[6] = bit == 0 ? lk0[6] : lk[6] - lk0[6];
                  lk[7] = bit == 0 ? lk0[7] : lk[7] - lk0[7];
                  }
               }
            }
            break;
         }
   }

void MerlinBitSet::Update(IntArray & current, IntArray & prev, double theta)
   {
   double omt = 1.0 - theta;

   ambiguous = 0;

   switch (type)
      {
      case BS_FOUNDER :
         {
         lk[0] = omt;
         lk[1] = theta;

         for (int i = 0; i < bits.Length(); i++)
            {
            int j = bits[i];

            if (current[j] != 2)
               {
               lk[0] *= prev[j] == current[j] ? omt : theta;
               lk[1] *= prev[j] != current[j] ? omt : theta;
               }
            else
               ambiguous++;
            }
         }
         break;
      case BS_COUPLE :
         {
         lk[0] = omt * omt * omt;
         lk[1] = lk[2] = lk[4] = theta * omt * omt;
         lk[3] = lk[5] = lk[6] = theta * theta * omt;
         lk[7] = theta * theta * theta;

         for (int i = 0; i < bits.Length(); i++)
            {
            int j = bits[i];

            if (current[j] != 2)
               {
               int bit = current[j];
               lk[0] *= prev[j] == bit ? omt : theta;
               lk[1] *= (prev[j] ^ p1[j]) == bit ? omt : theta;
               lk[2] *= (prev[j] ^ p2[j]) == bit ? omt : theta;
               lk[3] *= (prev[j] ^ p1[j] ^ p2[j]) == bit ? omt : theta;

               j = couple_flips[i];
               lk[4] *= (prev[j] ^ p3[j]) == bit ? omt : theta;
               lk[5] *= (prev[j] ^ p1[j] ^ p3[j]) == bit ? omt : theta;
               lk[6] *= (prev[j] ^ p2[j] ^ p3[j]) == bit ? omt : theta;
               lk[7] *= (prev[j] ^ p1[j] ^ p2[j] ^ p3[j]) == bit ? omt : theta;
               }
            else
               ambiguous++;
            }
         }
         break;
      case BS_DEFAULT :
         {
         lk[0] = 1.0;

         for (int i = 0; i < bits.Length(); i++)
            {
            int j = bits[i];

            if (current[j] != 2)
               lk[0] *= prev[j] == current[j] ? omt : theta;
            else
               ambiguous++;
            }
         }
         break;
      }
   }

void MerlinBitSet::CountRecombinants(IntArray & current, IntArray & prev)
   {
   switch (type)
      {
      case BS_FOUNDER :
         {
         rec[0] = 0;
         rec[1] = 1;

         for (int i = 0; i < bits.Length(); i++)
            {
            int j = bits[i];

            if (current[j] == 2)
               error("Internal Error -- CountRecombinants called with ambiguous meiosis\n");

            if (current[j] == prev[j])
               rec[1]++;
            else
               rec[0]++;
            }
         }
         break;
      case BS_COUPLE :
         {
         rec[0] = 0;
         rec[1] = rec[2] = rec[4] = 1;
         rec[3] = rec[5] = rec[6] = 2;
         rec[7] = 3;

         for (int i = 0; i < bits.Length(); i++)
            {
            int j = bits[i];

            if (current[j] == 2)
               error("Internal Error -- CountRecombinants called with ambiguous meiosis\n");

            int bit = current[j];
            rec[0] += (prev[j] != bit);
            rec[1] += (prev[j] ^ p1[j]) != bit;
            rec[2] += (prev[j] ^ p2[j]) != bit;
            rec[3] += (prev[j] ^ p1[j] ^ p2[j]) != bit;

            j = couple_flips[i];
            rec[4] += (prev[j] ^ p3[j]) != bit;
            rec[5] += (prev[j] ^ p1[j] ^ p3[j]) != bit;
            rec[6] += (prev[j] ^ p2[j] ^ p3[j]) != bit;
            rec[7] += (prev[j] ^ p1[j] ^ p2[j] ^ p3[j]) != bit;
            }
         }
         break;
      case BS_DEFAULT :
         {
         rec[0] = 0;

         for (int i = 0; i < bits.Length(); i++)
            {
            int j = bits[i];

            if (current[j] == 2)
               error("Internal Error -- CountRecombinants called with ambiguous meiosis\n");

            if (prev[j] != current[j])
               rec[0]++;
            }
         }
         break;
      }
   }

void MerlinBitSet::SetupCouple(Mantra & m, int couple)
   {
   p1.Dimension(m.bit_count);
   p2.Dimension(m.bit_count);
   p3.Dimension(m.bit_count);
   couple_flips.Clear();

   couple_founder1 = couple_founder2 = -1;

   // Retrieve the indexes for the founders in this couple
   for (int i = 0; i < m.f; i++)
      if (m.couple_index[i] == couple)
         if (m.founder_male[i * 2])
            couple_founder2 = i * 2;
         else
            couple_founder1 = i * 2;

   // This looks a bit messy because the actual bits vector is reversed!
   for (int i = 0, j = m.bit_count - 1; j >= 0; i++, j--)
      if (m.couple_bits[couple][i] == 2)
         {
         p1[j] = 1;
         p2[j] = 0;
         p3[j] = 0;
         couple_flips.Push(j - 1);
         // couple_founder1 = m.vector[m.bits[j]] & ~1;

         j--;
         i++;

         p1[j] = 0;
         p2[j] = 1;
         p3[j] = 0;
         couple_flips.Push(j + 1);
         // couple_founder2 = m.vector[m.bits[j]] & ~1;
         }
      else if (m.couple_bits[couple][i] == 1)
         {
         p1[j] = 0;
         p2[j] = 0;
         p3[j] = 1;
         couple_flips.Push(j);
         }
      else
         {
         p1[j] = p2[j] = p3[j] = 0;
         }
   }

void MerlinBitSet::MapTransition(Mantra & m, int path,
                   IntArray & current, IntArray & previous,
                   IntArray & crossover, IntArray & key)
   {
   switch (type)
      {
      case BS_FOUNDER :
         {
         // First set all meiosis to default state
         for (int i = m.two_f; i < m.two_n; i++)
            if ((m.vector[i] >> 1) == index)
               crossover[i] = path;

         // Then change non-hidden meiosis
         for (int i = 0; i < bits.Length(); i++)
            {
            // Because of ancient mistakes, mantra has bits reversed!
            int j = bits[i];

            crossover[m.bits[j]] = (current[j] != previous[j]) ^ path;
            }

         if (path) Swap(key[index * 2], key[index * 2 + 1]);
         }
         break;
      case BS_COUPLE  :
         {
         // Decode path
         int flip1 = (path & 1) ? 1 : 0;
         int flip2 = (path & 2) ? 1 : 0;
         int flip3 = (path & 4) ? 1 : 0;

         // Are the two grandparents flipped because of previous transitions
         int flipped = (key[couple_founder1] & ~1) != couple_founder1;

         // printf("flips: (%d, %d, %d)\n", flip1, flip2, flip3);

         // First set all meiosis to default state
         for (int i = m.two_f; i < m.two_n; i++)
            if (m.vector[i] < m.two_f && m.couple_index[m.vector[i] >> 1] == index)
               crossover[i] = ((i & 1) ^ flipped) ? flip2 : flip1;
            else if (m.vector[i] >= m.two_f && m.vector[m.vector[i]] < m.two_f &&
                     m.couple_index[m.vector[m.vector[i]] >> 1] == index)
               crossover[i] = flip3;

         int founder1 = flipped ? couple_founder2 : couple_founder1;
         int founder2 = flipped ? couple_founder1 : couple_founder2;

         if (flip3) Swap(key[key.Find(founder1)], key[key.Find(founder2)]);
         if (flip3) Swap(key[key.Find(founder1 + 1)], key[key.Find(founder2 + 1)]);
         if (flip1) Swap(key[founder1], key[founder1 + 1]);
         if (flip2) Swap(key[founder2], key[founder2 + 1]);

         if (flip3) for (int i = m.two_f; i < m.two_n; i++)
            if (m.vector[i] < m.two_f && m.couple_index[m.vector[i] >> 1] == index)
               { Swap(key[i], key[i + 1]); i++; }

         // Then change non-hidden meiosis
         for (int i = 0; i < bits.Length(); i++)
            {
            // Because of ancient mistakes, mantra has bits reversed!
            int j = bits[i];
            int k = flip3 ? couple_flips[i] : j;

            crossover[key[m.bits[j]]] = current[j] != (previous[k] ^
                          (p1[k] & flip1) ^ (p2[k] & flip2) ^ (p3[k] & flip3));
            }

         // crossover.Print("  crossover");
         }
         break;
      case BS_DEFAULT :
         {
         // These bits are easy to handle!
         for (int i = 0; i < bits.Length(); i++)
            {
            // Because of ancient mistakes, mantra has bits reversed!
            int j = bits[i];

            crossover[m.bits[j]] = current[j] != previous[j];
            }
         }
      }
   }

void MerlinBitSet::Rescale(double theta)
   {
   int meioses = bits.Length();

   if (type == BS_FOUNDER) meioses++;
   if (type == BS_COUPLE) meioses+=3;

   double scale = pow(1.0 - theta, -meioses);

   for (int i = 0; i < type; i++)
      lk[i] *= scale;
   }

// MerlinBitSets
//

void MerlinBitSets::BuildSets(Mantra & mantra)
   {
   // Allocate sets
   Dimension(mantra.f + 2);

   // First work out what set each meiosis maps to
   //
   // Founder and couple bits     sets 0 .. f - 1
   // All other female meiosis    set f
   // All other male meiosis      set f + 1
   //
   Index.Dimension(mantra.bit_count);

   // Assign all meiosis to the male / female meisoses groups by default
   for (int i = 0, j = mantra.bit_count - 1; i < mantra.bit_count; i++, j--)
      Index[i] = mantra.f + (mantra.bits[j] & 1);

   // Then extract out founder meiosis
   for (int i = 0; i < mantra.f; i++)
      {
      sets[i].type = BS_FOUNDER;
      sets[i].index = i;
      sets[i].bits.Clear();

      for (int j = 0; j < mantra.bit_count; j++)
         if (mantra.founder_bits[i << 1][j])
            Index[j] = i;
      }

   // And next extract founder couple meiosis
   for (int i = 0; i < mantra.couples; i++)
      {
      // Find out which founder this couple refers to
      int founder = 0;
      while (mantra.couple_index[founder] != i) founder++;

      sets[founder].type = BS_COUPLE;
      sets[founder].index = i;
      sets[founder].bits.Clear();

      for (int j = 0; j < mantra.bit_count; j++)
         if (mantra.couple_bits[i][j])
            Index[j] = founder;

      sets[founder].SetupCouple(mantra, i);
      }


   sets[mantra.f].type = BS_DEFAULT;
   sets[mantra.f].bits.Clear();

   sets[mantra.f + 1].type = BS_DEFAULT;
   sets[mantra.f + 1].bits.Clear();

   // Mantra bits are reversed!!!
   Index.Reverse();

   // Add bits to their appropriate bit_sets
   for (int i = mantra.bit_count - 1; i >= 0; i--)
      {
      sets[Index[i]].bits.Push(i);
      sets[Index[i]].last_bit = i;
      }

   // Set sex codes for each bit set
   for (int i = 0; i < count; i++)
      if (Index[sets[i].last_bit] == i)
         sets[i].maleMeiosis = mantra.bits[sets[i].last_bit] & 1;
   }

void MerlinBitSets::Dimension(int no_of_sets)
   {
   if (no_of_sets > size)
      {
      if (sets) delete [] sets;

      size = (no_of_sets | 0x0F) + 1;
      sets = new MerlinBitSet[size];
      }

   count = no_of_sets;

   for (int i = 0; i < count; i++)
      sets[i].last_bit = 0;
   }

void MerlinBitSets::UpdateLikelihood(int bit, double theta[2], double & lk,
                                    IntArray & current, IntArray & previous)
   {
   int index = Index[bit];

   if (sets[index].last_bit == bit)
      lk *= sets[index].Likelihood(current, previous, theta[sets[index].maleMeiosis]);
   }

void MerlinBitSets::UpdateMaxLikelihood(int bit, double theta[2], double & lk,
                                        IntArray & current, IntArray & previous)
   {
   int index = Index[bit];

   if (sets[index].last_bit == bit)
      lk *= sets[index].MaxLikelihood(current, previous, theta[sets[index].maleMeiosis]);
   }

void MerlinBitSets::SelectMostLikely(IntArray & current, IntArray & prev)
   {
   for (int i = 0; i < count; i++)
      if (Index[sets[i].last_bit] == i)
         sets[i].SelectMostLikely(current, prev);
   }

void MerlinBitSets::SelectBySampling(IntArray & current, IntArray & prev, double theta[2])
   {
   for (int i = 0; i < count; i++)
      if (Index[sets[i].last_bit] == i)
         sets[i].SelectBySampling(current, prev, theta[sets[i].maleMeiosis]);
   }

void MerlinBitSets::UpdateComponents(IntArray & current, IntArray & prev, double theta[2])
   {
   for (int i = 0; i < count; i++)
      if (Index[sets[i].last_bit] == i)
         sets[i].Update(current, prev, theta[sets[i].maleMeiosis]);
   }

void MerlinBitSets::CountRecombinants(IntArray & current, IntArray & prev)
   {
   for (int i = 0; i < count; i++)
      if (Index[sets[i].last_bit] == i)
         sets[i].CountRecombinants(current, prev);
   }

void MerlinBitSets::BestTransition(Mantra & m,
                   IntArray & previous, IntArray & current, double theta[2],
                   IntArray & crossover, IntArray & key)
   {
   UpdateComponents(current, previous, theta);

   for (int i = 0; i < count; i++)
      if (Index[sets[i].last_bit] == i)
         {
         int path = 0;

         for (int j = 1; j < sets[i].type; j++)
            if (sets[i].lk[path] < sets[i].lk[j])
               path = j;

         sets[i].MapTransition(m, path, current, previous, crossover, key);
         }

   // previous.Print("Original iv");
   // current.Print(" Updated iv");
   }

void MerlinBitSets::SampleTransition(Mantra & m,
                   IntArray & previous, IntArray & current, double theta[2],
                   IntArray & crossover, IntArray & key)
   {
   UpdateComponents(current, previous, theta);

   for (int i = 0; i < count; i++)
      if (Index[sets[i].last_bit] == i)
         {
         double sum = sets[i].lk[0];
         int path = 0;

         for (int j = 1; j < sets[i].type; j++)
            {
            sum += sets[i].lk[j];
            if (globalRandom < sets[i].lk[j] / sum)
               path = j;
            }

         sets[i].MapTransition(m, path, current, previous, crossover, key);
         }
   }

void MerlinBitSets::SampleTransition(Mantra & m,
                   IntArray & previous, IntArray & current, double theta[2],
                   IntArray & crossover, IntArray & key, int maxRecombinants)
   {
   UpdateComponents(current, previous, theta);
   CountRecombinants(current, previous);

   for (int i = 0; i < count; i++)
      if (Index[sets[i].last_bit] == i)
         sets[i].Rescale(theta[sets[i].maleMeiosis]);

   sampled_so_far = 0.0;
   if (!RecursiveSampler(m, current, previous, maxRecombinants))
      error("Sampling failed within MerlinBitSet\n");

   for (int i = 0; i < count; i++)
      if (Index[sets[i].last_bit] == i)
         sets[i].MapTransition(m, sets[i].best, current, previous, crossover, key);
   }

bool MerlinBitSets::RecursiveSampler(Mantra & m,
                               IntArray & current, IntArray & previous,
                               int maxRecombinants, double likelihood, int level)
   {
   if (level >= count)
      {
      sampled_so_far += likelihood;
      return (likelihood > 0) && (globalRandom < likelihood / sampled_so_far);
      }

   if (Index[sets[level].last_bit] != level)
      return RecursiveSampler(m, current, previous, maxRecombinants, likelihood, level+1);

   int result = false;
   for (int i = 0; i < sets[level].type; i++)
      if (maxRecombinants >= sets[level].rec[i] &&
          RecursiveSampler(m, current, previous, maxRecombinants - sets[level].rec[i],
          likelihood * sets[level].lk[i], level+1))
          {
          result = true;
          sets[level].best = i;
          }

   return result;
   }


 
