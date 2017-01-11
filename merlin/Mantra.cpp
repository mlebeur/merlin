////////////////////////////////////////////////////////////////////// 
// merlin/Mantra.cpp 
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
 
#include "Mantra.h"
#include "Error.h"

#include <string.h>

// For an X Chromosome version of this file, define the symbol
// __CHROMOSOME_X__ during compilation
//

bool Mantra::ignoreCoupleSymmetries = false;

Mantra::Mantra()
   {
   // no memory allocated at start
   founder_allocation = ibd_allocation = couple_allocation = 0;

   // initialize pointers
   founder_bits = couple_bits = NULL;
   ibd = NULL;
   pedigree = NULL;
   family = NULL;

   // no family selected
   n = two_n = f = two_f = hidden_bit_count = bit_count = aff_count = 0;
   quick_bit_count = male_bit_count = hidden_male_bit_count = 0;
   }

Mantra::~Mantra()
   {
   if (couple_bits != NULL) delete [] couple_bits;
   if (founder_bits != NULL) delete [] founder_bits;
   if (ibd != NULL) delete [] ibd;
   }

void Mantra::Dimension()
   {
   // Check if we need to reallocate founder bit array
   if (two_f > founder_allocation)
      {
      if (founder_bits != NULL) delete [] founder_bits;

      founder_allocation = two_f + 16;
      founder_bits = new IntArray[founder_allocation];
      }

   // Allocate the arrays that we are guaranteed to use

   // Arrays for basic pedigree information
   state.Dimension(two_n);
   vector.Dimension(two_n);
   desc.Dimension(two_n);

   // Arrays listing interesting meioses
   bits.Dimension(0);
   bit_mask.Dimension(two_n - two_f);
   bit_sex.Dimension(0);

   founder_male.Dimension(two_f);
   founder_bit_count.Dimension(two_f);
   for (int i = 0; i < two_f; i++)
      founder_bits[i].Dimension(two_n - two_f);

   // Arrays for genotype storage space
   genotype.Dimension(two_n);
   hetero.Dimension(two_n);
   branch_genotyped.Dimension(two_n);

   // Arrays for tracking founder couple symmetries
   if (family->generations < 3) return;

   partners.Dimension(f);
   couple_index.Dimension(f);
   }

void Mantra::Prepare(Pedigree & ped, Family & fam)
   {
   // Temporary list of monozygotic twins found in this pedigree
   IntArray mzTwins;
   IntArray aliases;

   // Initialize pointers
   pedigree = &ped;
   family   = &fam;

   // Total number of individuals
   n = fam.count;
   two_n = 2 * n;

   // Prepare aliases array, which lists co-twins for each individual
   aliases.Dimension(two_n);
   aliases.Zero();

   // Total number of founders
   f = fam.founders;
   two_f = 2 * f;

   Dimension();

#ifndef __CHROMOSOME_X__
   // Check if we have founder couple symmetries
   if (fam.generations >= 3)
      PrepareCouples(ped, fam);
   else
#endif
      couples = 0;

   // Setup founder vectors
   for (int i = 0; i < two_f; i++)
      {
      bool isMale = ped[fam.path[i >> 1]].sex == SEX_MALE;

#ifdef __CHROMOSOME_X__
      if ((i & 1) && isMale)
         {
         state[i]  = i - 1;
         vector[i] = i - 1;
         desc[i]   = 1;

         founder_male[i] = 1;
         founder_bit_count[i] = 0;
         founder_bits[i].Zero();
         }
      else
#endif
         {
         state[i]  = i;
         vector[i] = -1;
         desc[i]   = 1;    // Everyone is their own "descendant"

         founder_bit_count[i] = 0;
         founder_bits[i].Zero();
         founder_male[i] = isMale;
         }
      }

   quick_bit_count = bit_count = male_bit_count = 0;
   hidden_male_bit_count = hidden_bit_count = 0;

   // Setup offspring vectors
   //        even index -- maternal alleles
   //        odd index  -- paternal alleles
   for (int i = two_f, parent; i < two_n; i++)
      {
      Person & p = ped[fam.path[i >> 1]];

      if (i & 1)
#ifdef __CHROMOSOME_X__
         // For females, we draw the inheritance graph as usual,
         // for males, the inheritance graph forces homozygozity
         // (as if, males had 2 copies of the same chromosome)
         parent = p.sex == SEX_MALE ? i - 1 : p.father -> traverse << 1;
#else
         parent = p.father->traverse << 1;
#endif
      else
         parent = p.mother->traverse << 1;

      // Check if the parent has been aliased ...
      if (aliases[parent]) parent = aliases[parent];

      // Check if we are dealing with a monozygotic twin
      if (p.zygosity & 1)
         {
         bool newTwin = true;

         // If so, check if we have listed co-twin in meiosis tree
         for (int j = 0; j < mzTwins.Length(); j++)
            // and if we find a co-twin...
            if ((p.serial != mzTwins[j]) && p.isTwin(ped[mzTwins[j]]))
               {
               // Set this individuals genes to be the same as his co-twins
               parent = (ped[mzTwins[j]].traverse << 1) + (i & 1);

               // This affects the founder allele state and vector ...
               state[i]     = state[parent];
               vector[i]    = parent;
               desc[i]      = 1;
               desc[parent]++;

               // In addition, we will edit the inheritance graph so
               // all offspring point to cotwin
               aliases[i]   = parent;

               // This is not a new twinship!
               newTwin = false;

               break;
               }

         // We finish by adding a new twin to the list...
         if (newTwin)
            mzTwins.Push(p.serial);
         // Or proceeding to the next person
         else
            continue;
         }

      state[i]  = state[parent];
      vector[i] = parent;
      desc[i]   = 1;

#ifdef __CHROMOSOME_X__
      // Only female meiosis (with even indexes) are interesting for X
      if (i & 1)
         {
         if (parent != i - 1) desc[parent]++;
         continue;
         }
#endif

      // Pedigree complexity is reduced if we assume that
      // founders always transmit "maternal" allele to first-born
      // However, this requires some book-keeping to track each
      // founder's offspring...
      if (desc[parent] != 1 && parent < two_f)
         {
         founder_bit_count[parent]++;
         founder_bits[parent][bit_count++] = 1;
         bits.Push(i);
         }

      if (parent >= two_f)
         {
         bits.Push(i);
         bit_count++;
         }

      desc[parent]++;

      // Track founder couple symmetries
      if (couples)
         {
         // Symmetric founder couples have couple_index[] >= 0
         if (parent < two_f && couple_index[parent >> 1] >= 0)
            {
            int couple = couple_index[parent >> 1];

            // Offspring of symmetric founder couple imply fancy symmetries
            if (desc[parent] != 2)
               couple_bits[couple][bit_count - 1] = 2;

            if (!first_born[couple])
               first_born[couple] = i;
            }
         else if (parent >= two_f)
            {
            int grand_parent = vector[parent] >> 1;

            if (grand_parent < f && couple_index[grand_parent] >= 0)
               {
               int couple = couple_index[grand_parent];

               if (!first_grandchild[couple])
                  first_grandchild[couple] = i;

               if (parent == first_born[couple] && !first_of_first[couple])
                  first_of_first[couple] = i;

               couple_bits[couple][two_n - 1]++;
               couple_bits[couple][bit_count - 1] = 1;
               }
            }
         }
      }

   // Hide meiosis for one grandchild per symmetric founder couple
   for (int couple = 0; couple < couples; couple++)
      if (first_of_first[couple])
         // Hiding the first meiosis of the first born is preferred
         {
         HideMeiosis(first_of_first[couple]);
         bit_count--;
         }
      else if (first_grandchild[couple])
         {
         ShowMeiosis(first_born[couple]);
         ShowMeiosis(first_born[couple] + 1);
         HideMeiosis(vector[first_grandchild[couple]]);
         HideMeiosis(vector[first_grandchild[couple]] + 1);
         HideMeiosis(first_grandchild[couple]);
         bit_count--;
         }

   // Order bits so ancestral bits change slowest
   bits.Reverse();

   // Label male and female meiosis
   bit_sex.Dimension(bits.Length() + 1);
   bit_sex[0] = 0;
   for (int i = 0; i < bits.Length(); i++)
      bit_sex[i + 1] = (bits[i] & 1);
   male_bit_count = bit_sex.Sum();

   // Tag variable links in inheritance vector
   isBit.Dimension(two_n);
   isBit.Zero();
   for (int i = 0; i < bit_count; i++)
      isBit[bits[i]] = 1;

   // score the total number of descendants for each person
   desc.Set(1);
   for (int i = two_n - 1; i >= 0; i --)
      if (vector[i] != -1)
         {
         desc[vector[i]] += desc[i];

         if (isBit[i])
            desc[vector[i] ^ 1] += desc[i];
         }

   // Score the effective number of hidden bits
   // (that is, excluding trivial founder symmetries where a founder
   //  has only one descendant... this information is required for
   //  likelihood calculations)
   for (int i = 0; i < two_f; i += 2)
      {
      hidden_bit_count += founder_bit_count[i] != 0;
      quick_bit_count  += founder_bit_count[i] == 1;
      }

   for (int i = 0; i < two_f; i += 2)
      hidden_male_bit_count += founder_bit_count[i] != 0 && founder_male[i];

   for (int i = 0; i < couples; i ++)
      hidden_bit_count += couple_bits[i][two_n - 1] != 0;
   }

// PrepareIBD sets up a matrix of known, invariant IBD relationships

void Mantra::PrepareIBD()
   {
   // Check if we need to reallocate IBD matrix
   if (n > ibd_allocation)
      {
      if (ibd != NULL) delete [] ibd;

      ibd_allocation = n + 8;
      ibd = new IntArray[ibd_allocation];

      for (int i = 0; i < ibd_allocation; i++)
         ibd[i].Dimension(ibd_allocation);
      }

   // Founders are easy (IBD_ONE for self, IBD_ZERO for others)
   for (int i = 0; i < f; i++)
      {
#ifdef __CHROMOSOME_X__
      // Check whether we are dealing with male X chromosomes
      bool isMaleX = vector[(i << 1) + 1] != -1;

      // If so, take "counter-measures"
      ibd[i][i] = isMaleX ?  MANTRA_IBD_ONE_MALE : MANTRA_IBD_ONE;
#else
      ibd[i][i] = MANTRA_IBD_ONE;
#endif

      for (int j = 0; j < i; j++)
         ibd[j][i] = ibd[i][j] = MANTRA_IBD_ZERO;
      }

   // Non-founders are trickier
   // IBD_ONE for self if no inbreeding
   // IBD_HALF for parents if no inbreeding
   // IBD_ZERO for unrelated individuals
   // IBD_UNKNOWN for others
   for (int i = f; i < n; i++)
      {
      int mother = vector[(i << 1)] >> 1;
      int father = vector[(i << 1) + 1] >> 1;
#ifdef __CHROMOSOME_X__
      // Special handling for male X chromosomes
      if (father == i)
         {
         ibd[i][i] = MANTRA_IBD_ONE_MALE;

         for (int j = 0; j < i; j++)
            {
            int ibd_state;

            if (ibd[mother][j] == MANTRA_IBD_ZERO)
               ibd_state = MANTRA_IBD_ZERO;
            else if (mother != j)
               ibd_state = MANTRA_IBD_UNKNOWN;
            else
               ibd_state = MANTRA_IBD_HALF_MALE;

            ibd[j][i] = ibd[i][j] = ibd_state;
            }
         }
      else
         {
#endif
      int inbred = ibd[mother][father] != MANTRA_IBD_ZERO;

      ibd[i][i] = inbred ? MANTRA_IBD_UNKNOWN : MANTRA_IBD_ONE;

      for (int j = 0, ibd_state; j < i; j++)
         {
         if (ibd[mother][j] == MANTRA_IBD_ZERO &&
             ibd[father][j] == MANTRA_IBD_ZERO)
            ibd_state = MANTRA_IBD_ZERO;
         else if (inbred || (mother != j && father != j))
            ibd_state = MANTRA_IBD_UNKNOWN;
#ifdef __CHROMOSOME_X__
         // Special handling for father daughter IBD
         else if (father == j)
            ibd_state = MANTRA_IBD_HALF_MALE;
#endif
         else /* mother == j || father == j (AUTOSOME) or mother == j (X) */
            ibd_state = ibd[j][j] == MANTRA_IBD_ONE ?
                    MANTRA_IBD_HALF : MANTRA_IBD_UNKNOWN;

         ibd[j][i] = ibd[i][j] = ibd_state;
         }
#ifdef __CHROMOSOME_X__
         }
#endif
      }
   }

void Mantra::SelectMarker(int m)
   {
   Pedigree & ped = *pedigree;
   Family   & fam = *family;

   // This code assumes everyone inherits a grand-maternal alleles
   Reset();

   // Store genotypes
   for (int i = 0; i < two_n; i++)
      {
      genotype[i] = ped[fam.path[i >> 1]].markers[m].BinaryCoded();
      hetero[i] = ped[fam.path[i >> 1]].markers[m].isHeterozygous();
      branch_genotyped[i] = 0;
      }

   // Count genotype descendants on each branch
   for (int i = two_n - 1; i >= 0; i--)
      {
      if (genotype[i] != NOTZERO)
         branch_genotyped[i]++;

      if (vector[i] == -1)
         continue;

      branch_genotyped[vector[i]] += branch_genotyped[i];

      if (isBit[i])
         branch_genotyped[vector[i] ^ 1] += branch_genotyped[i];
      }

   // Mark completely uninformative meioses
   for (int i = 0; i < bit_count; i++)
      {
      int pivot = bits[i];
      int parent = vector[pivot];

      bit_mask[i] = branch_genotyped[pivot] == 0 ||
                    genotype[parent] != NOTZERO && !hetero[parent];
      }

   // Link to allele frequencies
   frequencies = ped.GetMarkerInfo(m)->freq;
   }

void Mantra::SelectAffection(int affection)
   {
   Pedigree & ped = *pedigree;
   Family   & fam = *family;

   // This code assumes everyone inherits a grand-maternal alleles
   Reset();

   // Dimension affection status storage arrays
   aff.Dimension(two_n);
   branch_affected.Dimension(two_n);

   // Setup affection status vector
   for (int i = 0; i < two_n; i++)
      {
      aff[i] = ped[fam.path[i >> 1]].affections[affection] == 2;
      branch_affected[i]  = 0;
      }

   // score the total number of affected descendants for each person
   // NOTE: The number of affected descendants includes self!
   for (int i = two_n - 1; i >= 0; i--)
      {
      if (aff[i])
         branch_affected[i]++;

      if (vector[i] == -1)
         continue;

      branch_affected[vector[i]] += branch_affected[i];

      if (isBit[i])
         branch_affected[vector[i] ^ 1] += branch_affected[i];
      }

   // list all affecteds in a affecteds[] array
   affecteds.Dimension(0);
   for (int i = 0; i < two_n; i+= 2)
      if (aff[i])
         affecteds.Push(i);
   aff_count = affecteds.Length();
   }

void Mantra::SelectBinaryTrait(int affection)
   {
   Pedigree & ped = *pedigree;
   Family   & fam = *family;

   // This code assumes everyone inherits a grand-maternal alleles
   Reset();

   // Dimension affection status storage arrays
   aff.Dimension(two_n);
   branch_affected.Dimension(two_n);

   // Setup affection status vector
   for (int i = 0; i < two_n; i++)
      {
      aff[i] = ped[fam.path[i >> 1]].affections[affection];
      branch_affected[i]  = 0;
      }

   // score the total number of affected descendants for each person
   // NOTE: The number of affected descendants includes self!
   for (int i = two_n - 1; i >= 0; i--)
      {
      if (aff[i])
         branch_affected[i]++;

      if (vector[i] == -1)
         continue;

      branch_affected[vector[i]] += branch_affected[i];

      if (isBit[i])
         branch_affected[vector[i] ^ 1] += branch_affected[i];
      }

   // list all affecteds in a affecteds[] array
   affecteds.Dimension(0);
   for (int i = 0; i < two_n; i+= 2)
      if (aff[i])
         affecteds.Push(i);
   aff_count = affecteds.Length();
   }

void Mantra::SelectTrait(int trait, double mean)
   {
   Pedigree & ped = *pedigree;
   Family   & fam = *family;

   // This code assumes everyone inherits a grand-maternal alleles
   Reset();

   // Arrays for storing phenotype information
   pheno.Dimension(two_n);
   branch_phenotyped.Dimension(two_n);

   // Setup phenotypic scores status
   for (int i = 0; i < two_n; i++)
      {
      pheno[i] = ped[fam.path[i >> 1]].traits[trait];
      branch_phenotyped[i]  = 0;
      }

   // score the total number of phenotyped descendants for each person
   // NOTE: The number of phenotyped descendants includes self!
   for (int i = two_n - 1; i >= 0; i--)
      {
      if (pheno[i] != _NAN_)
         {
         pheno[i] -= mean;
         branch_phenotyped[i]++;
         }
      else
         pheno[i] = 0.0;

      if (vector[i] == -1)
         continue;

      branch_phenotyped[vector[i]] += branch_phenotyped[i];

      if (isBit[i])
         branch_phenotyped[vector[i] ^ 1] += branch_phenotyped[i];
      }
   }


void Mantra::Reset()
   {
   for (int i = two_f; i < two_n; i ++)
      if (isBit[i])
         // Reset links in ancestral graph
         state[i] = state[vector[i] &= ~1];
      else
         // Some links are fixed (e.g. MZ twins, X in males)
         state[i] = state[vector[i]];
   }

Mantra & Mantra::operator = (const Mantra & rhs)
   {
   n = rhs.n;
   f = rhs.f;
   two_n = rhs.two_n;
   two_f = rhs.two_f;
   bit_count = rhs.bit_count;
   aff_count = rhs.aff_count;
   hidden_bit_count = rhs.hidden_bit_count;

   state = rhs.state;
   vector = rhs.vector;
   isBit = rhs.isBit;
   desc = rhs.desc;

   aff = rhs.aff;
   affecteds = rhs.affecteds;

   bits = rhs.bits;
   bit_mask = rhs.bit_mask;

   genotype = rhs.genotype;
   hetero = rhs.hetero;

   pheno = rhs.pheno;

   branch_genotyped = rhs.branch_genotyped;
   branch_affected = rhs.branch_affected;
   branch_phenotyped = rhs.branch_phenotyped;

   founder_bit_count = rhs.founder_bit_count;
   for (int i = 0; i < two_f; i++)
      founder_bits[i] = rhs.founder_bits[i];

   return *this;
   }

void Mantra::PrepareCouples(Pedigree & ped, Family & fam)
   {
   // No couples by default
   couples = 0;

   if (ignoreCoupleSymmetries) return;

   // First we find all unique partnerships
   partners.Set(PARTNER_NONE);

   for (int i = f; i < n; i++)
      {
      // For each couple
      int father = ped[fam.path[i]].father->traverse;
      int mother = ped[fam.path[i]].mother->traverse;

      // If already registered ... no problem
      if (father < f && partners[father] == mother)
         continue;

      // Track founder couples
      if (father < f && partners[father] == PARTNER_NONE &&
          mother < f && partners[mother] == PARTNER_NONE)
         {
         partners[father] = mother;
         partners[mother] = father;
         continue;
         }

      // No founder symmetries if multiple matings, or matings
      // with non-founders are involved.
      if (father < f)
         {
         if (partners[father] >= 0)
            partners[partners[father]] = PARTNER_ASSYMETRIC;
         partners[father] = PARTNER_ASSYMETRIC;
         }

      if (mother < f)
         {
         if (partners[mother] >= 0)
            partners[partners[mother]] = PARTNER_ASSYMETRIC;
         partners[mother] = PARTNER_ASSYMETRIC;
         }
      }

   // Then look for couples that are not phenotyped
   couple_index.Set(-1);
   for (int i = 0; i < f; i++)
      if (partners[i] > i &&
          EffectivelyIdentical(ped[fam.path[i]], ped[fam.path[partners[i]]]))
         couple_index[i] = couple_index[partners[i]] = couples++;

   // Finally... allocate memory as required
   if (couple_allocation < couples)
      {
      if (couple_allocation) delete [] couple_bits;
      couple_allocation = couples + 8;
      couple_bits = new IntArray[couple_allocation];
      }

   // Dimension temporary arrays
   first_born.Dimension(couples);
   first_grandchild.Dimension(couples);
   first_of_first.Dimension(couples);

   for (int i = 0; i < couples; i++)
      {
      // We only need two_n - two_f - couples bits, but
      // additional bits are used as indicators
      couple_bits[i].Dimension(two_n);
      couple_bits[i].Zero();
      }

   first_born.Zero();
   first_grandchild.Zero();
   first_of_first.Zero();
   }

void Mantra::HideMeiosis(int meiosis)
   {
   for (int i = 0; i < bit_count; i++)
      if (bits[i] == meiosis)
         {
         bits.Delete(i);

         for (int j = 0; j < two_f; j += 2)
            founder_bits[j].Delete(i);

         for (int j = 0; j < couples; j++)
            {
            couple_bits[j].Delete(i);
            couple_bits[j].InsertAt(two_n - 2, 0);
            }

         return;
         }
   }

void Mantra::ShowMeiosis(int meiosis)
   {
   int i = 0;

   while ((i < bit_count) && (bits[i] < meiosis))
      i++;

   bits.InsertAt(i, meiosis);

   for (int j = 0; j < two_f; j += 2)
      founder_bits[j].InsertAt(i, (vector[meiosis] & ~1) == j);

   for (int j = 0; j < couples; j++)
      {
      int parent = vector[meiosis];
      int flip = 0;

      if ((parent < two_f) && (couple_index[parent >> 1] == j))
         flip = 2;

      if ((parent >= two_f) && (vector[parent] < two_f) &&
          (couple_index[vector[parent >> 1]] == j))
         flip = 1;

      couple_bits[j].InsertAt(i, flip);
      couple_bits[j].Delete(two_n - 2);
      }
   }

bool Mantra::EffectivelyIdentical(Person & p1, Person & p2)
   {
   // This function is called to decide whether the founder
   // couple symmetry can be applied to a particular pair
   // of individuals
   if (p1.ngeno != p2.ngeno) return false;

   // If there are any differences in the genotypes of a
   // grandparental couple, then we can't use the founder
   // couple symmetry ...
   for (int m = 0; m < p1.markerCount; m++)
      if (p1.markers[m] != p2.markers[m])
         return false;

   // It must be disabled when we are carrying out any
   // quantitatitive trait analysis and the parental
   // phenotypes differ
   for (int t = 0; t < p1.traitCount; t++)
      if (p1.traits[t] != p2.traits[t])
         return false;

   // It must be disabled when we are carrying out any
   // affecteds only linkage analysis and the parental
   // phenotypes differ
   for (int a = 0; a < p1.affectionCount; a++)
      if ((p1.affections[a] == 2) ^ (p2.affections[a] == 2))
         return false;

   // Finally, it must be disabled when covariates are
   // enabled and they differ between the two individuals
   for (int c = 0; c < p1.covariateCount; c++)
      if (p1.covariates[c] != p2.covariates[c])
         return false;

   return true;
   }


 
