////////////////////////////////////////////////////////////////////// 
// merlin/Mantra.h 
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
 
#ifndef __MANTRA_H__
#define __MANTRA_H__

#include "Pedigree.h"
#include "LongArray.h"

// Information on IBD gleaned from pedigree
#define MANTRA_IBD_UNKNOWN     0
#define MANTRA_IBD_ZERO        1       // IBD = KINSHIP = 0
#define MANTRA_IBD_HALF        2       // IBD = 1/2, KINSHIP = 1/4
#define MANTRA_IBD_ONE         3       // IBD = 2/2, KINSHIP = 2/4

// Additional states that are possible for the X chromosome
#define MANTRA_IBD_ONE_MALE    4       // IBD = 1/1, KINSHIP = 1/1
#define MANTRA_IBD_HALF_MALE   5       // IBD = 1/2, KINSHIP = 1/2

#define PARTNER_NONE           -1
#define PARTNER_ASSYMETRIC     -2

class Mantra
   {
   public:
      // Family structure information
      IntArray    state;             // allelic state
      IntArray    vector;            // inheritance vector
      IntArray    isBit;             // variable elements of inheritance vector
      IntArray    desc;              // number of descendants

      IntArray    bits;              // index of informative bits in genealogy
      IntArray    bit_mask;          // mask for uninformative bits
      IntArray    bit_sex;           // 0 - female meiosis, 1 - male meiosis

      IntArray *  founder_bits;      // Gametes originating in this founder
      IntArray    founder_bit_count; // No. of gametes from this founder
      IntArray    founder_male;      // 0 - female, 1 - male

      // Information for tracking founder couple symmetries
      IntArray    partners;
      IntArray    couple_index;
      IntArray *  couple_bits;

      // Codomininant marker data
      LongArray   genotype;          // bit coded genotypes
      IntArray    hetero;            // heterozigosity indicators

      // Affection status information
      IntArray    aff;               // affection status
      IntArray    affecteds;         // index of affected individuals

      // Quantitative trait information
      Vector      pheno;             // phenotypic deviates

      // Other helpful structures
      IntArray *  ibd;               // matrix of known ibds

      IntArray branch_phenotyped;    // phenotyped descendants in branch
      IntArray branch_genotyped;     // genotyped descendants in branch
      IntArray branch_affected;      // affected descendants in branch

      int two_n, n, two_f, f, bit_count, male_bit_count, quick_bit_count;
      int aff_count, couples, hidden_bit_count, hidden_male_bit_count;

      Pedigree * pedigree;
      Family   * family;
      Vector     frequencies;        // Marker allele frequencies

      Mantra();
      ~Mantra();

      void Prepare(Pedigree & ped, Family & f);
      void PrepareCouples(Pedigree & ped, Family & f);
      void PrepareIBD();
      void SelectMarker(int m);
      void SelectAffection(int affection);
      void SelectBinaryTrait(int affection);
      void SelectTrait(int trait, double mean = 0.0);
      void Reset();

      Mantra & operator = (const Mantra & rhs);

      // Setting this flag disables checking for founder couple symmetries
      static bool ignoreCoupleSymmetries;

      Person & GetPerson(int i)
         { return (*pedigree)[family->path[i]]; }

   private:
      // Storage management
      int ibd_allocation, founder_allocation, couple_allocation;

      // Temporaries arrays used by the founder couple reduction
      IntArray first_born;
      IntArray first_grandchild;
      IntArray first_of_first;

      // Auxiliary routines for managing founder couple reduction
      void HideMeiosis(int meiosis);
      void ShowMeiosis(int meiosis);

      static bool EffectivelyIdentical(Person & p1, Person & p2);

      void Dimension();
   };

#endif

 
