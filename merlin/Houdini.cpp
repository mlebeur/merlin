////////////////////////////////////////////////////////////////////// 
// merlin/Houdini.cpp 
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
 
#include "Houdini.h"

#include <math.h>

double FuzzyInheritanceTree::perAlleleError = 0.0;
double FuzzyInheritanceTree::perGenotypeError = 0.0;

// This portion of the code handles generic memory allocation and freeing
//

FuzzyInheritanceTree::FuzzyInheritanceTree()
   {
   metaAllele = false;
   allocated_founders = 0;
   peelKey = 0;
   }

FuzzyInheritanceTree::~FuzzyInheritanceTree()
   {
   if (allocated_founders)
      {
      delete [] matrix;
      delete [] links;
      delete [] possibleAlleles;
      delete [] alleles;
      delete [] carriers;
      delete [] freq;
      delete [] condition;
      }
   }

void FuzzyInheritanceTree::AllocateMemory(Mantra & m)
   {
   if (m.two_f > allocated_founders)
      {
      int new_allocation = (m.two_f + 7) & ~7;

      ReallocateArray(matrix, allocated_founders, new_allocation);
      ReallocateArray(links, allocated_founders, new_allocation);
      ReallocateArray(possibleAlleles, allocated_founders, new_allocation);
      ReallocateArray(alleles, allocated_founders, new_allocation);
      ReallocateArray(carriers, allocated_founders, new_allocation);
      ReallocateVector(freq, allocated_founders, new_allocation);
      ReallocateVector(condition, allocated_founders, new_allocation);

      for (int i = 0; i < new_allocation; i++)
         {
         matrix[i].Dimension(new_allocation);

         if (i > allocated_founders)
            matrix[i].Zero();
         else
            for (int j = allocated_founders; j < new_allocation; j++)
               matrix[i][j] = 0;
         }

      allocated_founders = new_allocation;
      }

   int new_alleles = CountAlleles(m);
   for (int i = 0; i < m.two_f; i++)
      {
      int allocated_alleles = possibleAlleles[i].Length();

      if (allocated_alleles < new_alleles)
         {
         possibleAlleles[i].Dimension(new_alleles);

         for (int j = allocated_alleles; j < new_alleles; j++)
            possibleAlleles[i][j] = 0;
         }
      }

   // Vector for tracking likelihood of inheritance graph components
   graph_likelihood.Dimension(m.two_f);
   graph_likelihood.Set(1.0);

   // Vector for tracking frequency of unobserved alleles
   metaAlleles.Dimension(m.two_f);
   metaAlleles.Set(1.0);

   // Initialize founder allele graph
   graph.Dimension(m.two_f);
   graph.SetSequence(0, 1);

   weight.Dimension(m.two_f);
   weight.Set(1);

   // Array for tracking doneness in founder allele graph
   graph_descendants.Dimension(m.two_f);

   // Allocate memory for tracking allele states
   alleleState.Dimension(m.two_f);
   alleleStateIndex.Dimension(m.two_f);

   // Allocate memory for tracking genotype and founder allele processing
   if (phenoVisits.Length() < m.two_n)
      {
      peelKey = 0;
      phenoVisits.Dimension(m.two_n);
      }

   if (visited.Length() < m.two_f)
      {
      peelKey = 0;
      visited.Dimension(m.two_f);
      visitIndex.Dimension(m.two_f);
      }
   }

void FuzzyInheritanceTree::ReallocateArray(IntArray * & ptr, int old_size, int new_size)
   {
   IntArray * new_ptr = new IntArray [new_size];

   for (int i = 0; i < old_size; i++)
      new_ptr[i].Swap(ptr[i]);

   if (old_size) delete [] ptr;

   ptr = new_ptr;
   }

void FuzzyInheritanceTree::ReallocateVector(Vector * & ptr, int old_size, int new_size)
   {
   Vector * new_ptr = new Vector [new_size];

   for (int i = 0; i < old_size; i++)
      new_ptr[i].Swap(ptr[i]);

   if (old_size) delete [] ptr;

   ptr = new_ptr;
   }

// This portion of the code handles recursion and generic inheritance
// graph manipulation tasks

void FuzzyInheritanceTree::FuzzyScoreVectors(Mantra & m)
   {
   // Check that analysis should proceed
   if (ShortScoreVectors(m)) return;

   // Clear the tree and define its depth
   Clear();
   bit_count = m.bit_count;

   // Setup up basic graphs and allele states, using founder genotypes
   AllocateMemory(m);

   // Setup probabilities for genotype matches and mismatches
   FillPenetranceMatrix(m);

   // Summarize basic phenotype information
   PreparePhenotypes(m);

   // The work horse
   SetupAlleleList(m);
   ScoreRecursive(m, m.bit_count - 1, 0);
   CleanUpAlleleList(m);

   // The next step is crucial for parametric linkage analyses with
   // sex specific penetrances ...
   //
   // Deal with assymetric founder couples (!)
   UnravelAssymmetries(m);

#ifndef __PRINTOUT_STATS__
   // Compact form for calculations
   switch (mergingStrategy)
      {
      case MERGE_ALL :
         Trim(); break;
      case MERGE_ZEROS :
         TrimNoMerge(); break;
      case MERGE_BOOLEAN :
         MakeBooleanTree(); break;
      }
#endif
   }

void FuzzyInheritanceTree::UnravelAssymmetries(Mantra & m)
   {
   if (m.couples == 0)
      return;

   IntArray problemCouples;

   for (int i = 0; i < m.f; i++)
      if (m.couple_index[i] >= 0 &&                // Individual is part of founder couple pair
          m.partners[i] > i &&                     // And we haven't processed spouse
          !IdenticalPenetrances(i, m.partners[i])) // Penetrances are sex specific
         problemCouples.Push(i);

   if (problemCouples.Length() == 0)
      return;

   int flips = 1 << problemCouples.Length();

   Tree sum;
   sum.MakeMinimalTree(0.0, m.bit_count);

   for (int i = 1; i < flips; i++)
      {
      sum.Add(*this);

      // Swap penetrance matrices for assymetric founder couples
      for (int j = 1, k = 0; j < flips; k++, j *= 2)
         if (i & j)
            SwapPenetrances(problemCouples[k], m.partners[problemCouples[k]]);

      // Prepare Basic Phenotype Information
      PreparePhenotypes(m);

      // Repeat the calculations
      Clear();
      SetupAlleleList(m);
      ScoreRecursive(m, m.bit_count - 1, 0);
      CleanUpAlleleList(m);

      // Swap penetrance matrices for assymetric founder couples
      for (int j = 1, k = 0; j < flips; k++, j *= 2)
         if (i & j)
            SwapPenetrances(problemCouples[k], m.partners[problemCouples[k]]);
      }

   Add(sum);
   }


void FuzzyInheritanceTree::PreparePhenotypes(Mantra & m)
   {
   // Allocate arrays
   isTyped.Dimension(m.two_n);
   typedDescendants.Dimension(m.two_n);
   typedDescendants.Zero();

   // Mark typed individuals
   for (int i = 0; i < m.two_n; i+=2)
      isTyped[i] = isTyped[i + 1] = isPhenotyped(m, i);

   // Count genotype descendants on each branch
   for (int i = m.two_n - 1; i >= 0; i--)
      {
      if (isTyped[i])
         typedDescendants[i]++;

      if (m.vector[i] == -1)
         continue;

      typedDescendants[m.vector[i]] += typedDescendants[i];

      if (m.isBit[i])
         typedDescendants[m.vector[i] ^ 1] += typedDescendants[i];
      }

   // Initialize the number of genotyped descendants for each founder allele
   for (int i = 0; i < m.two_f; i++)
      graph_descendants[i] = typedDescendants[i];
   }

int FuzzyInheritanceTree::ScoreRecursive(Mantra & m, int bit, int last_pivot)
   {
   // Allocate a new tree node
   int node = NewNode();

   // Where in the inheritance vector are we now?
   int pivot = bit >= 0 ? m.bits[bit] : m.two_n;

   // Check whether there are no possible founder allele states
   bool impossible = false;

   // First we update the graph, by assigning an appropriate founder
   // allele to each individual and calculating the potential number
   // of genotyped descendants for each graph component
   int first_update = last_pivot >= m.two_f ? last_pivot : m.two_f;
   for (int i = first_update; i < pivot; i++)
      {
      m.state[i] = m.state[m.vector[i]];

      // Adjust number of genotyped descendants, taking care to avoid twins
      if (m.isBit[i])
         {
         // This section updates the number of potential genotyped descendants
         // for each graph component, by subtracting descendants for
         // non-transmitted alleles
         int component = FindComponent(m.state[m.vector[i] ^ 1]);

         // Gradually peel components with no additional phenotyped descendants
         if (((graph_descendants[component] -= typedDescendants[i]) == 0)
             && !impossible)
            if (!PeelComponent(m, component))
               impossible = true;
         }
      }

   // Next, update the list of possible alleles for each founder and
   // the graph of connected components in pedigree
   for (int i = last_pivot; i < pivot; i++)
      if (isTyped[i])
         {
         // Update the list of possible alleles
         UpdateAlleles(m.state[i], i);

         // Check if there are plausible founder alleles for this founder
         if (!metaAllele && alleles[m.state[i]].Length() == 0)
            impossible = true;

         // Update the founder allele graph
         if (i & 1)
            {
            int a = m.state[i - 1];
            int b = m.state[i];

            carriers[a].Push(i - 1);
            carriers[b].Push(i - 1);

            a = ProcessConnection(a, b);

            // Peel components with no additional phenotyped descendants
            if (((graph_descendants[a] -= 2) == 0) && !impossible)
               if (!PeelComponent(m, a))
                  impossible = true;
            }
         }

   // Is this configuration impossible?
   if (impossible)
      {
      nodes[node].type = TREE_NODE_ZERO;
      nodes[node].value = 0.0;
      }
   // Is this a leaf node?
   else if (bit < 0)
      // If so, evaluate likelihood
      {
      nodes[node].type = TREE_NODE_LEAF;
      nodes[node].value = CalculateLikelihood(m);
      }
   else
      // Find out if we can save time by using symmetries in the data
      switch (nodes[node].type = FindRedundancy(m, pivot))
         {
         case TREE_NODE_ZERO :
            nodes[node].value = 0.0;
            break;
         case TREE_NODE_ONE :
            m.state[pivot] = m.state[m.vector[pivot] &= ~1];
            SAFE_SET(nodes[node].child[0], ScoreRecursive(m, bit - 1, pivot));
            break;
         case TREE_NODE_TWO :
            m.state[pivot] = m.state[m.vector[pivot] &= ~1];
            SAFE_SET(nodes[node].child[0], ScoreRecursive(m, bit - 1, pivot));

            m.state[pivot] = m.state[m.vector[pivot] |= 1];
            SAFE_SET(nodes[node].child[1], ScoreRecursive(m, bit - 1, pivot));
            break;
         }

   // Undo changes in the graph
   for (int i = pivot - 1; i >= last_pivot; i--)
      if (isTyped[i])
         {
         // Undo updates to the list of possible alleles
         UndoAlleles(m.state[i], i);

         // Undo updates to the founder allele graph
         if (i & 1)
            {
            int a = m.state[i - 1];
            int b = m.state[i];

            carriers[a].Pop();
            carriers[b].Pop();

            graph_descendants[FindComponent(a)] += 2;

            RemoveConnection(a,b);
            }
         }

   // Undo changes in descendant counts
   for (int i = pivot - 1; i >= last_pivot; i--)
      // Adjust number of genotyped descendants, taking care to avoid twins
      if (m.isBit[i])
         {
         int component = FindComponent(m.state[m.vector[i] ^ 1]);

         graph_descendants[component] += typedDescendants[i];
         }

   return node;
   }

int FuzzyInheritanceTree::FindRedundancy(Mantra & m, int pivot)
   {
   // Are there are any phenotyped descendants?
   if (typedDescendants[pivot] == 0)
      return TREE_NODE_ONE;

   // Are parents homozygous?
   int parent = m.vector[pivot];
   int matallele = m.state[parent & ~1];
   int fatallele = m.state[parent | 1];

   if (matallele == fatallele)
      return TREE_NODE_ONE;

   // If all the phenotyped descendants of the two alternative
   // founder alleles are on this branch, it doesn't matter which
   // we pick -- one will go down, the other will be fixed (or not)
   if (typedDescendants[matallele] == typedDescendants[pivot] &&
       typedDescendants[fatallele] == typedDescendants[pivot] &&
       mergingStrategy == MERGE_ALL)
       return TREE_NODE_ONE;

   return TREE_NODE_TWO;
   }

// Likelihood calculations routines
//

double FuzzyInheritanceTree::CalculateLikelihood(Mantra & m)
   {
   double likelihood = priorLikelihood;

   /* Go through each component */
   for (int i = 0; i < m.two_f; i++)
      if (graph_descendants[i] == 0)
         likelihood *= graph_likelihood[i];

   return likelihood;
   }

// Founder allele graphing routines
//

int FuzzyInheritanceTree::JoinAlleles(int founder1, int founder2)
   {
   while (graph[founder1] != founder1) founder1 = graph[founder1];
   while (graph[founder2] != founder2) founder2 = graph[founder2];

   if (weight[founder1] > weight[founder2])
      {
      graph[founder2] = founder1;
      weight[founder1] += weight[founder2];
      graph_descendants[founder1] += graph_descendants[founder2];
      graph_journal.Push(founder2);
      return founder1;
      }
   else if (weight[founder2] > weight[founder1] || founder2 != founder1)
      {
      graph[founder1] = founder2;
      weight[founder2] += weight[founder1];
      graph_descendants[founder2] += graph_descendants[founder1];
      graph_journal.Push(founder1);
      return founder2;
      }

   graph_journal.Push(-1);
   return founder1;
   }

void FuzzyInheritanceTree::UndoGraph()
   {
   int founder1 = graph_journal.Pop();

   if (founder1 == -1) return;

   int founder2 = graph[founder1];

   graph[founder1] = founder1;
   weight[founder2] -= weight[founder1];
   graph_descendants[founder2] -= graph_descendants[founder1];
   }

void FuzzyInheritanceTree::ShowAllele(int founder, int allele)
   {
   if (possibleAlleles[founder][allele]++ == 0)
      {
      double frequency = GetAlleleFrequency(founder, allele);

      alleles[founder].Push(allele);
      freq[founder].Push(frequency);
      metaAlleles[founder] -= frequency;
      }
   }

void FuzzyInheritanceTree::HideAllele(int founder, int allele)
   {
   if (--possibleAlleles[founder][allele] == 0)
      for (int i = alleles[founder].Length() - 1; i >= 0; i--)
         if (alleles[founder][i] == allele)
            {
            double frequency = freq[founder][i];

            alleles[founder].Delete(i);
            freq[founder].DeleteDimension(i);
            metaAlleles[founder] += frequency;

            return;
            }
   }

int FuzzyInheritanceTree::ProcessConnection(int founder1, int founder2)
   {
   if (founder1 > founder2)
      {
      int temp = founder1;
      founder1 = founder2;
      founder2 = temp;
      }

   if (matrix[founder1][founder2]++ == 0)
      {
      links[founder1].Push(founder2);
      links[founder2].Push(founder1);

      return JoinAlleles(founder1, founder2);
      }

   while (founder1 != graph[founder1]) founder1 = graph[founder1];

   return founder1;
   }

void FuzzyInheritanceTree::RemoveConnection(int founder1, int founder2)
   {
   if (founder1 > founder2)
      {
      int temp = founder1;
      founder1 = founder2;
      founder2 = temp;
      }

   if (--matrix[founder1][founder2] == 0)
      {
      links[founder1].Pop();
      links[founder2].Pop();

      UndoGraph();
      }
   }

int FuzzyInheritanceTree::FindComponent(int allele)
   {
   while (graph[allele] != allele)
      allele = graph[allele];

   return allele;
   }

// This portion of the code implements genotype error modelling
//

void FuzzyInheritanceTree::FillPenetranceMatrix(Mantra & m)
   {
   // Allocate matrix for storing probability of observed genotypes
   // conditional on the number of matching alleles in the underlying
   // true genotype.
   penetrances.Dimension(6, m.n);

   // And arrays for storing high and low alleles for each genotype
   observedHi.Dimension(m.two_n);
   observedLo.Dimension(m.two_n);

   // Since most possible underlying genotypes will differ from
   // the true genotype, we scale all probabilities to reduce
   // the number of multiplications.
   priorLikelihood = 1.0;

   // Calculate probability of each observed genotype conditional
   // on true genotype, using either a per allele or per genotype
   // error rate.
   if (perAlleleError > 0.0)
      AlleleErrorProbs(m);
   else
      GenotypeErrorProbs(m);

   // Store allele frequency pointer for later use
   allele_frequencies = &m.frequencies;
   metaAllele = true;
   }

#define MATCH_LO_1 0
#define MATCH_LO_2 1
#define MATCH_HI_1 2
#define MATCH_FULL 3
#define MATCH_HI_2 5

void FuzzyInheritanceTree::AlleleErrorProbs(Mantra & m)
   {
   // Error rate and its complement
   double & e = perAlleleError;
   double  ce = 1.0 - e;

   for (int i = 0, allele1, allele2; i < m.two_n; i+=2)
      if (m.genotype[i] != NOTZERO)
         {
         // Find the two alleles in the current genotype
         ExtractAlleles(m.genotype[i], allele1, allele2);

         // Retrieve the corresponding allele frequencies for convenience
         double f1 = m.frequencies[allele1], f2 = m.frequencies[allele2];

         // Most potential true genotypes will not match the observed genotype
         // so define the base probability assuming both alleles are called
         // incorrectly.
         double base = f1 * f2 * e * e;

         // Adjustment factor for heterozygote genotypes
         if (allele1 != allele2) base *= 2;

         double bc = 1.0 / base;

         // Probability of observing a match between the two genotypes is ...
         //   (1-e)^2 + (f1 + f2) * (1-e) * e + base
         penetrances[MATCH_FULL][i/2] = (base + ((f1 + f2) * e + ce) * ce) * bc;

         // Probability that one allele matches the first or second allele is ...
         //  f2 * (1-e) * e + base
         penetrances[MATCH_LO_1][i/2] = (base + f2 * e * ce) * bc;
         penetrances[MATCH_HI_1][i/2] = (base + f1 * e * ce) * bc;

         // Probability that both alleles match the first or second allele is ...
         // 2 * f2 * (1-e) * e + base
         penetrances[MATCH_LO_2][i/2] = (base + 2 * f2 * e * ce) * bc;
         penetrances[MATCH_HI_2][i/2] = (base + 2 * f1 * e * ce) * bc;

         observedHi[i] = observedHi[i + 1] = allele1;
         observedLo[i] = observedLo[i + 1] = allele2;

         priorLikelihood *= base;
         }
   }

void FuzzyInheritanceTree::GenotypeErrorProbs(Mantra & m)
   {
   // Error rate and its complement
   double & e = perGenotypeError;
   double  ce = 1.0 - e;

   for (int i = 0, allele1, allele2; i < m.two_n; i+=2)
      if (m.genotype[i] != NOTZERO)
         {
         // Find the two alleles in the current genotype
         ExtractAlleles(m.genotype[i], allele1, allele2);

         // Most potential true genotypes will not match the observed genotype
         // so define the base probability assuming both alleles are called
         // incorrectly.
         double base = m.frequencies[allele1] * m.frequencies[allele2] * e;

         // Adjustement factor for heterozygotes
         if (allele1 != allele2) base *= 2;

         // Only need to track penetrances for matching genotypes, the default
         // prior likelihood assumes all genotypes mismatch ...
         penetrances[MATCH_FULL][i/2] = (base + ce) / base;

         observedHi[i] = observedHi[i + 1] = allele1;
         observedLo[i] = observedLo[i + 1] = allele2;

         priorLikelihood *= base;
         }

   // Set all the other penetrances to 1.0 (e.g. the same as the base)
   for (int i = 0; i < penetrances.rows; i++)
      if (i != MATCH_FULL)
         penetrances[i].Set(1.0);
   }

void FuzzyInheritanceTree::ExtractAlleles(longint geno, int & allele1, int & allele2)
   {
   longint mask = 1;

   // Identify first allele
   allele1 = 1;
   while ((mask & geno) == 0)
      mask <<= 1, allele1++;

   // Identify second allele
   allele2 = allele1;
   if (mask != geno)
      do {
         mask <<= 1;
         allele2 ++;
      } while ((mask & geno) == 0);
   }

bool FuzzyInheritanceTree::ShortScoreVectors(Mantra & m)
   {
   // If there is no genotyping error in the model, use faster routines
   if (perAlleleError == 0.0 && perGenotypeError == 0.0)
      {
      ScoreVectors(m);
      return true;
      }

   return false;
   }

double FuzzyInheritanceTree::ShortEvaluateVector(
        Mantra & m, int * vector, int * graph,
        longint * alleles1, longint * alleles2, int * fixed)
   {
   // If there is no genotyping error in the model, use faster routines
   if (perAlleleError == 0.0 && perGenotypeError == 0.0)
      return EvaluateInheritanceVector(m, vector, graph, alleles1, alleles2, fixed);

   return _NAN_;
   }

double FuzzyInheritanceTree::GetAlleleFrequency(int founder, int allele)
   {
   return (*allele_frequencies)[allele];
   }

void FuzzyInheritanceTree::UpdateAlleles(int founder, int carrier)
   {
   ShowAllele(founder, observedHi[carrier]);
   ShowAllele(founder, observedLo[carrier]);
   }

void FuzzyInheritanceTree::UndoAlleles(int founder, int carrier)
   {
   HideAllele(founder, observedLo[carrier]);
   HideAllele(founder, observedHi[carrier]);
   }

int FuzzyInheritanceTree::CountAlleles(Mantra & m)
   {
   return m.frequencies.Length();
   }

bool FuzzyInheritanceTree::isPhenotyped(Mantra & m, int who)
   {
   return m.genotype[who] != NOTZERO;
   }

double FuzzyInheritanceTree::JointProbability(Mantra & m)
   {
   // Retrieve allele frequency portion of likelihood for current state
   double lk = product.Peek();

   // Adjust likelihood according to the probability of observing each phenotype
   for (int i = 0; i < phenoList.Length(); i++)
      {
      int who = phenoList[i];

      int al1 = alleleState[m.state[who]];
      int al2 = alleleState[m.state[who + 1]];

      int lo = observedHi[who];
      int hi = observedLo[who];

      int match1 = (al1 == lo) ? 1 : (al1 == hi ? 3 : 0);
      int match2 = (al2 == hi) ? 3 : (al2 == lo ? 1 : 0);
      int match = match1 + match2;

      if (match > 0)
         lk *= penetrances[match - 1][who / 2];
      }

   return lk;
   }

// Peeling routines
//

bool FuzzyInheritanceTree::PeelComponent(Mantra & m, int component)
   {
   // Update a unique index for each round of peeling
   if (peelKey++ == 0)
      {
      visited.Zero();
      phenoVisits.Zero();
      }

   // Print debugging information on current component
   // m.state.Print("Allele States");
   // observedHi.Print("Allele 1");
   // observedLo.Print("Allele 2");
   // printf("Component %d\n", component);
   // for (int i = 0; i < m.two_f; i++)
   //  if (FindComponent(i) == component)
   //     {
   //     printf("Founder allele %d, %d ", i, alleles[i].Length());
   //     alleles[i].Print("states");
   //     printf("  Frequencies: %.3f, ", metaAlleles[i]);
   //     freq[i].Print();
   //     }

   // Peel components
   PeelFounder(m, component, 0);
   PeelBranch(m, 0, false);
   peelSet.Clear();

   double lk = metaAllele ? metaAlleles[component] * condition[component].Last() : 0.0;

   // Update likelihood
   for (int i = 0; i < alleles[component].Length(); i++)
      lk += freq[component][i] * condition[component][i];

   // Update likelihood
   return (graph_likelihood[component] = lk) > 0.0;
   }

int FuzzyInheritanceTree::PeelFounder(Mantra & m, int founder, int order)
   {
   int peelOrder = order;

   peelSet.Push(founder);
   visited[founder] = peelKey;

   condition[founder].Dimension(0);
   visitIndex[founder] = order;

   for (int i = 0; i < links[founder].Length(); i++)
      {
      int founder2 = links[founder][i];

      if (visited[founder2] != peelKey)
         {
         int branchStart = peelSet.Length() - 1;
         int branchOrder = PeelFounder(m, founder2, order + 1);

         if (branchOrder >= order)
            PeelBranch(m, branchStart);
         else
            {
            // Keep current allele (where conditioning occurs)
            // at the top of peel stack
            peelSet.Swap(peelSet.Length() - 1, branchStart);

            // Update peel order, if appropriate
            if (branchOrder < peelOrder)
               peelOrder = branchOrder;
            }
         }
      else
         if (visitIndex[founder2] < peelOrder)
            peelOrder = visitIndex[founder2];
      }

   return peelOrder;
   }

void FuzzyInheritanceTree::PeelBranch(Mantra & m, int branch, bool partial)
   {
   int founder = peelSet[branch];
   int lastAllele = peelSet.Length();

   // Have we allocated conditioning probabilities for this founder yet?
   if (condition[founder].Length() == 0)
      {
      condition[founder].Dimension(alleles[founder].Length() + 1);
      condition[founder].Set(1.0);
      }

   // List genotypes to be peeled ...
   phenoList.Clear();

   // The partial flag is set to "true" (1) for the final round of
   // peeling for each founder allele. This setting ensures that
   // phenotyped individuals who carry no other alleles are processed
   // properly.
   for (int i = branch + partial; i < lastAllele; i++)
      {
      int founder2 = peelSet[i];

      for (int j = 0; j < carriers[founder2].Length(); j++)
         {
         int carrier = carriers[founder2][j];

         if (phenoVisits[carrier] != peelKey)
            {
            phenoVisits[carrier] = peelKey;
            phenoList.Push(carrier);
            }
         }
      }

   InitializeAlleleSet(branch);

   // Partial likelihood conditional on current conditioning state
   double lkSum = 0.0;

   // Iterate over possible states ...
   while (true)
      {
      // Update the overall probability
      lkSum += JointProbability(m);

      // Update allele state
      if (!IncrementAlleleSet(branch + 1))
         {
         int state = alleleStateIndex[founder];

         if (state == -1)
            condition[founder].Last() *= lkSum;
         else
            condition[founder][state] *= lkSum;

         if (state + 1 >= alleles[founder].Length())
            break;

         alleleState[founder] = alleles[founder][++alleleStateIndex[founder]];

         ResetAlleleStates(branch + 1);

         lkSum = 0.0;
         }
      }

   // Remove alleles from graph
   peelSet.Dimension(branch + 1);
   }

void FuzzyInheritanceTree::InitializeAlleleSet(int branch)
   {
   // Separate initialization for the first allele, which is not
   // included in allele frequency product ...
   int founder = peelSet[branch];

   if (metaAllele && metaAlleles[founder] > 1e-10)
      {
      alleleStateIndex[founder] = -1;
      alleleState[founder] = -1;
      }
   else
      {
      alleleStateIndex[founder] = 0;
      alleleState[founder] = alleles[founder][0];
      }

   product.Clear();
   product.Push(1.0);

   ResetAlleleStates(branch + 1);
   }

void FuzzyInheritanceTree::ResetAlleleStates(int branch)
   {
   double partial = product.Peek();

   // Reset allele states and update allele frequency product
   for (int i = branch; i < peelSet.Length(); i++)
      {
      int founder = peelSet[i];

      if (metaAllele && metaAlleles[founder] > 1e-10)
         {
         alleleStateIndex[founder] = -1;
         alleleState[founder] = -1;

         partial *= metaAlleles[founder];

         if (condition[founder].Length())
            partial *= condition[founder].Last();
         }
      else
         {
         int state = alleles[founder][0];

         alleleStateIndex[founder] = 0;
         alleleState[founder] = state;

         partial *= freq[founder][0];

         if (condition[founder].Length())
            partial *= condition[founder][0];
         }

      product.Push(partial);
      }
   }

bool FuzzyInheritanceTree::IncrementAlleleSet(int branch)
   {
   int update = peelSet.Length() - 1;

   if (update < branch)
      return false;

   int founder2 = peelSet[update];
   product.Pop();

   // Find the first allele that can be incremented
   while (alleleStateIndex[founder2] + 1 == alleles[founder2].Length())
      // Check if all allele states have been cycled through
      if (update == branch)
         return false;
      else
         {
         founder2 = peelSet[--update];
         product.Pop();
         }

   int newIndex = ++alleleStateIndex[founder2];
   int newState = alleles[founder2][newIndex];

   alleleState[founder2] = newState;

   product.Push(condition[founder2].Length() ?
      product.Peek() * freq[founder2][newIndex] * condition[founder2][newIndex] :
      product.Peek() * freq[founder2][newIndex]);

   ResetAlleleStates(update + 1);

   return true;
   }

double  FuzzyInheritanceTree::FuzzyEvaluateInheritanceVector
        (Mantra & m, int * vector, int * graph, longint * alleles1, longint * alleles2, int * fixed)
   {
   // Check whether an alternative algorithm exists to solve the problem
   double quickSolution = ShortEvaluateVector(m, vector, graph, alleles1, alleles2, fixed);

   if (quickSolution != _NAN_)
      return quickSolution;

   // Clear the tree and define its depth
   Clear();
   bit_count = m.bit_count;

   // Setup up basic graphs and allele states, using founder genotypes
   AllocateMemory(m);

   // Setup probabilities for genotype matches and mismatches
   FillPenetranceMatrix(m);

   // Summarize basic phenotype information
   PreparePhenotypes(m);

   // Setup the inheritance vector
   for (int i = 0; i < m.bit_count; i++)
      {
      int bit = m.bits[i];

      m.vector[bit] = vector[i] == 0 ? m.vector[bit] & ~1 : m.vector[bit] | 1;
      }

   // Initialize founder allele states to missing
   for (int i = 0; i < m.two_f; i++)
      alleles1[i] = -1;

   // First we update the graph, by assigning an appropriate founder
   // allele to each individual and calculating the potential number
   // of genotyped descendants for each graph component
   for (int i = m.two_f; i < m.two_n; i++)
      {
      m.state[i] = m.state[m.vector[i]];

      // Adjust number of genotyped descendants, taking care to avoid twins
      if (m.isBit[i])
         {
         // This section updates the number of potential genotyped descendants
         // for each graph component, by subtracting descendants for
         // non-transmitted alleles
         int component = FindComponent(m.state[m.vector[i] ^ 1]);

         // Gradually peel components with no additional phenotyped descendants
         if ((graph_descendants[component] -= typedDescendants[i]) == 0)
            FindModalState(m, component, alleles1);
         }
      }

   // Next, update the list of possible alleles for each founder and
   // the graph of connected components in pedigree
   for (int i = 0; i < m.two_n; i++)
      if (isTyped[i])
         {
         // Update the list of possible alleles
         UpdateAlleles(m.state[i], i);

         // Update the founder allele graph
         if (i & 1)
            {
            int a = m.state[i - 1];
            int b = m.state[i];

            carriers[a].Push(i - 1);
            carriers[b].Push(i - 1);

            a = ProcessConnection(a, b);

            // Peel components with no additional phenotyped descendants
            if ((graph_descendants[a] -= 2) == 0)
               FindModalState(m, a, alleles1);
            }
         }

   // Reset data structures
   for (int i = 0; i < m.two_f; i++)
      {
      // Clear list of possible alleles
      for (int j = alleles[i].Length() - 1; j >= 0; j--)
         possibleAlleles[i][alleles[i].Pop()] = 0;

      // Discard allele frequency information
      freq[i].Clear();

      // Clear founder allele graph
      for (int j = links[i].Length() - 1; j >= 0; j--)
         matrix[i][links[i].Pop()] = 0;

      // Discard list of carriers
      carriers[i].Clear();
      }

   // Mark all allele states as part of singleton
   // graphs and as known
   for (int i = 0; i < m.two_f; i++)
      {
      alleles1[i] = alleles1[i] == NOTZERO ? NOTZERO : longint(1) << (int) (alleles1[i] - 1);
      fixed[i] = alleles1[i] != NOTZERO;
      graph[i] = i;
      }

   return CalculateLikelihood(m);
   }

void FuzzyInheritanceTree::FindModalState(Mantra & m, int component, longint * bestAlleles)
   {
   // Update a unique index for each round of peeling
   if (peelKey++ == 0)
      {
      visited.Zero();
      phenoVisits.Zero();
      }

   // Find founder alleles in this component
   for (int i = 0; i < m.two_f; i++)
      if (m.vector[i] == -1 && FindComponent(i) == component)
         peelSet.Push(i);

   // List genotypes to be peeled ...
   phenoList.Clear();
   for (int i = 0; i < peelSet.Length(); i++)
      {
      int founder2 = peelSet[i];

      for (int j = 0; j < carriers[founder2].Length(); j++)
         {
         int carrier = carriers[founder2][j];

         if (phenoVisits[carrier] != peelKey)
            {
            phenoVisits[carrier] = peelKey;
            phenoList.Push(carrier);
            }
         }

      condition[founder2].Clear();
      }

   // The first allele in the current peel set
   int founder = peelSet[0];

   // Initialize allele set and allele frequency product
   InitializeAlleleSet(0);

   // Update allele frequency product to include the first allele
   // (Usually left out to allow for conditioning)
   product *= alleleStateIndex[founder] == -1 ?
              metaAlleles[founder] :
              freq[founder][alleleStateIndex[founder]];

   // Partial likelihood conditional on current conditioning state
   double lkSum = 0.0, lkMax = 0.0;

   // Iterate over possible states ...
   while (true)
      {
      // Probability of observed data for current state
      double lk = JointProbability(m);

      // Update conditional probability
      lkSum += lk;

      // Update modal founder allele state
      if (lk > lkMax)
         {
         for (int i = 0; i < peelSet.Length(); i++)
            bestAlleles[peelSet[i]] = alleleState[peelSet[i]];
         lkMax = lk;
         }

      // Update allele state
      if (!IncrementAlleleSet(1))
         {
         int state = ++alleleStateIndex[founder];

         if (state >= alleles[founder].Length())
            break;

         alleleState[founder] = alleles[founder][state];
         product[0] = freq[founder][state];

         ResetAlleleStates(1);
         }
      }

   // Store component likelihood
   graph_likelihood[component] = lkSum;

   // Remove alleles from graph
   peelSet.Clear();
   }

void FuzzyInheritanceTree::SetupAlleleList(Mantra & )
   {
   }

void FuzzyInheritanceTree::CleanUpAlleleList(Mantra & )
   {
   }

bool FuzzyInheritanceTree::IdenticalPenetrances(int founder1, int founder2)
   {
   for (int i = 0; i < 6; i++)
      if (penetrances[i][founder1] != penetrances[i][founder2])
         return false;

   if (observedHi[founder1 * 2] != observedHi[founder2 * 2] ||
       observedHi[founder1 * 2 + 1] != observedHi[founder2 * 2 + 1] ||
       observedLo[founder1 * 2] != observedLo[founder2 * 2] ||
       observedLo[founder1 * 2 + 1] != observedLo[founder2 * 2 + 1])
       return false;

   return true;
   }

#define SWAPINT(a,b)    {int tmp = (a); (a) = (b); (b) = tmp; }

void FuzzyInheritanceTree::SwapPenetrances(int founder1, int founder2)
   {
   penetrances.SwapColumns(founder1, founder2);

   SWAPINT(observedHi[founder1 * 2], observedHi[founder2 * 2]);
   SWAPINT(observedHi[founder1 * 2 + 1], observedHi[founder2 * 2 + 1]);
   SWAPINT(observedLo[founder1 * 2], observedLo[founder2 * 2]);
   SWAPINT(observedLo[founder1 * 2 + 1], observedLo[founder2 * 2 + 1]);
   }


 
