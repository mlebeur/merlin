////////////////////////////////////////////////////////////////////// 
// merlin/MerlinHaplotype.cpp 
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
 
#include "MerlinHaplotype.h"
#include "MerlinCluster.h"
#include "MerlinFamily.h"
#include "MathConstant.h"
#include "MapFunction.h"
#include "Houdini.h"
#include "Random.h"

#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>

// Maximum number of haplotype vectors to list
#define MAXIMUM_HAPLOTYPES    50

// Constructor and destructor
//

bool MerlinHaplotype::outputFounders = false;
bool MerlinHaplotype::outputGraph = true;
bool MerlinHaplotype::outputHorizontal = false;

MerlinHaplotype::MerlinHaplotype()
   {
   output = NULL;
   foutput = NULL;
   goutput = NULL;
   maximum_markers = 0;
   }

MerlinHaplotype::~MerlinHaplotype()
   {
   if (maximum_markers)
      delete [] inheritanceVector;
   }

// Routine for enumerating all zero-recombinant haplotypes
//

void MerlinHaplotype::All(FamilyAnalysis & which)
   {
   // Keep a pointer to FamilyAnalysis structure for later use
   family = &which;

   // Only try to list all haplotypes if there is no recombination
   if (family->analysisPositions.dim != 1)
      error("Cannot list all haplotypes in maps with recombination");

   // The current inheritance vector
   int bits = family->mantra.bit_count;
   current.Dimension(bits);
   current.Zero();

   // Initialize founder key
   flipKey.Dimension(family->mantra.two_n);
   flipKey.SetSequence(0, 1);

   // Add all haplotypes to a chain
   double vectors = pow(2.0, bits);
   HaplotypeChain chain;
   family->right[0].UnPack();
   ListAllRecursively(family->right[0], 0, bits - 1, vectors, chain);
   family->right[0].RePack();

   // Dummy recombination strings
   StringArray      recombination(family->markerCount);
   HaplotypeChain * current = &chain;

   int sets = chain.Length();

   // Output all haplotypes to a file
   do {
      OutputHaplotypes(current->haplo, recombination,
                       "[Set %d/%d including %.0f vectors]",
                       current->index, sets, current->count);
      OutputFounders(current->haplo, "[Zero-Recombinant]", (int) current->count);
      current = current->next;
      }
   while (current != NULL && current->index < MAXIMUM_HAPLOTYPES);

   if (current != NULL)
      fprintf(output, "FAMILY %s Additional %.0f vectors ignored\n\n",
              (const char *) family->family->famid, current->count);
   }

// Output a random inheritance vector for uninformative families
// (eg. where no one is genotyped or founders are homozygous)
//

void MerlinHaplotype::UninformativeFamily(FamilyAnalysis & which, bool sample)
   {
   // Keep a pointer to FamilyAnalysis structure for later use
   family = &which;

   // Allocate inheritance vector
   AllocateVectors(family->markerCount, family->mantra.f);

   // Choose a random inheritance vector
   inheritanceVector[0].Dimension(family->mantra.bit_count);
   for (int i = 0; i < family->mantra.bit_count; i++)
      inheritanceVector[0][i] = globalRandom.Binary();

   // User readable (ie, text) haplotypes
   StringArray * haploString = new StringArray[family->markerCount];
   StringArray recombinantString(family->markerCount);

   // Label all haplotypes
   for (int m = 0; m < family->markerCount; m++)
      {
      if (clusters.IsClustered(family->markers[m]))
         m += LabelCluster(m, inheritanceVector[0], haploString + m, sample);
      else
         LabelChromosomes(m, inheritanceVector[0], haploString[m]);
      }

   OutputHaplotypes(haploString, recombinantString, "[Uninformative]");
   OutputFounders(haploString, "[Uninformative]");

   // Free temporary storage
   delete [] haploString;
   }

void MerlinHaplotype::HaplotypeFamily(FamilyAnalysis & which, bool sample)
   {
   // Initialize pointers which much be freed if memory runs out ...
   Tree * right = NULL;
   StringArray * haploString = NULL;

   try {
      // Keep a pointer to FamilyAnalysis structure for later use
      family = &which;

      // Allocate inheritance vector
      AllocateVectors(family->markerCount, family->mantra.f);

      // Perform multipoint calculation
      right = (sample || family->zeroRecombination) ? family->right :
              new Tree[family->informativeCount];
      Vector scale(family->informativeCount);

      if (!sample && !family->zeroRecombination)
         if (!ScoreConditional(right, scale))
            return;

      // Find the best vector at the first marker
      family->ProgressReport("Choosing Haplotype", 0, family->informativeCount);

      if (right == family->right)
         family->RetrieveMultipoint(0);
      else
         right[0].UnPack();

      if (!sample)
         FindBest(right[0], inheritanceVector[0]);
      else
         Sample(right[0], inheritanceVector[0]);

      if (right == family->right)
         family->ReStoreMultipoint(0);
      else
         right[0].RePack();

      // User readable (ie, text) haplotypes
      int marker = 0;
      haploString = new StringArray[family->markerCount];
      StringArray recombinantString(family->markerCount);

      int positions = family->zeroRecombination ? 1 : family->informativeCount;

      if (positions > 1) bitSet.BuildSets(family->mantra);

      // Find the best vector at other markers
      if (!family->zeroRecombination)
         for (int m = 1; m < positions; m++)
            {
            family->CalculateThetas();

            theta[0] = family->keyFemaleTheta[m - 1];
            theta[1] = family->keyMaleTheta[m - 1];

            family->ProgressReport("Choosing Haplotype", m, positions);

            if (right == family->right)
               family->RetrieveMultipoint(m);
            else
               right[m].UnPack();

            if (!sample)
               FindBest(right[m], inheritanceVector[m], inheritanceVector[m-1]);
            else
               Sample(right[m], inheritanceVector[m], inheritanceVector[m-1]);

            if (right == family->right)
               family->ReStoreMultipoint(m);
            else
               right[m].RePack();

            // Haplotype based on previous marker
            double current = family->informativePositions[m];

            while (family->markerPositions[marker] < current)
               {
               if (clusters.IsClustered(family->markers[marker]))
                  marker += LabelCluster(marker, inheritanceVector[m-1], haploString + marker, sample);
               else
                  LabelChromosomes(marker, inheritanceVector[m-1], haploString[marker]);
               marker++;
               }

            // Label recombinants
            LabelRecombinants(inheritanceVector[m-1], inheritanceVector[m],
                              recombinantString[marker], sample);
            }

      while (marker < family->markerCount)
         {
         if (clusters.IsClustered(family->markers[marker]))
            marker += LabelCluster(marker, inheritanceVector[positions-1], haploString + marker, sample);
         else
            LabelChromosomes(marker, inheritanceVector[positions-1], haploString[marker]);
         marker++;
         }

      OutputHaplotypes(haploString, recombinantString,
                       sample ? "[Sampled]" : "[Most Likely]");
      OutputFounders(haploString, sample ? "[Sampled]" : "[Most Likely]");

      delete [] haploString;
      if (!sample && !family->zeroRecombination) delete [] right;
      }
   catch (TreesTooBig & problem)
      {
      // Free temporary storage
      if (haploString != NULL)
         delete [] haploString;
      if (!sample && !family->zeroRecombination && right != NULL)
         delete [] right;
      throw;
      }
   }

void MerlinHaplotype::OutputHaplotypes(
     StringArray * haploString, StringArray & recombString,
     const char * header_format, ...)
   {
   // Construct a header for each family
   // TODO -- should include likelihood, additional info
   String header;

   header.printf("FAMILY %s ", (const char *) family->family->famid);

   va_list ap;
   va_start(ap, header_format);
   header.vcatprintf(header_format, ap);
   va_end(ap);

   // Output header
   fprintf(output, "%s\n", (const char *) header);
   if (goutput != NULL) fprintf(goutput, "%s\n", (const char *) header);

   if (outputHorizontal)
      {
      HorizontalOutput(haploString);
      return;
      }

   fprintf(output, "\n");
   if (goutput != NULL) fprintf(goutput, "\n");

   String label;

   int colwidth = 20;

   // Work out how wide individual columns should be
   for (int i = 0; i < family->mantra.n; i++)
      {
      Person & person = family->ped[family->family->path[i]];

      int width = person.pid.Length() + 5;

      if (!person.isFounder())
         width += person.father->pid.Length() + person.mother->pid.Length();

      colwidth = min(max(width, colwidth), 80);
      }

   // Adjust width so we have neat use of space on 80-column screens
   if (colwidth > 40) colwidth = 79;
   else if (colwidth > 26) colwidth = 40;
   else if (colwidth > 20) colwidth = 26;

   int columns = 80 / colwidth;

   // Output individuals in groups
   for (int i = 0; i < family->mantra.n; i += columns)
      {
      int end_of_line = min(i + columns, family->mantra.n);

      // Output row of individual labels
      for (int j = i; j < end_of_line; j++)
         {
         Person & person = family->ped[family->family->path[j]];

         if (person.isFounder())
            label = person.pid +  " (F)";
         else
            label = person.pid + " (" + person.motid + "," + person.fatid + ")";

         int width = min(colwidth, label.Length());
         int left_pad = (colwidth - width) >> 1;
         int right_pad = colwidth - width - left_pad;

         fprintf(output, "%*s%*.*s%*s", left_pad, "",
                 width, width, (const char *) label, right_pad, "" );

         if (goutput != NULL)
            fprintf(goutput, "%*s%*.*s%*s", left_pad, "",
                    width, width, (const char *) label, right_pad, "" );
         }
      fprintf(output, "\n");
      if (goutput != NULL) fprintf(goutput, "\n");

      int right = colwidth / 2 - 2;
      int left = colwidth - right - 4;

      // Output haplotype at each marker
      for (int m = 0; m < family->markerCount; m++)
         for (int j = i; j < end_of_line; j++)
            {
            String & maternal  = haploString[m][j << 1];
            String & paternal  = haploString[m][(j << 1) + 1];
            char recombination = recombString[m].IsEmpty() ? ':' :
                                 recombString[m][j];

#ifdef __CHROMOSOME_X__
            Person & person = family->ped[family->family->path[j]];

            if (person.sex == SEX_MALE)
               fprintf(output, "%*s %c %-*s%c",
                       left, (const char *) maternal,
                       ':',
                       right, ".",
                       j == end_of_line - 1 ? '\n' : ' ');
            else
#endif
            fprintf(output, "%*s %c %-*s%c",
                    left, (const char *) maternal,
                    recombination,
                    right, (const char *) paternal,
                    j == end_of_line - 1 ? '\n' : ' ');

            if (goutput != NULL)
               {
               String & maternal  = haploString[m][(j << 1) + family->mantra.two_n];
               String & paternal  = haploString[m][(j << 1) + 1 + family->mantra.two_n];

#ifdef __CHROMOSOME_X__
               if (person.sex == SEX_MALE)
                  fprintf(goutput, "%*s %c %-*s%c",
                          left, (const char *) maternal,
                          ':',
                          right, ".",
                          j == end_of_line - 1 ? '\n' : ' ');
               else
#endif
                  fprintf(goutput, "%*s %c %-*s%c",
                          left, (const char *) maternal,
                          recombination,
                          right, (const char *) paternal,
                          j == end_of_line - 1 ? '\n' : ' ');
               }
            }

      fprintf(output, "\n");
      if (goutput != NULL) fprintf(goutput, "\n");
      }

   // Separators
   fprintf(output, "\n\n\n");
   if (goutput != NULL) fprintf(goutput, "\n\n\n");
   }

void MerlinHaplotype::HorizontalOutput(StringArray * haplo)
   {
   String buffer;

   // Output two haplotypes for each individual
   for (int i = 0; i < family->mantra.two_n; i ++)
      {
      Person & person = family->ped[family->family->path[i >> 1]];

#ifdef __CHROMOSOME_X__
      if ((person.sex == SEX_MALE) && (i & 1)) continue;
#endif

      buffer.printf("%12s %10s ", (const char *) person.pid,
         person.isFounder() ? "(FOUNDER)" : (i&1) ? "(PATERNAL)" : "(MATERNAL)");

      fprintf(output, "%s", (const char *) buffer);
      if (goutput)
         fprintf(goutput, "%s", (const char *) buffer);

      String prefix, suffix;

      for (int m = 0; m < family->markerCount; m++)
         {
         prefix.Clear();
         suffix.Clear();

         for (int j = 0; j < haplo[m][i].Length(); j++)
            if (haplo[m][i][j] == ' ')
               continue;
            else if (isalpha(haplo[m][i][j]))
               prefix += haplo[m][i][j];
            else
               suffix += haplo[m][i][j];

         fprintf(output, "%s%s ", (const char *) prefix, (const char *) suffix);
         if (goutput != NULL) fprintf(goutput, "%2s ", (const char *) haplo[m][i + family->mantra.two_n]);
         }
      fprintf(output, "\n");
      if (goutput != NULL) fprintf(goutput, "\n");
      }
   fprintf(output, "\n");
   if (goutput != NULL) fprintf(goutput, "\n");
   }

void MerlinHaplotype::OutputFounders(StringArray * haplo, const char *, int weight)
   {
   if (foutput == NULL) return;

#ifdef __CHROMOSOME_X__
   int founderHaplotypes = family->mantra.two_f;

   for (int i = 0; i < family->mantra.f; i++)
      if (family->ped[family->family->path[i]].sex == SEX_MALE)
         founderHaplotypes--;
#else
   int founderHaplotypes = family->mantra.two_f;
#endif

   fprintf(foutput, "FAMILY %s, HAPLOTYPES %d, WEIGHT %d\n",
          (const char *) family->family->famid, founderHaplotypes, weight);

   String prefix, suffix;

   for (int i = 0; i < family->mantra.two_f; i ++)
      {
      Person & person = family->ped[family->family->path[i >> 1]];

#ifdef __CHROMOSOME_X__
      if ((i & 1) && person.sex == SEX_MALE) continue;
#endif

      fprintf(foutput, "%s%c ", (const char *) person.pid, i & 1 ? 'B' : 'A');

      for (int m = 0; m < family->markerCount; m++)
         {
         prefix.Clear();
         suffix.Clear();

         for (int j = 0; j < haplo[m][i].Length(); j++)
            if (haplo[m][i][j] == ' ')
               continue;
            else if (isalpha(haplo[m][i][j]))
               prefix += haplo[m][i][j];
            else
               suffix += haplo[m][i][j];

         fprintf(foutput, "%s%s ", (const char *) prefix, (const char *) suffix);
         }

      fprintf(foutput, "\n");
      }
   }

int MerlinHaplotype::LabelCluster
   (int firstMarker, IntArray & vector, StringArray * states, bool sample)
   {
   Mantra & mantra = family->mantra;

   MarkerCluster * cluster = clusters.GetCluster(firstMarker);
   int markers = cluster->markerIds.Length();

   IntArray haplotypes(mantra.two_f);

   cluster->GetHaplotypes(mantra, vector, haplotypes, sample);

   StringArray  founder_alleles;
   for (int m = 0; m < markers; m++)
      {
      MarkerInfo * markerInfo = Pedigree::GetMarkerInfo(cluster->markerIds[m]);

      // First label founder alleles
      founder_alleles.Clear();
      for (int i = 0; i < mantra.two_f; i++)
         if (haplotypes[i] == -1)
            founder_alleles.Add('?');
         else
            founder_alleles.Add(markerInfo->GetAlleleLabel(
                                cluster->RetrieveAllele(haplotypes[i], m)));

      // Now label chromosomes, taking into account symmetries and flips
      for (int i = 0; i < mantra.two_n; i++)
         if (i & 1)
            {
            states[m].Add(" ");
            states[m].Last() += founder_alleles[mantra.state[flipKey[i]]];
            }
         else
            {
            states[m].Add(founder_alleles[mantra.state[flipKey[i]]]);
            states[m].Last() += ' ';
            }

      if (goutput != NULL)
         LabelDescent(vector, states[m]);
      }

   return markers - 1;
   }

void MerlinHaplotype::LabelDescent(IntArray & vector, StringArray & states)
   {
   Mantra & m = family->mantra;

   // Descent graph labels will be added to the end of allele state array
   states.Dimension(m.two_n);

   String label('A');

   // Labeling founder chromosomes is straight-forward
   for (int i = 0; i < m.two_f; i++)
      {
      states.Add(label);

      if (label.Last() != 'Z')
         label.Last()++;
      else
         if (label.Length() == 1)
            label = "AA";
         else
            label[0]++, label[1] = 'A';
      }

   // m.vector.Print("inheritance");
   // m.state.Print("fndr states");

   // Now label descendant chromosomes, taking into account symmetries and flips
   for (int i = m.two_f; i < m.two_n; i++)
      {
      states.Add(states[flipKey.Find(m.state[flipKey[i]]) + m.two_n]);

      // printf("%d maps to %d -> %s\n",
      //       i, m.state[flipKey[i]], (const char *) states[flipKey.Find(m.state[flipKey[i]]) + m.two_n]);
      }
   }

void MerlinHaplotype::LabelChromosomes
   (int marker, IntArray & vector, StringArray & states)
   {
   Mantra & m = family->mantra;
   IntArray graph(m.two_f), fixed(m.two_f);
   LongArray alleles(m.two_f), alleles2(m.two_f);

   m.SelectMarker(family->markers[marker]);

   FuzzyInheritanceTree engine;

   // Clear previous allele labeles
   states.Clear();

   // Approportion recombination
   double likelihood = engine.FuzzyEvaluateInheritanceVector(family->mantra,
                              vector, graph, alleles, alleles2, fixed);

   // The selected inheritance vector is incompatible with this marker
   // (either we made a mistake, or there is bad inheritance)
   if (likelihood == 0.0)
      {
      family->PrintMessage("  WARNING: Marker %s can't be haplotyped",
              (const char *) family->ped.markerNames[family->markers[marker]]);

      for (int i = 0; i < m.two_n; i++)
         states.Add(i & 1 ? " ?" : "? ");

      return;
      }

   // Labels for sets of haplotypes of uncertain phase
   String groups(' ', m.two_f);
   char   group = 'A';

   // Label each group with unique label A, B, ...
   for (int i = 0; i < m.two_f; i++)
      if (!fixed[i] && graph[i] == i && alleles[i] != NOTZERO)
         groups[i] = group++;

   // Assign the same label to all alleles in a group
   for (int i = 0; i < m.two_f; i++)
      groups[i] = groups[graph[i]];

   // Allele states for each founder haplotype
   StringArray founder_alleles;

   // Need to get markerInfo to translate internal allele codes into allele names
   MarkerInfo * markerInfo = Pedigree::GetMarkerInfo(family->markers[marker]);

   // Label each founder haplotype with possible allele states
   for (int i = 0; i < m.two_f; i++)
      if (fixed[i])
         founder_alleles.Add(NameAllele(markerInfo, alleles[i]));
      else if (alleles[i] == NOTZERO)
         founder_alleles.Add("?");
      else
         founder_alleles.Add(NameAllele(markerInfo, alleles[i] ^ alleles2[i]) + "," +
                             NameAllele(markerInfo, alleles2[i]));

   String state;

   // Now label chromosomes, taking into account symmetries and flips
   for (int i = 0; i < m.two_n; i++)
      {
      int allele = m.state[flipKey[i]];
      state = (i & 1) ? groups[allele] : (char) 0;
      state += founder_alleles[allele];
      state += (i & 1) ? (char) 0 : groups[allele];

      // Flag genotypes where the most likely true genotype
      // does not match the observed genotype
      if ((i & 1) && (m.genotype[i] != NOTZERO))
         {
         int allele2 = m.state[flipKey[i - 1]];

         if ((alleles[allele] | alleles[allele2]) != m.genotype[i])
            state += " **";
         }

      states.Add(state);
      }

   if (goutput != NULL)
      LabelDescent(vector, states);
   }

void MerlinHaplotype::LabelRecombinants
   (IntArray & vector1, IntArray & vector2, String & recomb, bool sample)
   {
   Mantra & m = family->mantra;

   IntArray & meiosis = current;
   meiosis.Dimension(m.two_n);
   meiosis.Zero();

   // flipKey.Print("Before");

   if (!sample)
      bitSet.BestTransition(m, vector1, vector2, theta, meiosis, flipKey);
   else if (Multipoint::maximum_recombinants == 0)
      bitSet.SampleTransition(m, vector1, vector2, theta, meiosis, flipKey);
   else
      bitSet.SampleTransition(m, vector1, vector2, theta, meiosis, flipKey,
                              Multipoint::maximum_recombinants);

   // flipKey.Print(" After");

   // Reset recombination vector
   recomb = ':';
   recomb *= m.f;

   // Mark recombinations
   for (int i = m.two_f; i < m.two_n; i += 2)
      switch (meiosis[i] + meiosis[i + 1] * 2)
         {
         case 0 : recomb += '|'; break;
         case 1 : recomb += '/'; break;
         case 2 : recomb += '\\'; break;
         case 3 : recomb += '+'; break;
         }
   }

String MerlinHaplotype::NameAllele(MarkerInfo * markerInfo, longint binary_code)
   {
   String result;
   longint bit = 1;
   int allele = 1;

   while ((binary_code & bit) == (unsigned int) 0)
      {
      bit <<= 1;
      allele++;
      }

   result = markerInfo->GetAlleleLabel(allele);

   return result;
   }

void MerlinHaplotype::AllocateVectors(int markers, int founders)
   {
   // Memory for tracking of alleles for every marker
   if (markers > maximum_markers)
      {
      if (maximum_markers)
         delete [] inheritanceVector;

      maximum_markers = (markers + 4) & ~3;
      inheritanceVector = new IntArray [maximum_markers];
      }

   // Initialize founder key
   flipKey.Dimension(family->mantra.two_n);
   flipKey.SetSequence(0, 1);
   }

bool MerlinHaplotype::ScoreConditional(Tree * right, Vector & scale)
   {
   MultipointHaplotyping engine;

   if (family->useSparse) engine.SetupConditioning(family->mantra);

   int last = family->informativeCount - 1;

   scale[last] = 1.0;

   family->RetrieveSinglepoint(last);
   right[last].Copy(family->singlepoint[last]);
   family->ReStoreSinglepoint(last);

   for (int m = last - 1; m >= 0; m--)
      {
      family->CalculateThetas();
      family->ProgressReport("Haplotyping", last - m, family->informativeCount);

      theta[0] = family->keyFemaleTheta[m];
      theta[1] = family->keyMaleTheta[m];

      // Retrieve haplotypes at m+1 given m + 1, m + 2, ...
      engine.Clear();
      engine.Copy(right[m+1]);
      right[m+1].Pack();

      if (family->useSparse &&
         (family->informationBounds[m+1] + family->informationBounds[m] > 1.0))
         {
         // Try fast conditioning for informative locations?
         family->RetrieveSinglepoint(m);
         right[m].Copy(family->singlepoint[m]);
         family->ReStoreSinglepoint(m);

         // Calculate conditional probability without intermediates
         engine.Condition(right[m], theta, scale[m+1]);
         }
      else
         {
         // Calculate inheritance at m given m + 1, m + 2, ...
         engine.MoveAlong(family->mantra, theta, scale[m+1]);

         // Calculate inheritance at m given *** m ***, m + 1, ...
         family->RetrieveSinglepoint(m);
         engine.Multiply(family->singlepoint[m]);
         family->ReStoreSinglepoint(m);

         // Store it for future use
         right[m].Copy(engine);
         }

      if (stats.GetMean(right[m]) == 0.0)
         {
         family->PrintMessage
                ("  SKIPPED: Requires impossible recombination pattern");
         return false;
         }

      // If the tree is very unlikely, multiply by a constant to avoid underflow
      scale[m] = stats.mean < family->lowBound ? family->rescale : 1.0;
      }

   right[0].Pack();
   return true;
   }

void MerlinHaplotype::FindBest(Tree & tree, IntArray & best_vector)
   {
   // The current inheritance vector
   current.Dimension(tree.bit_count);
   current.Zero();

   // Likelihood for the most likely vector so far
   sum = 0.0;

   FindBestRecursively(tree, 0, tree.bit_count - 1, best_vector);
   }

void MerlinHaplotype::FindBestRecursively(Tree & tree, int node, int bit,
                                          IntArray & best)
   {
   switch (tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         if (tree.nodes[node].value > sum )
            {
            while (bit >= 0)
               {
               current[bit] = 0;
               bit--;
               }
            sum = tree.nodes[node].value;
            best = current;
            }
         break;
      case TREE_NODE_ONE :
         current[bit] = 0;
         FindBestRecursively(tree, tree.nodes[node].child[0], bit - 1, best);
         break;
      case TREE_NODE_TWO :
         current[bit] = 0;
         FindBestRecursively(tree, tree.nodes[node].child[0], bit - 1, best);
         current[bit] = 1;
         FindBestRecursively(tree, tree.nodes[node].child[1], bit - 1, best);
         break;
      }
   }

void MerlinHaplotype::FindBest(Tree & tree, IntArray & best, IntArray & flanker)
   {
   // If there is no recombination
   if (theta[0] == 0.0 && theta[1] == 0.0)
      {
      best = flanker;
      return;
      }

   // if the loci are unlinked
   if (theta[0] >= 0.499999 && theta[1] >= 0.499999)
      {
      FindBest(tree, best);
      return;
      }

   // If there is reduced recombination ...
   if (Multipoint::maximum_recombinants)
      {
      MultipointHaplotyping new_tree;

      new_tree.MakeTree(flanker, tree.bit_count, 1.0);
      new_tree.MoveAlong(family->mantra, theta);
      new_tree.Multiply(tree);

      FindBest(new_tree, best);
      return;
      }

   // The current inheritance vector
   current.Dimension(tree.bit_count);
   current.Zero();

   // Likelihood for the most likely vector so far
   sum = 0.0;

   // Perform sampling
   FindBestRecursively(tree, 0, tree.bit_count - 1, 1.0, best, flanker);
   }

void MerlinHaplotype::FindBestRecursively(Tree & tree, int node, int bit, double lk,
                                          IntArray & best, IntArray & flanker)
   {
   switch (tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         while (bit >= 0)
            {
            current[bit] = 2;

            bitSet.UpdateMaxLikelihood(bit, theta, lk, current, flanker);
            bit--;
            }
         if (lk * tree.nodes[node].value > sum)
            {
            bitSet.SelectMostLikely(best = current, flanker);
            sum = lk * tree.nodes[node].value;
            }
         break;
      case TREE_NODE_ONE :
         {
         current[bit] = 2;

         bitSet.UpdateMaxLikelihood(bit, theta, lk, current, flanker);
         FindBestRecursively(tree, tree.nodes[node].child[0], bit - 1, lk,
                             best, flanker);
         }
         break;
      case TREE_NODE_TWO :
         {
         double lk2 = lk;

         // First search the more likely path (with no recombinants)
         current[bit] = flanker[bit];
         bitSet.UpdateMaxLikelihood(bit, theta, lk, current, flanker);
         FindBestRecursively(tree, tree.nodes[node].child[flanker[bit]],
                             bit - 1, lk, best, flanker);

         // Now search the other path
         current[bit] = !flanker[bit];
         bitSet.UpdateMaxLikelihood(bit, theta, lk2, current, flanker);
         FindBestRecursively(tree, tree.nodes[node].child[!flanker[bit]],
                             bit - 1, lk2, best, flanker);
         }
         break;
      }
   }

void MerlinHaplotype::Sample(Tree & tree, IntArray & best_vector)
   {
   // The current inheritance vector
   current.Dimension(tree.bit_count);
   current.Zero();

   // Likelihood for the most likely vector so far
   sum = 0.0;

   SampleRecursively(tree, 0, tree.bit_count - 1, best_vector);
   }

void MerlinHaplotype::SampleRecursively(Tree & tree, int node, int bit,
                                        IntArray & best, double weight)
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
         SampleRecursively(tree, tree.nodes[node].child[0],
                           bit - 1, best, weight);
         break;
      case TREE_NODE_TWO :
         current[bit] = 0;
         SampleRecursively(tree, tree.nodes[node].child[0],
                           bit - 1, best, weight * 0.5);
         current[bit] = 1;
         SampleRecursively(tree, tree.nodes[node].child[1],
                           bit - 1, best, weight * 0.5);
         break;
      }
   }

void MerlinHaplotype::SampleRecursively(Tree & tree, int node, int bit, double lk,
                                          IntArray & best, IntArray & flanker)
   {
   switch (tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         {
         while (bit >= 0)
            {
            current[bit] = 2;
            bitSet.UpdateLikelihood(bit, theta, lk, current, flanker);
            bit--;
            }
         lk *= tree.nodes[node].value;
         sum += lk;
         if (lk > 0 && (globalRandom < lk / sum))
            bitSet.SelectBySampling(best = current, flanker, theta);
         sum += lk * tree.nodes[node].value;
         }
         break;
      case TREE_NODE_ONE :
         {
         current[bit] = 2;

         bitSet.UpdateLikelihood(bit, theta, lk, current, flanker);
         SampleRecursively(tree, tree.nodes[node].child[0], bit - 1, lk,
                           best, flanker);
         }
         break;
      case TREE_NODE_TWO :
         {
         double lk2 = lk;

         // First search the more likely path (with no recombinants)
         current[bit] = flanker[bit];
         bitSet.UpdateLikelihood(bit, theta, lk, current, flanker);
         SampleRecursively(tree, tree.nodes[node].child[flanker[bit]],
                           bit - 1, lk, best, flanker);

         // Now search the other path
         current[bit] = !flanker[bit];
         bitSet.UpdateLikelihood(bit, theta, lk2, current, flanker);
         SampleRecursively(tree, tree.nodes[node].child[!flanker[bit]],
                           bit - 1, lk2, best, flanker);
         }
         break;
      }
   }

void MerlinHaplotype::Sample(Tree & tree,
                             IntArray & best_vector, IntArray & flanker)
   {
   // If there is no recombination
   if (theta[0] == 0.0 && theta[1] == 0.0)
      {
      best_vector = flanker;
      return;
      }

   // if the loci are unlinked
   if (theta[0] >= 0.499999 && theta[1] >= 0.499999)
      {
      Sample(tree, best_vector);
      return;
      }

   // If there is reduced recombination ...
   if (Multipoint::maximum_recombinants)
      {
      MultipointHaplotyping new_tree;

      new_tree.MakeTree(flanker, tree.bit_count, 1.0);
      new_tree.MoveAlong(family->mantra, theta);
      new_tree.Multiply(tree);

      Sample(new_tree, best_vector);
      return;
      }

   // The current inheritance vector
   current.Dimension(tree.bit_count);
   current.Zero();

   // Likelihood for the most likely vector so far
   sum = 0.0;

   // Perform sampling
   SampleRecursively(tree, 0, tree.bit_count - 1, 1.0, best_vector, flanker);
   }

void MerlinHaplotype::ListAllRecursively(Tree & tree, int node, int bit,
                                         double n, HaplotypeChain & head)
   {
   switch (tree.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         // TODO -- Haplotype and add to chain...
         {
         StringArray * haplo = new StringArray[family->markerCount];
         for (int m = 0; m < family->markerCount; m++)
            LabelChromosomes(m, current, haplo[m]);
         head.Append(haplo, n, family->markerCount);
         }
         break;
      case TREE_NODE_ONE :
         current[bit] = 0;
         ListAllRecursively(tree, tree.nodes[node].child[0], bit-1, n, head);
         break;
      case TREE_NODE_TWO :
         current[bit] = 0;
         ListAllRecursively(tree, tree.nodes[node].child[0], bit-1, n*.5, head);
         current[bit] = 1;
         ListAllRecursively(tree, tree.nodes[node].child[1], bit-1, n*.5, head);
         break;
      }
   }

void HaplotypeChain::Append(StringArray * states, double vectors, int markers)
   {
   HaplotypeChain * current = this;

   if (count == 0.0)
      {
      count = vectors;
      haplo = states;
      return;
      }

   do {
      bool match = true;

      for (int m = 0; m < markers; m++)
         if (current->haplo[m] != states[m])
            {
            match = false;
            break;
            }

      if (match || current->index == MAXIMUM_HAPLOTYPES)
         {
         current->count += vectors;
         delete [] states;
         return;
         }

      if (current->next == NULL)
         break;

      current = current->next;
   } while (true);

   current->next = new HaplotypeChain(current->index + 1);
   current->next->Append(states, vectors, markers);
   }

int HaplotypeChain::Length()
   {
   HaplotypeChain * current = this;

   int len = 1;

   while (current->next)
      {
      current = current->next;
      len++;
      }

   return len;
   }

FILE * MerlinHaplotype::OpenFile(const char * extension)
   {
   String filename(MerlinCore::filePrefix);
   filename += extension;

   FILE * f = fopen(filename, "wt");

   if (f == NULL)
      error("Error opening file [%s]", (const char *) filename);

   return f;
   }

void MerlinHaplotype::OpenFile()
   {
   output = OpenFile(".chr");

   if (outputFounders)
      foutput = OpenFile(".hap");

   if (outputGraph)
      goutput = OpenFile(".flow");
   }

void MerlinHaplotype::CloseFile()
   {
   fclose(output);
   printf("Haplotyping results in file [%s.chr]\n", (const char *) MerlinCore::filePrefix);

   if (foutput != NULL)
      {
      fclose(foutput);
      printf("Founder haplotypes in file [%s.hap]\n", (const char *) MerlinCore::filePrefix);
      }

   if (goutput != NULL)
      {
      fclose(goutput);
      printf("Gene flow graphs in file [%s.flow]\n", (const char *) MerlinCore::filePrefix);
      }
   }




 
