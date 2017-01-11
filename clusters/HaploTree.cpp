////////////////////////////////////////////////////////////////////// 
// clusters/HaploTree.cpp 
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
 
#include "HaploTree.h"

HaploTree::HaploTree()
   {
   depth = 0;
   size  = 0;
   count = 0;
   leafs = 0;
   branches = NULL;
   }

HaploTree::~HaploTree()
   {
   for (int i = 0; i < count; i++)
      delete branches[i];

   if (branches != NULL)
      delete [] branches;
   }

int HaploTree::AllocateFrequency()
   {
   return leafs++;
   }

int HaploTree::AllocateBranch(int level)
   {
   if (count + 1 >= size)
      Grow(count + 1);

   branches[count] = new IntArray;
   branches[count]->Dimension(alleleCounts[level]);
   branches[count]->Set(-1);

   levels[count] = level;
   entries[count] = 0;

   return count++;
   }

void HaploTree::Grow(int minsize)
   {
   int newsize = size ? size : 256;

   while (newsize <= minsize)
      newsize *= 2;

   if (newsize > size)
      {
      IntArray **newbranches = new IntArray * [newsize];

      for (int i = 0; i < count; i++)
         newbranches[i] = branches[i];

      if (branches != NULL)
         delete [] branches;

      branches = newbranches;
      size = newsize;

      levels.Dimension(newsize);
      entries.Dimension(newsize);
      }
   }

int HaploTree::PeekBranch(int branch, int allele) const
   {
   if (count == 0)
      return -1;

   return (*branches[branch])[allele];
   }

int HaploTree::GetBranch(int branch, int allele)
   {
   if (branch == 0 && count == 0)
      AllocateBranch(0);

   if ((*branches[branch])[allele] != -1)
      return (*branches[branch])[allele];

   int result = levels[branch] < depth - 1 ?
                AllocateBranch(levels[branch] + 1) :
                AllocateFrequency();

   (*branches[branch])[allele] = result;

   entries[branch]++;

   return result;
   }

int HaploTree::AddHaplotype(const IntArray & state)
   {
   int branch = 0;

   for (int i = 0; i < depth; i++)
      branch = GetBranch(branch, state[i]);
   
   return branch;
   }

void HaploTree::SetupTraversal(IntArray & pointer, IntArray & state) const
   {
   pointer.Clear();
   pointer.Push(0);

   state.Dimension(depth);
   state.Set(-1);
   }

int HaploTree::Traverse(IntArray & pointer, IntArray & state) const
   {
   while (pointer.Length())
      {
      int branch = pointer.Peek();
      int & allele = state[levels[branch]];
      int next = -1;

      while (++allele < alleleCounts[levels[branch]])
         if ((next = PeekBranch(branch, allele)) != -1)
            break;

      if (next == -1)
         {
         allele = -1;
         pointer.Pop();
         continue;
         }

      if (levels[branch] == depth - 1)
         return next;

      pointer.Push(next);
      }

   return -1;
   }

void HaploTree::MakeEmptyTree(const IntArray & counts)
   {
   alleleCounts = counts;
   depth = counts.Length();
   leafs = 0;
   }

void HaploTree::MakeShrub(int alleles)
   {
   alleleCounts.Dimension(1);
   alleleCounts[0] = alleles;

   freqs.Dimension(alleles);
   new_freqs.Dimension(alleles);

   depth = 1;
   leafs = alleles;

   AllocateBranch(0);

   entries[0] = alleles;

   for (int i = 0; i < alleles; i++)
      (*branches[0])[i] = i;
   }

int HaploTree::Traverse(IntArray & pointer, IntArray & state, double minfreq) const
   {
   while (true)
      {
      int result = Traverse(pointer, state);

      if (result == -1)
         return -1;

      if (freqs[result] > minfreq)
         return result;
      }
   }

void HaploTree::Print(double minfreq) const
   {
   IntArray pointer, state;

   SetupTraversal(pointer, state);

   int    leaf;
   int    count = 0;
   double total = 0.0;

   while ((leaf = Traverse(pointer, state, minfreq)) != -1)
      {
      printf("%6.2f%% ", freqs[leaf] * 100);

      for (int i = 0; i < depth; i++)
         printf("%d", state[i] + 1);

      printf("\n");

      count++;
      total += freqs[leaf];
      }

   printf("\nThese %d haplotypes account for %.2f%% of total probability\n",
          count, total * 100.0);
   }

void HaploTree::Print(StringArray * labels, double minfreq) const
   {
   IntArray pointer, state;

   SetupTraversal(pointer, state);

   int    leaf;
   int    count = 0;
   double total = 0.0;

   while ((leaf = Traverse(pointer, state, minfreq)) != -1)
      {
      printf("%6.2f%% ", freqs[leaf] * 100);

      for (int i = 0; i < depth; i++)
         if (labels[i].Length())
            printf("%s", (const char *) labels[i][state[i] + 1]);
         else
            printf("%d", state[i] + 1);

      printf("\n");

      count++;
      total += freqs[leaf];
      }

   printf("\nThese %d haplotypes account for %.2f%% of total probability\n",
          count, total * 100.0);
   }

void HaploTree::Copy(const HaploTree & rhs)
   {
   Clear();

   alleleCounts = rhs.alleleCounts;
   depth = rhs.depth;

   IntArray pointer, state;

   rhs.SetupTraversal(pointer, state);

   while (rhs.Traverse(pointer, state, 1e-5) != -1)
      AddHaplotype(state);

   SetupFrequencies();

   freqs = rhs.freqs;
   }

void HaploTree::Merge(const HaploTree & rhs, const HaploTree & lhs)
   {
   Clear();

   alleleCounts = rhs.alleleCounts;
   alleleCounts.Stack(lhs.alleleCounts);

   depth = rhs.depth + lhs.depth;

   IntArray rPointer, rState;

   rhs.SetupTraversal(rPointer, rState);

   IntArray lPointer, lState;
   while (rhs.Traverse(rPointer, rState, 1e-5) != -1)
      {
      lhs.SetupTraversal(lPointer, lState);

      while (lhs.Traverse(lPointer, lState, 1e-5) != -1)
         {
         rState.Stack(lState);
         AddHaplotype(rState);
         rState.Dimension(rhs.depth);
         }
      }

   SetupFrequencies();
   }

void HaploTree::Clear()
   {
   for (int i = 0; i < count; i++)
      delete branches[i];

   leafs = 0;
   count = 0;
   depth = 0;
   }

void HaploTree::SetupFrequencies()
   {
   freqs.Dimension(leafs);
   new_freqs.Dimension(leafs);
   }
 
