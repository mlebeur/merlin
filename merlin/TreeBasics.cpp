////////////////////////////////////////////////////////////////////// 
// merlin/TreeBasics.cpp 
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
 
#include "TreeBasics.h"
#include "MathConstant.h"

#include <string.h>
#include <stdio.h>

int BasicTree::maxNodes = 0;
int BasicTree::totalNodes = 0;

void BasicTree::Grow()
   {
   totalNodes -= count;
   count = count ? count * 2 : 256;
   totalNodes += count;

   if (maxNodes && totalNodes > maxNodes) MemoryCeiling();
   TreeNode * _nodes = (TreeNode *) realloc(nodes, count * sizeof(TreeNode));
   if (_nodes == NULL) OutOfMemory();
   nodes = _nodes;
   }

void BasicTree::TrimNoMerge(int node)
   {
   int child, child2;

   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return;
      case TREE_NODE_TWO :
         TrimNoMerge(child = nodes[node].child[0]);
         TrimNoMerge(child2 = nodes[node].child[1]);
         if (nodes[child].type == nodes[child2].type &&
             nodes[child].type == TREE_NODE_ZERO)
             nodes[node].type = TREE_NODE_ZERO;
         return;
      case TREE_NODE_ONE :
         TrimNoMerge(child = nodes[node].child[0]);
         if (nodes[child].type == TREE_NODE_ZERO)
            nodes[node].type = TREE_NODE_ZERO;
         else if (nodes[child].type == TREE_NODE_LEAF)
            nodes[node].type = TREE_NODE_LEAF,
            nodes[node].value = nodes[child].value;
         return;
      case TREE_NODE_LEAF :
         if (nodes[node].value == 0.0)
            nodes[node].type = TREE_NODE_ZERO;
         return;
      }
   }

void BasicTree::Trim(int node)
   {
   int child, child2;

   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return;
      case TREE_NODE_TWO :
         Trim(child = nodes[node].child[0]);
         Trim(child2 = nodes[node].child[1]);
         if (nodes[child].type == nodes[child2].type)
            if (nodes[child].type == TREE_NODE_ZERO)
               nodes[node].type = TREE_NODE_ZERO;
            else if (nodes[child].type == TREE_NODE_LEAF &&
                     nodes[child].value == nodes[child2].value)
               {
               nodes[node].type = TREE_NODE_LEAF;
               nodes[node].value = nodes[child].value;
               }
         return;
      case TREE_NODE_ONE :
         Trim(child = nodes[node].child[0]);
         if (nodes[child].type == TREE_NODE_ZERO)
            nodes[node].type = TREE_NODE_ZERO;
         else if (nodes[child].type == TREE_NODE_LEAF)
            {
            nodes[node].type = TREE_NODE_LEAF;
            nodes[node].value = nodes[child].value;
            }
         return;
      case TREE_NODE_LEAF :
         if (nodes[node].value == 0.0)
            nodes[node].type = TREE_NODE_ZERO;
         return;
      }
   }

void BasicTree::MakeBooleanTree(int node)
   {
   int child, child2;

   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return;
      case TREE_NODE_TWO :
         MakeBooleanTree(child = nodes[node].child[0]);
         MakeBooleanTree(child2 = nodes[node].child[1]);
         if (nodes[child].type == nodes[child2].type)
            if (nodes[child].type == TREE_NODE_ZERO)
               nodes[node].type = TREE_NODE_ZERO;
         return;
      case TREE_NODE_ONE :
         MakeBooleanTree(child = nodes[node].child[0]);
         if (nodes[child].type == TREE_NODE_ZERO)
            nodes[node].type = TREE_NODE_ZERO;
         return;
      case TREE_NODE_LEAF :
         if (nodes[node].value == 0.0)
            nodes[node].type = TREE_NODE_ZERO;
         else
            nodes[node].value = 1.0;
         return;
      }
   }

int BasicTree::Copy(int node)
   {
   int new_node = NewNode();

   switch (nodes[new_node].type = nodes[node].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         nodes[new_node].value = nodes[node].value;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0], Copy(nodes[node].child[0]));
         SAFE_SET(nodes[new_node].child[1], Copy(nodes[node].child[1]));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0], Copy(nodes[node].child[0]));
         break;
      }

   return new_node;
   }

int BasicTree::Copy(const BasicTree & source, int node)
   {
   int new_node = NewNode();

   switch (nodes[new_node].type = source.nodes[node].type)
      {
      case TREE_NODE_ZERO :
         break;
      case TREE_NODE_LEAF :
         nodes[new_node].value = source.nodes[node].value;
         break;
      case TREE_NODE_TWO :
         SAFE_SET(nodes[new_node].child[0], Copy(source, source.nodes[node].child[0]));
         SAFE_SET(nodes[new_node].child[1], Copy(source, source.nodes[node].child[1]));
         break;
      case TREE_NODE_ONE :
         SAFE_SET(nodes[new_node].child[0], Copy(source, source.nodes[node].child[0]));
         break;
      }

   return new_node;
   }

void BasicTree::Copy(const BasicTree & rhs)
   {
   Clear();
   logOffset = rhs.logOffset;
   bit_count = rhs.bit_count;
   Copy(rhs, 0);
   }

void BasicTree::Print(int node, int pad)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         printf("%*s ZERO\n", pad, "");
         break;
      case TREE_NODE_LEAF :
         printf("%*s %.5g [%d]\n", pad, "", nodes[node].value, node);
         break;
      case TREE_NODE_TWO :
         printf("%*s %s\n", pad, "", "Left");
         Print(nodes[node].child[0], pad + 5);
         printf("%*s %s\n", pad, "", "Right");
         Print(nodes[node].child[1], pad + 5);
         break;
      case TREE_NODE_ONE :
         printf("%*s %s\n", pad, "", "Left");
         Print(nodes[node].child[0], pad + 5);
         break;
      }
   }

void BasicTree::MakeMinimalTree(double value, int bits)
   {
   Clear();

   int node = NewNode();

   nodes[node].type = TREE_NODE_LEAF;
   nodes[node].value = value;

   bit_count = bits;
   logOffset = 0;
   }

void BasicTree::MakeEmptyTree(int bits)
   {
   Clear();

   int node = NewNode();

   nodes[node].type = TREE_NODE_ZERO;

   bit_count = bits;
   logOffset = 0;
   }

void BasicTree::CarbonCopy(const BasicTree & rhs)
   {
   if (count < rhs.count)
      {
      totalNodes -= count;
      count = rhs.count;
      totalNodes += count;

      if (maxNodes && totalNodes > maxNodes) MemoryCeiling();
      TreeNode * _nodes = (TreeNode *) realloc(nodes, count * sizeof(TreeNode));
      if (_nodes == NULL) OutOfMemory();
      nodes = _nodes;
      }

   nextFree = rhs.nextFree;
   bit_count = rhs.bit_count;

   memcpy(nodes, rhs.nodes, sizeof(TreeNode) * nextFree);
   }

// This section handles tree swapping
//

// Trees any smaller than this will not be swapped ...
//
int  BasicTree::SWAP_MIN = 16 * 1024;

void BasicTree::UpdateSwapThreshold(int nTrees)
   {
   // Update swap tresholds so that at least nTrees can be stored
   // in memory (it is assumed that about 256 Mb are available)

   int availableMemory = maxNodes ? maxNodes : 512 * 1024 * 1024 / sizeof(TreeNode);

   int threshold = availableMemory / nTrees;

   if (threshold > 16 * 1024 && tmpfileInfo.skeleton != NULL)
      SWAP_MIN = 16 * 1024;
   else if (threshold >= 8)
      SWAP_MIN = threshold;
   else
      SWAP_MIN = 8;
   }

#define SWAP_PACKET    (128 * 1024)
#define SWAP_MASK      (SWAP_PACKET - 1)

char *      BasicTree::skeleton = NULL;
double *    BasicTree::leaves = NULL;
MiniDeflate BasicTree::zip;

TreeManager BasicTree::tmpfileInfo;

void BasicTree::AllocateBuffers()
   {
   if (skeleton == NULL)
      {
      skeleton = new char [SWAP_PACKET];
      leaves = new double [SWAP_PACKET];
      }
   }

void BasicTree::FreeBuffers()
   {
   if (skeleton != NULL)
      {
      delete [] skeleton;
      delete [] leaves;

      skeleton = NULL;
      leaves = NULL;
      }
   }

void BasicTree::SetupSwap()
   {
   if (tmpfileInfo.OpenFiles())
      AllocateBuffers();
   }

void BasicTree::CloseSwap(bool quiet)
   {
   if (tmpfileInfo.skeleton != NULL)
      {
      if (!quiet)
         {
         double fileSize = tmpfileInfo.GetFileSize();

         if (fileSize > 0.0)
            printf("Swap file usage peaked at < %.0fM \n", fileSize * 1e-6 + 1);
         else if (fileSize < 0.0)
            printf("Swap file usage peaked at > 2048 MB.\n");
         }

      tmpfileInfo.CloseFiles();

      FreeBuffers();
      }
   }

void BasicTree::FreeSwap()
   {
   tmpfileInfo.Free();
   }

void BasicTree::RePack()
   {
   if (tmpfileInfo.skeleton == NULL || count < SWAP_MIN) return;

   totalNodes -= count;
   free(nodes);
   nodes = NULL;
   }

void BasicTree::Discard()
   {
//   if (count < SWAP_MIN) return;

   if (nodes != NULL)
      {
      totalNodes -= count;
      free(nodes);
      nodes = NULL;
      }

   count = 0;
   }

void BasicTree::Pack()
   {
   Pack(tmpfileInfo);
   }

void BasicTree::PackOnly()
   {
   PackOnly(tmpfileInfo);
   }

void BasicTree::PackOrDiscard(TreeManager & output)
   {
   if (output.skeleton == NULL)
      Discard();
   else
      Pack(output);
   }

void BasicTree::PackOnly(TreeManager & output)
   {
   if (output.skeleton == NULL || count < SWAP_MIN) return;

   int nextSkel = 0;
   int nextLeaf = 0;

   tmpfileOffset = output.position;
   output.SetPosition(output.position);

   leafCount = 0;

   stack.Clear();
   stack.Push(0);
   nextFree = 0;

   while (stack.Length())
     {
     int i = stack.Pop();

     skeleton[nextSkel++] = (char) nodes[i].type;

     if ((nextSkel & SWAP_MASK) == 0)
       zip.Deflate(output.skeleton, skeleton, SWAP_PACKET),
       nextSkel = 0;

     switch(nodes[i].type)
       {
       case TREE_NODE_LEAF :
         leafCount++;
         leaves[nextLeaf++] = nodes[i].value;
         if ((nextLeaf & SWAP_MASK) == 0)
            zip.Deflate(output.leaves, leaves, SWAP_PACKET * sizeof(double)),
            nextLeaf = 0;
         break;
       case TREE_NODE_TWO :
         stack.Push(nodes[i].child[1]);
       case TREE_NODE_ONE :
         stack.Push(nodes[i].child[0]);
         break;
       }

     nextFree++;
     }

   if (nextSkel) zip.Deflate(output.skeleton, skeleton, nextSkel);
   if (nextLeaf) zip.Deflate(output.leaves, leaves, nextLeaf * sizeof(double));

   output.UpdateFileSize();
   output.UpdatePosition();
   }

void BasicTree::Pack(TreeManager & output)
   {
   if (output.skeleton == NULL || count < SWAP_MIN) return;

   PackOnly(output);

   totalNodes -= count;
   free(nodes);
   nodes = NULL;
   }

void BasicTree::UnPack()
   {
   UnPack(tmpfileInfo);
   }

void BasicTree::UnPack(TreeManager & input)
   {
   if (input.skeleton == NULL || count < SWAP_MIN) return;

   totalNodes += count;

   if (maxNodes && totalNodes > maxNodes) MemoryCeiling();
   TreeNode * _nodes = (TreeNode *) realloc(nodes, count * sizeof(TreeNode));
   if (_nodes == NULL) OutOfMemory();
   nodes = _nodes;

   input.SetPosition(tmpfileOffset);

   int nextSkel = 0;
   int nextLeaf = 0;
   int leavesToGo = leafCount;

   stack.Clear();

   for (int i = 0; i < nextFree; i++)
      {
      if ((nextSkel & SWAP_MASK) == 0)
         zip.Inflate(input.skeleton, skeleton, min(nextFree - i, SWAP_PACKET)),
         nextSkel = 0;

      nodes[i].type = skeleton[nextSkel++];

      if (nodes[i].type == TREE_NODE_LEAF)
         {
         if ((nextLeaf & SWAP_MASK) == 0)
            zip.Inflate(input.leaves, leaves, min(leavesToGo, SWAP_PACKET) * sizeof(double)),
            nextLeaf = 0, leavesToGo -= SWAP_PACKET;

         nodes[i].value = leaves[nextLeaf++];
         }

      switch (nodes[i].type)
         {
         case TREE_NODE_ZERO :
         case TREE_NODE_LEAF :
            if (stack.Length())
               nodes[stack.Pop()].child[1] = i + 1;
            break;
         case TREE_NODE_TWO :
            stack.Push(i);
         case TREE_NODE_ONE :
            nodes[i].child[0] = i + 1;
            break;
         }
      }
   }

void BasicTree::MemoryCeiling()
   {
   // Calculated requested memory
   int request = totalNodes / 1024 * sizeof(TreeNode) / 1024 + 1;

   // Free Memory
   free(nodes);
   totalNodes -= count;
   nodes = NULL;

   // Throw Exception
   throw TreesTooBig(request);
   }

void BasicTree::OutOfMemory()
   {
   // Calculated requested memory
   int request = totalNodes / 1024 * sizeof(TreeNode) / 1024 + 1;

   // Update memory usage
   totalNodes -= count;

   // Reset count to zero to avoid double counting
   count = 0;

   // Throw Exception
   throw TreesTooBig(request);
   }

void BasicTree::ReadFromFile(FILE * input)
   {
   AllocateBuffers();
   Free();

   fread(&bit_count, sizeof(bit_count), 1, input);
   fread(&count, sizeof(count), 1, input);

   totalNodes += count;

   if (maxNodes && totalNodes > maxNodes) MemoryCeiling();
   nodes = (TreeNode *) malloc(count * sizeof(TreeNode));
   if (nodes == NULL) OutOfMemory();

   nextFree = count;

   int nextSkel = 0;
   int nextLeaf = 0;
   int packetSize = 0;

   stack.Clear();

   for (int i = 0; i < nextFree; i++)
      {
      if ((nextSkel & SWAP_MASK) == 0)
         {
         fread(&packetSize, sizeof(packetSize), 1, input);
         if (packetSize) zip.Inflate(input, skeleton, packetSize);
         nextSkel = 0;

         fread(&packetSize, sizeof(packetSize), 1, input);
         if (packetSize) zip.Inflate(input, leaves, packetSize * sizeof(double));
         nextLeaf = 0;
         }

      nodes[i].type = skeleton[nextSkel++];

      if (nodes[i].type == TREE_NODE_LEAF)
         nodes[i].value = leaves[nextLeaf++];

      switch (nodes[i].type)
         {
         case TREE_NODE_ZERO :
         case TREE_NODE_LEAF :
            if (stack.Length())
               nodes[stack.Pop()].child[1] = i + 1;
            break;
         case TREE_NODE_TWO :
            stack.Push(i);
         case TREE_NODE_ONE :
            nodes[i].child[0] = i + 1;
            break;
         }
      }
   }

void BasicTree::WriteToFile(FILE * output)
   {
   AllocateBuffers();

   fwrite(&bit_count, sizeof(bit_count), 1, output);

   int placeholder = ftell(output), nodeCount = 0;

   fwrite(&nodeCount, sizeof(nodeCount), 1, output);

   int nextSkel = 0;
   int nextLeaf = 0;

   stack.Clear();
   stack.Push(0);

   while (stack.Length())
      {
      int i = stack.Pop();

      skeleton[nextSkel++] = (char) nodes[i].type;

      switch(nodes[i].type)
         {
         case TREE_NODE_LEAF :
            leaves[nextLeaf++] = nodes[i].value;
            break;
         case TREE_NODE_TWO :
            stack.Push(nodes[i].child[1]);
         case TREE_NODE_ONE :
            stack.Push(nodes[i].child[0]);
            break;
         }

      if ((nextSkel & SWAP_MASK) == 0)
         {
         fwrite(&nextSkel, sizeof(nextSkel), 1, output);
         if (nextSkel) zip.Deflate(output, skeleton, nextSkel);
         nextSkel = 0;

         fwrite(&nextLeaf, sizeof(nextLeaf), 1, output);
         if (nextLeaf) zip.Deflate(output, leaves, nextLeaf * sizeof(double));
         nextLeaf = 0;
         }

      nodeCount++;
      }

   fwrite(&nextSkel, sizeof(nextSkel), 1, output);
   if (nextSkel) zip.Deflate(output, skeleton, nextSkel),

   fwrite(&nextLeaf, sizeof(nextLeaf), 1, output);
   if (nextLeaf) zip.Deflate(output, leaves, nextLeaf * sizeof(double));

   fseek(output, placeholder, SEEK_SET);
   fwrite(&nodeCount, sizeof(nodeCount), 1, output);
   }

bool BasicTree::FindNonZero(int node)
   {
   switch (nodes[node].type)
      {
      case TREE_NODE_ZERO :
         return false;
      case TREE_NODE_LEAF :
         return nodes[node].value != 0.0;
      case TREE_NODE_ONE :
         return FindNonZero(nodes[node].child[0]);
      case TREE_NODE_TWO :
         return FindNonZero(nodes[node].child[0]) || FindNonZero(nodes[node].child[1]);
      }
   return false;
   }


 
