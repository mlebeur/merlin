////////////////////////////////////////////////////////////////////// 
// merlin/TreeManager.cpp 
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
 
#include "TreeManager.h"
#include "TreeBasics.h"

void TreeManager::GetPosition(TreePosition & pos)
   {
   fgetpos(skeleton, &pos.skeleton);
   fgetpos(leaves, &pos.leaves);
   }

void TreeManager::UpdateFileSize()
   {
   if (skeletonSize == -1)
      return;

   int newSkeletonUsage = ftell(skeleton);
   int newLeafUsage = ftell(leaves);

   if (newSkeletonUsage < skeletonUsage || newLeafUsage < leafUsage)
      {
      skeletonSize = leafSize = skeletonUsage = leafUsage = -1;
      return;
      }

   skeletonUsage = newSkeletonUsage;
   leafUsage = newLeafUsage;

   if (skeletonUsage > skeletonSize)
      skeletonSize = skeletonUsage;

   if (leafUsage > leafSize)
      leafSize = leafUsage;
   }


void TreeManager::SetPosition(const TreePosition & pos)
   {
   fsetpos(skeleton, &pos.skeleton);
   fsetpos(leaves, &pos.leaves);

//   printf("/** Moved pointer to [%d, %d] **/\n", ftell(skeleton), ftell(leaves));
   }

bool TreeManager::OpenFiles()
   {
   if (skeleton == NULL)
      {
      skeleton = tmpfile();
      leaves = tmpfile();

      return true;
      }

   return false;
   }

void TreeManager::CloseFiles()
   {
   if (skeleton != NULL)
      {
      fclose(skeleton);
      fclose(leaves);
      }

   skeleton = leaves = NULL;
   }

void TreeManager::Free()
   {
   if (skeleton == NULL) return;

   fseek(skeleton, 0, SEEK_SET);
   fseek(leaves, 0, SEEK_SET);

   fgetpos(skeleton, &position.skeleton);
   fgetpos(leaves, &position.leaves);

   skeletonUsage = leafUsage = 0;

//   printf("/** Freed Swap File **/\n");
   }

double TreeManager::GetFileUsage()
   {
   if (skeletonUsage == -1 || leafUsage == -1)
      return -1.0;

   return (double) skeletonUsage + (double) leafUsage;
   }

double TreeManager::GetFileSize()
   {
   if (skeletonSize == -1 || leafSize == -1)
      return -1.0;

   return (double) skeletonSize + (double) leafSize;
   }

TreeManager::TreeManager()
   {
   skeleton = NULL;
   leaves = NULL;

   skeletonSize = leafSize = 0;
   skeletonUsage = leafUsage = 0;
   }
 
