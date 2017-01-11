////////////////////////////////////////////////////////////////////// 
// merlin/TreeManager.h 
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
 
#ifndef __TREEMANAGER_H__
#define __TREEMANAGER_H__

#include "MiniDeflate.h"

#define   FILEINDEX_t fpos_t

class BasicTree;

class TreePosition
   {
   public:
      FILEINDEX_t leaves;
      FILEINDEX_t skeleton;
   };


class TreeManager
   {
   public:
      // Constructor
      //
      TreeManager();

      // Swap file management functions
      //
      bool OpenFiles();       // Open file pointers
      void CloseFiles();         // Close files
      void Free();               // Discard contents and reuse swap file

      // These function returns the estimated file size / current usage
      // or a negative number when files get too big
      //
      double GetFileSize();
      double GetFileUsage();

      // This function should be called after the files are written to,
      // and is used to keep track of file sizes in a reasonably portable
      // manner ...
      void UpdateFileSize();

      // Routines for swapping tree in and out of memory
      FILEINDEX_t Pack(BasicTree & tree);
      FILEINDEX_t UnPack(BasicTree & tree);

      // Swap file pointers
      //
      FILE * skeleton;
      FILE * leaves;

      TreePosition  position;

      void GetPosition(TreePosition & position);
      void SetPosition(const TreePosition & position);

      void UpdatePosition()
         {
         GetPosition(position);
         }

   private:
      // These variables are used to track file size
      int leafSize;
      int skeletonSize;
      int leafUsage;
      int skeletonUsage;
   };

#endif


 
