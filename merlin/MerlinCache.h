////////////////////////////////////////////////////////////////////// 
// merlin/MerlinCache.h 
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
 
#ifndef __MERLINCACHE_H__
#define __MERLINCACHE_H__

#include "Mantra.h"
#include "TreeBasics.h"

class MerlinCache
   {
   public:
      MerlinCache();
      ~MerlinCache();

      void OpenCache(Mantra & m);
      void CloseCache();

      bool RetrieveFromCache(Mantra & m, BasicTree & tree, String & marker);
      void SaveToCache(Mantra & m, BasicTree & tree, String & marker);

      static String directory;

      static void FreeBuffers() { BasicTree::FreeBuffers(); }

   private:
      FILE * indexFile;
      FILE * dataFile;

      StringIntHash markers;
      String tag;

      bool isActive;

      bool IndexMatches(Mantra & m);
   };


#endif
 
