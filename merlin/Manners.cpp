////////////////////////////////////////////////////////////////////// 
// merlin/Manners.cpp 
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
 
#include "Manners.h"

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

void SetupCrashHandlers()
   {
   signal(SIGINT, (signal_handler) UserBreak);
   signal(SIGSEGV, (signal_handler) OutOfMemory);
   signal(SIGABRT, (signal_handler) OutOfMemory);
   }

void OutOfMemory(int)
   {
   printf("\n\nMERLIN HAS CRASHED\n\n"
          "The operating system has decided to terminate Merlin,\n"
          "probably due to a memory access problem.\n\n"
          "If you set a memory usage limit (--megabytes:9999),\n"
          "Merlin may be able to recover gracefully from out of\n"
          "memory errors.\n\n"
          "You can reduce Merlin's memory requirements with the\n"
          "--swap option, which uses some of your free disk space\n"
          "as memory. The --smallSwap option uses even less memory\n"
          "and disk space.\n\n"
          "Alternatively, you could try a --singlepoint analysis,\n"
          "or use one of the options --zero, --one, --two, --three,\n"
          "to find an approximate solution if you are analysing a\n"
          "dense map.\n\n"
          "To help improve this program, please report bugs by\n"
          "e-mailing goncalo@umich.edu\n\n");

   exit(EXIT_FAILURE);
   }

void UserBreak(int)
   {
   printf("\n\nMERLIN STOPPED BY USER\n\n");

   exit(EXIT_FAILURE);
   }

 
