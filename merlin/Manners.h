////////////////////////////////////////////////////////////////////// 
// merlin/Manners.h 
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
 
#ifndef __MANNERS_H__
#define __MANNERS_H__

// This function sets up signal handlers to ensure
// merlin terminates gracefully when something goes
// wrong

void SetupCrashHandlers();

typedef void (*signal_handler)(int);

void OutOfMemory(int);
void UserBreak(int);

#endif
 
