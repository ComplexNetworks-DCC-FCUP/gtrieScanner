/* -------------------------------------------------
      _       _     ___                            
 __ _| |_ _ _(_)___/ __| __ __ _ _ _  _ _  ___ _ _ 
/ _` |  _| '_| / -_)__ \/ _/ _` | ' \| ' \/ -_) '_|
\__, |\__|_| |_\___|___/\__\__,_|_||_|_||_\___|_|  
|___/                                          
    
gtrieScanner: quick discovery of network motifs
Released under Artistic License 2.0
(see README and LICENSE)

Pedro Ribeiro - CRACS & INESC-TEC, DCC/FCUP

----------------------------------------------------
Timers for time measurement

Last Update: 11/02/2012
---------------------------------------------------- */

#ifndef _TIMER_
#define _TIMER_

#define MAX_TIMERS 1000 // Maximum number of timers

#include "Common.h"

#include <sys/time.h>

class Timer {
 private:
  static struct timeval cstart[MAX_TIMERS]; // Start times
  static struct timeval   cend[MAX_TIMERS]; // End times

 public:
  static void     start(int n);             // Start the clock of a timer
  static void      stop(int n);             // Stop the clock of a timer
  static double elapsed(int n);             // Elapsed time of a timer
};

#endif
