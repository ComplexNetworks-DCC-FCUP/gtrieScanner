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

#include "Timer.h"

// Static variables
struct timeval Timer::cstart[MAX_TIMERS];
struct timeval Timer::cend[MAX_TIMERS];

// Stop the clock of a timer
void Timer::start(int n) {
  if (n>=0 && n<MAX_TIMERS)
    gettimeofday(&cstart[n], NULL);
}

// Elapsed time of a timer
void Timer::stop(int n) {
  if (n>=0 && n<MAX_TIMERS)
    gettimeofday(&cend[n], NULL);
}

double Timer::elapsed(int n) {
  if (n>=0 && n<MAX_TIMERS)
    return (cend[n].tv_sec - cstart[n].tv_sec) + (double)(cend[n].tv_usec - cstart[n].tv_usec)/1000000.0;
  else
    return 0;
}
