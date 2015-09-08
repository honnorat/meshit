#ifndef FILE_GLOBAL
#define FILE_GLOBAL


/**************************************************************************/
/* File:   global.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

/*
  global functions and variables
*/
#include <string>
#include "../general/array.hpp"
#include "meshtype.hpp"

namespace meshit
{
  extern double GetTime ();
  extern void ResetTime ();

  extern int testmode;

  extern Array<int> tets_in_qualclass;

  extern DebugParameters debugparam;
}

#endif
