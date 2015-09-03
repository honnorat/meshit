/**************************************************************************/
/* File:   symbolta.cc                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   Abstract data type Symbol Table
*/

#include <meshit.hpp>
#include "symbolta.hpp"

#ifndef FILE_SYMBOLTABLECC
#define FILE_SYMBOLTABLECC
// necessary for SGI ????


namespace meshit
{
  //using namespace netgen;

  BASE_SYMBOLTABLE :: BASE_SYMBOLTABLE ()
  {
    ;
  }


  BASE_SYMBOLTABLE :: ~BASE_SYMBOLTABLE()
  {
    DelNames();
  }


  void BASE_SYMBOLTABLE :: DelNames()
  {
    for (int i = 0; i < names.size(); i++)
      delete [] names[i];
    names.resize (0);
  }

  int BASE_SYMBOLTABLE :: Index (const char * name) const
  {
    if (!name) return 0;
    for (int i = 0; i < names.size(); i++)
      if (strcmp (names[i], name) == 0) return i+1;
    return 0;
  }
}

#endif
