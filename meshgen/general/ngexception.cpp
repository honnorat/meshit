/**************************************************************************/
/* File:   ngexception.cpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   16. Jan. 02                                                    */
/**************************************************************************/

//#include <myadt.hpp>
#include <meshgen.hpp>
#include <string>
#include "ngexception.hpp"

namespace netgen
{
  //using namespace netgen;



  NgException :: NgException (const std::string & s) 
    : what(s)
  {
    ; 
  }


  NgException :: ~NgException () 
  {
    ;
  }

  /// append string to description
  void NgException :: Append (const std::string & s)
  { 
    what += s; 
  }

}
