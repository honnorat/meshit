#include <meshgen.hpp>
#include "global.hpp"


namespace netgen
{
  //  Flags parameters;

  int silentflag = 0;
  int testmode = 0;

  volatile multithreadt multithread;

  Array<int> tets_in_qualclass;

  int h_argc = 0;
  char ** h_argv = NULL;

  multithreadt :: multithreadt()
  {
    pause =0;
    testmode = 0;
    redraw = 0;
    drawing = 0;
    terminate = 0;
    running = 0;
    percent = 0;
    task = "";
  }

  DebugParameters debugparam;
  bool verbose = 0;

  int timestamp = 0;
  int GetTimeStamp() 
  { 
    return timestamp; 
  }

  int NextTimeStamp()
  {
    timestamp++;
    return timestamp;
  }
}
