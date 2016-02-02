#include "../meshit.hpp"
#include "msghandler.hpp"
#include "../general/array.hpp"
#include "global.hpp"

namespace meshit {

    bool printdots = true;

    static int _meshit_logLevel = MESHIT_INFO_LOG_LEVEL;

    int GetLogLevel()
    {
        return _meshit_logLevel;
    }

    void SetLogLevel(int logLevel)
    {
        _meshit_logLevel = logLevel;
    }

    void PrintDot(char ch)
    {
        //the dots for progression of program
        if (printdots) {
            MESHIT_LOG_INFO_CONT(ch << '\0' << std::flush);
        }
    }
}
