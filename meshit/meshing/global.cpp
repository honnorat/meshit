#include "global.hpp"

namespace meshit {

    DebugParameters debugparam;

    int timestamp = 0;

    int NextTimeStamp()
    {
        timestamp++;
        return timestamp;
    }
}  // namespace meshit
