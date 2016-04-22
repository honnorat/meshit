#include "logging.hpp"

namespace meshit {

static int _meshit_logLevel = MESHIT_INFO_LOG_LEVEL;

int GetLogLevel()
{
    return _meshit_logLevel;
}

void SetLogLevel(int logLevel)
{
    _meshit_logLevel = logLevel;
}

}  // namespace meshit
