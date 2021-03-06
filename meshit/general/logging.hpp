#ifndef MESHIT_LOGGING_HPP
#define MESHIT_LOGGING_HPP
/**
 * meshit - a 2d mesh generator
 *
 * Copyright © 1995-2015 Joachim Schoeberl <joachim.schoeberl@tuwien.ac.at>
 * Copyright © 2015-2016 Marc Honnorat <marc.honnorat@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library in the file LICENSE.LGPL; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 */

#include <iostream>

#define MESHIT_FATAL_LOG_LEVEL 500
#define MESHIT_ERROR_LOG_LEVEL 400
#define MESHIT_WARN_LOG_LEVEL 300
#define MESHIT_INFO_LOG_LEVEL 200
#define MESHIT_DEBUG_LOG_LEVEL 100
#define MESHIT_TRACE_LOG_LEVEL 0

#if defined(__GNUC__) || (defined(__ICC) && (__ICC >= 600))
#define _MESHIT_CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__FUNCSIG__)
#define _MESHIT_CURRENT_FUNCTION __FUNCSIG__
#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600))
#define _MESHIT_CURRENT_FUNCTION __FUNCTION__
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
#define _MESHIT_CURRENT_FUNCTION __func__
#else
#define _MESHIT_CURRENT_FUNCTION "(unknown)"
#endif

#define _MESHIT_LOG_POSITION(file, line, function)                                                           \
    do {                                                                                                     \
        std::cerr << " (at " << file << ":" << line << " in '" << function << "'";                           \
    } while (0);

#define _MESHIT_LOG_MACRO_BODY(logEvent, logLevel, retline)                                                  \
    do {                                                                                                     \
        if (logLevel >= meshit::GetLogLevel()) {                                                             \
            if (logLevel >= MESHIT_FATAL_LOG_LEVEL) {                                                        \
                std::cerr << "FATAL: " << logEvent;                                                          \
                _MESHIT_LOG_POSITION(__FILE__, __LINE__, _MESHIT_CURRENT_FUNCTION);                          \
                std::cerr << std::endl;                                                                      \
            } else if (logLevel >= MESHIT_ERROR_LOG_LEVEL) {                                                 \
                if (retline) std::cerr << "ERROR: ";                                                         \
                std::cerr << logEvent;                                                                       \
                if (retline) {                                                                               \
                    _MESHIT_LOG_POSITION(__FILE__, __LINE__, _MESHIT_CURRENT_FUNCTION);                      \
                    std::cerr << std::endl;                                                                  \
                }                                                                                            \
            } else if (logLevel >= MESHIT_WARN_LOG_LEVEL) {                                                  \
                std::cout << "WARNING: " << logEvent << std::endl;                                           \
            } else {                                                                                         \
                std::cout << logEvent;                                                                       \
                if (retline) std::cout << std::endl;                                                         \
            }                                                                                                \
        }                                                                                                    \
    } while (0)

#define MESHIT_LOG_FATAL(logEvent) _MESHIT_LOG_MACRO_BODY(logEvent, MESHIT_FATAL_LOG_LEVEL, true)
#define MESHIT_LOG_ERROR(logEvent) _MESHIT_LOG_MACRO_BODY(logEvent, MESHIT_ERROR_LOG_LEVEL, true)
#define MESHIT_LOG_WARNING(logEvent) _MESHIT_LOG_MACRO_BODY(logEvent, MESHIT_WARN_LOG_LEVEL, true)
#define MESHIT_LOG_INFO(logEvent) _MESHIT_LOG_MACRO_BODY(logEvent, MESHIT_INFO_LOG_LEVEL, true)
#define MESHIT_LOG_DEBUG(logEvent) _MESHIT_LOG_MACRO_BODY(logEvent, MESHIT_DEBUG_LOG_LEVEL, true)

#define MESHIT_LOG_ERROR_CONT(logEvent) _MESHIT_LOG_MACRO_BODY(logEvent, MESHIT_ERROR_LOG_LEVEL, false)
#define MESHIT_LOG_INFO_CONT(logEvent) _MESHIT_LOG_MACRO_BODY(logEvent, MESHIT_INFO_LOG_LEVEL, false)
#define MESHIT_LOG_DEBUG_CONT(logEvent) _MESHIT_LOG_MACRO_BODY(logEvent, MESHIT_DEBUG_LOG_LEVEL, false)

namespace meshit {

int GetLogLevel();
void SetLogLevel(int logLevel);

}

#endif
