#ifndef FILE_MSGHANDLER
#define FILE_MSGHANDLER

/**************************************************************************/
/* File:   msghandler.hh                                                  */
/* Author: Johannes Gerstmayr                                             */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/

#include "../general/mystring.hpp"

#define FATAL_LOG_LEVEL 500
#define ERROR_LOG_LEVEL 400
#define WARN_LOG_LEVEL 300
#define INFO_LOG_LEVEL 200
#define DEBUG_LOG_LEVEL 100
#define TRACE_LOG_LEVEL 0

#if defined(__GNUC__) || (defined(__ICC) && (__ICC >= 600))
#define MESHIT_CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__FUNCSIG__)
#define MESHIT_CURRENT_FUNCTION __FUNCSIG__
#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600))
#define MESHIT_CURRENT_FUNCTION __FUNCTION__
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
#define MESHIT_CURRENT_FUNCTION __func__
#else
#define MESHIT_CURRENT_FUNCTION "(unknown)"
#endif

#define LOG_POSITION(file, line, function)      \
    do {                                        \
        std::cerr << " (at " << file            \
                  << ":"     << line            \
                  << " in '" << function        \
                  << "'";      \
    }                                           \
    while(0);

#define LOG_MACRO_BODY(logEvent, logLevel, retline)                         \
    do {                                                                    \
        if ( logLevel >= meshit::GetLogLevel() ) {                          \
            if ( logLevel >= FATAL_LOG_LEVEL ) {                            \
                std::cerr << "FATAL: " << logEvent;                         \
                LOG_POSITION(__FILE__, __LINE__, MESHIT_CURRENT_FUNCTION);  \
                std::cerr << std::endl;                                     \
            } else if ( logLevel >= ERROR_LOG_LEVEL ) {                     \
                std::cerr << "ERROR: " << logEvent;                         \
                LOG_POSITION(__FILE__, __LINE__, MESHIT_CURRENT_FUNCTION);  \
                std::cerr << std::endl;                                     \
            } else if ( logLevel >= WARN_LOG_LEVEL ) {                      \
                std::cout << "WARNING: " << logEvent << std::endl;          \
            } else {                                                        \
                std::cout << logEvent;                                      \
                if ( retline ) std::cout << std::endl;                      \
            }                                                               \
        }                                                                   \
    } while (0)    

#define LOG_FATAL(logEvent)     LOG_MACRO_BODY (logEvent, FATAL_LOG_LEVEL, true)
#define LOG_ERROR(logEvent)     LOG_MACRO_BODY (logEvent, ERROR_LOG_LEVEL, true)
#define LOG_WARNING(logEvent)   LOG_MACRO_BODY (logEvent, WARN_LOG_LEVEL, true)
#define LOG_INFO(logEvent)      LOG_MACRO_BODY (logEvent, INFO_LOG_LEVEL, true)
#define LOG_DEBUG(logEvent)     LOG_MACRO_BODY (logEvent, DEBUG_LOG_LEVEL, true)

#define LOG_INFO_CONT(logEvent)      LOG_MACRO_BODY (logEvent, INFO_LOG_LEVEL, false)
#define LOG_DEBUG_CONT(logEvent)     LOG_MACRO_BODY (logEvent, DEBUG_LOG_LEVEL, false)

namespace meshit {

    void PrintDot(char ch = '.');
    extern int printmessage_importance;
    extern int printdots;

    int GetLogLevel();
    void SetLogLevel(int logLevel);

    //Message Pipeline:

    //importance: importance of message: 1=very important, 3=middle, 5=low, 7=unimportant
    //    void PrintMessage(int importance, const std::stringstream& s);
    //
    void PrintMessage(int importance,
            const MyStr& s1, const MyStr& s2 = MyStr());
    void PrintMessage(int importance,
            const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4 = MyStr());
    void PrintMessage(int importance,
            const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4,
            const MyStr& s5, const MyStr& s6 = MyStr(), const MyStr& s7 = MyStr(), const MyStr& s8 = MyStr());

    // CR without line-feed
    void PrintWarning(const MyStr& s1, const MyStr& s2 = "", const MyStr& s3 = "", const MyStr& s4 = "",
            const MyStr& s5 = "", const MyStr& s6 = "", const MyStr& s7 = "", const MyStr& s8 = "");
    void PrintError(const MyStr& s1, const MyStr& s2 = "", const MyStr& s3 = "", const MyStr& s4 = "",
            const MyStr& s5 = "", const MyStr& s6 = "", const MyStr& s7 = "", const MyStr& s8 = "");
    void PrintSysError(const MyStr& s1, const MyStr& s2 = "", const MyStr& s3 = "", const MyStr& s4 = "",
            const MyStr& s5 = "", const MyStr& s6 = "", const MyStr& s7 = "", const MyStr& s8 = "");
    void SetStatMsg(const MyStr& s);

    void PushStatus(const MyStr& s);
    void PopStatus();
    void GetStatus(MyStr & s, double & percentage);
}

#endif

