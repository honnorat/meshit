//File for handling warnings, errors, messages
#include <meshit.hpp>
#include "msghandler.hpp"
#include "../general/array.hpp"
#include "global.hpp"

namespace meshit {

    int printmessage_importance = 5;
    int printwarnings = 1;
    int printerrors = 1;
    int printdots = 1;
    int printfnstart = 0;

    static int _meshit_logLevel = INFO_LOG_LEVEL;

    int GetLogLevel()
    {
        return _meshit_logLevel;
    }

    void SetLogLevel(int logLevel)
    {
        _meshit_logLevel = logLevel;
    }

    std::ostream & MESHIT_COUT = std::cerr;

    //the dots for progression of program

    void PrintDot(char ch)
    {
        if (printdots) {
            MESHIT_COUT << ch << '\0' << std::flush;
        }
    }

    void PrintMessage(int importance,
            const MyStr& s1, const MyStr& s2)
    {
        if (importance <= printmessage_importance) {
            MESHIT_COUT << " " << s1 << s2 << std::endl;
        }
    }

    void PrintMessage(int importance,
            const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4)
    {
        if (importance <= printmessage_importance) {
            MESHIT_COUT << " " << s1 << s2 << s3 << s4 << std::endl;
        }
    }

    void PrintMessage(int importance,
            const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4,
            const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
    {
        if (importance <= printmessage_importance) {
            MESHIT_COUT << " " << s1 << s2 << s3 << s4 << s5 << s6 << s7 << s8 << std::endl;
        }
    }

    void PrintWarning(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4,
            const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
    {
        if (printwarnings)
            MESHIT_COUT << " WARNING: " << s1 << s2 << s3 << s4 << s5 << s6 << s7 << s8 << std::endl;
    }

    void PrintError(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4,
            const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
    {
        if (printerrors)
            MESHIT_COUT << " ERROR: " << s1 << s2 << s3 << s4 << s5 << s6 << s7 << s8 << std::endl;
    }

    void PrintSysError(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4,
            const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
    {
        if (printerrors)
            MESHIT_COUT << " SYSTEM ERROR: " << s1 << s2 << s3 << s4 << s5 << s6 << s7 << s8 << std::endl;
    }

    static Array<MyStr*> msgstatus_stack(0);
    static Array<double> threadpercent_stack(0);
    static MyStr msgstatus = "";

    void PushStatus(const MyStr& s)
    {
        msgstatus_stack.push_back(new MyStr(s));
        SetStatMsg(s);
        threadpercent_stack.push_back(0);
    }

    void PopStatus()
    {
        if (msgstatus_stack.size()) {
            if (msgstatus_stack.size() > 1)
                // SetStatMsg (*msgstatus_stack.Last());
                SetStatMsg(*msgstatus_stack[msgstatus_stack.size() - 2]);
            else
                SetStatMsg("");
            delete msgstatus_stack.Last();
            msgstatus_stack.DeleteLast();
        }
        else {
            PrintSysError("PopStatus failed");
        }
    }

    void SetStatMsg(const MyStr& s)
    {
        msgstatus = s;
    }
}
