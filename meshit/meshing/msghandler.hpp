#ifndef FILE_MSGHANDLER
#define FILE_MSGHANDLER

/**************************************************************************/
/* File:   msghandler.hh                                                  */
/* Author: Johannes Gerstmayr                                             */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/

#include "../general/mystring.hpp"

namespace meshit {

    void PrintDot(char ch = '.');
    extern int printmessage_importance;
    extern int printdots;

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

