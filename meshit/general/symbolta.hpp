#ifndef FILE_SYMBOLTA_HPP
#define FILE_SYMBOLTA_HPP

#include <vector>
#include <cstring>

#include "index.hpp"

/**************************************************************************/
/* File:   symbolta.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

namespace meshit
{
/**
   Base class for the generic SYMBOLTABLE.
   An array of identifiers is maintained.
*/
    class BASE_SYMBOLTABLE
    {
     protected:
        /// identifiers
        std::vector<char*> names;

     public:
        /// Constructor
        BASE_SYMBOLTABLE() { }

        ~BASE_SYMBOLTABLE();

        void DelNames();

        /// Index of symbol name, returns 0 if not used.
        size_t Index(const char* name) const;
    };


/** 
    Abstract data type Symbol Table.
   
    To a string an value of the generic type T is associated.
    The string is not copied into the symbol table class!
*/
    template<class T>
    class SYMBOLTABLE : public BASE_SYMBOLTABLE
    {
     private:
        /// Associated data
        std::vector<T> data;

     public:
        /// Creates a symboltable
        SYMBOLTABLE() { }

        /// Returns size of symboltable
        size_t Size() const
        {
            return data.size();
        }

        /// Returns element, error if not used
        const T& Get(const char* name) const
        {
            size_t i = Index(name);
            return data[std::max(0UL, i - 1)];
        }

        /// Associates el to the string name, overrides if name is used
        inline void Set(const char* name, const T& el);

        /// Checks whether name is used
        bool Used(const char* name) const
        {
            return (Index(name) > 0);
        }

        /// Deletes symboltable
        inline void DeleteAll();

        inline T& operator[](size_t i) { return data[i]; }

        inline const T& operator[](size_t i) const { return data[i]; }

     private:
        /// Prevents from copying symboltable by pointer assignment
        SYMBOLTABLE<T>& operator=(SYMBOLTABLE<T>&) { return *this; }
    };

    template<class T>
    inline void SYMBOLTABLE<T>::Set(const char* name, const T& el)
    {
        size_t i = Index(name);
        if (i) {
            data[i - 1] = el;
        } else {
            data.push_back(el);
            char* hname = new char[strlen(name) + 1];
            strcpy(hname, name);
            names.push_back(hname);
        }
    }

    template<class T>
    inline void SYMBOLTABLE<T>::DeleteAll()
    {
        DelNames();
        data.resize(0);
    }

}  // namespace meshit

#endif
