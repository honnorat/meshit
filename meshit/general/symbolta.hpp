#ifndef FILE_SYMBOLTA_HPP
#define FILE_SYMBOLTA_HPP

#include "array.hpp"
#include "template.hpp"

/**************************************************************************/
/* File:   symbolta.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

namespace meshit {

/**
   Base class for the generic SYMBOLTABLE.
   An array of identifiers is maintained.
*/
    class BASE_SYMBOLTABLE
    {
     protected:
        /// identifiers
        Array<char*> names;

     public:
        /// Constructor
        BASE_SYMBOLTABLE() { }

        ~BASE_SYMBOLTABLE();

        void DelNames();
        /// Index of symbol name, returns 0 if not used.
        int Index(const char* name) const;
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
        Array<T> data;

     public:
        /// Creates a symboltable
        SYMBOLTABLE() { }

        /// Returns size of symboltable
        INDEX Size() const
        {
            return data.size();
        }

        /// Returns element, error if not used
        inline const T& Get(const char* name) const;
        /// Returns i-th element
        inline const T& Get(int i) const;
        /// Associates el to the string name, overrides if name is used
        inline void Set(const char* name, const T& el);
        /// Checks whether name is used
        inline bool Used(const char* name) const;
        /// Deletes symboltable
        inline void DeleteAll();

        inline T& operator[](int i) { return data[i]; }

        inline const T& operator[](int i) const { return data[i]; }

     private:
        /// Prevents from copying symboltable by pointer assignment
        SYMBOLTABLE<T>& operator=(SYMBOLTABLE<T>&) { return *this; }
    };

    template<class T>
    inline const T& SYMBOLTABLE<T>::Get(const char* name) const
    {
        int i = Index(name);
        return data[std::max(0, i - 1)];
    }

    template<class T>
    inline const T& SYMBOLTABLE<T>::Get(int i) const
    {
        return data[i - 1];
    }

    template<class T>
    inline void SYMBOLTABLE<T>::Set(const char* name, const T& el)
    {
        int i = Index(name);
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
    inline bool SYMBOLTABLE<T>::Used(const char* name) const
    {
        return static_cast<bool>(Index(name));
    }

    template<class T>
    inline void SYMBOLTABLE<T>::DeleteAll()
    {
        DelNames();
        data.DeleteAll();
    }

}  // namespace meshit

#endif
