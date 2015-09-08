/* 
 * File:   meshit.hpp
 * Author: honnorat
 *
 * Created on 31 août 2015, 16:21
 */

#ifndef MESHIT_HPP
#define	MESHIT_HPP

#ifdef WIN32
   #if NGINTERFACE_EXPORTS || NGLIB_EXPORTS || nglib_EXPORTS
      #define DLL_HEADER   __declspec(dllexport)
   #else
      #define DLL_HEADER   __declspec(dllimport)
   #endif
#else
   #define DLL_HEADER 
#endif

#include "../meshing/msghandler.hpp"

#endif	/* MESHIT_HPP */
