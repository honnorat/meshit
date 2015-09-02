#include <cstring>
#include <iostream>
#include <iomanip>
#include <meshgen.hpp>
#include "dynamicmem.hpp"

using namespace std;

namespace netgen
{

  BaseDynamicMem * BaseDynamicMem::first = 0;
  BaseDynamicMem * BaseDynamicMem::last = 0;


  BaseDynamicMem :: BaseDynamicMem ()
  {
    prev = last;
    next = 0;

    if (last) last->next = this;
    last = this;
    if (!first) first = this;

    size = 0;
    ptr = 0;
    name = 0;
  }
 
  BaseDynamicMem :: ~BaseDynamicMem ()
  {
    Free();

    if (next) next->prev = prev;
    else last = prev;
    if (prev) prev->next = next;
    else first = next;

    delete [] name;
  }

  void BaseDynamicMem :: SetName (const char * aname)
  {
    delete [] name;
    name = NULL;
    if (aname)
      {
	name = new char[strlen(aname)+1];
	strcpy (name, aname);
      }
  }


  void BaseDynamicMem :: Alloc (size_t s)
  {
    size = s;
    ptr = new char[s];

    if (!ptr)
      {
	std::cerr << "BaseynamicMem, cannot allocate " << s << " bytes" << std::endl;
	Print ();
	throw ("BaseDynamicMem::Alloc: out of memory");
      }
    // ptr = (char*)malloc (s);
    // ptr = (char*) _mm_malloc (s,16);
  }

  void BaseDynamicMem :: ReAlloc (size_t s)
  {
    if (size == s) return;

    char * old = ptr;
    ptr = new char[s];

    if (!ptr)
      {
	std::cerr << "BaseynamicMem, cannot Reallocate " << s << " bytes" << std::endl;
	Print ();
	throw ("BaseDynamicMem::Alloc: out of memory");
      }


    // ptr = (char*)malloc(s);
    // ptr = (char*) _mm_malloc (s,16);
    memmove (ptr, old, (s < size) ? s : size);
    delete [] old;
    // free (old);
    // _mm_free (old);
    size = s;
  }

  void BaseDynamicMem :: Free ()
  {
    delete [] ptr;
    // free (ptr);
    // _mm_free (ptr);
    ptr = 0;
  }

  void BaseDynamicMem :: Swap (BaseDynamicMem & m2)
  {
    size_t hi;
    char * cp;
    hi = size; size  = m2.size; m2.size = hi;
    cp = ptr; ptr = m2.ptr; m2.ptr = cp;
    cp = name; name = m2.name; m2.name = cp;
  }


  void BaseDynamicMem :: Print ()
  {
    std::cout << "****************** Dynamic Mem Report ****************" << std::endl;
    BaseDynamicMem * p = first;
    size_t mem = 0;
    int cnt = 0;
    while (p)
      {
	mem += p->size;
	cnt++;

	std::cout << std::setw(10) << p->size << " Bytes";
	std::cout << ", addr = " << (void*)p->ptr;
	if (p->name)
	  std::cout << " in block " << p->name;
	std::cout << std::endl;

	p = p->next;
      }

    if (mem > 100000000)
      std::cout << "memory in dynamic memory: " << mem/1048576 << " MB" << std::endl;
    else if (mem > 100000)
      std::cout << "memory in dynamic memory: " << mem/1024 << " kB" << std::endl;
    else
      std::cout << "memory in dynamic memory: " << mem << " Bytes" << std::endl;
    std::cout << "number of blocks:         " << cnt << std::endl;
    //  std::cout << "******************************************************" <<std::endl;
  }


#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:1684)
#endif

  void BaseDynamicMem :: GetUsed (int nr, char * ch)
  {
    BaseDynamicMem * p = first;

    for (int i = 0; i < nr; i++)
      ch[i] = '0';

    while (p)
      {
        long unsigned hptr = (long unsigned) (p->ptr);
	// uintptr_t hptr = reinterpret_cast<uintptr_t>(p->ptr); //??

	hptr /= (1024*1024);
	hptr /= (4096/nr);

	size_t blocks = p->size / (1024*1024);
	blocks /= (4096/nr);
	
	// std::cout << "ptr = " << (void*)(p->ptr) << ", size = " << p->size << ", hptr = " << hptr << " blocks = " << blocks <<std::endl;

	for (size_t i = 0; i <= blocks; i++)
	  ch[hptr+i] = '1';

	p = p->next;
      }
    
    {

      /*
    BaseMoveableMem * pm = BaseMoveableMem::first;
    while (pm)
      {
        long unsigned hptr = (long unsigned) pm->ptr;
        // uintptr_t hptr = reinterpret_cast<uintptr_t>(pm->ptr);

	hptr /= (1024*1024);
	hptr /= (4096/nr);

	size_t blocks = pm->size / (1024*1024);
	blocks /= (4096/nr);
	
	// std::cout << "moveable, ptr = " << (void*)(pm->ptr) << ", size = " << pm->size << ", hptr = " << hptr << " blocks = " << blocks <<std::endl;

	for (size_t i = 0; i <= blocks; i++)
	  ch[hptr+i] = '1';

	pm = pm->next;
      }
      */
    }



  }

#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

}
