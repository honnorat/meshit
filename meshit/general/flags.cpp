/**************************************************************************/
/* File:   flags.cc                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                    */
/**************************************************************************/

/* 
   Datatype Flags
*/

#include <meshit.hpp>
#include <sstream>
#include <fstream>
#include "flags.hpp"

namespace meshit
{
  //using namespace netgen;

  Flags :: Flags ()
  {
    ;
  }
  
  Flags :: ~Flags ()
  {
    DeleteFlags ();
  }
  
  void Flags :: DeleteFlags ()
  {
    for (int i = 0; i < strflags.Size(); i++)
      delete [] strflags[i];
    for (int i = 0; i < numlistflags.Size(); i++)
      delete numlistflags[i];
    strflags.DeleteAll();
    numflags.DeleteAll();
    defflags.DeleteAll();
    strlistflags.DeleteAll();
    numlistflags.DeleteAll();
  }
  
  void Flags :: SetFlag (const char * name, const char * val)
  {
    char * hval = new char[strlen (val) + 1];
    strcpy (hval, val);
    strflags.Set (name, hval);
  }
  
  void Flags :: SetFlag (const char * name, double val)
  {
    numflags.Set (name, val);
  }
  
  void Flags :: SetFlag (const char * name)
  {
    defflags.Set (name, 1);
  }


  void Flags :: SetFlag (const char * name, const Array<char*> & val)
  {
    Array<char*> * strarray = new Array<char*>;
    for (int i = 1; i <= val.size(); i++)
      {
	strarray->push_back (new char[strlen(val.Get(i))+1]);
	strcpy (strarray->Last(), val.Get(i));
      }
    strlistflags.Set (name, strarray);
  }

  void Flags :: SetFlag (const char * name, const Array<double> & val)
  {
    Array<double> * numarray = new Array<double>;
    for (int i = 1; i <= val.size(); i++)
      numarray->push_back (val.Get(i));
    numlistflags.Set (name, numarray);
  }




  
  const char * 
  Flags :: GetStringFlag (const char * name, const char * def) const
  {
    if (strflags.Used (name))
      return strflags.Get(name);
    else
      return def;
  }

  double Flags :: GetNumFlag (const char * name, double def) const
  {
    if (numflags.Used (name))
      return numflags.Get(name);
    else
      return def;
  }
  
  const double * Flags :: GetNumFlagPtr (const char * name) const
  {
    if (numflags.Used (name))
      return & ((SYMBOLTABLE<double>&)numflags).Elem(name);
    else
      return NULL;
  }
  
  double * Flags :: GetNumFlagPtr (const char * name) 
  {
    if (numflags.Used (name))
      return & ((SYMBOLTABLE<double>&)numflags).Elem(name);
    else
      return NULL;
  }
  
  bool Flags :: GetDefineFlag (const char * name) const
  {
    return defflags.Used (name);
  }


  const Array<char*> & 
  Flags :: GetStringListFlag (const char * name) const
  {
    if (strlistflags.Used (name))
      return *strlistflags.Get(name);
    else
      {
	static Array<char*> dummy_array(0);
	return dummy_array;
      }
  }

  const Array<double> & 
  Flags ::GetNumListFlag (const char * name) const
  {
    if (numlistflags.Used (name))
      return *numlistflags.Get(name);
    else
      {
	static Array<double> dummy_array(0);
	return dummy_array;
      }
  }


  bool Flags :: StringFlagDefined (const char * name) const
  {
    return strflags.Used (name);
  }

  bool Flags :: NumFlagDefined (const char * name) const
  {
    return numflags.Used (name);
  }

  bool Flags :: StringListFlagDefined (const char * name) const
  {
    return strlistflags.Used (name);
  }

  bool Flags :: NumListFlagDefined (const char * name) const
  {
    return numlistflags.Used (name);
  }


  void Flags :: SaveFlags (const char * filename) const 
  {
    int i;
    std::ofstream outfile (filename);
  
    for (i = 1; i <= strflags.Size(); i++)
      outfile << strflags.GetName(i) << " = " << strflags.Get(i) << std::endl;
    for (i = 1; i <= numflags.Size(); i++)
      outfile << numflags.GetName(i) << " = " << numflags.Get(i) << std::endl;
    for (i = 1; i <= defflags.Size(); i++)
      outfile << defflags.GetName(i) << std::endl;
  }
 


  void Flags :: PrintFlags (std::ostream & ost) const 
  {
    int i;
  
    for (i = 1; i <= strflags.Size(); i++)
      ost << strflags.GetName(i) << " = " << strflags.Get(i) << std::endl;
    for (i = 1; i <= numflags.Size(); i++)
      ost << numflags.GetName(i) << " = " << numflags.Get(i) << std::endl;
    for (i = 1; i <= defflags.Size(); i++)
      ost << defflags.GetName(i) << std::endl;
  }
 

  void Flags :: LoadFlags (const char * filename) 
  {
    char name[100], str[100];
    char ch;
    double val;
    std::ifstream infile(filename);

    //  (*logout) << "Load flags from " << filename <<std::endl <<std::endl;
    while (infile.good())
      {
	infile >> name;
	if (strlen (name) == 0) break;

	if (name[0] == '/' && name[1] == '/')
	  {
	    //	  (*logout) << "comment: ";
	    ch = 0;
	    while (ch != '\n' && infile.good())
	      {
		ch = infile.get();
		//	      (*logout) << ch;
	      }
	    continue;
	  }

	//      (*logout)  << name;
	ch = 0;
	infile >> ch;
	if (ch != '=')
	  {
	    //	  (*logout) <<std::endl;
	    infile.putback (ch);
	    SetFlag (name);
	  }
	else
	  {
	    infile >> val;
	    if (!infile.good())
	      {
		infile.clear();
		infile >> str;
		SetFlag (name, str);
		//	      (*logout) << " = " << str <<std::endl;
	      }
	    else
	      {
		SetFlag (name, val);
		//	      (*logout) << " = " << val <<std::endl;
	      }
	  }
      }
    //  (*logout) <<std::endl;
  }


  void Flags :: SetCommandLineFlag (const char * st)
  {
    //  std::cout << "clflag = " << st <<std::endl;
    std::istringstream inst( (char *)st);
    // istrstream defined with char *  (not const char *  ?????)

    char name[100];
    double val;


    if (st[0] != '-')
      {
	std::cerr << "flag must start with '-'" << std::endl;
	return;
      }
  
    const char * pos = strchr (st, '=');
  
    if (!pos)
      {
	//      (std::cout) << "Add def flag: " << st+1 <<std::endl;
	SetFlag (st+1);
      }
    else
      {
	//      std::cout << "pos = " << pos <<std::endl;

	strncpy (name, st+1, (pos-st)-1);
	name[pos-st-1] = 0;

	//      std::cout << "name = " << name <<std::endl;

	pos++;
	char * endptr = NULL;

	val = strtod (pos, &endptr);

	//      std::cout << "val = " << val <<std::endl;

	if (endptr == pos)
	  {
	    //	  (std::cout) << "Add String Flag: " << name << " = " << pos <<std::endl;
	    SetFlag (name, pos);
	  }
	else
	  {
	    //	  (std::cout) << "Add Num Flag: " << name << " = " << val <<std::endl;
	    SetFlag (name, val);
	  }
      }


    /*
      inst >> name;
      (*mystd::cout) << "name = " << name <<std::endl;

      ch = 0;
      inst >> ch;
      if (ch != '=')
      {
      SetFlag (name);
      }
      else
      {
      inst >> val;
      if (!inst.good())
      {
      inst.clear();
      inst >> str;
      SetFlag (name, str);
      (*mystd::cout) << "str = " << str <<std::endl;
      }
      else
      {
      SetFlag (name, val);
      (*mystd::cout) << "val = " << val <<std::endl;
      }
      }
    */
  }
}
