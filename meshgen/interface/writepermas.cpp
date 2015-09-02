//
// Write Permas file
// for Intes GmbH, Stuttgart
//

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

#include <string>

using namespace std;

namespace netgen
{
#include "writeuser.hpp"
    // Forward declarations (don't know, where to define them, sorry)
    int addComponent(string &strComp, string &strSitu, ofstream &out);


    // This should be the new function to export a PERMAS file
    void WritePermasFormat (const Mesh &mesh, const string &filename, 
                            string &strComp, string &strSitu) 
    {
        ofstream outfile (filename.c_str());
        addComponent(strComp, strSitu, outfile);
        WritePermasFormat ( mesh, filename);
    }

    void WritePermasFormat (const Mesh &mesh, const string &filename)
    {
        string strComp, strSitu;
        ofstream outfile (filename.c_str());
        
        outfile.precision(8);

        strSitu = strComp = "";
        if (addComponent(strComp, strSitu, outfile) == 1) {
            printf("Error while exporting PERMAS dat!\n");
            return;
        }
    
        int np = mesh.GetNP();
        int ne = mesh.GetNE();
        int nse = mesh.GetNSE();
        int i, j, k;
    
        if (ne == 0)
        {
            // pure surface mesh
            std::cout << "\nWrite Permas Surface Mesh" <<std::endl;
            
            int elnr = 0;
            for (j = 1; j <= 2; j++)
            {
                int nelp(0);
                switch (j)
                {
                case 1:
                    nelp = 3;
                    outfile << "$ELEMENT TYPE = TRIA3  ESET = ALLQUAD" <<std::endl;		  
                    break;
                case 2:
                    nelp = 4;
                    outfile << "$ELEMENT TYPE = QUAD4  ESET = ALLQUAD" <<std::endl;		  
                    break;
                }
                
                for (i = 1; i <= nse; i++)
                {
                    const Element2d & el = mesh.SurfaceElement(i);
                    if (el.GetNP() != nelp)
                        continue;
                    
                    elnr++;
                    outfile << elnr << "  ";
                    for (k = 1; k <= nelp; k++)
                        outfile << " " << el.PNum(k);
                    outfile <<std::endl;
                    
                }
            }
        }
        else
        {
            std::cout << "\nWrite Permas Volume Mesh" <<std::endl;
            
            int secondorder = (mesh.VolumeElement(1).GetNP() == 10);
            
            if (!secondorder)
            {
                outfile << "$ELEMENT TYPE = TET4  ESET = ALLTET" <<std::endl;
                for (i = 1; i <= ne; i++)
                {
                    const Element & el = mesh.VolumeElement(i);
                    outfile << i 
                            << " " << el.PNum(1) 
                            << " " << el.PNum(2) 
                            << " " << el.PNum(3) 
                            << " " << el.PNum(4) <<std::endl;
                }
            }
            else
            {
                outfile << "$ELEMENT TYPE = TET10  ESET = ALLTET" <<std::endl;
                for (i = 1; i <= ne; i++)
                {
                    const Element & el = mesh.VolumeElement(i);
                    outfile << i 
                            << " " << el.PNum(1) 
                            << " " << el.PNum(5) 
                            << " " << el.PNum(2) 
                            << " " << el.PNum(8) 
                            << " " << el.PNum(3) 
                            << " " << el.PNum(6) <<std::endl << "& "
                            << " " << el.PNum(7) 
                            << " " << el.PNum(9) 
                            << " " << el.PNum(10) 
                            << " " << el.PNum(4) <<std::endl;
                }
            }
            
            outfile <<std::endl <<std::endl;
            
            
            outfile << "$SURFACE GEO  SURFID = 1  SFSET = ALLSUR" <<std::endl;
            for (i = 1; i <= nse; i++)
            {
                const Element2d & el = mesh.SurfaceElement(i);
                if (el.GetNP() == 3)
                    outfile << "STRIA3"
                            << " " << el.PNum(1) 
                            << " " << el.PNum(2) 
                            << " " << el.PNum(3) <<std::endl;
            }    
            
            for (i = 1; i <= nse; i++)
            {
                const Element2d & el = mesh.SurfaceElement(i);
                if (el.GetNP() == 4)
                    outfile << "SQUAD4"
                            << " " << el.PNum(1) 
                            << " " << el.PNum(2) 
                            << " " << el.PNum(3) 
                            << " " << el.PNum(4) <<std::endl;
            }      
            
            for (i = 1; i <= nse; i++)
            {
                const Element2d & el = mesh.SurfaceElement(i);
                if (el.GetNP() == 6)
                    outfile << "STRIA6"
                            << " " << el.PNum(1) 
                            << " " << el.PNum(4) 
                            << " " << el.PNum(2) 
                            << " " << el.PNum(5) 
                            << " " << el.PNum(3) 
                            << " " << el.PNum(6) <<std::endl;
            }      
        }
        
        
        outfile <<std::endl <<std::endl;
        
        outfile << "$COOR  NSET = ALLNODES" <<std::endl;
        
        outfile.precision(6);
        outfile.setf (ios::fixed, ios::floatfield);
        outfile.setf (ios::showpoint);
        
        for (i = 1; i <= np; i++)
        {
            outfile << i << " ";
            outfile << mesh.Point(i)(0) << " ";
            outfile << mesh.Point(i)(1) << " ";
            outfile << mesh.Point(i)(2) << "\n";
        }
    }
    ////////////////////////////////////////////////////////////////////////////////// 
    // \brief Writes PERMAS configuration header into export file
    //        Returns >0 in case of errors
    // \par string &strComp  : Reference to component description
    // \par string &strComp  : Reference to situation description
    ////////////////////////////////////////////////////////////////////////////////// 
    int addComponent(string &strComp, string &strSitu, ofstream &out)
    {
        if (strComp.size() > 12 || strSitu > 12) 
            return 1;

        if (0 == strComp.size()) 
            strComp = "KOMPO1";
        
        if (0 == strSitu.size())
            strSitu = "SIT1";

        // Writing description header of configuration
        out << "$ENTER COMPONENT  NAME = " << strComp << "  DOFTYPE = DISP MATH" <<std::endl <<std::endl;
        out << "   $SITUATION  NAME = " << strSitu <<std::endl;
        out << "   $END SITUATION" <<std::endl <<std::endl;
        out << "   $STRUCTURE" <<std::endl;
        
        return 0;
    }
    
}
