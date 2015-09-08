/*

2d Spline curve for Mesh generator

 */

#include <meshit.hpp>
#include "geometry2d.hpp"
#include "../general/flags.hpp"
#include "../general/ngexception.hpp"

namespace meshit {

    SplineGeometry2d::~SplineGeometry2d()
    {
        for (int i = 0; i < bcnames.size(); i++)
            delete bcnames[i];
        for (int i = 0; i < materials.size(); i++)
            delete [] materials[i];
    }

    void SplineGeometry2d::Load(const char * filename)
    {
        std::ifstream infile;
        Point<2> x;
        char buf[50];

        infile.open(filename);

        if (!infile.good())
            throw NgException(std::string("Input file '") +
                std::string(filename) +
                std::string("' not available!"));

        TestComment(infile);

        infile >> buf; // file recognition

        tensormeshing.resize(0);
        quadmeshing.resize(0);

        TestComment(infile);
        if (strcmp(buf, "splinecurves2dnew") == 0) {
            LoadDataNew(infile);
        }
        else if (strcmp(buf, "splinecurves2dv2") == 0) {
            LoadDataV2(infile);
        }
        else {
            LoadData(infile);
        }
        infile.close();
    }

    void SplineGeometry2d::TestComment(std::istream & infile)
    {
        bool comment = true;
        char ch;
        while (comment == true && !infile.eof()) {
            infile.get(ch);
            if (ch == '#') { // skip comments
                while (ch != '\n' && !infile.eof()) {
                    infile.get(ch);
                }
            }
            else if (ch == '\n') { // skip empty lines
                ;
            }
            else if (isspace(ch)) { // skip whitespaces
                ;
            }
            else { // end of comment
                infile.putback(ch);
                comment = false;
            }
        }
        return;
    }

    void SplineGeometry2d::LoadData(std::istream & infile)
    {
        int nump, numseg, leftdom, rightdom;
        Point<2> x;
        int hi1, hi2, hi3;
        double hd;
        char buf[50], ch;

        materials.resize(0);
        maxh.resize(0);
        infile >> elto0;

        TestComment(infile);

        infile >> nump;
        for (int i = 0; i < nump; i++) {
            TestComment(infile);
            for (int j = 0; j < 2; j++)
                infile >> x(j);
            infile >> hd;

            Flags flags;

            ch = 'a';
            // infile >> ch;
            do {
                infile.get(ch);
            } while (isspace(ch) && ch != '\n');
            while (ch == '-') {
                char flag[100];
                flag[0] = '-';
                infile >> (flag + 1);
                flags.SetCommandLineFlag(flag);
                ch = 'a';
                do {
                    infile.get(ch);
                } while (isspace(ch) && ch != '\n');
            }

            if (infile.good())
                infile.putback(ch);

            geompoints.push_back(GeomPoint<2>(x, hd));
            geompoints.Last().hpref = flags.GetDefineFlag("hpref");
            geompoints.Last().hmax = 1e99;
        }

        PrintMessage(3, nump, " points loaded");
        TestComment(infile);

        infile >> numseg;
        bcnames.resize(numseg);
        for (int i = 0; i < numseg; i++)
            bcnames[i] = 0; // "default";

        SplineSeg<2> * spline = 0;

        PrintMessage(3, numseg, " segments loaded");
        for (int i = 0; i < numseg; i++) {
            TestComment(infile);

            infile >> leftdom >> rightdom;
            infile >> buf;
            // type of spline segement
            if (strcmp(buf, "2") == 0) { // a line
                infile >> hi1 >> hi2;
                spline = new LineSeg<2>(geompoints[hi1 - 1],
                        geompoints[hi2 - 1]);
            }
            else if (strcmp(buf, "3") == 0) { // a rational spline
                infile >> hi1 >> hi2 >> hi3;
                spline = new SplineSeg3<2> (geompoints[hi1 - 1],
                        geompoints[hi2 - 1],
                        geompoints[hi3 - 1]);
            }
            else if (strcmp(buf, "4") == 0) { // an arc
                infile >> hi1 >> hi2 >> hi3;
                spline = new CircleSeg<2> (geompoints[hi1 - 1],
                        geompoints[hi2 - 1],
                        geompoints[hi3 - 1]);
                // 	  break;
            }
            else if (strcmp(buf, "discretepoints") == 0) {
                int npts;
                infile >> npts;
                Array< Point<2> > pts(npts);
                for (int j = 0; j < npts; j++)
                    for (int k = 0; k < 2; k++)
                        infile >> pts[j](k);

                spline = new DiscretePointsSeg<2> (pts);
            }

            SplineSegExt * spex = new SplineSegExt(*spline);

            infile >> spex->reffak;
            spex -> leftdom = leftdom;
            spex -> rightdom = rightdom;
            spex -> hmax = 1e99;
            splines.push_back(spex);

            Flags flags;
            ch = 'a';
            infile >> ch;
            while (ch == '-') {
                char flag[100];
                flag[0] = '-';
                infile >> (flag + 1);
                flags.SetCommandLineFlag(flag);
                ch = 'a';
                infile >> ch;
            }

            if (infile.good())
                infile.putback(ch);

            spex->bc = int (flags.GetNumFlag("bc", i + 1));
            spex->hpref_left = int (flags.GetDefineFlag("hpref")) ||
                    int (flags.GetDefineFlag("hprefleft"));
            spex->hpref_right = int (flags.GetDefineFlag("hpref")) ||
                    int (flags.GetDefineFlag("hprefright"));
            spex->copyfrom = int (flags.GetNumFlag("copy", -1));
            if (flags.StringFlagDefined("bcname")) {
                int mybc = spex->bc - 1;
                delete bcnames[mybc];
                bcnames[mybc] = new std::string(flags.GetStringFlag("bcname", ""));
            }
        }
    }

    void SplineGeometry2d::LoadDataNew(std::istream & infile)
    {
        int nump, numseg, leftdom, rightdom;
        Point<2> x;
        int hi1, hi2, hi3;
        double hd;
        char buf[50], ch;
        int pointnr;


        TestComment(infile);
        infile >> elto0;
        TestComment(infile);

        infile >> nump;
        geompoints.resize(nump);

        for (int i = 0; i < nump; i++) {
            TestComment(infile);
            infile >> pointnr;
            if (pointnr > nump) {
                throw NgException(std::string("Point number greater than total number of points"));
            }
            for (int j = 0; j < 2; j++)
                infile >> x(j);


            // hd is now optional, default 1
            //  infile >> hd;
            hd = 1;

            Flags flags;


            // get flags, 
            ch = 'a';
            // infile >> ch;
            do {

                infile.get(ch);
                // if another int-value, set refinement flag to this value
                // (corresponding to old files)
                if (int (ch) >= 48 && int(ch) <= 57) {
                    infile.putback(ch);
                    infile >> hd;
                    infile.get(ch);
                }
            } while (isspace(ch) && ch != '\n');
            while (ch == '-') {
                char flag[100];
                flag[0] = '-';
                infile >> (flag + 1);
                flags.SetCommandLineFlag(flag);
                ch = 'a';
                do {
                    infile.get(ch);
                } while (isspace(ch) && ch != '\n');
            }

            if (infile.good())
                infile.putback(ch);

            if (hd == 1)
                hd = flags.GetNumFlag("ref", 1.0);
            //       geompoints.Append (GeomPoint<D>(x, hd));
            geompoints[pointnr - 1] = GeomPoint<2>(x, hd);
            geompoints[pointnr - 1].hpref = flags.GetDefineFlag("hpref");
        }

        TestComment(infile);

        infile >> numseg;
        bcnames.resize(numseg);
        for (int i = 0; i < numseg; i++)
            bcnames[i] = 0; //new"default";

        SplineSeg<2> * spline = 0;
        for (int i = 0; i < numseg; i++) {
            TestComment(infile);

            infile >> leftdom >> rightdom;

            // std::cout << "add spline " << i << ", left = " << leftdom <<std::endl;

            infile >> buf;
            // type of spline segement
            if (strcmp(buf, "2") == 0) { // a line
                infile >> hi1 >> hi2;
                spline = new LineSeg<2> (geompoints[hi1 - 1],
                        geompoints[hi2 - 1]);
            }
            else if (strcmp(buf, "3") == 0) { // a rational spline
                infile >> hi1 >> hi2 >> hi3;
                spline = new SplineSeg3<2> (geompoints[hi1 - 1],
                        geompoints[hi2 - 1],
                        geompoints[hi3 - 1]);
            }
            else if (strcmp(buf, "4") == 0) { // an arc
                infile >> hi1 >> hi2 >> hi3;
                spline = new CircleSeg<2> (geompoints[hi1 - 1],
                        geompoints[hi2 - 1],
                        geompoints[hi3 - 1]);
                // 	  break;
            }
            else if (strcmp(buf, "discretepoints") == 0) {
                int npts;
                infile >> npts;
                Array< Point<2> > pts(npts);
                for (int j = 0; j < npts; j++)
                    for (int k = 0; k < 2; k++)
                        infile >> pts[j](k);

                spline = new DiscretePointsSeg<2> (pts);
            }

            //      infile >> spline->reffak;

            SplineSegExt * spex = new SplineSegExt(*spline);

            spex -> leftdom = leftdom;
            spex -> rightdom = rightdom;
            splines.push_back(spex);

            // hd is now optional, default 1
            //  infile >> hd;
            hd = 1;
            infile >> ch;

            // get refinement parameter, if it is there
            // infile.get (ch);
            // if another int-value, set refinement flag to this value
            // (corresponding to old files)
            if (int (ch) >= 48 && int(ch) <= 57) {
                infile.putback(ch);
                infile >> hd;
                infile >> ch;
            }

            Flags flags;
            while (ch == '-') {
                char flag[100];
                flag[0] = '-';
                infile >> (flag + 1);
                flags.SetCommandLineFlag(flag);
                ch = 'a';
                infile >> ch;
            }

            if (infile.good())
                infile.putback(ch);

            spex->bc = int (flags.GetNumFlag("bc", i + 1));
            spex->hpref_left = int (flags.GetDefineFlag("hpref")) ||
                    int (flags.GetDefineFlag("hprefleft"));
            spex->hpref_right = int (flags.GetDefineFlag("hpref")) ||
                    int (flags.GetDefineFlag("hprefright"));
            spex->copyfrom = int (flags.GetNumFlag("copy", -1));
            spex->reffak = flags.GetNumFlag("ref", 1);
            spex->hmax = flags.GetNumFlag("maxh", 1e99);

            if (flags.StringFlagDefined("bcname")) {
                int mybc = spex->bc - 1;
                if (bcnames[mybc]) delete bcnames[mybc];
                bcnames[mybc] = new std::string(flags.GetStringFlag("bcname", ""));
            }

            if (hd != 1)
                spex->reffak = hd;
        }
        if (!infile.good())
            return;
        TestComment(infile);
        int numdomains;
        int domainnr;
        char material[100];

        if (!infile.good())
            return;

        infile >> numdomains;
        materials.resize(numdomains);
        maxh.resize(numdomains);
        maxh = 1e99;

        TestComment(infile);

        for (int i = 0; i < numdomains; i++)
            materials [ i ] = new char (100);

        for (int i = 0; i < numdomains && infile.good(); i++) {
            TestComment(infile);
            infile >> domainnr;
            infile >> material;
            strcpy(materials[domainnr - 1], material);

            Flags flags;
            ch = 'a';
            infile >> ch;
            while (ch == '-') {
                char flag[100];
                flag[0] = '-';
                infile >> (flag + 1);
                flags.SetCommandLineFlag(flag);
                ch = 'a';
                infile >> ch;
            }

            if (infile.good())
                infile.putback(ch);

            maxh[domainnr - 1] = flags.GetNumFlag("maxh", 1e99);
        }
        return;
    }

    void SplineGeometry2d::LoadDataV2(std::istream & infile)
    {

        enum
        {
            D = 2
        };
        // new parser by Astrid Sinwel

        PrintMessage(1, "Load 2D Geometry V2");
        int nump, leftdom, rightdom;
        Point<D> x;
        int hi1, hi2, hi3;
        double hd;
        char buf[50], ch;
        int pointnr;

        std::string keyword;

        Array < GeomPoint<D> > infilepoints(0);
        Array <int> pointnrs(0);
        nump = 0;
        int numdomains = 0;


        TestComment(infile);
        // refinement factor
        infile >> elto0;
        TestComment(infile);


        // test if next ch is a letter, i.e. new keyword starts
        bool ischar = false;

        while (infile.good()) {
            infile >> keyword;

            ischar = false;

            if (keyword == "points") {
                PrintMessage(3, "load points");
                infile.get(ch);
                infile.putback(ch);

                // test if ch is a letter
                if (int(ch) >= 65 && int(ch) <= 90)
                    ischar = true;
                if (int(ch) >= 97 && int(ch) <= 122)
                    ischar = true;

                while (!ischar) {
                    TestComment(infile);
                    infile >> pointnr;
                    // pointnrs 1-based
                    if (pointnr > nump) nump = pointnr;
                    pointnrs.push_back(pointnr);

                    for (int j = 0; j < D; j++)
                        infile >> x(j);
                    // hd is now optional, default 1
                    //  infile >> hd;
                    hd = 1;

                    Flags flags;


                    // get flags, 
                    ch = 'a';
                    // infile >> ch;
                    do {
                        infile.get(ch);
                        // if another int-value, set refinement flag to this value
                        // (corresponding to old files)
                        if (int (ch) >= 48 && int(ch) <= 57) {
                            infile.putback(ch);
                            infile >> hd;
                            infile.get(ch);
                        }
                    } while (isspace(ch) && ch != '\n');
                    while (ch == '-') {
                        char flag[100];
                        flag[0] = '-';
                        infile >> (flag + 1);
                        flags.SetCommandLineFlag(flag);
                        ch = 'a';
                        do {
                            infile.get(ch);
                        } while (isspace(ch) && ch != '\n');
                    }
                    if (infile.good())
                        infile.putback(ch);

                    if (hd == 1)
                        hd = flags.GetNumFlag("ref", 1.0);
                    //       geompoints.Append (GeomPoint<D>(x, hd));

                    infilepoints.push_back(GeomPoint<D>(x, hd));
                    infilepoints.Last().hpref = flags.GetDefineFlag("hpref");
                    infilepoints.Last().hmax = flags.GetNumFlag("maxh", 1e99);

                    TestComment(infile);
                    infile.get(ch);
                    infile.putback(ch);

                    // test if letter
                    if (int(ch) >= 65 && int(ch) <= 90)
                        ischar = true;
                    if (int(ch) >= 97 && int(ch) <= 122)
                        ischar = true;
                }

                //	  infile.putback (ch);

                geompoints.resize(nump);
                for (int i = 0; i < nump; i++) {
                    geompoints[pointnrs[i] - 1] = infilepoints[i];
                    geompoints[pointnrs[i] - 1].hpref = infilepoints[i].hpref;
                }
                TestComment(infile);
            }

            else if (keyword == "segments") {
                PrintMessage(3, "load segments");

                bcnames.resize(0);
                infile.get(ch);
                infile.putback(ch);
                int i = 0;

                // test if ch is a letter
                if (int(ch) >= 65 && int(ch) <= 90)
                    ischar = true;
                if (int(ch) >= 97 && int(ch) <= 122)
                    ischar = true;

                while (!ischar) //ch != 'p' && ch != 'm' )
                {
                    i++;
                    TestComment(infile);

                    SplineSeg<D> * spline = 0;
                    TestComment(infile);

                    infile >> leftdom >> rightdom;

                    if (leftdom > numdomains) numdomains = leftdom;
                    if (rightdom > numdomains) numdomains = rightdom;


                    infile >> buf;
                    // type of spline segement
                    if (strcmp(buf, "2") == 0) { // a line
                        infile >> hi1 >> hi2;
                        spline = new LineSeg<D>(geompoints[hi1 - 1],
                                geompoints[hi2 - 1]);
                    }
                    else if (strcmp(buf, "3") == 0) { // a rational spline
                        infile >> hi1 >> hi2 >> hi3;
                        spline = new SplineSeg3<D> (geompoints[hi1 - 1],
                                geompoints[hi2 - 1],
                                geompoints[hi3 - 1]);
                    }
                    else if (strcmp(buf, "4") == 0) { // an arc
                        infile >> hi1 >> hi2 >> hi3;
                        spline = new CircleSeg<D> (geompoints[hi1 - 1],
                                geompoints[hi2 - 1],
                                geompoints[hi3 - 1]);
                    }
                    else if (strcmp(buf, "discretepoints") == 0) {
                        int npts;
                        infile >> npts;
                        Array< Point<D> > pts(npts);
                        for (int j = 0; j < npts; j++)
                            for (int k = 0; k < D; k++)
                                infile >> pts[j](k);

                        spline = new DiscretePointsSeg<D> (pts);
                    }
                    else if (strcmp(buf, "bsplinepoints") == 0) {
                        int npts, order;
                        infile >> npts;
                        infile >> order;
                        Array< Point<D> > pts(npts);
                        for (int j = 0; j < npts; j++)
                            for (int k = 0; k < D; k++)
                                infile >> pts[j](k);
                        if (order < 2)
                            std::cerr << "Minimum order of 2 is required!!" << std::endl;
                        else if (order == 2)
                            spline = new BSplineSeg<D, 2> (pts);
                        else if (order == 3)
                            spline = new BSplineSeg<D, 3> (pts);
                        else if (order == 4)
                            spline = new BSplineSeg<D, 4> (pts);
                        else if (order > 4)
                            std::cerr << "Maximum allowed order is 4!!" << std::endl;
                    }

                    //      infile >> spline->reffak;
                    SplineSegExt * spex = new SplineSegExt(*spline);

                    spex -> leftdom = leftdom;
                    spex -> rightdom = rightdom;
                    splines.push_back(spex);


                    // hd is now optional, default 1
                    hd = 1;

                    infile >> ch;

                    // get flags, 
                    Flags flags;
                    while (ch == '-') {
                        char flag[100];
                        flag[0] = '-';
                        infile >> (flag + 1);
                        flags.SetCommandLineFlag(flag);
                        ch = 'a';
                        infile >> ch;
                    }

                    if (infile.good())
                        infile.putback(ch);

                    spex->bc = int (flags.GetNumFlag("bc", i + 1));
                    spex->hpref_left = int (flags.GetDefineFlag("hpref")) ||
                            int (flags.GetDefineFlag("hprefleft"));
                    spex->hpref_right = int (flags.GetDefineFlag("hpref")) ||
                            int (flags.GetDefineFlag("hprefright"));
                    spex->copyfrom = int (flags.GetNumFlag("copy", -1));
                    spex->reffak = flags.GetNumFlag("ref", 1);
                    spex->hmax = flags.GetNumFlag("maxh", 1e99);
                    if (hd != 1) spex->reffak = hd;

                    if (flags.StringFlagDefined("bcname")) {
                        int mybc = spex->bc - 1;
                        for (int ii = bcnames.size(); ii <= mybc; ii++)
                            bcnames.push_back(new std::string("default"));
                        if (bcnames[mybc]) delete bcnames[mybc];
                        bcnames[mybc] = new std::string(flags.GetStringFlag("bcname", ""));
                    }

                    TestComment(infile);
                    infile.get(ch);
                    infile.putback(ch);

                    // test if ch is a letter
                    if (int(ch) >= 65 && int(ch) <= 90)
                        ischar = true;
                    if (int(ch) >= 97 && int(ch) <= 122)
                        ischar = true;

                }

                infile.get(ch);
                infile.putback(ch);


            }
            else if (keyword == "materials") {
                TestComment(infile);
                int domainnr;
                char material[100];

                if (!infile.good())
                    return;

                materials.resize(numdomains);
                maxh.resize(numdomains);
                for (int i = 0; i < numdomains; i++)
                    maxh[i] = 1000;
                quadmeshing.resize(numdomains);
                quadmeshing = false;
                tensormeshing.resize(numdomains);
                tensormeshing = false;
                layer.resize(numdomains);
                layer = 1;


                TestComment(infile);

                for (int i = 0; i < numdomains; i++)
                    materials [ i ] = new char[100];

                for (int i = 0; i < numdomains && infile.good(); i++) {
                    TestComment(infile);
                    infile >> domainnr;
                    infile >> material;

                    strcpy(materials[domainnr - 1], material);

                    Flags flags;
                    ch = 'a';
                    infile >> ch;
                    while (ch == '-') {
                        char flag[100];
                        flag[0] = '-';
                        infile >> (flag + 1);
                        flags.SetCommandLineFlag(flag);
                        ch = 'a';
                        infile >> ch;
                    }

                    if (infile.good())
                        infile.putback(ch);

                    maxh[domainnr - 1] = flags.GetNumFlag("maxh", 1000);
                    if (flags.GetDefineFlag("quad")) quadmeshing[domainnr - 1] = true;
                    if (flags.GetDefineFlag("tensor")) tensormeshing[domainnr - 1] = true;
                    layer[domainnr - 1] = int(flags.GetNumFlag("layer", 1));
                }
            }
        }
        return;
    }

    void SplineGeometry2d::AddLine(
            const std::vector<Point2d>& point_list,
            double hmax,
            bool hole,
            int bc)
    {

        size_t nold_points = geompoints.size();
        size_t nnew_points = point_list.size();

        std::vector<GeomPoint<2> > gpts;

        gpts.reserve(nnew_points);
        geompoints.reserve(nold_points + nnew_points);

        for (size_t i = 0; i < nnew_points; i++) {
            gpts.push_back(GeomPoint<2>(point_list[i]));
            geompoints.push_back(gpts[i]);
        }
        for (size_t i = 0; i < nnew_points; i++) {
            size_t i0 = (hole) ? nnew_points - i - 1 : i;
            size_t i1 = (hole) ? (i0 == 0) ? nnew_points - 1 : i0 - 1
                    : (i0 == nnew_points - 1) ? 0 : i0 + 1;
            SplineSeg<2> * spline = new LineSeg<2>(gpts[i0], gpts[i1]);
            SplineSegExt * seg = new SplineSegExt(*spline);
            if (hole) {
                seg->leftdom = 0;
                seg->rightdom = 1;
            }
            else {
                seg->leftdom = 1;
                seg->rightdom = 0;
            }
            seg->bc = bc;
            seg->hpref_left = false;
            seg->hpref_right = false;
            seg->reffak = 1; // Refinement factor
            seg->hmax = hmax;
            splines.push_back(seg);
        }
    }

    void SplineGeometry2d::FakeData()
    {
        int numdomains = 1;
        materials.resize(numdomains);
        maxh.assign(numdomains, 1000);
        quadmeshing.assign(numdomains, false);
        layer.assign(numdomains, 1);

        for (int i = 0; i < numdomains; i++) {
            materials[i] = new char[1];
            materials[i][0] = '\0';
        }
    }

    std::string SplineGeometry2d::GetBCName(const int bcnr) const
    {
        if (bcnames.size() >= bcnr)
            if (bcnames[bcnr - 1])
                return *bcnames[bcnr - 1];
        return "default";
    }

    std::string * SplineGeometry2d::BCNamePtr(const int bcnr)
    {
        if (bcnr > bcnames.size())
            return 0;
        else
            return bcnames[bcnr - 1];
    }

    void SplineGeometry2d::GetMaterial(const int domnr, char* & material)
    {
        if (materials.size() >= domnr)
            material = materials[domnr - 1];
        else
            material = 0;
    }

    double SplineGeometry2d::GetDomainMaxh(const int domnr)
    {
        if (maxh.size() >= domnr && domnr > 0)
            return maxh[domnr - 1];
        else
            return -1;
    }

    int SplineGeometry2d::GenerateMesh(Mesh*& mesh, MeshingParameters & mp,
            int perfstepsstart, int perfstepsend)
    {
        mesh->BuildFromSpline2D(*this, mp);
        return 0;
    }

    Refinement & SplineGeometry2d::GetRefinement() const
    {
        return * new Refinement2d(*this);
    }

}