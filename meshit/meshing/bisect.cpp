#include <sstream>
#include <stdexcept>

#include <meshit.hpp>
#include "bisect.hpp"
#include "../gprim/geomfuncs.hpp"
#include "validate.hpp"
#include "../general/mystring.hpp"

namespace meshit {

    class MarkedTet;
    class MarkedPrism;
    class MarkedIdentification;
    class MarkedTri;
    class MarkedQuad;

    typedef Array<MarkedTet> T_MTETS;
    typedef Array<MarkedPrism> T_MPRISMS;
    typedef Array<MarkedIdentification> T_MIDS;
    typedef Array<MarkedTri> T_MTRIS;
    typedef Array<MarkedQuad> T_MQUADS;

    class MarkedTet
    {
      public:
        /// pnums of tet
        PointIndex pnums[4];
        /// material number
        int matindex;
        /// element marked for refinement
        /// marked = 1: marked by element marker, marked = 2 due to closure
        unsigned int marked : 2;
        /// flag of Arnold-Mukherjee algorithm
        unsigned int flagged : 1;
        /// tetedge (local coordinates 0..3)
        unsigned int tetedge1 : 3;
        unsigned int tetedge2 : 3;
        // marked edge of faces
        // face_j : face without node j,
        // mark_k : edge without node k

        char faceedges[4];
        // unsigned char faceedges[4];
        bool incorder;
        unsigned int order : 6;

        MarkedTet()
        {
            for (int i = 0; i < 4; i++) {
                faceedges[i] = 255;
            }
        }
    };

    std::ostream & operator<<(std::ostream & ost, const MarkedTet & mt)
    {
        for (int i = 0; i < 4; i++)
            ost << mt.pnums[i] << " ";

        ost << mt.matindex << " " << int(mt.marked) << " " << int(mt.flagged) << " " << int(mt.tetedge1) << " " << int(mt.tetedge2) << " ";

        ost << "faceedges = ";
        for (int i = 0; i < 4; i++)
            ost << int(mt.faceedges[i]) << " ";

        ost << " order = ";
        ost << mt.incorder << " " << int(mt.order) << "\n";
        return ost;
    }

    std::istream & operator>>(std::istream & ost, MarkedTet & mt)
    {
        for (int i = 0; i < 4; i++)
            ost >> mt.pnums[i];

        ost >> mt.matindex;

        int auxint;
        ost >> auxint;
        mt.marked = auxint;
        ost >> auxint;
        mt.flagged = auxint;
        ost >> auxint;
        mt.tetedge1 = auxint;
        ost >> auxint;
        mt.tetedge2 = auxint;

        char auxchar;

        for (int i = 0; i < 4; i++) {
            ost >> auxchar;
            mt.faceedges[i] = auxchar;
        }

        ost >> mt.incorder;
        ost >> auxint;
        mt.order = auxint;
        return ost;
    }

    class MarkedPrism
    {
      public:
        /// 6 point numbers
        PointIndex pnums[6];
        /// material number
        int matindex;
        /// marked for refinement
        int marked;
        /// edge without node k (0,1,2)
        int markededge;

        bool incorder;
        unsigned int order : 6;
    };

    std::ostream & operator<<(std::ostream & ost, const MarkedPrism & mp)
    {
        for (int i = 0; i < 6; i++)
            ost << mp.pnums[i] << " ";

        ost << mp.matindex << " " << mp.marked << " " << mp.markededge << " " << mp.incorder << " " << int(mp.order) << "\n";
        return ost;
    }

    std::istream & operator>>(std::istream & ist, MarkedPrism & mp)
    {
        for (int i = 0; i < 6; i++)
            ist >> mp.pnums[i];

        ist >> mp.matindex >> mp.marked >> mp.markededge >> mp.incorder;
        int auxint;
        ist >> auxint;
        mp.order = auxint;
        return ist;
    }

    class MarkedIdentification
    {
      public:
        // number of points of one face (3 or 4)
        int np;
        /// 6 or 8 point numbers
        PointIndex pnums[8];
        /// marked for refinement
        int marked;
        /// edge starting with node k (0,1,2, or 3)
        int markededge;

        bool incorder;
        unsigned int order : 6;
    };

    std::ostream & operator<<(std::ostream & ost, const MarkedIdentification & mi)
    {
        ost << mi.np << " ";
        for (int i = 0; i < 2 * mi.np; i++)
            ost << mi.pnums[i] << " ";
        ost << mi.marked << " " << mi.markededge << " " << mi.incorder << " " << int(mi.order) << "\n";
        return ost;
    }

    std::istream & operator>>(std::istream & ist, MarkedIdentification & mi)
    {
        ist >> mi.np;
        for (int i = 0; i < 2 * mi.np; i++)
            ist >> mi.pnums[i];
        ist >> mi.marked >> mi.markededge >> mi.incorder;
        int auxint;
        ist >> auxint;
        mi.order = auxint;
        return ist;
    }

    class MarkedTri
    {
      public:
        /// three point numbers
        PointIndex pnums[3];
        /// three geominfos
        PointGeomInfo pgeominfo[3];
        /// marked for refinement
        int marked;
        /// edge without node k
        int markededge;
        /// surface id
        int surfid;

        bool incorder;
        unsigned int order : 6;
    };

    std::ostream & operator<<(std::ostream & ost, const MarkedTri & mt)
    {
        for (int i = 0; i < 3; i++)
            ost << mt.pnums[i] << " ";
        for (int i = 0; i < 3; i++)
            ost << mt.pgeominfo[i] << " ";
        ost << mt.marked << " " << mt.markededge << " " << mt.surfid << " " << mt.incorder << " " << int(mt.order) << "\n";
        return ost;
    }

    std::istream & operator>>(std::istream & ist, MarkedTri & mt)
    {
        for (int i = 0; i < 3; i++)
            ist >> mt.pnums[i];
        for (int i = 0; i < 3; i++)
            ist >> mt.pgeominfo[i];
        ist >> mt.marked >> mt.markededge >> mt.surfid >> mt.incorder;
        int auxint;
        ist >> auxint;
        mt.order = auxint;
        return ist;
    }

    class MarkedQuad
    {
      public:
        /// point numbers
        PointIndex pnums[4];
        PointGeomInfo pgeominfo[4];
        /// marked for refinement
        int marked;
        /// marked edge: 0/2 = vertical, 1/3 = horizontal
        int markededge;
        /// surface id
        int surfid;

        bool incorder;
        unsigned int order : 6;
    };

    std::ostream & operator<<(std::ostream & ost, const MarkedQuad & mt)
    {
        for (int i = 0; i < 4; i++)
            ost << mt.pnums[i] << " ";
        for (int i = 0; i < 4; i++)
            ost << mt.pgeominfo[i] << " ";
        ost << mt.marked << " " << mt.markededge << " " << mt.surfid << " " << mt.incorder << " " << int(mt.order) << "\n";
        return ost;
    }

    std::istream & operator>>(std::istream & ist, MarkedQuad & mt)
    {
        for (int i = 0; i < 4; i++)
            ist >> mt.pnums[i];
        for (int i = 0; i < 4; i++)
            ist >> mt.pgeominfo[i];
        ist >> mt.marked >> mt.markededge >> mt.surfid >> mt.incorder;
        int auxint;
        ist >> auxint;
        mt.order = auxint;
        return ist;
    }

    void PrettyPrint(std::ostream & ost, const MarkedTet & mt)
    {
        int te1 = mt.tetedge1;
        int te2 = mt.tetedge2;
        int order = mt.order;

        ost << "MT: " << mt.pnums[0] << " - " << mt.pnums[1] << " - "
                << mt.pnums[2] << " - " << mt.pnums[3] << std::endl
                << "marked edge: " << te1 << " - " << te2
                << ", order = " << order << std::endl;
        for (int k = 0; k < 4; k++) {
            ost << "face";
            for (int j = 0; j < 4; j++)
                if (j != k)
                    ost << " " << mt.pnums[j];
            for (int i = 0; i < 3; i++)
                for (int j = i + 1; j < 4; j++)
                    if (i != k && j != k && int(mt.faceedges[k]) == 6 - k - i - j)
                        ost << " marked edge " << mt.pnums[i] << " " << mt.pnums[j] << std::endl;
        }
        ost << std::endl;
    }

    int BTSortEdges(const Mesh & mesh,
            const Array< Array<int, PointIndex::BASE>* > & idmaps,
            INDEX_2_CLOSED_HASHTABLE<int> & edgenumber)
    {
        LOG_DEBUG("sorting ... ");

        //  if (mesh.PureTetMesh())
        if (1) {
            // new, fast version

            Array<INDEX_2> edges;
            Array<int> eclasses;

            int i, j, k;
            int cntedges = 0;
            int go_on;
            int ned(0);

            // enumerate edges:
            for (i = 1; i <= mesh.GetNE(); i++) {
                const Element & el = mesh.VolumeElement(i);
                static int tetedges[6][2] = {
                    { 1, 2},
                    { 1, 3},
                    { 1, 4},
                    { 2, 3},
                    { 2, 4},
                    { 3, 4}
                };
                static int prismedges[9][2] = {
                    { 1, 2},
                    { 1, 3},
                    { 2, 3},
                    { 4, 5},
                    { 4, 6},
                    { 5, 6},
                    { 1, 4},
                    { 2, 5},
                    { 3, 6}
                };
                int pyramidedges[6][2] = {
                    { 1, 2},
                    { 3, 4},
                    { 1, 5},
                    { 2, 5},
                    { 3, 5},
                    { 4, 5}
                };

                int (*tip)[2] = NULL;

                switch (el.GetType()) {
                    case TET:
                    case TET10:
                    {
                        tip = tetedges;
                        ned = 6;
                        break;
                    }
                    case PRISM:
                    case PRISM12:
                    {
                        tip = prismedges;
                        ned = 6;
                        break;
                    }
                    case PYRAMID:
                    {
                        tip = pyramidedges;
                        ned = 6;
                        break;
                    }
                    default:
                        throw std::runtime_error("Bisect, element type not handled in switch");
                }

                for (j = 0; j < ned; j++) {
                    INDEX_2 i2(el.PNum(tip[j][0]), el.PNum(tip[j][1]));
                    i2.Sort();
                    //std::cerr << "edge " << i2 <<std::endl;
                    if (!edgenumber.Used(i2)) {
                        cntedges++;
                        edges.push_back(i2);
                        edgenumber.Set(i2, cntedges);
                    }
                }
            }

            // additional surface edges:
            for (i = 1; i <= mesh.GetNSE(); i++) {
                const Element2d & el = mesh.SurfaceElement(i);
                static int trigedges[3][2] = {
                    { 1, 2},
                    { 2, 3},
                    { 3, 1}
                };

                static int quadedges[4][2] = {
                    { 1, 2},
                    { 2, 3},
                    { 3, 4},
                    { 4, 1}
                };

                int (*tip)[2] = NULL;

                switch (el.GetType()) {
                    case TRIG:
                    case TRIG6:
                    {
                        tip = trigedges;
                        ned = 3;
                        break;
                    }
                    case QUAD:
                    case QUAD6:
                    {
                        tip = quadedges;
                        ned = 4;
                        break;
                    }
                    default:
                    {
                        std::cerr << "Error: Sort for Bisect, SE has " << el.GetNP() << " points" << std::endl;
                        ned = 0;
                    }
                }

                for (j = 0; j < ned; j++) {
                    INDEX_2 i2(el.PNum(tip[j][0]), el.PNum(tip[j][1]));
                    i2.Sort();
                    if (!edgenumber.Used(i2)) {
                        cntedges++;
                        edges.push_back(i2);
                        edgenumber.Set(i2, cntedges);
                    }
                }
            }

            eclasses.resize(cntedges);
            for (i = 1; i <= cntedges; i++)
                eclasses.Elem(i) = i;

            // identify edges in element stack
            do {
                go_on = 0;
                for (i = 1; i <= mesh.GetNE(); i++) {
                    const Element & el = mesh.VolumeElement(i);

                    if (el.GetType() != PRISM &&
                            el.GetType() != PRISM12 &&
                            el.GetType() != PYRAMID)
                        continue;

                    int prismpairs[3][4] = {
                        { 1, 2, 4, 5},
                        { 2, 3, 5, 6},
                        { 1, 3, 4, 6}
                    };

                    int pyramidpairs[3][4] = {
                        { 1, 2, 4, 3},
                        { 1, 5, 4, 5},
                        { 2, 5, 3, 5}
                    };

                    int (*pairs)[4] = NULL;
                    switch (el.GetType()) {
                        case PRISM:
                        case PRISM12:
                        {
                            pairs = prismpairs;
                            break;
                        }
                        case PYRAMID:
                        {
                            pairs = pyramidpairs;
                            break;
                        }
                        default:
                            throw std::runtime_error("Bisect, element type not handled in switch, 2");
                    }

                    for (j = 0; j < 3; j++) {
                        INDEX_2 e1(el.PNum(pairs[j][0]),
                                el.PNum(pairs[j][1]));
                        INDEX_2 e2(el.PNum(pairs[j][2]),
                                el.PNum(pairs[j][3]));
                        e1.Sort();
                        e2.Sort();

                        int eclass1 = edgenumber.Get(e1);
                        int eclass2 = edgenumber.Get(e2);

                        if (eclasses.Get(eclass1) >
                                eclasses.Get(eclass2)) {
                            eclasses.Elem(eclass1) =
                                    eclasses.Get(eclass2);
                            go_on = 1;
                        }
                        else if (eclasses.Get(eclass2) >
                                eclasses.Get(eclass1)) {
                            eclasses.Elem(eclass2) =
                                    eclasses.Get(eclass1);
                            go_on = 1;
                        }
                    }
                }

                for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++) {
                    const Element2d & el2d = mesh.SurfaceElement(sei);

                    for (i = 0; i < el2d.GetNP(); i++) {
                        INDEX_2 e1(el2d[i], el2d[(i + 1) % el2d.GetNP()]);
                        e1.Sort();
                        INDEX_2 e2;

                        for (k = 0; k < idmaps.size(); k++) {
                            e2.I1() = (*idmaps[k])[e1.I1()];
                            e2.I2() = (*idmaps[k])[e1.I2()];

                            if (e2.I1() == 0 || e2.I2() == 0 ||
                                    e1.I1() == e2.I1() || e1.I2() == e2.I2())
                                continue;

                            e2.Sort();
                            if (!edgenumber.Used(e2))
                                continue;

                            int eclass1 = edgenumber.Get(e1);
                            int eclass2 = edgenumber.Get(e2);

                            if (eclasses.Get(eclass1) >
                                    eclasses.Get(eclass2)) {
                                eclasses.Elem(eclass1) =
                                        eclasses.Get(eclass2);

                                go_on = 1;
                            }
                            else if (eclasses.Get(eclass2) >
                                    eclasses.Get(eclass1)) {
                                eclasses.Elem(eclass2) =
                                        eclasses.Get(eclass1);
                                go_on = 1;
                            }
                        }
                    }

                }

            } while (go_on);

            // compute classlength:
            Array<double> edgelength(cntedges);

            for (i = 1; i <= cntedges; i++) {
                INDEX_2 edge = edges.Get(i);
                double elen = Dist(mesh.Point(edge.I1()),
                        mesh.Point(edge.I2()));
                edgelength.Elem(i) = elen;
            }

            for (i = 1; i <= cntedges; i++) {
                if (eclasses.Get(i) != i) {
                    if (edgelength.Get(i) < edgelength.Get(eclasses.Get(i)))
                        edgelength.Elem(eclasses.Get(i)) = edgelength.Get(i);
                    edgelength.Elem(i) = 1e20;
                }
            }

            TABLE<int> eclasstab(cntedges);
            for (i = 1; i <= cntedges; i++)
                eclasstab.Add1(eclasses.Get(i), i);

            // sort edges:
            Array<int> sorted(cntedges);

            QuickSort(edgelength, sorted);

            int cnt = 0;
            for (i = 1; i <= cntedges; i++) {
                int ii = sorted.Get(i);
                for (j = 1; j <= eclasstab.EntrySize(ii); j++) {
                    cnt++;
                    edgenumber.Set(edges.Get(eclasstab.Get(ii, j)), cnt);
                }
            }
            return cnt;
        }

        else {
            // old version

            int i, j;
            int cnt = 0;
            int found;
            double len2, maxlen2;
            INDEX_2 ep;

            // sort edges by length, parallel edges (on prisms)
            // are added in blocks

            do {
                found = 0;
                maxlen2 = 1e30;

                for (i = 1; i <= mesh.GetNE(); i++) {
                    const Element & el = mesh.VolumeElement(i);
                    int ned;
                    int tetedges[6][2] = {
                        { 1, 2},
                        { 1, 3},
                        { 1, 4},
                        { 2, 3},
                        { 2, 4},
                        { 3, 4}
                    };
                    int prismedges[6][2] = {
                        { 1, 2},
                        { 1, 3},
                        { 2, 4},
                        { 4, 5},
                        { 4, 6},
                        { 5, 6}
                    };
                    int pyramidedges[6][2] = {
                        { 1, 2},
                        { 3, 4},
                        { 1, 5},
                        { 2, 5},
                        { 3, 5},
                        { 4, 5}
                    };

                    int (*tip)[2];

                    switch (el.GetType()) {
                        case TET:
                        {
                            tip = tetedges;
                            ned = 6;
                            break;
                        }
                        case PRISM:
                        {
                            tip = prismedges;
                            ned = 6;
                            break;
                        }
                        case PYRAMID:
                        {
                            tip = pyramidedges;
                            ned = 6;
                            break;
                        }
                        default:
                            throw std::runtime_error("Bisect, element type not handled in switch, 3");
                    }

                    for (j = 0; j < ned; j++) {
                        INDEX_2 i2(el.PNum(tip[j][0]), el.PNum(tip[j][1]));
                        i2.Sort();
                        if (!edgenumber.Used(i2)) {
                            len2 = Dist(mesh.Point(i2.I1()),
                                    mesh.Point(i2.I2()));
                            if (len2 < maxlen2) {
                                maxlen2 = len2;
                                ep = i2;
                                found = 1;
                            }
                        }
                    }
                }
                if (found) {
                    cnt++;
                    edgenumber.Set(ep, cnt);

                    // find connected edges:
                    int go_on = 0;
                    do {
                        go_on = 0;
                        for (i = 1; i <= mesh.GetNE(); i++) {
                            const Element & el = mesh.VolumeElement(i);
                            if (el.GetNP() != 6) continue;

                            int prismpairs[3][4] = {
                                { 1, 2, 4, 5},
                                { 2, 3, 5, 6},
                                { 1, 3, 4, 6}
                            };

                            int pyramidpairs[3][4] = {
                                { 1, 2, 4, 3},
                                { 1, 5, 4, 5},
                                { 2, 5, 3, 5}
                            };

                            int (*pairs)[4];
                            switch (el.GetType()) {
                                case PRISM:
                                {
                                    pairs = prismpairs;
                                    break;
                                }
                                case PYRAMID:
                                {
                                    pairs = pyramidpairs;
                                    break;
                                }
                                default:
                                    throw std::runtime_error("Bisect, element type not handled in switch, 3a");
                            }

                            for (j = 0; j < 3; j++) {
                                INDEX_2 e1(el.PNum(pairs[j][0]),
                                        el.PNum(pairs[j][1]));
                                INDEX_2 e2(el.PNum(pairs[j][2]),
                                        el.PNum(pairs[j][3]));
                                e1.Sort();
                                e2.Sort();

                                int used1 = edgenumber.Used(e1);
                                int used2 = edgenumber.Used(e2);

                                if (used1 && !used2) {
                                    cnt++;
                                    edgenumber.Set(e2, cnt);
                                    go_on = 1;
                                }
                                if (used2 && !used1) {
                                    cnt++;
                                    edgenumber.Set(e1, cnt);
                                    go_on = 1;
                                }
                            }
                        }
                    } while (go_on);
                }
            } while (found);

            return cnt;
        }
    }

    void BTDefineMarkedTet(const Element & el,
            INDEX_2_CLOSED_HASHTABLE<int> & edgenumber,
            MarkedTet & mt)
    {
        for (int i = 0; i < 4; i++)
            mt.pnums[i] = el[i];

        mt.marked = 0;
        mt.flagged = 0;

        mt.incorder = 0;
        mt.order = 1;

        int val = 0;
        // find marked edge of tet:
        for (int i = 0; i < 3; i++)
            for (int j = i + 1; j < 4; j++) {
                INDEX_2 i2(mt.pnums[i], mt.pnums[j]);
                i2.Sort();
                int hval = edgenumber.Get(i2);
                if (hval > val) {
                    val = hval;
                    mt.tetedge1 = i;
                    mt.tetedge2 = j;
                }
            }

        // find marked edges of faces:
        for (int k = 0; k < 4; k++) {
            val = 0;
            for (int i = 0; i < 3; i++)
                for (int j = i + 1; j < 4; j++)
                    if (i != k && j != k) {
                        INDEX_2 i2(mt.pnums[i], mt.pnums[j]);
                        i2.Sort();
                        int hval = edgenumber.Get(i2);
                        if (hval > val) {
                            val = hval;
                            int hi = 6 - k - i - j;
                            mt.faceedges[k] = char(hi);
                        }
                    }
        }
    }

    void BTDefineMarkedPrism(const Element & el,
            INDEX_2_CLOSED_HASHTABLE<int> & edgenumber,
            MarkedPrism & mp)
    {
        if (el.GetType() == PRISM ||
                el.GetType() == PRISM12) {
            for (int i = 0; i < 6; i++)
                mp.pnums[i] = el[i];
        }
        else if (el.GetType() == PYRAMID) {
            static int map[6] = {1, 2, 5, 4, 3, 5};
            for (int i = 0; i < 6; i++)
                mp.pnums[i] = el.PNum(map[i]);
        }
        else if (el.GetType() == TET ||
                el.GetType() == TET10) {
            static int map[6] = {1, 4, 3, 2, 4, 3};
            for (int i = 0; i < 6; i++)
                mp.pnums[i] = el.PNum(map[i]);

        }
        else {
            LOG_ERROR("Define marked prism called for non-prism and non-pyramid");
        }

        mp.marked = 0;

        mp.incorder = 0;
        mp.order = 1;

        int val = 0;
        for (int i = 0; i < 2; i++)
            for (int j = i + 1; j < 3; j++) {
                INDEX_2 i2(mp.pnums[i], mp.pnums[j]);
                i2.Sort();
                int hval = edgenumber.Get(i2);
                if (hval > val) {
                    val = hval;
                    mp.markededge = 3 - i - j;
                }
            }
    }

    bool BTDefineMarkedId(const Element2d & el,
            INDEX_2_CLOSED_HASHTABLE<int> & edgenumber,
            const Array<int, PointIndex::BASE> & idmap,
            MarkedIdentification & mi)
    {

        bool identified = true;
        mi.np = el.GetNP();
        int min1(0), min2(0);
        for (int j = 0; identified && j < mi.np; j++) {
            mi.pnums[j] = el[j];
            mi.pnums[j + mi.np] = idmap[el[j]];

            if (j == 0 || el[j] < min1)
                min1 = el[j];
            if (j == 0 || mi.pnums[j + mi.np] < min2)
                min2 = mi.pnums[j + mi.np];

            identified = (mi.pnums[j + mi.np] != 0 && mi.pnums[j + mi.np] != mi.pnums[j]);
        }

        identified = identified && (min1 < min2);

        if (identified) {
            mi.marked = 0;

            mi.incorder = 0;
            mi.order = 1;

            int val = 0;
            for (int i = 0; i < mi.np; i++) {
                INDEX_2 i2(mi.pnums[i], mi.pnums[(i + 1) % mi.np]);
                i2.Sort();
                int hval = edgenumber.Get(i2);
                if (hval > val) {
                    val = hval;
                    mi.markededge = i;
                }
            }
        }

        return identified;
    }

    void BTDefineMarkedTri(const Element2d & el,
            INDEX_2_CLOSED_HASHTABLE<int> & edgenumber,
            MarkedTri & mt)
    {
        for (int i = 0; i < 3; i++) {
            mt.pnums[i] = el[i];
            mt.pgeominfo[i] = el.GeomInfoPi(i + 1);
        }

        mt.marked = 0;
        mt.surfid = el.GetIndex();

        mt.incorder = 0;
        mt.order = 1;

        int val = 0;
        for (int i = 0; i < 2; i++)
            for (int j = i + 1; j < 3; j++) {
                INDEX_2 i2(mt.pnums[i], mt.pnums[j]);
                i2.Sort();
                int hval = edgenumber.Get(i2);
                if (hval > val) {
                    val = hval;
                    mt.markededge = 3 - i - j;
                }
            }
    }

    void PrettyPrint(std::ostream & ost, const MarkedTri & mt)
    {
        ost << "MarkedTrig: " << std::endl;
        ost << "  pnums = ";
        for (int i = 0; i < 3; i++) ost << mt.pnums[i] << " ";
        ost << std::endl;
        ost << "  marked = " << mt.marked << ", markededge=" << mt.markededge << std::endl;
        for (int i = 0; i < 2; i++)
            for (int j = i + 1; j < 3; j++)
                if (mt.markededge == 3 - i - j)
                    ost << "  marked edge pnums = " << mt.pnums[i] << " " << mt.pnums[j] << std::endl;
    }

    void PrettyPrint(std::ostream & ost, const MarkedQuad & mq)
    {
        ost << "MarkedQuad: " << std::endl;
        ost << "  pnums = ";
        for (int i = 0; i < 4; i++) ost << mq.pnums[i] << " ";
        ost << std::endl;
        ost << "  marked = " << mq.marked << ", markededge=" << mq.markededge << std::endl;
    }

    void BTDefineMarkedQuad(const Element2d & el,
            INDEX_2_CLOSED_HASHTABLE<int> & edgenumber,
            MarkedQuad & mq)
    {
        for (int i = 0; i < 4; i++)
            mq.pnums[i] = el[i];
        std::swap(mq.pnums[2], mq.pnums[3]);

        mq.marked = 0;
        mq.markededge = 0;
        mq.surfid = el.GetIndex();
    }

    // mark elements due to local h

    int BTMarkTets(T_MTETS & mtets,
            T_MPRISMS & mprisms,
            const Mesh & mesh)
    {
        int marked = 0;

        int np = mesh.GetNP();
        Vector hv(np);
        for (int i = 0; i < np; i++)
            hv(i) = mesh.GetH(mesh.Point(i + 1));

        double hfac = 1;

        for (int step = 1; step <= 2; step++) {
            for (int i = 1; i <= mtets.size(); i++) {
                double h = 0;

                for (int j = 0; j < 3; j++)
                    for (int k = j + 1; k < 4; k++) {
                        const Point<3> & p1 = mesh.Point(mtets.Get(i).pnums[j]);
                        const Point<3> & p2 = mesh.Point(mtets.Get(i).pnums[k]);
                        double hh = Dist2(p1, p2);
                        if (hh > h) h = hh;
                    }
                h = sqrt(h);

                double hshould = 1e10;
                for (int j = 0; j < 4; j++) {
                    double hi = hv(mtets.Get(i).pnums[j] - 1);
                    if (hi < hshould)
                        hshould = hi;
                }

                if (step == 1) {
                    if (h / hshould > hfac)
                        hfac = h / hshould;
                }
                else {
                    if (h > hshould * hfac) {
                        mtets.Elem(i).marked = 1;
                        marked = 1;
                    }
                    else
                        mtets.Elem(i).marked = 0;
                }

            }
            for (int i = 1; i <= mprisms.size(); i++) {
                double h = 0;

                for (int j = 0; j < 2; j++)
                    for (int k = j + 1; k < 3; k++) {
                        const Point<3> & p1 = mesh.Point(mprisms.Get(i).pnums[j]);
                        const Point<3> & p2 = mesh.Point(mprisms.Get(i).pnums[k]);
                        double hh = Dist2(p1, p2);
                        if (hh > h) h = hh;
                    }
                h = sqrt(h);

                double hshould = 1e10;
                for (int j = 0; j < 6; j++) {
                    double hi = hv(mprisms.Get(i).pnums[j] - 1);
                    if (hi < hshould)
                        hshould = hi;
                }

                if (step == 1) {
                    if (h / hshould > hfac)
                        hfac = h / hshould;
                }
                else {
                    if (h > hshould * hfac) {
                        mprisms.Elem(i).marked = 1;
                        marked = 1;
                    }
                    else
                        mprisms.Elem(i).marked = 0;
                }

            }

            if (step == 1) {
                if (hfac > 2)
                    hfac /= 2;
                else
                    hfac = 1;
            }
        }
        return marked;
    }

    void BTBisectTet(const MarkedTet & oldtet, int newp,
            MarkedTet & newtet1, MarkedTet & newtet2)
    {
        // points vis a vis from tet-edge
        int vis1, vis2;
        vis1 = 0;
        while (vis1 == oldtet.tetedge1 || vis1 == oldtet.tetedge2)
            vis1++;
        vis2 = 6 - vis1 - oldtet.tetedge1 - oldtet.tetedge2;

        // is tet of type P ?
        int istypep = 0;
        for (int i = 0; i < 4; i++) {
            int cnt = 0;
            for (int j = 0; j < 4; j++)
                if (oldtet.faceedges[j] == i)
                    cnt++;
            if (cnt == 3)
                istypep = 1;
        }

        for (int i = 0; i < 4; i++) {
            newtet1.pnums[i] = oldtet.pnums[i];
            newtet2.pnums[i] = oldtet.pnums[i];
        }
        newtet1.flagged = istypep && !oldtet.flagged;
        newtet2.flagged = istypep && !oldtet.flagged;

        int nm = oldtet.marked - 1;
        if (nm < 0) nm = 0;
        newtet1.marked = nm;
        newtet2.marked = nm;

        for (int i = 0; i < 4; i++) {
            if (i == oldtet.tetedge1) {
                newtet2.pnums[i] = newp;
                newtet2.faceedges[i] = oldtet.faceedges[i]; // inherited face
                newtet2.faceedges[vis1] = i; // cut faces
                newtet2.faceedges[vis2] = i;

                int j = 0;
                while (j == i || j == oldtet.faceedges[i])
                    j++;
                int k = 6 - i - oldtet.faceedges[i] - j;
                newtet2.tetedge1 = j; // tet-edge
                newtet2.tetedge2 = k;

                // new face:
                if (istypep && oldtet.flagged) {
                    int hi = 6 - oldtet.tetedge1 - j - k;
                    newtet2.faceedges[oldtet.tetedge2] = char(hi);
                }
                else
                    newtet2.faceedges[oldtet.tetedge2] = oldtet.tetedge1;
            }

            if (i == oldtet.tetedge2) {
                newtet1.pnums[i] = newp;
                newtet1.faceedges[i] = oldtet.faceedges[i]; // inherited face
                newtet1.faceedges[vis1] = i;
                newtet1.faceedges[vis2] = i;
                int j = 0;
                while (j == i || j == oldtet.faceedges[i])
                    j++;
                int k = 6 - i - oldtet.faceedges[i] - j;
                newtet1.tetedge1 = j;
                newtet1.tetedge2 = k;

                // new face:
                if (istypep && oldtet.flagged) {
                    int hi = 6 - oldtet.tetedge2 - j - k;
                    newtet1.faceedges[oldtet.tetedge1] = char(hi);
                }
                else
                    newtet1.faceedges[oldtet.tetedge1] = oldtet.tetedge2;
            }
        }

        newtet1.matindex = oldtet.matindex;
        newtet2.matindex = oldtet.matindex;
        newtet1.incorder = 0;
        newtet1.order = oldtet.order;
        newtet2.incorder = 0;
        newtet2.order = oldtet.order;

        std::cerr << "newtet1 =  " << newtet1 << std::endl;
        std::cerr << "newtet2 =  " << newtet2 << std::endl;
    }

    void BTBisectPrism(const MarkedPrism & oldprism, int newp1, int newp2,
            MarkedPrism & newprism1, MarkedPrism & newprism2)
    {
        for (int i = 0; i < 6; i++) {
            newprism1.pnums[i] = oldprism.pnums[i];
            newprism2.pnums[i] = oldprism.pnums[i];
        }

        int pe1 = 0;
        if (pe1 == oldprism.markededge)
            pe1++;
        int pe2 = 3 - oldprism.markededge - pe1;

        newprism1.pnums[pe2] = newp1;
        newprism1.pnums[pe2 + 3] = newp2;
        newprism1.markededge = pe2;
        newprism2.pnums[pe1] = newp1;
        newprism2.pnums[pe1 + 3] = newp2;
        newprism2.markededge = pe1;

        newprism1.matindex = oldprism.matindex;
        newprism2.matindex = oldprism.matindex;

        int nm = oldprism.marked - 1;
        if (nm < 0) nm = 0;
        newprism1.marked = nm;
        newprism2.marked = nm;

        newprism1.incorder = 0;
        newprism1.order = oldprism.order;
        newprism2.incorder = 0;
        newprism2.order = oldprism.order;
    }

    void BTBisectIdentification(const MarkedIdentification & oldid,
            Array<PointIndex> & newp,
            MarkedIdentification & newid1,
            MarkedIdentification & newid2)
    {
        for (int i = 0; i < 2 * oldid.np; i++) {
            newid1.pnums[i] = oldid.pnums[i];
            newid2.pnums[i] = oldid.pnums[i];
        }
        newid1.np = newid2.np = oldid.np;

        if (oldid.np == 3) {
            newid1.pnums[(oldid.markededge + 1) % 3] = newp[0];
            newid1.pnums[(oldid.markededge + 1) % 3 + 3] = newp[1];
            newid1.markededge = (oldid.markededge + 2) % 3;

            newid2.pnums[oldid.markededge] = newp[0];
            newid2.pnums[oldid.markededge + 3] = newp[1];
            newid2.markededge = (oldid.markededge + 1) % 3;
        }
        else if (oldid.np == 4) {
            newid1.pnums[(oldid.markededge + 1) % 4] = newp[0];
            newid1.pnums[(oldid.markededge + 2) % 4] = newp[2];
            newid1.pnums[(oldid.markededge + 1) % 4 + 4] = newp[1];
            newid1.pnums[(oldid.markededge + 2) % 4 + 4] = newp[3];
            newid1.markededge = (oldid.markededge + 3) % 4;

            newid2.pnums[oldid.markededge] = newp[0];
            newid2.pnums[(oldid.markededge + 3) % 4] = newp[2];
            newid2.pnums[oldid.markededge + 4] = newp[1];
            newid2.pnums[(oldid.markededge + 3) % 4 + 4] = newp[3];
            newid2.markededge = (oldid.markededge + 1) % 4;
        }

        int nm = oldid.marked - 1;
        if (nm < 0) nm = 0;
        newid1.marked = newid2.marked = nm;

        newid1.incorder = newid2.incorder = 0;
        newid1.order = newid2.order = oldid.order;
    }

    void BTBisectTri(const MarkedTri & oldtri, int newp, const PointGeomInfo & newpgi,
            MarkedTri & newtri1, MarkedTri & newtri2)
    {
        for (int i = 0; i < 3; i++) {
            newtri1.pnums[i] = oldtri.pnums[i];
            newtri1.pgeominfo[i] = oldtri.pgeominfo[i];
            newtri2.pnums[i] = oldtri.pnums[i];
            newtri2.pgeominfo[i] = oldtri.pgeominfo[i];
        }

        int pe1 = 0;
        if (pe1 == oldtri.markededge)
            pe1++;
        int pe2 = 3 - oldtri.markededge - pe1;

        newtri1.pnums[pe2] = newp;
        newtri1.pgeominfo[pe2] = newpgi;
        newtri1.markededge = pe2;

        newtri2.pnums[pe1] = newp;
        newtri2.pgeominfo[pe1] = newpgi;
        newtri2.markededge = pe1;

        newtri1.surfid = oldtri.surfid;
        newtri2.surfid = oldtri.surfid;

        int nm = oldtri.marked - 1;
        if (nm < 0) nm = 0;
        newtri1.marked = nm;
        newtri2.marked = nm;

        newtri1.incorder = 0;
        newtri1.order = oldtri.order;
        newtri2.incorder = 0;
        newtri2.order = oldtri.order;
    }

    void BTBisectQuad(const MarkedQuad & oldquad,
            int newp1, const PointGeomInfo & npgi1,
            int newp2, const PointGeomInfo & npgi2,
            MarkedQuad & newquad1, MarkedQuad & newquad2)
    {
        for (int i = 0; i < 4; i++) {
            newquad1.pnums[i] = oldquad.pnums[i];
            newquad1.pgeominfo[i] = oldquad.pgeominfo[i];
            newquad2.pnums[i] = oldquad.pnums[i];
            newquad2.pgeominfo[i] = oldquad.pgeominfo[i];
        }

        if (oldquad.markededge == 0 || oldquad.markededge == 2) {
            newquad1.pnums[1] = newp1;
            newquad1.pgeominfo[1] = npgi1;
            newquad1.pnums[3] = newp2;
            newquad1.pgeominfo[3] = npgi2;

            newquad2.pnums[0] = newp1;
            newquad2.pgeominfo[0] = npgi1;
            newquad2.pnums[2] = newp2;
            newquad2.pgeominfo[2] = npgi2;
        }
        else // 1 || 3 
        {
            newquad1.pnums[2] = newp1;
            newquad1.pgeominfo[2] = npgi1;
            newquad1.pnums[3] = newp2;
            newquad1.pgeominfo[3] = npgi2;

            newquad2.pnums[0] = newp1;
            newquad2.pgeominfo[0] = npgi1;
            newquad2.pnums[1] = newp2;
            newquad2.pgeominfo[1] = npgi2;
        }
        newquad1.surfid = oldquad.surfid;
        newquad2.surfid = oldquad.surfid;

        int nm = oldquad.marked - 1;
        if (nm < 0) nm = 0;

        newquad1.marked = nm;
        newquad2.marked = nm;

        if (nm == 1) {
            newquad1.markededge = 1;
            newquad2.markededge = 1;
        }
        else {
            newquad1.markededge = 0;
            newquad2.markededge = 0;
        }

    }

    int MarkHangingIdentifications(T_MIDS & mids,
            const INDEX_2_CLOSED_HASHTABLE<PointIndex> & cutedges)
    {
        int hanging = 0;
        for (int i = 1; i <= mids.size(); i++) {
            if (mids.Elem(i).marked) {
                hanging = 1;
                continue;
            }

            const int np = mids.Get(i).np;
            for (int j = 0; j < np; j++) {
                INDEX_2 edge1(mids.Get(i).pnums[j],
                        mids.Get(i).pnums[(j + 1) % np]);
                INDEX_2 edge2(mids.Get(i).pnums[j + np],
                        mids.Get(i).pnums[((j + 1) % np) + np]);

                edge1.Sort();
                edge2.Sort();
                if (cutedges.Used(edge1) ||
                        cutedges.Used(edge2)) {
                    mids.Elem(i).marked = 1;
                    hanging = 1;
                }
            }
        }

        return hanging;
    }

    int MarkHangingTets(T_MTETS & mtets,
            const INDEX_2_CLOSED_HASHTABLE<PointIndex> & cutedges)
    {
        int hanging = 0;
        for (int i = 1; i <= mtets.size(); i++) {
            MarkedTet & teti = mtets.Elem(i);

            if (teti.marked) {
                hanging = 1;
                continue;
            }

            for (int j = 0; j < 3; j++)
                for (int k = j + 1; k < 4; k++) {
                    INDEX_2 edge(teti.pnums[j],
                            teti.pnums[k]);
                    edge.Sort();
                    if (cutedges.Used(edge)) {
                        teti.marked = 1;
                        hanging = 1;
                    }
                }
        }

        return hanging;
    }

    int MarkHangingPrisms(T_MPRISMS & mprisms,
            const INDEX_2_CLOSED_HASHTABLE<PointIndex> & cutedges)
    {
        int hanging = 0;
        for (int i = 1; i <= mprisms.size(); i++) {
            if (mprisms.Elem(i).marked) {
                hanging = 1;
                continue;
            }

            for (int j = 0; j < 2; j++)
                for (int k = j + 1; k < 3; k++) {
                    INDEX_2 edge1(mprisms.Get(i).pnums[j],
                            mprisms.Get(i).pnums[k]);
                    INDEX_2 edge2(mprisms.Get(i).pnums[j + 3],
                            mprisms.Get(i).pnums[k + 3]);
                    edge1.Sort();
                    edge2.Sort();
                    if (cutedges.Used(edge1) ||
                            cutedges.Used(edge2)) {
                        mprisms.Elem(i).marked = 1;
                        hanging = 1;
                    }
                }
        }
        return hanging;
    }

    int MarkHangingTris(T_MTRIS & mtris,
            const INDEX_2_CLOSED_HASHTABLE<PointIndex> & cutedges)
    {
        int hanging = 0;
        for (int i = 1; i <= mtris.size(); i++) {
            if (mtris.Get(i).marked) {
                hanging = 1;
                continue;
            }
            for (int j = 0; j < 2; j++)
                for (int k = j + 1; k < 3; k++) {
                    INDEX_2 edge(mtris.Get(i).pnums[j],
                            mtris.Get(i).pnums[k]);
                    edge.Sort();
                    if (cutedges.Used(edge)) {
                        mtris.Elem(i).marked = 1;
                        hanging = 1;
                    }
                }
        }
        return hanging;
    }

    int MarkHangingQuads(T_MQUADS & mquads,
            const INDEX_2_CLOSED_HASHTABLE<PointIndex> & cutedges)
    {
        int hanging = 0;
        for (int i = 1; i <= mquads.size(); i++) {
            if (mquads.Elem(i).marked) {
                hanging = 1;
                continue;
            }

            INDEX_2 edge1(mquads.Get(i).pnums[0],
                    mquads.Get(i).pnums[1]);
            INDEX_2 edge2(mquads.Get(i).pnums[2],
                    mquads.Get(i).pnums[3]);
            edge1.Sort();
            edge2.Sort();
            if (cutedges.Used(edge1) ||
                    cutedges.Used(edge2)) {
                mquads.Elem(i).marked = 1;
                mquads.Elem(i).markededge = 0;
                hanging = 1;
                continue;
            }

            // he/sz: second case: split horizontally
            INDEX_2 edge3(mquads.Get(i).pnums[1],
                    mquads.Get(i).pnums[3]);
            INDEX_2 edge4(mquads.Get(i).pnums[2],
                    mquads.Get(i).pnums[0]);

            edge3.Sort();
            edge4.Sort();
            if (cutedges.Used(edge3) ||
                    cutedges.Used(edge4)) {
                mquads.Elem(i).marked = 1;
                mquads.Elem(i).markededge = 1;
                hanging = 1;
                continue;
            }

        }
        return hanging;
    }

    void ConnectToNodeRec(int node, int tonode,
            const TABLE<int> & conto, Array<int> & connecttonode)
    {
        //  std::cerr << "connect " << node << " to " << tonode <<std::endl;
        for (int i = 1; i <= conto.EntrySize(node); i++) {
            int n2 = conto.Get(node, i);
            if (!connecttonode.Get(n2)) {
                connecttonode.Elem(n2) = tonode;
                ConnectToNodeRec(n2, tonode, conto, connecttonode);
            }
        }
    }

    T_MTETS mtets;
    T_MPRISMS mprisms;
    T_MIDS mids;
    T_MTRIS mtris;
    T_MQUADS mquads;

    void WriteMarkedElements(std::ostream & ost)
    {
        ost << "Marked Elements\n";

        ost << mtets.size() << "\n";
        for (int i = 0; i < mtets.size(); i++)
            ost << mtets[i];

        ost << mprisms.size() << "\n";
        for (int i = 0; i < mprisms.size(); i++)
            ost << mprisms[i];

        ost << mids.size() << "\n";
        for (int i = 0; i < mids.size(); i++)
            ost << mids[i];

        ost << mtris.size() << "\n";
        for (int i = 0; i < mtris.size(); i++)
            ost << mtris[i];

        ost << mquads.size() << "\n";
        for (int i = 0; i < mquads.size(); i++)
            ost << mquads[i];
        ost << std::endl;
    }

    bool ReadMarkedElements(std::istream & ist, const Mesh & mesh)
    {
        std::string auxstring("");
        if (ist)
            ist >> auxstring;

        if (auxstring != "Marked")
            return false;

        if (ist)
            ist >> auxstring;

        if (auxstring != "Elements")
            return false;

        int size;

        ist >> size;
        mtets.resize(size);
        for (int i = 0; i < size; i++) {
            ist >> mtets[i];
            if (mtets[i].pnums[0] > mesh.GetNV() ||
                    mtets[i].pnums[1] > mesh.GetNV() ||
                    mtets[i].pnums[2] > mesh.GetNV() ||
                    mtets[i].pnums[3] > mesh.GetNV())
                return false;
        }

        ist >> size;
        mprisms.resize(size);
        for (int i = 0; i < size; i++)
            ist >> mprisms[i];

        ist >> size;
        mids.resize(size);
        for (int i = 0; i < size; i++)
            ist >> mids[i];

        ist >> size;
        mtris.resize(size);
        for (int i = 0; i < size; i++)
            ist >> mtris[i];

        ist >> size;
        mquads.resize(size);
        for (int i = 0; i < size; i++)
            ist >> mquads[i];

        return true;
    }

    void BisectTetsCopyMesh(
            Mesh & mesh, const class CSGeometry *,
            BisectionOptions & opt,
            const Array< Array<int, PointIndex::BASE>* > & idmaps,
            const std::string & refinfofile)
    {
        //int np = mesh.GetNP();
        int ne = mesh.GetNE();
        int nse = mesh.GetNSE();

        bool readok = false;

        if (refinfofile != "") {
            LOG_DEBUG("Reading marked-element information from \"" << refinfofile << "\"");
            std::ifstream ist(refinfofile.c_str());

            readok = ReadMarkedElements(ist, mesh);

            ist.close();
        }

        if (!readok) {
            LOG_DEBUG("resetting marked-element information");
            mtets.resize(0);
            mprisms.resize(0);
            mids.resize(0);
            mtris.resize(0);
            mquads.resize(0);

            INDEX_2_HASHTABLE<int> shortedges(100);
            for (int i = 1; i <= ne; i++) {
                const Element & el = mesh.VolumeElement(i);
                if (el.GetType() == PRISM ||
                        el.GetType() == PRISM12) {
                    for (int j = 1; j <= 3; j++) {
                        INDEX_2 se(el.PNum(j), el.PNum(j + 3));
                        se.Sort();
                        shortedges.Set(se, 1);
                    }
                }
            }

            // INDEX_2_HASHTABLE<int> edgenumber(np);
            INDEX_2_CLOSED_HASHTABLE<int> edgenumber(9 * ne + 4 * nse);

            BTSortEdges(mesh, idmaps, edgenumber);

            for (int i = 1; i <= ne; i++) {
                const Element & el = mesh.VolumeElement(i);

                switch (el.GetType()) {
                    case TET:
                    case TET10:
                    {
                        // if tet has short edge, it is handled as degenerated prism

                        int foundse = 0;
                        for (int j = 1; j <= 3; j++)
                            for (int k = j + 1; k <= 4; k++) {
                                INDEX_2 se(el.PNum(j), el.PNum(k));
                                se.Sort();
                                if (shortedges.Used(se)) {
                                    //		      std::cout << "tet converted to prism" <<std::endl;

                                    foundse = 1;
                                    int p3 = 1;
                                    while (p3 == j || p3 == k)
                                        p3++;
                                    int p4 = 10 - j - k - p3;

                                    // even permutation ?
                                    int pi[4];
                                    pi[0] = j;
                                    pi[1] = k;
                                    pi[2] = p3;
                                    pi[3] = p4;
                                    int cnt = 0;
                                    for (int l = 1; l <= 4; l++)
                                        for (int m = 0; m < 3; m++)
                                            if (pi[m] > pi[m + 1]) {
                                                std::swap(pi[m], pi[m + 1]);
                                                cnt++;
                                            }
                                    if (cnt % 2)
                                        std::swap(p3, p4);

                                    Element hel = el;
                                    hel.PNum(1) = el.PNum(j);
                                    hel.PNum(2) = el.PNum(k);
                                    hel.PNum(3) = el.PNum(p3);
                                    hel.PNum(4) = el.PNum(p4);

                                    MarkedPrism mp;
                                    BTDefineMarkedPrism(hel, edgenumber, mp);
                                    mp.matindex = el.GetIndex();
                                    mprisms.push_back(mp);
                                }
                            }
                        if (!foundse) {
                            MarkedTet mt;
                            BTDefineMarkedTet(el, edgenumber, mt);
                            mt.matindex = el.GetIndex();
                            mtets.push_back(mt);
                        }
                        break;
                    }
                    case PYRAMID:
                    {
                        // eventually rotate
                        MarkedPrism mp;

                        INDEX_2 se(el.PNum(1), el.PNum(2));
                        se.Sort();
                        if (shortedges.Used(se)) {
                            Element hel = el;
                            hel.PNum(1) = el.PNum(2);
                            hel.PNum(2) = el.PNum(3);
                            hel.PNum(3) = el.PNum(4);
                            hel.PNum(4) = el.PNum(1);
                            BTDefineMarkedPrism(hel, edgenumber, mp);
                        }
                        else {
                            BTDefineMarkedPrism(el, edgenumber, mp);
                        }

                        mp.matindex = el.GetIndex();
                        mprisms.push_back(mp);
                        break;
                    }
                    case PRISM:
                    case PRISM12:
                    {
                        MarkedPrism mp;
                        BTDefineMarkedPrism(el, edgenumber, mp);
                        mp.matindex = el.GetIndex();
                        mprisms.push_back(mp);
                        break;
                    }
                    default:
                        throw std::runtime_error("Bisect, element type not handled in switch, 4");
                }
            }

            for (int i = 1; i <= nse; i++) {
                const Element2d & el = mesh.SurfaceElement(i);
                if (el.GetType() == TRIG ||
                        el.GetType() == TRIG6) {
                    MarkedTri mt;
                    BTDefineMarkedTri(el, edgenumber, mt);
                    mtris.push_back(mt);
                }
                else {
                    MarkedQuad mq;
                    BTDefineMarkedQuad(el, edgenumber, mq);
                    mquads.push_back(mq);
                }

                MarkedIdentification mi;
                for (int j = 0; j < idmaps.size(); j++)
                    if (BTDefineMarkedId(el, edgenumber, *idmaps[j], mi))
                        mids.push_back(mi);
            }
        }

        mesh.mlparentelement.resize(ne);
        for (int i = 1; i <= ne; i++) {
            mesh.mlparentelement.Elem(i) = 0;
        }
        mesh.mlparentsurfaceelement.resize(nse);
        for (int i = 1; i <= nse; i++) {
            mesh.mlparentsurfaceelement.Elem(i) = 0;
        }

        LOG_DEBUG("copied " << mtets.size() << " tets, " << mprisms.size() << " prisms");
        LOG_DEBUG("       " << mtris.size() << " trigs, " << mquads.size() << " quads");
    }

    void UpdateEdgeMarks(Mesh & mesh,
            const Array< Array<int, PointIndex::BASE>* > & idmaps)
    {
        T_MTETS mtets_old(mtets);
        T_MPRISMS mprisms_old(mprisms);
        T_MIDS mids_old(mids);
        T_MTRIS mtris_old(mtris);
        T_MQUADS mquads_old(mquads);

        mtets.resize(0);
        mprisms.resize(0);
        mids.resize(0);
        mtris.resize(0);
        mquads.resize(0);

        INDEX_2_CLOSED_HASHTABLE<int> edgenumber(9 * mesh.GetNE() + 4 * mesh.GetNSE());

        int maxnum = BTSortEdges(mesh, idmaps, edgenumber);

        for (int m = 0; m < mtets_old.size(); m++) {
            MarkedTet & mt = mtets_old[m];

            INDEX_2 edge(mt.pnums[mt.tetedge1], mt.pnums[mt.tetedge2]);
            edge.Sort();
            if (edgenumber.Used(edge)) {
                int val = edgenumber.Get(edge);
                if (val <= maxnum) {
                    val += 2 * maxnum;
                    edgenumber.Set(edge, val);
                }
                else if (val <= 2 * maxnum) {
                    val += maxnum;
                    edgenumber.Set(edge, val);
                }
            }

            for (int k = 0; k < 4; k++)
                for (int i = 0; i < 3; i++)
                    for (int j = i + 1; i != k && j < 4; j++)
                        if (j != k && int(mt.faceedges[k]) == 6 - k - i - j) {
                            edge[0] = mt.pnums[i];
                            edge[1] = mt.pnums[j];
                            edge.Sort();
                            if (edgenumber.Used(edge)) {
                                int val = edgenumber.Get(edge);
                                if (val <= maxnum) {
                                    val += maxnum;
                                    edgenumber.Set(edge, val);
                                }
                            }
                        }
        }

        for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++) {

            const Element & el = mesh[ei];

            switch (el.GetType()) {
                case TET:
                case TET10:
                {
                    MarkedTet mt;
                    BTDefineMarkedTet(el, edgenumber, mt);
                    mt.matindex = el.GetIndex();
                    mtets.push_back(mt);

                    break;
                }
                case PYRAMID:
                {
                    std::cerr << "Refinement :: UpdateEdgeMarks not yet implemented for pyramids"
                            << std::endl;
                    break;
                }

                case PRISM:
                case PRISM12:
                {
                    std::cerr << "Refinement :: UpdateEdgeMarks not yet implemented for prisms"
                            << std::endl;
                    break;
                }
                default:
                    throw std::runtime_error("Bisect, element type not handled in switch, 6");
            }

        }

        for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++) {
            const Element2d & el = mesh.SurfaceElement(sei);

            switch (el.GetType()) {
                case TRIG:
                case TRIG6:
                {
                    MarkedTri mt;
                    BTDefineMarkedTri(el, edgenumber, mt);
                    mtris.push_back(mt);
                    break;
                }

                case QUAD:
                case QUAD6:
                {
                    MarkedQuad mt;
                    BTDefineMarkedQuad(el, edgenumber, mt);
                    mquads.push_back(mt);
                    break;
                }
                default:
                    throw std::runtime_error("Bisect, element type not handled in switch, 5");
            }

            MarkedIdentification mi;
            for (int j = 0; j < idmaps.size(); j++)
                if (BTDefineMarkedId(el, edgenumber, *idmaps[j], mi))
                    mids.push_back(mi);

        }

    }
    void ReadEnclString(std::istream & in, std::string & str, const char encl)
    {
        char currchar;
        str = "";

        in.get(currchar);
        while (in && (currchar == ' ' || currchar == '\t' || currchar == '\n')) {
            in.get(currchar);
        }
        if (currchar == encl) {
            in.get(currchar);
            while (in && currchar != encl) {
                str += currchar;
                in.get(currchar);
            }
        }
        else {
            in.putback(currchar);
            in >> str;
        }
    }
    
    void Refinement::Bisect(Mesh & mesh, BisectionOptions & opt, Array<double> * quality_loss) const
    {
        LOG_DEBUG("Mesh bisection");

        LocalizeEdgePoints(mesh);

        Array< Array<int, PointIndex::BASE>* > idmaps;
        for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++) {
            if (mesh.GetIdentifications().GetType(i) == Identifications::PERIODIC) {
                idmaps.push_back(new Array<int, PointIndex::BASE>);
                mesh.GetIdentifications().GetMap(i, *idmaps.Last(), true);
            }
        }

        std::string refelementinfofileread = "";
        std::string refelementinfofilewrite = "";

        if (opt.refinementfilename) {
            std::ifstream inf(opt.refinementfilename);
            std::string st;
            inf >> st;
            if (st == "refinementinfo") {
                while (inf) {
                    while (inf && st != "markedelementsfile")
                        inf >> st;

                    if (inf)
                        inf >> st;

                    if (st == "read" && inf)
                        ReadEnclString(inf, refelementinfofileread, '\"');
                    else if (st == "write" && inf)
                        ReadEnclString(inf, refelementinfofilewrite, '\"');
                }
            }
            inf.close();
        }

        if (mesh.mglevels == 1 || idmaps.size() > 0)
            BisectTetsCopyMesh(mesh, NULL, opt, idmaps, refelementinfofileread);

        mesh.ComputeNVertices();

        int np = mesh.GetNV();
        mesh.SetNP(np);
        mesh.mglevels++;

        if (opt.refine_p) {
            int ne = mesh.GetNE();
            int nse = mesh.GetNSE();
            int ox, oy, oz;
            for (ElementIndex ei = 0; ei < ne; ei++)
                if (mesh[ei].TestRefinementFlag()) {
                    mesh[ei].GetOrder(ox, oy, oz);
                    mesh[ei].SetOrder(ox + 1, oy + 1, oz + 1);
                    if (mesh[ei].TestStrongRefinementFlag())
                        mesh[ei].SetOrder(ox + 2, oy + 2, oz + 2);
                }
            for (SurfaceElementIndex sei = 0; sei < nse; sei++) {
                Element2d & el = mesh.SurfaceElement(sei);
                if (el.TestRefinementFlag()) {
                    el.GetOrder(ox, oy);
                    el.SetOrder(ox + 1, oy + 1);
                    if (el.TestStrongRefinementFlag()) {
                        el.SetOrder(ox + 2, oy + 2);
                    }
                }
            }

#ifndef SABINE //Nachbarelemente mit ordx,ordy,ordz 

            Array<int, PointIndex::BASE> v_order(mesh.GetNP());
            v_order = 0;

            for (ElementIndex ei = 0; ei < ne; ei++)
                for (int j = 0; j < mesh[ei].GetNP(); j++)
                    if (mesh[ei].GetOrder() > v_order[mesh[ei][j]])
                        v_order[mesh[ei][j]] = mesh[ei].GetOrder();

            for (SurfaceElementIndex sei = 0; sei < nse; sei++) {
                const Element2d & el = mesh.SurfaceElement(sei);
                for (int j = 0; j < el.GetNP(); j++) {
                    if (el.GetOrder() > v_order[el[j]]) {
                        v_order[el[j]] = el.GetOrder();
                    }
                }
            }
            for (ElementIndex ei = 0; ei < ne; ei++) {
                for (int j = 0; j < mesh[ei].GetNP(); j++) {
                    if (mesh[ei].GetOrder() < v_order[mesh[ei][j]] - 1) {
                        mesh[ei].SetOrder(v_order[mesh[ei][j]] - 1);
                    }
                }
            }

            for (SurfaceElementIndex sei = 0; sei < nse; sei++) {
                Element2d & el = mesh.SurfaceElement(sei);
                for (int j = 0; j < el.GetNP(); j++) {
                    if (el.GetOrder() < v_order[el[j]] - 1) {
                        el.SetOrder(v_order[el[j]] - 1);
                    }
                }
            }
#endif
            return;
        }

        INDEX_2_CLOSED_HASHTABLE<PointIndex> cutedges(10 + 9 * (mtets.size() + mprisms.size() + mtris.size() + mquads.size()));

        bool noprojection = false;

        for (int l = 1; l <= 1; l++) {
            int marked = 0;
            if (opt.refinementfilename) {
                std::ifstream inf(opt.refinementfilename);
                LOG_DEBUG("load refinementinfo from file " << opt.refinementfilename);

                std::string st;
                inf >> st;
                if (st == "refinementinfo")
                    // new version
                {
                    for (int i = 1; i <= mtets.size(); i++)
                        mtets.Elem(i).marked = 0;
                    for (int i = 1; i <= mprisms.size(); i++)
                        mprisms.Elem(i).marked = 0;
                    for (int i = 1; i <= mtris.size(); i++)
                        mtris.Elem(i).marked = 0;
                    for (int i = 1; i <= mquads.size(); i++)
                        mquads.Elem(i).marked = 0;
                    for (int i = 1; i <= mprisms.size(); i++)
                        mids.Elem(i).marked = 0;

                    inf >> st;
                    while (inf) {
                        if (st[0] == '#') {
                            inf.ignore(10000, '\n');
                            inf >> st;
                        }
                        else if (st == "markedelementsfile") {
                            inf >> st;
                            ReadEnclString(inf, st, '\"');
                            inf >> st;
                        }
                        else if (st == "noprojection") {
                            noprojection = true;
                            inf >> st;
                        }
                        else if (st == "refine") {
                            inf >> st;
                            if (st == "elements") {
                                inf >> st;
                                bool isint = true;
                                for (std::string::size_type ii = 0; isint && ii < st.size(); ii++)
                                    isint = (isdigit(st[ii]) != 0);

                                while (inf && isint) {
                                    mtets.Elem(atoi(st.c_str())).marked = 3;
                                    marked = 1;

                                    inf >> st;
                                    isint = true;
                                    for (std::string::size_type ii = 0; isint && ii < st.size(); ii++)
                                        isint = (isdigit(st[ii]) != 0);
                                }
                            }
                            else if (st == "orthobrick") {
                                double bounds[6];
                                for (int i = 0; i < 6; i++)
                                    inf >> bounds[i];

                                int cnt = 0;

                                for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++) {
                                    const Element & el = mesh[ei];

                                    //
                                    Point<3> center(0, 0, 0);
                                    for (int i = 0; i < el.GetNP(); i++) {
                                        const MeshPoint & point = mesh[el[i]];
                                        center(0) += point(0);
                                        center(1) += point(1);
                                        center(2) += point(2);
                                    }
                                    for (int i = 0; i < 3; i++)
                                        center(i) *= 1. / double(el.GetNP());
                                    if (bounds[0] <= center(0) && center(0) <= bounds[3] &&
                                            bounds[1] <= center(1) && center(1) <= bounds[4] &&
                                            bounds[2] <= center(2) && center(2) <= bounds[5]) {
                                        mtets[ei].marked = 3;
                                        cnt++;
                                    }
                                }

                                std::ostringstream strstr;
                                strstr.precision(2);
                                strstr << "marked " << float(cnt) / float(mesh.GetNE())*100.
#ifdef WIN32
                                        << "%%"
#else
                                        << "%"
#endif
                                        << " of the elements";
                                LOG_DEBUG(strstr.str());

                                if (cnt > 0)
                                    marked = 1;

                                inf >> st;
                            }
                            else {
                                throw std::runtime_error("something wrong with refinementinfo file");
                            }
                        }
                    }
                }
                else {
                    inf.close();
                    inf.open(opt.refinementfilename);

                    char ch;
                    for (int i = 1; i <= mtets.size(); i++) {
                        inf >> ch;
                        if (!inf)
                            throw std::runtime_error("something wrong with refinementinfo file (old format)");
                        mtets.Elem(i).marked = (ch == '1');
                    }
                    marked = 1;
                }
                inf.close();
            }

            else if (opt.usemarkedelements) {
                int cntm = 0;

                // all in one !
                if (mprisms.size()) {
                    int cnttet = 0;
                    int cntprism = 0;
                    for (int i = 1; i <= mesh.GetNE(); i++) {
                        if (mesh.VolumeElement(i).GetType() == TET ||
                                mesh.VolumeElement(i).GetType() == TET10) {
                            cnttet++;
                            mtets.Elem(cnttet).marked =
                                    3 * mesh.VolumeElement(i).TestRefinementFlag();
                            if (mtets.Elem(cnttet).marked)
                                cntm++;
                        }
                        else {
                            cntprism++;
                            mprisms.Elem(cntprism).marked =
                                    2 * mesh.VolumeElement(i).TestRefinementFlag();
                            if (mprisms.Elem(cntprism).marked)
                                cntm++;
                        }

                    }
                }
                else
                    for (int i = 1; i <= mtets.size(); i++) {
                        mtets.Elem(i).marked = 3 * mesh.VolumeElement(i).TestRefinementFlag();
                        if (mtets.Elem(i).marked)
                            cntm++;
                    }

                LOG_DEBUG("marked elements: " << cntm);

                int cnttrig = 0;
                int cntquad = 0;
                for (int i = 1; i <= mesh.GetNSE(); i++) {
                    if (mesh.SurfaceElement(i).GetType() == TRIG ||
                            mesh.SurfaceElement(i).GetType() == TRIG6) {
                        cnttrig++;
                        mtris.Elem(cnttrig).marked =
                                mesh.SurfaceElement(i).TestRefinementFlag() ? 2 : 0;
                        // mtris.Elem(cnttrig).marked = 0;
                        if (mtris.Elem(cnttrig).marked)
                            cntm++;
                    }
                    else {
                        cntquad++;
                        // 2d: marked=2, 3d prisms: marked=1
                        mquads.Elem(cntquad).marked =
                                mesh.SurfaceElement(i).TestRefinementFlag() ? 4 - mesh.GetDimension() : 0;
                        // mquads.Elem(cntquad).marked = 0;
                        if (mquads.Elem(cntquad).marked)
                            cntm++;
                    }
                }

                LOG_DEBUG("with surface-elements: " << cntm);

                marked = (cntm > 0);
            }
            else {
                marked = BTMarkTets(mtets, mprisms, mesh);
            }

            if (!marked) break;

            if (opt.refine_p) {
                LOG_DEBUG("refine p");

                for (int i = 1; i <= mtets.size(); i++)
                    mtets.Elem(i).incorder = mtets.Elem(i).marked ? 1 : 0;

                for (int i = 1; i <= mtets.size(); i++)
                    if (mtets.Elem(i).incorder)
                        mtets.Elem(i).marked = 0;

                for (int i = 1; i <= mprisms.size(); i++)
                    mprisms.Elem(i).incorder = mprisms.Elem(i).marked ? 1 : 0;

                for (int i = 1; i <= mprisms.size(); i++)
                    if (mprisms.Elem(i).incorder)
                        mprisms.Elem(i).marked = 0;

                for (int i = 1; i <= mtris.size(); i++)
                    mtris.Elem(i).incorder = mtris.Elem(i).marked ? 1 : 0;

                for (int i = 1; i <= mtris.size(); i++) {
                    if (mtris.Elem(i).incorder)
                        mtris.Elem(i).marked = 0;
                }
            }

            if (opt.refine_hp) {
                LOG_DEBUG("refine hp");
                BitArray singv(np);
                singv.Clear();

                if (mesh.GetDimension() == 3) {
                    for (int i = 1; i <= mesh.GetNSeg(); i++) {
                        const Segment & seg = mesh.LineSegment(i);
                        singv.Set(seg[0]);
                        singv.Set(seg[1]);
                    }
                }
                else {
                    // vertices with 2 different bnds
                    Array<int> bndind(np);
                    bndind = 0;
                    for (int i = 1; i <= mesh.GetNSeg(); i++) {
                        const Segment & seg = mesh.LineSegment(i);
                        for (int j = 0; j < 2; j++) {
                            int pi = (j == 0) ? seg[0] : seg[1];
                            if (bndind.Elem(pi) == 0)
                                bndind.Elem(pi) = seg.edgenr;
                            else if (bndind.Elem(pi) != seg.edgenr)
                                singv.Set(pi);
                        }
                    }
                }

                for (int i = 1; i <= mtets.size(); i++)
                    mtets.Elem(i).incorder = 1;
                for (int i = 1; i <= mtets.size(); i++) {
                    if (!mtets.Elem(i).marked)
                        mtets.Elem(i).incorder = 0;
                    for (int j = 0; j < 4; j++)
                        if (singv.Test(mtets.Elem(i).pnums[j]))
                            mtets.Elem(i).incorder = 0;
                }
                for (int i = 1; i <= mtets.size(); i++)
                    if (mtets.Elem(i).incorder)
                        mtets.Elem(i).marked = 0;

                for (int i = 1; i <= mprisms.size(); i++)
                    mprisms.Elem(i).incorder = 1;
                for (int i = 1; i <= mprisms.size(); i++) {
                    if (!mprisms.Elem(i).marked)
                        mprisms.Elem(i).incorder = 0;
                    for (int j = 0; j < 6; j++)
                        if (singv.Test(mprisms.Elem(i).pnums[j]))
                            mprisms.Elem(i).incorder = 0;
                }
                for (int i = 1; i <= mprisms.size(); i++)
                    if (mprisms.Elem(i).incorder)
                        mprisms.Elem(i).marked = 0;

                for (int i = 1; i <= mtris.size(); i++)
                    mtris.Elem(i).incorder = 1;
                for (int i = 1; i <= mtris.size(); i++) {
                    if (!mtris.Elem(i).marked)
                        mtris.Elem(i).incorder = 0;
                    for (int j = 0; j < 3; j++)
                        if (singv.Test(mtris.Elem(i).pnums[j]))
                            mtris.Elem(i).incorder = 0;
                }
                for (int i = 1; i <= mtris.size(); i++) {
                    if (mtris.Elem(i).incorder)
                        mtris.Elem(i).marked = 0;
                }
            }

            int hangingvol, hangingsurf, hangingedge;

            std::cerr << "refine volume elements" << std::endl;
            do {
                // refine volume elements

                int nel = mtets.size();
                for (int i = 1; i <= nel; i++)
                    if (mtets.Elem(i).marked) {
                        MarkedTet oldtet;
                        MarkedTet newtet1, newtet2;
                        PointIndex newp;

                        oldtet = mtets.Get(i);

                        INDEX_2 edge(oldtet.pnums[oldtet.tetedge1],
                                oldtet.pnums[oldtet.tetedge2]);
                        edge.Sort();
                        if (cutedges.Used(edge)) {
                            newp = cutedges.Get(edge);
                        }
                        else {
                            Point<3> npt = Center(mesh.Point(edge.I1()),
                                    mesh.Point(edge.I2()));
                            newp = mesh.AddPoint(npt);
                            cutedges.Set(edge, newp);
                        }

                        BTBisectTet(oldtet, newp, newtet1, newtet2);

                        mtets.Elem(i) = newtet1;
                        mtets.push_back(newtet2);

                        mesh.mlparentelement.push_back(i);
                    }

                int npr = mprisms.size();
                for (int i = 1; i <= npr; i++)
                    if (mprisms.Elem(i).marked) {
                        MarkedPrism oldprism;
                        MarkedPrism newprism1, newprism2;
                        PointIndex newp1, newp2;

                        oldprism = mprisms.Get(i);
                        int pi1 = 0;
                        if (pi1 == oldprism.markededge)
                            pi1++;
                        int pi2 = 3 - pi1 - oldprism.markededge;

                        INDEX_2 edge1(oldprism.pnums[pi1],
                                oldprism.pnums[pi2]);
                        INDEX_2 edge2(oldprism.pnums[pi1 + 3],
                                oldprism.pnums[pi2 + 3]);
                        edge1.Sort();
                        edge2.Sort();

                        if (cutedges.Used(edge1))
                            newp1 = cutedges.Get(edge1);
                        else {
                            Point<3> npt = Center(mesh.Point(edge1.I1()),
                                    mesh.Point(edge1.I2()));
                            newp1 = mesh.AddPoint(npt);
                            cutedges.Set(edge1, newp1);
                        }
                        if (cutedges.Used(edge2))
                            newp2 = cutedges.Get(edge2);
                        else {
                            Point<3> npt = Center(mesh.Point(edge2.I1()),
                                    mesh.Point(edge2.I2()));
                            newp2 = mesh.AddPoint(npt);
                            cutedges.Set(edge2, newp2);
                        }

                        BTBisectPrism(oldprism, newp1, newp2, newprism1, newprism2);
                        mprisms.Elem(i) = newprism1;
                        mprisms.push_back(newprism2);
                    }

                int nid = mids.size();
                for (int i = 1; i <= nid; i++)
                    if (mids.Elem(i).marked) {
                        MarkedIdentification oldid, newid1, newid2;
                        Array<PointIndex> newp;

                        oldid = mids.Get(i);

                        Array<INDEX_2> edges;
                        edges.push_back(INDEX_2(oldid.pnums[oldid.markededge],
                                oldid.pnums[(oldid.markededge + 1) % oldid.np]));
                        edges.push_back(INDEX_2(oldid.pnums[oldid.markededge + oldid.np],
                                oldid.pnums[(oldid.markededge + 1) % oldid.np + oldid.np]));

                        if (oldid.np == 4) {
                            edges.push_back(INDEX_2(oldid.pnums[(oldid.markededge + 2) % oldid.np],
                                    oldid.pnums[(oldid.markededge + 3) % oldid.np]));
                            edges.push_back(INDEX_2(oldid.pnums[(oldid.markededge + 2) % oldid.np + oldid.np],
                                    oldid.pnums[(oldid.markededge + 3) % oldid.np + oldid.np]));
                        }
                        for (int j = 0; j < edges.size(); j++) {
                            edges[j].Sort();

                            if (cutedges.Used(edges[j]))
                                newp.push_back(cutedges.Get(edges[j]));
                            else {
                                Point<3> npt = Center(mesh.Point(edges[j].I1()),
                                        mesh.Point(edges[j].I2()));
                                newp.push_back(mesh.AddPoint(npt));
                                cutedges.Set(edges[j], newp[j]);
                            }
                        }

                        BTBisectIdentification(oldid, newp, newid1, newid2);
                        mids.Elem(i) = newid1;
                        mids.push_back(newid2);
                    }

                hangingvol =
                        MarkHangingTets(mtets, cutedges) +
                        MarkHangingPrisms(mprisms, cutedges) +
                        MarkHangingIdentifications(mids, cutedges);

                int nsel = mtris.size();

                for (int i = 1; i <= nsel; i++)
                    if (mtris.Elem(i).marked) {
                        MarkedTri oldtri;
                        MarkedTri newtri1, newtri2;
                        PointIndex newp;

                        oldtri = mtris.Get(i);
                        int oldpi1 = oldtri.pnums[(oldtri.markededge + 1) % 3];
                        int oldpi2 = oldtri.pnums[(oldtri.markededge + 2) % 3];
                        INDEX_2 edge(oldpi1, oldpi2);
                        edge.Sort();

                        if (cutedges.Used(edge)) {
                            newp = cutedges.Get(edge);
                        }
                        else {
                            Point<3> npt = Center(mesh.Point(edge.I1()),
                                    mesh.Point(edge.I2()));
                            newp = mesh.AddPoint(npt);
                            cutedges.Set(edge, newp);
                        }
                        int si = mesh.GetFaceDescriptor(oldtri.surfid).SurfNr();
                        PointGeomInfo npgi;

                        if (mesh[newp].Type() != EDGEPOINT)
                            PointBetween(mesh.Point(oldpi1), mesh.Point(oldpi2),
                                0.5, si,
                                oldtri.pgeominfo[(oldtri.markededge + 1) % 3],
                                oldtri.pgeominfo[(oldtri.markededge + 2) % 3],
                                mesh.Point(newp), npgi);

                        BTBisectTri(oldtri, newp, npgi, newtri1, newtri2);

                        mtris.Elem(i) = newtri1;
                        mtris.push_back(newtri2);
                        mesh.mlparentsurfaceelement.push_back(i);
                    }

                int nquad = mquads.size();
                for (int i = 1; i <= nquad; i++)
                    if (mquads.Elem(i).marked) {
                        MarkedQuad oldquad;
                        MarkedQuad newquad1, newquad2;
                        PointIndex newp1, newp2;

                        oldquad = mquads.Get(i);

                        INDEX_2 edge1, edge2;
                        PointGeomInfo pgi11, pgi12, pgi21, pgi22;
                        if (oldquad.markededge == 0 || oldquad.markededge == 2) {
                            edge1.I1() = oldquad.pnums[0];
                            pgi11 = oldquad.pgeominfo[0];
                            edge1.I2() = oldquad.pnums[1];
                            pgi12 = oldquad.pgeominfo[1];
                            edge2.I1() = oldquad.pnums[2];
                            pgi21 = oldquad.pgeominfo[2];
                            edge2.I2() = oldquad.pnums[3];
                            pgi22 = oldquad.pgeominfo[3];
                        }
                        else // 3 || 1
                        {
                            edge1.I1() = oldquad.pnums[0];
                            pgi11 = oldquad.pgeominfo[0];
                            edge1.I2() = oldquad.pnums[2];
                            pgi12 = oldquad.pgeominfo[2];
                            edge2.I1() = oldquad.pnums[1];
                            pgi21 = oldquad.pgeominfo[1];
                            edge2.I2() = oldquad.pnums[3];
                            pgi22 = oldquad.pgeominfo[3];
                        }

                        edge1.Sort();
                        edge2.Sort();

                        if (cutedges.Used(edge1)) {
                            newp1 = cutedges.Get(edge1);
                        }
                        else {
                            Point<3> np1 = Center(mesh.Point(edge1.I1()),
                                    mesh.Point(edge1.I2()));
                            newp1 = mesh.AddPoint(np1);
                            cutedges.Set(edge1, newp1);
                        }

                        if (cutedges.Used(edge2)) {
                            newp2 = cutedges.Get(edge2);
                        }
                        else {
                            Point<3> np2 = Center(mesh.Point(edge2.I1()),
                                    mesh.Point(edge2.I2()));
                            newp2 = mesh.AddPoint(np2);
                            cutedges.Set(edge2, newp2);
                        }

                        PointGeomInfo npgi1, npgi2;

                        int si = mesh.GetFaceDescriptor(oldquad.surfid).SurfNr();

                        PointBetween(
                                mesh.Point(edge1.I1()), mesh.Point(edge1.I2()),
                                0.5, si,
                                pgi11,
                                pgi12,
                                mesh.Point(newp1), npgi1);

                        PointBetween(
                                mesh.Point(edge2.I1()), mesh.Point(edge2.I2()),
                                0.5, si,
                                pgi21,
                                pgi22,
                                mesh.Point(newp2), npgi2);

                        BTBisectQuad(oldquad, newp1, npgi1, newp2, npgi2,
                                newquad1, newquad2);

                        mquads.Elem(i) = newquad1;
                        mquads.push_back(newquad2);
                    }

                hangingsurf =
                        MarkHangingTris(mtris, cutedges) +
                        MarkHangingQuads(mquads, cutedges);

                hangingedge = 0;

                int nseg = mesh.GetNSeg();
                for (int i = 1; i <= nseg; i++) {
                    Segment & seg = mesh.LineSegment(i);
                    INDEX_2 edge(seg[0], seg[1]);
                    edge.Sort();
                    if (cutedges.Used(edge)) {
                        hangingedge = 1;
                        Segment nseg1 = seg;
                        Segment nseg2 = seg;

                        int newpi = cutedges.Get(edge);

                        nseg1[1] = newpi;
                        nseg2[0] = newpi;

                        EdgePointGeomInfo newepgi;

                        PointBetween(
                                mesh.Point(seg[0]), mesh.Point(seg[1]),
                                0.5, seg.surfnr1, seg.surfnr2,
                                seg.epgeominfo[0], seg.epgeominfo[1],
                                mesh.Point(newpi), newepgi);

                        nseg1.epgeominfo[1] = newepgi;
                        nseg2.epgeominfo[0] = newepgi;

                        mesh.LineSegment(i) = nseg1;
                        mesh.AddSegment(nseg2);
                    }
                }
            } while (hangingvol || hangingsurf || hangingedge);

            LOG_DEBUG(mtets.size() << " tets");
            LOG_DEBUG(mtris.size() << " trigs");
            if (mprisms.size()) {
                LOG_DEBUG(mprisms.size() << " prisms");
                LOG_DEBUG(mquads.size() << " quads");
            }
            LOG_DEBUG(mesh.GetNP() << " points");
        }

        if (opt.refine_hp) {
            //
            Array<int> v_order(mesh.GetNP());
            v_order = 0;
            if (mesh.GetDimension() == 3) {
                for (int i = 1; i <= mtets.size(); i++)
                    if (mtets.Elem(i).incorder)
                        mtets.Elem(i).order++;

                for (int i = 0; i < mtets.size(); i++)
                    for (int j = 0; j < 4; j++)
                        if (int(mtets[i].order) > v_order.Elem(mtets[i].pnums[j]))
                            v_order.Elem(mtets[i].pnums[j]) = mtets[i].order;
                for (int i = 0; i < mtets.size(); i++)
                    for (int j = 0; j < 4; j++)
                        if (int(mtets[i].order) < v_order.Elem(mtets[i].pnums[j]) - 1)
                            mtets[i].order = v_order.Elem(mtets[i].pnums[j]) - 1;
            }
            else {
                for (int i = 1; i <= mtris.size(); i++)
                    if (mtris.Elem(i).incorder) {
                        mtris.Elem(i).order++;
                    }

                for (int i = 0; i < mtris.size(); i++)
                    for (int j = 0; j < 3; j++)
                        if (int(mtris[i].order) > v_order.Elem(mtris[i].pnums[j]))
                            v_order.Elem(mtris[i].pnums[j]) = mtris[i].order;
                for (int i = 0; i < mtris.size(); i++) {
                    for (int j = 0; j < 3; j++)
                        if (int(mtris[i].order) < v_order.Elem(mtris[i].pnums[j]) - 1)
                            mtris[i].order = v_order.Elem(mtris[i].pnums[j]) - 1;
                }
            }
        }

        mtets.reserve(mtets.size());
        mprisms.reserve(mprisms.size());
        mids.reserve(mids.size());
        mtris.reserve(mtris.size());
        mquads.reserve(mquads.size());

        mesh.ClearVolumeElements();
        mesh.VolumeElements().reserve(mtets.size() + mprisms.size());
        for (int i = 1; i <= mtets.size(); i++) {
            Element el(TET);
            el.SetIndex(mtets.Get(i).matindex);
            for (int j = 1; j <= 4; j++)
                el.PNum(j) = mtets.Get(i).pnums[j - 1];
            el.SetOrder(mtets.Get(i).order);
            mesh.AddVolumeElement(el);
        }
        for (int i = 1; i <= mprisms.size(); i++) {
            Element el(PRISM);
            el.SetIndex(mprisms.Get(i).matindex);
            for (int j = 1; j <= 6; j++)
                el.PNum(j) = mprisms.Get(i).pnums[j - 1];
            el.SetOrder(mprisms.Get(i).order);

            // degenerated prism ?
            static const int map1[] = {3, 2, 5, 6, 1};
            static const int map2[] = {1, 3, 6, 4, 2};
            static const int map3[] = {2, 1, 4, 5, 3};

            const int * map = NULL;
            int deg1 = 0, deg2 = 0, deg3 = 0;
            // int deg = 0;
            if (el.PNum(1) == el.PNum(4)) {
                map = map1;
                deg1 = 1;
            }
            if (el.PNum(2) == el.PNum(5)) {
                map = map2;
                deg2 = 1;
            }
            if (el.PNum(3) == el.PNum(6)) {
                map = map3;
                deg3 = 1;
            }

            switch (deg1 + deg2 + deg3) {
                case 1:
                {
                    for (int j = 1; j <= 5; j++)
                        el.PNum(j) = mprisms.Get(i).pnums[map[j - 1] - 1];

                    el.SetType(PYRAMID);
                    break;
                }
                case 2:
                {
                    static const int tetmap1[] = {1, 2, 3, 4};
                    static const int tetmap2[] = {2, 3, 1, 5};
                    static const int tetmap3[] = {3, 1, 2, 6};
                    if (!deg1) map = tetmap1;
                    if (!deg2) map = tetmap2;
                    if (!deg3) map = tetmap3;
                    for (int j = 1; j <= 4; j++)
                        el.PNum(j) = mprisms.Get(i).pnums[map[j - 1] - 1];

                    el.SetType(TET);
                    break;
                }
                default:
                    ;
            }
            mesh.AddVolumeElement(el);
        }

        mesh.ClearSurfaceElements();
        for (int i = 1; i <= mtris.size(); i++) {
            Element2d el(TRIG);
            el.SetIndex(mtris.Get(i).surfid);
            el.SetOrder(mtris.Get(i).order);
            for (int j = 1; j <= 3; j++) {
                el.PNum(j) = mtris.Get(i).pnums[j - 1];
                el.GeomInfoPi(j) = mtris.Get(i).pgeominfo[j - 1];
            }
            mesh.AddSurfaceElement(el);
        }
        for (int i = 1; i <= mquads.size(); i++) {
            Element2d el(QUAD);
            el.SetIndex(mquads.Get(i).surfid);
            for (int j = 1; j <= 4; j++)
                el.PNum(j) = mquads.Get(i).pnums[j - 1];
            std::swap(el.PNum(3), el.PNum(4));
            mesh.AddSurfaceElement(el);
        }

        // write multilevel hierarchy to mesh:
        np = mesh.GetNP();
        mesh.mlbetweennodes.resize(np);
        if (mesh.mglevels <= 2) {
            LOG_DEBUG("RESETTING mlbetweennodes");
            for (int i = 1; i <= np; i++) {
                mesh.mlbetweennodes.Elem(i).I1() = 0;
                mesh.mlbetweennodes.Elem(i).I2() = 0;
            }
        }

        BitArray isnewpoint(np);
        isnewpoint.Clear();

        for (int i = 1; i <= cutedges.Size(); i++)
            if (cutedges.UsedPos(i)) {
                INDEX_2 edge;
                PointIndex newpi;
                cutedges.GetData(i, edge, newpi);
                isnewpoint.Set(newpi);
                mesh.mlbetweennodes.Elem(newpi) = edge;
            }

        //  mesh.BuildConnectedNodes();

        mesh.ComputeNVertices();
        mesh.RebuildSurfaceElementLists();

        // update identification tables
        for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++) {
            Array<int, PointIndex::BASE> identmap;

            mesh.GetIdentifications().GetMap(i, identmap);

            for (int j = 1; j <= cutedges.Size(); j++)
                if (cutedges.UsedPos(j)) {
                    INDEX_2 i2;
                    PointIndex newpi;
                    cutedges.GetData(j, i2, newpi);
                    INDEX_2 oi2(identmap.Get(i2.I1()),
                            identmap.Get(i2.I2()));
                    oi2.Sort();
                    if (cutedges.Used(oi2)) {
                        PointIndex onewpi = cutedges.Get(oi2);
                        mesh.GetIdentifications().Add(newpi, onewpi, i);
                    }
                }
        }

        // Repair works only for tets!
        bool do_repair = mesh.PureTetMesh();

        do_repair = false; // JS, March 2009: multigrid crashes

        //noprojection = true;

        if (noprojection) {
            do_repair = false;
            for (int ii = 1; ii <= mesh.GetNP(); ii++) {
                if (isnewpoint.Test(ii) && mesh.mlbetweennodes[ii][0] > 0) {
                    mesh.Point(ii) = Center(mesh.Point(mesh.mlbetweennodes[ii][0]),
                            mesh.Point(mesh.mlbetweennodes[ii][1]));
                }
            }
        }

        // Check/Repair

        static bool repaired_once;
        if (mesh.mglevels == 1)
            repaired_once = false;

        //mesh.Save("before.vol");

        Array<ElementIndex> bad_elts;
        Array<double> pure_badness;

        if (do_repair || quality_loss != NULL) {
            pure_badness.resize(mesh.GetNP() + 2);
            GetPureBadness(mesh, pure_badness, isnewpoint);
        }

        if (do_repair) // by Markus W
        {
            const double max_worsening = 1;

            const bool uselocalworsening = false;

            // bool repaired = false;

            Validate(mesh, bad_elts, pure_badness, max_worsening, uselocalworsening);

            for (int ii = 0; ii < bad_elts.size(); ii++) {
                LOG_DEBUG("bad element " << bad_elts[ii]);
            }

            if (repaired_once || bad_elts.size() > 0) {
                clock_t t1(clock());

                // update id-maps
                int j = 0;
                for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++) {
                    if (mesh.GetIdentifications().GetType(i) == Identifications::PERIODIC) {
                        mesh.GetIdentifications().GetMap(i, *idmaps[j], true);
                        j++;
                    }
                }

                // do the repair
                try {
                    RepairBisection(mesh, bad_elts, isnewpoint, *this,
                            pure_badness,
                            max_worsening, uselocalworsening,
                            idmaps);
                    // repaired = true;
                    repaired_once = true;
                }
                catch (std::runtime_error & ex) {
                    LOG_DEBUG(std::string("Problem: ") + ex.what());
                }

                LOG_DEBUG("Time for Repair: " << double(clock() - t1) / double(CLOCKS_PER_SEC));
                LOG_DEBUG("bad elements after repair: " << bad_elts);

                if (quality_loss != NULL)
                    Validate(mesh, bad_elts, pure_badness, 1e100, uselocalworsening, quality_loss);

                if (idmaps.size() == 0)
                    UpdateEdgeMarks(mesh, idmaps);

                //mesh.ImproveMesh();

            }
        }

        for (int i = 0; i < idmaps.size(); i++)
            delete idmaps[i];
        idmaps.DeleteAll();

        mesh.UpdateTopology();

        if (refelementinfofilewrite != "") {
            LOG_DEBUG("writing marked-elements information to \"" << refelementinfofilewrite << "\"");
            std::ofstream ofst(refelementinfofilewrite.c_str());

            WriteMarkedElements(ofst);

            ofst.close();
        }

        mesh.CalcSurfacesOfNode();

        LOG_DEBUG("Bisection done");
    }

    BisectionOptions::BisectionOptions()
    {
        outfilename = NULL;
        mlfilename = NULL;
        refinementfilename = NULL;
        femcode = NULL;
        maxlevel = 50;
        usemarkedelements = 0;
        refine_hp = 0;
        refine_p = 0;
    }

    Refinement::Refinement()
    {
        optimizer2d = NULL;
    }

    void Refinement::PointBetween(const Point<3> & p1, const Point<3> & p2, double secpoint,
            int surfi,
            const PointGeomInfo & gi1,
            const PointGeomInfo & gi2,
            Point<3> & newp, PointGeomInfo & newgi) const
    {
        newp = p1 + secpoint * (p2 - p1);
    }

    void Refinement::PointBetween(const Point<3> & p1, const Point<3> & p2, double secpoint,
            int surfi1, int surfi2,
            const EdgePointGeomInfo & ap1,
            const EdgePointGeomInfo & ap2,
            Point<3> & newp, EdgePointGeomInfo & newgi) const
    {
        std::cout << "base class edge point between" << std::endl;
        newp = p1 + secpoint * (p2 - p1);
    }

    Vec<3> Refinement::GetTangent(const Point<3> & p, int surfi1, int surfi2,
            const EdgePointGeomInfo & ap1) const
    {
        std::cerr << "Refinement::GetTangent not overloaded" << std::endl;
        return Vec<3> (0, 0, 0);
    }

    Vec<3> Refinement::GetNormal(const Point<3> & p, int surfi1,
            const PointGeomInfo & gi) const
    {
        std::cerr << "Refinement::GetNormal not overloaded" << std::endl;
        return Vec<3> (0, 0, 0);
    }

    void Refinement::ProjectToSurface(Point<3> & p, int surfi) const
    {
        LOG_ERROR("Refinement::ProjectToSurface: no geometry set");
    };

    void Refinement::ProjectToEdge(Point<3> & p, int surfi1, int surfi2, const EdgePointGeomInfo & egi) const
    {
        LOG_ERROR("Refinement::ProjectToEdge not overloaded");
    }
}
