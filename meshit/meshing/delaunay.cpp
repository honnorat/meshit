#include <meshit.hpp>
#include "meshclass.hpp"
#include "../general/seti.hpp"
#include "../gprim/adtree.hpp"
#include "../gprim/geomtest3d.hpp"
#include "../gprim/geomfuncs.hpp"
#include "adfront3.hpp"
#include "meshing3.hpp"
#include "improve3.hpp"
#include "meshtool.hpp"
#include "global.hpp"

namespace meshit {

    static const int deltetfaces[][3] = {
        { 1, 2, 3},
        { 2, 0, 3},
        { 0, 1, 3},
        { 1, 0, 2}
    };

    class DelaunayTet
    {
        PointIndex pnums[4];
        int nb[4];

      public:

        DelaunayTet() { }

        DelaunayTet(const DelaunayTet & el)
        {
            for (int i = 0; i < 4; i++)
                pnums[i] = el[i];
        }

        DelaunayTet(const Element & el)
        {
            for (int i = 0; i < 4; i++)
                pnums[i] = el[i];
        }

        PointIndex & operator[](int i)
        {
            return pnums[i];
        }

        PointIndex operator[](int i) const
        {
            return pnums[i];
        }

        int & NB1(int i)
        {
            return nb[i - 1];
        }

        int NB1(int i) const
        {
            return nb[i - 1];
        }

        int & NB(int i)
        {
            return nb[i];
        }

        int NB(int i) const
        {
            return nb[i];
        }

        int FaceNr(INDEX_3 & face) const // which face nr is it ?
        {
            for (int i = 0; i < 3; i++)
                if (pnums[i] != face.I1() &&
                        pnums[i] != face.I2() &&
                        pnums[i] != face.I3())
                    return i;
            return 3;
        }

        void GetFace1(int i, INDEX_3 & face) const
        {
            face.I(1) = pnums[deltetfaces[i - 1][0]];
            face.I(2) = pnums[deltetfaces[i - 1][1]];
            face.I(3) = pnums[deltetfaces[i - 1][2]];
        }

        void GetFace(int i, INDEX_3 & face) const
        {
            face.I(1) = pnums[deltetfaces[i][0]];
            face.I(2) = pnums[deltetfaces[i][1]];
            face.I(3) = pnums[deltetfaces[i][2]];
        }

        INDEX_3 GetFace1(int i) const
        {
            return INDEX_3(pnums[deltetfaces[i - 1][0]],
                    pnums[deltetfaces[i - 1][1]],
                    pnums[deltetfaces[i - 1][2]]);
        }

        INDEX_3 GetFace(int i) const
        {
            return INDEX_3(pnums[deltetfaces[i][0]],
                    pnums[deltetfaces[i][1]],
                    pnums[deltetfaces[i][2]]);
        }

        void GetFace1(int i, Element2d & face) const
        {
            // face.SetType(TRIG);
            face[0] = pnums[deltetfaces[i - 1][0]];
            face[1] = pnums[deltetfaces[i - 1][1]];
            face[2] = pnums[deltetfaces[i - 1][2]];
        }
    };

    /*
      Table to maintain neighbour elements
     */
    class MeshNB
    {
        // face nodes -> one element
        INDEX_3_CLOSED_HASHTABLE<int> faces;

        // 
        Array<DelaunayTet> & tets;

      public:

        // estimated number of points

        MeshNB(Array<DelaunayTet> & atets, int np)
            : faces(200), tets(atets)
        {
            ;
        }

        // add element with 4 nodes
        void Add(int elnr);

        // delete element with 4 nodes

        void Delete(int elnr)
        {
            DelaunayTet & el = tets.Elem(elnr);
            for (int i = 0; i < 4; i++)
                faces.Set(el.GetFace(i).Sort(), el.NB(i));
        }

        // get neighbour of element elnr in direction fnr 

        int GetNB(int elnr, int fnr)
        {
            return tets.Get(elnr).NB1(fnr);
        }

        //

        void ResetFaceHT(int size)
        {
            faces.SetSize(size);
        }
    };

    void MeshNB::Add(int elnr)
    {
        DelaunayTet & el = tets.Elem(elnr);

        for (int i = 0; i < 4; i++) {
            INDEX_3 i3 = INDEX_3::Sort(el.GetFace(i));

            int posnr;

            if (!faces.PositionCreate(i3, posnr)) {
                // face already in use
                int othertet = faces.GetData(posnr);

                el.NB(i) = othertet;
                if (othertet) {
                    int fnr = tets.Get(othertet).FaceNr(i3);
                    tets.Elem(othertet).NB(fnr) = elnr;
                }
            }
            else {
                faces.SetData(posnr, elnr);
                el.NB(i) = 0;
            }
        }
    }

    /*
      connected lists of cosphereical elements
     */
    class SphereList
    {
        Array<int> links;
      public:

        SphereList()
        {
            ;
        }

        void AddElement(int elnr)
        {
            if (elnr > links.size())
                links.push_back(1);
            links.Elem(elnr) = elnr;
        }

        void DeleteElement(int elnr)
        {
            links.Elem(elnr) = 0;
        }

        void ConnectElement(int eli, int toi)
        {
            links.Elem(eli) = links.Get(toi);
            links.Elem(toi) = eli;
        }

        void GetList(int eli, Array<int> & linked) const;
    };

    void SphereList::GetList(int eli, Array<int> & linked) const
    {
        linked.resize(0);
        int pi = eli;

        do {
            if (pi <= 0 || pi > links.size()) {
                std::cerr << "link, error " << std::endl;
                std::cerr << "pi = " << pi << " linked.s = " << linked.size() << std::endl;
                exit(1);
            }
            if (linked.size() > links.size()) {
                std::cerr << "links have loop" << std::endl;
                exit(1);
            }

            linked.push_back(pi);
            pi = links.Get(pi);
        } while (pi != eli);
    }

    void AddDelaunayPoint(PointIndex newpi, const Point3d & newp,
            Array<DelaunayTet> & tempels,
            Mesh & mesh,
            Box3dTree & tettree,
            MeshNB & meshnb,
            Array<Point<3> > & centers, Array<double> & radi2,
            Array<int> & connected, Array<int> & treesearch,
            Array<int> & freelist, SphereList & list,
            IndexSet & insphere, IndexSet & closesphere)
    {
        /*
          find any sphere, such that newp is contained in
         */

        DelaunayTet el;
        int cfelind = -1;

        const Point<3> * pp[4];
        Point<3> pc;
        double r2;
        Point3d tpmin, tpmax;

        tettree.GetIntersecting(newp, newp, treesearch);
        double quot, minquot(1e20);

        for (int j = 0; j < treesearch.size(); j++) {
            int jjj = treesearch[j];
            quot = Dist2(centers.Get(jjj), newp) / radi2.Get(jjj);

            if ((cfelind == -1 || quot < 0.99 * minquot) && quot < 1) {
                minquot = quot;
                el = tempels.Get(jjj);
                cfelind = jjj;
                if (minquot < 0.917632)
                    break;
            }
        }

        if (cfelind == -1) {
            LOG_WARNING("Delaunay, point not in any sphere");
            return;
        }


        /*
          insphere:     point is in sphere -> delete element
          closesphere:  point is close to sphere -> considered for same center
         */

        // save overestimate
        insphere.SetMaxIndex(2 * tempels.size() + 5 * mesh.GetNP());
        closesphere.SetMaxIndex(2 * tempels.size() + 5 * mesh.GetNP());

        insphere.Clear();
        closesphere.Clear();


        insphere.Add(cfelind);

        int changed = 1;
        int nstarti = 1, starti;

        while (changed) {
            changed = 0;
            starti = nstarti;
            nstarti = insphere.GetArray().size() + 1;


            // if point in sphere, then it is also closesphere
            for (int j = starti; j < nstarti; j++) {
                int helind = insphere.GetArray().Get(j);
                if (!closesphere.IsIn(helind))
                    closesphere.Add(helind);
            }

            // add connected spheres to insphere - list
            for (int j = starti; j < nstarti; j++) {
                list.GetList(insphere.GetArray().Get(j), connected);
                for (int k = 0; k < connected.size(); k++) {
                    int celind = connected[k];

                    if (tempels.Get(celind)[0] != -1 &&
                            !insphere.IsIn(celind)) {
                        changed = 1;
                        insphere.Add(celind);
                    }
                }
            }

            // check neighbour-tets
            for (int j = starti; j < nstarti; j++)
                for (int k = 1; k <= 4; k++) {
                    int helind = insphere.GetArray().Get(j);
                    int nbind = meshnb.GetNB(helind, k);

                    if (nbind && !insphere.IsIn(nbind)) {
                        //changed
                        //int prec = testout->precision();
                        //testout->precision(12);
                        //std::cerr << "val1 " << Dist2 (centers.Get(nbind), newp)
                        //	     << " val2 " << radi2.Get(nbind) * (1+1e-8)
                        //	     << " val3 " << radi2.Get(nbind)
                        //	     << " val1 / val3 " << Dist2 (centers.Get(nbind), newp)/radi2.Get(nbind) <<std::endl;
                        //testout->precision(prec);
                        if (Dist2(centers.Get(nbind), newp)
                                < radi2.Get(nbind) * (1 + 1e-8))
                            closesphere.Add(nbind);

                        if (Dist2(centers.Get(nbind), newp)
                                < radi2.Get(nbind) * (1 + 1e-12)) {
                            // point is in sphere -> remove tet
                            insphere.Add(nbind);
                            changed = 1;
                        }
                        else {
                            /*
                            Element2d face;
                            tempels.Get(helind).GetFace (k, face);

                            const Point3d & p1 = mesh.Point (face.PNum(1));
                            const Point3d & p2 = mesh.Point (face[1]);
                            const Point3d & p3 = mesh.Point (face[2]);
                             */

                            INDEX_3 i3 = tempels.Get(helind).GetFace(k - 1);

                            const Point3d & p1 = mesh.Point(i3.I1());
                            const Point3d & p2 = mesh.Point(i3.I2());
                            const Point3d & p3 = mesh.Point(i3.I3());


                            Vec3d v1(p1, p2);
                            Vec3d v2(p1, p3);
                            Vec3d n = Cross(v1, v2);
                            n /= n.Length();

                            if (n * Vec3d(p1, mesh.Point(tempels.Get(helind)[k - 1])) > 0)
                                n *= -1;

                            double dist = n * Vec3d(p1, newp);


                            if (dist > -1e-10) // 1e-10
                            {
                                insphere.Add(nbind);
                                changed = 1;
                            }


                        }
                    }
                }
        } // while (changed)
        //      std::cerr << "newels: " <<std::endl;
        Array<Element> newels;

        Element2d face(TRIG);

        for (int j = 1; j <= insphere.GetArray().size(); j++)
            for (int k = 1; k <= 4; k++) {
                //	    int elind = insphere.GetArray().Get(j);
                int celind = insphere.GetArray().Get(j);
                int nbind = meshnb.GetNB(celind, k);

                if (!nbind || !insphere.IsIn(nbind)) {
                    tempels.Get(celind).GetFace1(k, face);

                    Element newel(TET);
                    for (int l = 0; l < 3; l++)
                        newel[l] = face[l];
                    newel[3] = newpi;

                    newels.push_back(newel);

                    Vec<3> v1 = mesh[face[1]] - mesh[face[0]];
                    Vec<3> v2 = mesh[face[2]] - mesh[face[0]];
                    Vec<3> n = Cross(v1, v2);

                    n.Normalize();
                    if (n * Vec3d(mesh.Point(face[0]),
                            mesh.Point(tempels.Get(insphere.GetArray().Get(j))[k - 1]))
                            > 0)
                        n *= -1;

                    double hval = n * (newp - mesh[face[0]]);

                    if (hval > -1e-12) {
                        std::cerr << "vec to outer" << std::endl;
                        std::cerr << "vec to outer, hval = " << hval << std::endl;
                        std::cerr << "v1 x v2 = " << Cross(v1, v2) << std::endl;
                        std::cerr << "facep: "
                                << mesh.Point(face[0]) << " "
                                << mesh.Point(face[1]) << " "
                                << mesh.Point(face[2]) << std::endl;
                    }
                }
            }

        meshnb.ResetFaceHT(10 * insphere.GetArray().size() + 1);

        for (int j = 1; j <= insphere.GetArray().size(); j++) {
            //	  int elind = 
            int celind = insphere.GetArray().Get(j);

            meshnb.Delete(celind);
            list.DeleteElement(celind);

            for (int k = 0; k < 4; k++)
                tempels.Elem(celind)[k] = -1;

            ((ADTree6&) tettree.Tree()).DeleteElement(celind);
            freelist.push_back(celind);
        }

        int hasclose = 0;
        for (int j = 1; j <= closesphere.GetArray().size(); j++) {
            int ind = closesphere.GetArray().Get(j);
            if (!insphere.IsIn(ind) &&
                    fabs(Dist2(centers.Get(ind), newp) - radi2.Get(ind)) < 1e-8)
                hasclose = 1;
        }

        for (int j = 1; j <= newels.size(); j++) {
            int nelind;

            if (!freelist.size()) {
                tempels.push_back(newels.Get(j));
                nelind = tempels.size();
            }
            else {
                nelind = freelist.Last();
                freelist.DeleteLast();

                tempels.Elem(nelind) = newels.Get(j);
            }

            meshnb.Add(nelind);
            list.AddElement(nelind);

            for (int k = 0; k < 4; k++)
                pp[k] = &mesh.Point(newels.Get(j)[k]);

            if (CalcSphereCenter(&pp[0], pc)) {
                LOG_ERROR("Delaunay: New tet is flat");

                for (int k = 1; k <= 4; k++)
                    std::cerr << newels.Get(j).PNum(k) << " ";
                std::cerr << std::endl;
                for (int k = 1; k <= 4; k++)
                    std::cerr << *pp[k - 1] << " ";
                std::cerr << std::endl;
            }

            r2 = Dist2(*pp[0], pc);
            if (hasclose)
                for (int k = 1; k <= closesphere.GetArray().size(); k++) {
                    int csameind = closesphere.GetArray().Get(k);
                    if (!insphere.IsIn(csameind) &&
                            fabs(r2 - radi2.Get(csameind)) < 1e-10 &&
                            Dist(pc, centers.Get(csameind)) < 1e-10) {
                        pc = centers.Get(csameind);
                        r2 = radi2.Get(csameind);
                        list.ConnectElement(nelind, csameind);
                        break;
                    }
                }

            if (centers.size() < nelind) {
                centers.push_back(pc);
                radi2.push_back(r2);
            }
            else {
                centers.Elem(nelind) = pc;
                radi2.Elem(nelind) = r2;
            }

            closesphere.Add(nelind);

            tpmax = tpmin = *pp[0];
            for (int k = 1; k <= 3; k++) {
                tpmin.SetToMin(*pp[k]);
                tpmax.SetToMax(*pp[k]);
            }
            tpmax = tpmax + 0.01 * (tpmax - tpmin);
            tettree.Insert(tpmin, tpmax, nelind);
        }
    }

    void Delaunay1(
            Mesh & mesh, const MeshingParameters & mp,
            AdFront3 * adfront, Array<DelaunayTet> & tempels,
            int oldnp, DelaunayTet & startel, Point3d & pmin, Point3d & pmax)
    {
        Array<Point<3> > centers;
        Array<double> radi2;

        Point3d tpmin, tpmax;

        // new: local box
        mesh.GetBox(pmax, pmin); // lower bound for pmax, upper for pmin
        for (int i = 1; i <= adfront->GetNF(); i++) {
            const MiniElement2d & face = adfront->GetFace(i);
            for (int j = 0; j < face.GetNP(); j++) {
                pmin.SetToMin(mesh.Point(face[j]));
                pmax.SetToMax(mesh.Point(face[j]));
            }
        }

        for (int i = 0; i < mesh.LockedPoints().size(); i++) {
            pmin.SetToMin(mesh.Point(mesh.LockedPoints()[i]));
            pmax.SetToMax(mesh.Point(mesh.LockedPoints()[i]));
        }



        Vec3d vdiag(pmin, pmax);
        // double r1 = vdiag.Length();
        double r1 = sqrt(3.0) * max3(vdiag.X(), vdiag.Y(), vdiag.Z());
        vdiag = Vec3d(r1, r1, r1);
        //double r2;

        Point3d pmin2 = pmin - 8 * vdiag;
        Point3d pmax2 = pmax + 8 * vdiag;

        Point3d cp1(pmin2), cp2(pmax2), cp3(pmax2), cp4(pmax2);
        cp2.X() = pmin2.X();
        cp3.Y() = pmin2.Y();
        cp4.Z() = pmin2.Z();




        int np = mesh.GetNP();

        startel[0] = mesh.AddPoint(cp1);
        startel[1] = mesh.AddPoint(cp2);
        startel[2] = mesh.AddPoint(cp3);
        startel[3] = mesh.AddPoint(cp4);

        // flag points to use for Delaunay:
        BitArrayChar<PointIndex::BASE> usep(np);
        usep.Clear();
        for (int i = 1; i <= adfront->GetNF(); i++) {
            const MiniElement2d & face = adfront->GetFace(i);
            for (int j = 0; j < face.GetNP(); j++)
                usep.Set(face[j]);
        }

        for (int i = oldnp + PointIndex::BASE;
                i < np + PointIndex::BASE; i++)
            usep.Set(i);

        for (int i = 0; i < mesh.LockedPoints().size(); i++)
            usep.Set(mesh.LockedPoints()[i]);


        Array<int> freelist;


        int cntp = 0;

        MeshNB meshnb(tempels, mesh.GetNP() + 5);
        SphereList list;

        pmin2 = pmin2 + 0.1 * (pmin2 - pmax2);
        pmax2 = pmax2 + 0.1 * (pmax2 - pmin2);

        Box3dTree tettree(pmin2, pmax2);


        tempels.push_back(startel);
        meshnb.Add(1);
        list.AddElement(1);
        Array<int> connected, treesearch;


        tpmin = tpmax = mesh.Point(startel[0]);
        for (int k = 1; k < 4; k++) {
            tpmin.SetToMin(mesh.Point(startel[k]));
            tpmax.SetToMax(mesh.Point(startel[k]));
        }
        tpmax = tpmax + 0.01 * (tpmax - tpmin);
        tettree.Insert(tpmin, tpmax, 1);

        Point<3> pc;

        const Point<3> * pp[4];
        for (int k = 0; k < 4; k++)
            pp[k] = &mesh.Point(startel[k]);
        CalcSphereCenter(&pp[0], pc);

        centers.push_back(pc);
        radi2.push_back(Dist2(*pp[0], pc));


        IndexSet insphere(mesh.GetNP());
        IndexSet closesphere(mesh.GetNP());

        // "random" reordering of points  (speeds a factor 3 - 5 !!!)
        Array<PointIndex, PointIndex::BASE, PointIndex> mixed(np);
        int prims[] = {11, 13, 17, 19, 23, 29, 31, 37};
        int prim;

        {
            int i = 0;
            while (np % prims[i] == 0) i++;
            prim = prims[i];
        }

        for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End() - 4; pi++)
            mixed[pi] = PointIndex((prim * pi) % np + PointIndex::BASE);

        for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End() - 4; pi++) {
            if (pi % 1000 == 0) {
                if (pi % 10000 == 0)
                    PrintDot('+');
                else
                    PrintDot('.');
            }

            PointIndex newpi = mixed[pi];

            if (!usep.Test(newpi))
                continue;

            cntp++;

            const MeshPoint & newp = mesh[newpi];

            AddDelaunayPoint(newpi, newp, tempels, mesh,
                    tettree, meshnb, centers, radi2,
                    connected, treesearch, freelist, list, insphere, closesphere);

        }

        for (int i = tempels.size(); i >= 1; i--)
            if (tempels.Get(i)[0] <= 0)
                tempels.DeleteElement(i);

        PrintDot('\n');

        LOG_DEBUG("Points: " << cntp);
        LOG_DEBUG("Elements: " << tempels.size());

    }

    void Meshing3::Delaunay(Mesh & mesh, int domainnr, const MeshingParameters & mp)
    {
        int np, ne;

        LOG_DEBUG("Delaunay meshing");
        LOG_DEBUG("  number of points: " << mesh.GetNP());

        Array<DelaunayTet> tempels;
        Point3d pmin, pmax;

        DelaunayTet startel;

        int oldnp = mesh.GetNP();
        if (mp.blockfill) {
            BlockFillLocalH(mesh, mp);
            LOG_DEBUG("number of points: " << mesh.GetNP());
        }

        np = mesh.GetNP();

        Delaunay1(mesh, mp, adfront, tempels, oldnp, startel, pmin, pmax);

        {
            // improve delaunay - mesh by std::swapping !!!!

            Mesh tempmesh;
            for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End(); pi++)
                tempmesh.AddPoint(mesh[pi]);

            for (int i = 1; i <= tempels.size(); i++) {
                Element el(4);
                for (int j = 0; j < 4; j++)
                    el[j] = tempels.Elem(i)[j];

                el.SetIndex(1);

                const Point3d & lp1 = mesh.Point(el[0]);
                const Point3d & lp2 = mesh.Point(el[1]);
                const Point3d & lp3 = mesh.Point(el[2]);
                const Point3d & lp4 = mesh.Point(el[3]);
                Vec3d v1(lp1, lp2);
                Vec3d v2(lp1, lp3);
                Vec3d v3(lp1, lp4);

                Vec3d n = Cross(v1, v2);
                double vol = n * v3;
                if (vol > 0) std::swap(el[2], el[3]);

                tempmesh.AddVolumeElement(el);
            }


            MeshQuality3d(tempmesh);

            tempmesh.AddFaceDescriptor(FaceDescriptor(1, 1, 0, 0));
            tempmesh.AddFaceDescriptor(FaceDescriptor(2, 1, 0, 0));



            for (int i = 1; i <= mesh.GetNOpenElements(); i++) {
                Element2d sel = mesh.OpenElement(i);
                sel.SetIndex(1);
                tempmesh.AddSurfaceElement(sel);
                std::swap(sel[1], sel[2]);
                tempmesh.AddSurfaceElement(sel);
            }


            for (int i = 1; i <= 4; i++) {
                Element2d self(TRIG);
                self.SetIndex(1);
                startel.GetFace1(i, self);
                tempmesh.AddSurfaceElement(self);
            }


            //  for (i = mesh.GetNP() - 3; i <= mesh.GetNP(); i++)
            //    tempmesh.AddLockedPoint (i);
            for (PointIndex pi = tempmesh.Points().Begin();
                    pi < tempmesh.Points().End(); pi++)
                tempmesh.AddLockedPoint(pi);

            //    tempmesh.PrintMemInfo(std::cout);
            // tempmesh.Save ("tempmesh.vol");

            for (int i = 1; i <= 2; i++) {
                tempmesh.FindOpenElements();

                LOG_DEBUG("Num open: " << tempmesh.GetNOpenElements());
                tempmesh.CalcSurfacesOfNode();

                tempmesh.FreeOpenElementsEnvironment(1);

                MeshOptimize3d meshopt(mp);
                // tempmesh.CalcSurfacesOfNode();
                meshopt.SwapImprove(tempmesh, OPT_CONFORM);
            }

            MeshQuality3d(tempmesh);

            tempels.resize(0);
            for (int i = 1; i <= tempmesh.GetNE(); i++)
                tempels.push_back(tempmesh.VolumeElement(i));
        }



        // remove degenerated

        BitArray badnode(mesh.GetNP());
        badnode.Clear();
        int ndeg = 0;
        for (int i = 1; i <= tempels.size(); i++) {
            Element el(4);
            for (int j = 0; j < 4; j++)
                el[j] = tempels.Elem(i)[j];
            //      Element & el = tempels.Elem(i);
            const Point3d & lp1 = mesh.Point(el[0]);
            const Point3d & lp2 = mesh.Point(el[1]);
            const Point3d & lp3 = mesh.Point(el[2]);
            const Point3d & lp4 = mesh.Point(el[3]);
            Vec3d v1(lp1, lp2);
            Vec3d v2(lp1, lp3);
            Vec3d v3(lp1, lp4);
            Vec3d n = Cross(v1, v2);
            double vol = n * v3;

            double h = v1.Length() + v2.Length() + v3.Length();
            if (fabs(vol) < 1e-8 * (h * h * h) &&
                    (el[0] <= np && el[1] <= np &&
                    el[2] <= np && el[3] <= np)) // old: 1e-12
            {
                badnode.Set(el[0]);
                badnode.Set(el[1]);
                badnode.Set(el[2]);
                badnode.Set(el[3]);
                ndeg++;
                std::cerr << "vol = " << vol << " h = " << h << std::endl;
            }

            if (vol > 0)
                std::swap(el[2], el[3]);
        }

        ne = tempels.size();
        for (int i = ne; i >= 1; i--) {
            const DelaunayTet & el = tempels.Get(i);
            if (badnode.Test(el[0]) ||
                    badnode.Test(el[1]) ||
                    badnode.Test(el[2]) ||
                    badnode.Test(el[3]))
                tempels.DeleteElement(i);
        }


        LOG_DEBUG(ndeg << " degenerated elements removed");

        // find surface triangles which are no face of any tet

        INDEX_3_HASHTABLE<int> openeltab(mesh.GetNOpenElements() + 3);
        Array<int> openels;
        for (int i = 1; i <= mesh.GetNOpenElements(); i++) {
            const Element2d & tri = mesh.OpenElement(i);
            INDEX_3 i3(tri[0], tri[1], tri[2]);
            i3.Sort();
            openeltab.Set(i3, i);
        }

        for (int i = 1; i <= tempels.size(); i++) {
            for (int j = 0; j < 4; j++) {
                INDEX_3 i3 = tempels.Get(i).GetFace(j);
                i3.Sort();
                if (openeltab.Used(i3))
                    openeltab.Set(i3, 0);
            }
        }

        // and store them in openels
        for (int i = 1; i <= openeltab.GetNBags(); i++)
            for (int j = 1; j <= openeltab.GetBagSize(i); j++) {
                INDEX_3 i3;
                int fnr;
                openeltab.GetData(i, j, i3, fnr);
                if (fnr)
                    openels.push_back(fnr);
            }





        // find open triangle with close edge (from halfening of surface squares)

        INDEX_2_HASHTABLE<INDEX_2> twotrias(mesh.GetNOpenElements() + 5);
        //  for (i = 1; i <= mesh.GetNOpenElements(); i++)
        for (int ii = 1; ii <= openels.size(); ii++) {
            int i = openels.Get(ii);
            const Element2d & el = mesh.OpenElement(i);
            for (int j = 1; j <= 3; j++) {
                INDEX_2 hi2(el.PNumMod(j), el.PNumMod(j + 1));
                hi2.Sort();
                if (twotrias.Used(hi2)) {
                    INDEX_2 hi3;
                    hi3 = twotrias.Get(hi2);
                    hi3.I2() = el.PNumMod(j + 2);
                    twotrias.Set(hi2, hi3);
                }
                else {
                    INDEX_2 hi3(el.PNumMod(j + 2), 0);
                    twotrias.Set(hi2, hi3);
                }
            }
        }

        INDEX_2_HASHTABLE<int> tetedges(tempels.size() + 5);
        for (int i = 1; i <= tempels.size(); i++) {
            const DelaunayTet & el = tempels.Get(i);
            INDEX_2 i2;
            for (int j = 1; j <= 6; j++) {
                switch (j) {
                    case 1: i2.I1() = el[0];
                        i2.I2() = el[1];
                        break;
                    case 2: i2.I1() = el[0];
                        i2.I2() = el[2];
                        break;
                    case 3: i2.I1() = el[0];
                        i2.I2() = el[3];
                        break;
                    case 4: i2.I1() = el[1];
                        i2.I2() = el[2];
                        break;
                    case 5: i2.I1() = el[1];
                        i2.I2() = el[3];
                        break;
                    case 6: i2.I1() = el[2];
                        i2.I2() = el[3];
                        break;
                    default: i2.I1() = i2.I2() = 0;
                        break;
                }
                i2.Sort();
                tetedges.Set(i2, 1);
            }
        }
        //  std::cout << "tetedges:";
        //  tetedges.PrintMemInfo (std::cout);


        for (INDEX_2_HASHTABLE<INDEX_2>::Iterator it = twotrias.Begin();
                it != twotrias.End(); it++) {
            INDEX_2 hi2, hi3;
            twotrias.GetData(it, hi2, hi3);
            hi3.Sort();
            if (tetedges.Used(hi3)) {
                const Point3d & p1 = mesh.Point(PointIndex(hi2.I1()));
                const Point3d & p2 = mesh.Point(PointIndex(hi2.I2()));
                const Point3d & p3 = mesh.Point(PointIndex(hi3.I1()));
                const Point3d & p4 = mesh.Point(PointIndex(hi3.I2()));
                Vec3d v1(p1, p2);
                Vec3d v2(p1, p3);
                Vec3d v3(p1, p4);
                Vec3d n = Cross(v1, v2);
                double vol = n * v3;

                double h = v1.Length() + v2.Length() + v3.Length();
                if (fabs(vol) < 1e-4 * (h * h * h)) // old: 1e-12
                {
                    badnode.Set(hi3.I1());
                    badnode.Set(hi3.I2());
                }
            }
        }

        ne = tempels.size();
        for (int i = ne; i >= 1; i--) {
            const DelaunayTet & el = tempels.Get(i);
            if (badnode.Test(el[0]) ||
                    badnode.Test(el[1]) ||
                    badnode.Test(el[2]) ||
                    badnode.Test(el[3]))
                tempels.DeleteElement(i);
        }

        // find intersecting:
        LOG_DEBUG("Remove intersecting");
        if (openels.size()) {
            Box3dTree setree(pmin, pmax);

            for (int i = 1; i <= openels.size(); i++) {
                int fnr;
                fnr = openels.Get(i);
                if (fnr) {
                    const Element2d & tri = mesh.OpenElement(fnr);

                    Point3d ltpmin(mesh.Point(tri[0]));
                    Point3d ltpmax(ltpmin);

                    for (int k = 2; k <= 3; k++) {
                        ltpmin.SetToMin(mesh.Point(tri.PNum(k)));
                        ltpmax.SetToMax(mesh.Point(tri.PNum(k)));
                    }
                    setree.Insert(ltpmin, ltpmax, fnr);
                }
            }

            Array<int> neartrias;
            for (int i = 1; i <= tempels.size(); i++) {
                const Point<3> *pp[4];
                int tetpi[4];
                DelaunayTet & el = tempels.Elem(i);

                int intersect = 0;

                for (int j = 0; j < 4; j++) {
                    pp[j] = &mesh.Point(el[j]);
                    tetpi[j] = el[j];
                }

                Point3d tetpmin(*pp[0]);
                Point3d tetpmax(tetpmin);
                for (int j = 1; j < 4; j++) {
                    tetpmin.SetToMin(*pp[j]);
                    tetpmax.SetToMax(*pp[j]);
                }
                tetpmin = tetpmin + 0.01 * (tetpmin - tetpmax);
                tetpmax = tetpmax + 0.01 * (tetpmax - tetpmin);

                setree.GetIntersecting(tetpmin, tetpmax, neartrias);


                //      for (j = 1; j <= mesh.GetNSE(); j++)
                //	{
                for (int jj = 1; jj <= neartrias.size(); jj++) {
                    int j = neartrias.Get(jj);

                    const Element2d & tri = mesh.OpenElement(j);
                    const Point<3> *tripp[3];
                    int tripi[3];

                    for (int k = 1; k <= 3; k++) {
                        tripp[k - 1] = &mesh.Point(tri.PNum(k));
                        tripi[k - 1] = tri.PNum(k);
                    }

                    if (IntersectTetTriangle(&pp[0], &tripp[0], tetpi, tripi)) {
                        /*
                        int il1, il2;
                        std::cerr << "intersect !" <<std::endl;
                        std::cerr << "triind: ";
                        for (il1 = 0; il1 < 3; il1++)
                          std::cerr << " " << tripi[il1];
                        std::cerr <<std::endl;
                        std::cerr << "tetind: ";
                        for (il2 = 0; il2 < 4; il2++)
                          std::cerr << " " << tetpi[il2];
                        std::cerr <<std::endl;
		  
                        std::cerr << "trip: ";
                        for (il1 = 0; il1 < 3; il1++)
                          std::cerr << " " << *tripp[il1];
                        std::cerr <<std::endl;
                        std::cerr << "tetp: ";
                        for (il2 = 0; il2 < 4; il2++)
                          std::cerr << " " << *pp[il2];
                        std::cerr <<std::endl;
                         */


                        intersect = 1;
                        break;
                    }
                }


                if (intersect) {
                    tempels.DeleteElement(i);
                    i--;
                }
            }
        }

        LOG_DEBUG("Remove outer");

        // find connected tets (with no face between, and no hole due
        // to removed intersecting tets.
        //  INDEX_3_HASHTABLE<INDEX_2> innerfaces(np);


        INDEX_3_HASHTABLE<int> boundaryfaces(mesh.GetNOpenElements() / 3 + 1);
        for (int i = 1; i <= mesh.GetNOpenElements(); i++) {
            const Element2d & tri = mesh.OpenElement(i);
            INDEX_3 i3(tri[0], tri[1], tri[2]);
            i3.Sort();
            boundaryfaces.PrepareSet(i3);
        }
        boundaryfaces.AllocateElements();
        for (int i = 1; i <= mesh.GetNOpenElements(); i++) {
            const Element2d & tri = mesh.OpenElement(i);
            INDEX_3 i3(tri[0], tri[1], tri[2]);
            i3.Sort();
            boundaryfaces.Set(i3, 1);
        }

        for (int i = 0; i < tempels.size(); i++)
            for (int j = 0; j < 4; j++)
                tempels[i].NB(j) = 0;

        TABLE<int, PointIndex::BASE> elsonpoint(mesh.GetNP());
        for (int i = 0; i < tempels.size(); i++) {
            const DelaunayTet & el = tempels[i];
            INDEX_4 i4(el[0], el[1], el[2], el[3]);
            i4.Sort();
            elsonpoint.IncSizePrepare(i4.I1());
            elsonpoint.IncSizePrepare(i4.I2());
        }

        elsonpoint.AllocateElementsOneBlock();

        for (int i = 0; i < tempels.size(); i++) {
            const DelaunayTet & el = tempels[i];
            INDEX_4 i4(el[0], el[1], el[2], el[3]);
            i4.Sort();
            elsonpoint.Add(i4.I1(), i + 1);
            elsonpoint.Add(i4.I2(), i + 1);
        }

        INDEX_3_CLOSED_HASHTABLE<INDEX_2> faceht(100);

        Element2d hel(TRIG);
        for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End(); pi++) {
            faceht.SetSize(4 * elsonpoint[pi].size());
            for (int ii = 0; ii < elsonpoint[pi].size(); ii++) {
                int i = elsonpoint[pi][ii];
                const DelaunayTet & el = tempels.Get(i);

                for (int j = 1; j <= 4; j++) {
                    el.GetFace1(j, hel);
                    hel.Invert();
                    hel.NormalizeNumbering();

                    if (hel[0] == pi) {
                        INDEX_3 i3(hel[0], hel[1], hel[2]);

                        if (!boundaryfaces.Used(i3)) {
                            if (faceht.Used(i3)) {
                                INDEX_2 i2 = faceht.Get(i3);

                                tempels.Elem(i).NB1(j) = i2.I1();
                                tempels.Elem(i2.I1()).NB1(i2.I2()) = i;
                            }
                            else {
                                hel.Invert();
                                hel.NormalizeNumbering();
                                INDEX_3 i3i(hel[0], hel[1], hel[2]);
                                INDEX_2 i2(i, j);
                                faceht.Set(i3i, i2);
                            }
                        }
                    }
                }
            }
        }

        /*
          for (i = 1; i <= tempels.Size(); i++)
          {
          const DelaunayTet & el = tempels.Get(i);
          for (j = 1; j <= 4; j++)
          {
          INDEX_3 i3;
          Element2d face;
          el.GetFace1 (j, face);
          for (int kk = 1; kk <= 3; kk++)
          i3.I(kk) = face.PNum(kk);

          i3.Sort();
          if (!boundaryfaces.Used (i3))
          {
          if (innerfaces.Used(i3))
          {
          INDEX_2 i2;
          i2 = innerfaces.Get(i3);
          i2.I2() = i;
          innerfaces.Set (i3, i2);
          }
          else
          {
          INDEX_2 i2;
          i2.I1() = i;
          i2.I2() = 0;
          innerfaces.Set (i3, i2);
          }
          }
          }
          }
         */

        /*
          std::cerr << "nb elements:" <<std::endl;
          for (i = 1; i <= tempels.Size(); i++)
          {
          std::cerr << i << " ";
          for (j = 1; j <= 4; j++)
          std::cerr << tempels.Get(i).NB1(j) << " ";
          std::cerr <<std::endl;
          }
  
          std::cerr << "pairs:" <<std::endl;
          for (i = 1; i <= innerfaces.GetNBags(); i++)
          for (j = 1; j <= innerfaces.GetBagSize(i); j++)
          {
          INDEX_3 i3;
          INDEX_2 i2;
          innerfaces.GetData (i, j, i3, i2);
          std::cerr << i2 <<std::endl;
          }
         */







        /*
          std::cout << "innerfaces: ";
          innerfaces.PrintMemInfo (std::cout);
         */

        //  std::cout << "boundaryfaces: ";
        //  boundaryfaces.PrintMemInfo (std::cout);


        LOG_DEBUG("tables filled");


        ne = tempels.size();
        BitArray inner(ne), outer(ne);
        inner.Clear();
        outer.Clear();
        Array<int> elstack;

        /*
          int starti = 0;
          for (i = 1; i <= ne; i++)
          {
          const Element & el = tempels.Get(i);
          for (j = 1; j <= 4; j++)
          for (k = 1; k <= 4; k++)
          if (el.PNum(j) == startel.PNum(k))
          {
          outer.Set(i);
          starti = i;
          }
          }
         */

        while (1) {
            int inside;
            bool done = 1;

            int i;
            for (i = 1; i <= ne; i++)
                if (!inner.Test(i) && !outer.Test(i)) {
                    done = 0;
                    break;
                }

            if (done) break;

            const DelaunayTet & el = tempels.Get(i);
            const Point3d & p1 = mesh.Point(el[0]);
            const Point3d & p2 = mesh.Point(el[1]);
            const Point3d & p3 = mesh.Point(el[2]);
            const Point3d & p4 = mesh.Point(el[3]);

            Point3d ci = Center(p1, p2, p3, p4);

            inside = adfront->Inside(ci);

            /*
              std::cout << "startel: " << i <<std::endl;
              std::cout << "inside = " << inside <<std::endl;
              std::cout << "ins2 = " << adfront->Inside (Center (ci, p1)) <<std::endl;
              std::cout << "ins3 = " << adfront->Inside (Center (ci, p2)) <<std::endl;
             */

            elstack.resize(0);
            elstack.push_back(i);

            while (elstack.size()) {
                int ei = elstack.Last();
                elstack.DeleteLast();

                if (!inner.Test(ei) && !outer.Test(ei)) {
                    if (inside)
                        inner.Set(ei);
                    else
                        outer.Set(ei);


                    for (int j = 1; j <= 4; j++) {
                        INDEX_3 i3 = tempels.Get(ei).GetFace1(j);
                        /*
                        Element2d face;
                        tempels.Get(ei).GetFace(j, face);
                        for (int kk = 1; kk <= 3; kk++)
                          i3.I(kk) = face.PNum(kk);
                         */
                        i3.Sort();


                        if (tempels.Get(ei).NB1(j))
                            elstack.push_back(tempels.Get(ei).NB1(j));

                        /*
                          if (innerfaces.Used(i3))
                          {
                          INDEX_2 i2 = innerfaces.Get(i3);
                          int other = i2.I1() + i2.I2() - ei;

                          if (other != tempels.Get(ei).NB1(j))
                          std::cerr << "different1 !!" <<std::endl;

                          if (other)
                          {
                          elstack.Append (other);
                          }
                          }
                          else
                          if (tempels.Get(ei).NB1(j))
                          std::cerr << "different2 !!" <<std::endl;
                         */

                    }
                }
            }
        }



        // check outer elements
        if (debugparam.slowchecks) {
            for (int i = 1; i <= ne; i++) {
                const DelaunayTet & el = tempels.Get(i);
                const Point3d & p1 = mesh.Point(el[0]);
                const Point3d & p2 = mesh.Point(el[1]);
                const Point3d & p3 = mesh.Point(el[2]);
                const Point3d & p4 = mesh.Point(el[3]);

                Point3d ci = Center(p1, p2, p3, p4);

                //       if (adfront->Inside (ci) != adfront->Inside (Center (ci, p1)))
                // 	std::cout << "ERROR: outer test unclear !!!" <<std::endl;	

                if (inner.Test(i) != adfront->Inside(ci)) {
                    /*
                      std::cout << "ERROR: outer test wrong !!!" 
                      << "inner = " << int(inner.Test(i))
                      << "outer = " << int(outer.Test(i))
                      <<std::endl;
	      
                      std::cout << "Vol = " << Determinant(Vec3d(p1, p2),
                      Vec3d(p1, p3),
                      Vec3d(p1, p4)) <<std::endl;
	      
                     */
                    for (int j = 1; j <= 4; j++) {
                        Point3d hp;
                        switch (j) {
                            case 1: hp = Center(ci, p1);
                                break;
                            case 2: hp = Center(ci, p2);
                                break;
                            case 3: hp = Center(ci, p3);
                                break;
                            case 4: hp = Center(ci, p4);
                                break;
                        }
                        //		  std::cout << "inside(" << hp << ") = " << adfront->Inside(hp) <<std::endl;
                    }

                }

                if (adfront->Inside(ci))
                    outer.Clear(i);
                else
                    outer.Set(i);
            }
        }


        /*

        // find bug in innerfaces

        tempmesh.DeleteVolumeElements();

        for (i = 1; i <= innerfaces.GetNBags(); i++)
        for (j = 1; j <= innerfaces.GetBagSize(i); j++)
        {
        INDEX_3 i3;
        INDEX_2 i2;
        innerfaces.GetData (i, j, i3, i2);
        if (i2.I2())
        {
        if (outer.Test(i2.I1()) != outer.Test(i2.I2()))
        {
        tempmesh.AddVolumeElement (tempels.Get(i2.I1()));
        tempmesh.AddVolumeElement (tempels.Get(i2.I2()));
        std::cerr << "outer flag different for connected els" <<std::endl;
        }
        }
        }


        std::cout << "Check intersectiong once more" <<std::endl;

        for (i = 1; i <= openels.Size(); i++)
        {
        tempmesh.SurfaceElement(2*openels.Get(i)).SetIndex(2);
        tempmesh.SurfaceElement(2*openels.Get(i)-1).SetIndex(2);
        }

        //  for (i = 1; i <= tempmesh.GetNE(); i++)
        //    for (j = 1; j <= tempmesh.GetNSE(); j++)
        i = 6; j = 403;
        if (i <= tempmesh.GetNE() && j <= tempmesh.GetNSE())
        if (tempmesh.SurfaceElement(j).GetIndex()==2)
        {
        const Element & el = tempmesh.VolumeElement(i);
        const Element2d & sel = tempmesh.SurfaceElement(j);

        const Point3d *tripp[3];
        const Point3d *pp[4];
        int tetpi[4], tripi[3];

        for (k = 1; k <= 4; k++)
        {
        pp[k-1] = &tempmesh.Point(el.PNum(k));
        tetpi[k-1] = el.PNum(k);
        }

        for (k = 1; k <= 3; k++)
        {
        tripp[k-1] = &tempmesh.Point (sel.PNum(k));
        tripi[k-1] = sel.PNum(k);
        }

        std::cerr << "Check Triangle " << j << ":";
        for (k = 1; k <= 3; k++)
        std::cerr << " " << sel.PNum(k);
        for (k = 1; k <= 3; k++)
        std::cerr << " " << tempmesh.Point(sel.PNum(k));
        std::cerr <<std::endl;

        std::cerr << "Check Tet " << i << ":";
        for (k = 1; k <= 4; k++)
        std::cerr << " " << el.PNum(k);
        for (k = 1; k <= 4; k++)
        std::cerr << " " << tempmesh.Point(el.PNum(k));
        std::cerr <<std::endl;

        if (IntersectTetTriangle (&pp[0], &tripp[0], tetpi, tripi))
        {
        std::cout << "Intesection detected !!" <<std::endl;
        }
        }

        tempmesh.Save ("temp.vol");

        // end bug search
         */


        for (int i = ne; i >= 1; i--) {
            if (outer.Test(i))
                tempels.DeleteElement(i);
        }


        // mesh.points.SetSize(mesh.points.Size()-4);

        for (int i = 0; i < tempels.size(); i++) {
            Element el(4);
            for (int j = 0; j < 4; j++)
                el[j] = tempels[i][j];
            mesh.AddVolumeElement(el);
        }

        LOG_DEBUG("outer removed");

        mesh.FindOpenElements(domainnr);
        mesh.Compress();
    }
}
