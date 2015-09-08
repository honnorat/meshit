#include <meshit.hpp>
#include <stdexcept>
#include "meshfunc.hpp"
#include "meshtool.hpp"
#include "global.hpp"
#include "improve3.hpp"

namespace meshit {
    extern const char * tetrules[];
    // extern const char * tetrules2[];
    extern const char * prismrules2[];
    extern const char * pyramidrules[];
    extern const char * pyramidrules2[];

    // extern double teterrpow; 

    MESHING3_RESULT MeshVolume(MeshingParameters & mp, Mesh& mesh3d)
    {
        int oldne;
        int meshed;

        Array<INDEX_2> connectednodes;

        if (&mesh3d.LocalHFunction() == NULL) mesh3d.CalcLocalH(mp.grading);

        mesh3d.Compress();

        if (mp.checkoverlappingboundary)
            if (mesh3d.CheckOverlappingBoundary())
                throw std::runtime_error("Stop meshing since boundary mesh is overlapping");

        int nonconsist = 0;
        for (int k = 1; k <= mesh3d.GetNDomains(); k++) {
            LOG_DEBUG("Check subdomain " << k << " / " << mesh3d.GetNDomains());

            mesh3d.FindOpenElements(k);

            bool res = (mesh3d.CheckConsistentBoundary() != 0);
            if (res) {
                LOG_WARNING("Surface mesh not consistent");
                nonconsist = 1;
            }
        }

        if (nonconsist) {
            LOG_ERROR("Stop meshing since surface mesh not consistent");
            throw std::runtime_error("Stop meshing since surface mesh not consistent");
        }

        double globmaxh = mp.maxh;

        for (int k = 1; k <= mesh3d.GetNDomains(); k++) {
            if (multithread.terminate)
                break;

            LOG_DEBUG("Meshing subdomain " << k << " of " << mesh3d.GetNDomains());

            mp.maxh = min2(globmaxh, mesh3d.MaxHDomain(k));

            mesh3d.CalcSurfacesOfNode();
            mesh3d.FindOpenElements(k);

            if (!mesh3d.GetNOpenElements())
                continue;

            Box<3> domain_bbox(Box<3>::EMPTY_BOX);

            for (SurfaceElementIndex sei = 0; sei < mesh3d.GetNSE(); sei++) {
                const Element2d & el = mesh3d.SurfaceElement(sei);
                if (el.IsDeleted()) continue;

                if (mesh3d.GetFaceDescriptor(el.GetIndex()).DomainIn() == k ||
                        mesh3d.GetFaceDescriptor(el.GetIndex()).DomainOut() == k)

                    for (int j = 0; j < el.GetNP(); j++)
                        domain_bbox.Add(mesh3d[el[j]]);
            }
            domain_bbox.Increase(0.01 * domain_bbox.Diam());


            for (int qstep = 1; qstep <= 3; qstep++) {
                // std::cout << "openquads = " << mesh3d.HasOpenQuads() <<std::endl;
                if (mesh3d.HasOpenQuads()) {
                    const char ** rulep = NULL;
                    switch (qstep) {
                        case 1:
                            rulep = prismrules2;
                            break;
                        case 2: // connect pyramid to triangle
                            rulep = pyramidrules2;
                            break;
                        case 3: // connect to vis-a-vis point
                            rulep = pyramidrules;
                            break;
                    }

                    Meshing3 meshing(rulep);

                    MeshingParameters mpquad = mp;

                    mpquad.giveuptol = 15;
                    mpquad.baseelnp = 4;
                    mpquad.starshapeclass = 1000;
                    mpquad.check_impossible = qstep == 1; // for prisms only (air domain in trafo)


                    for (PointIndex pi = mesh3d.Points().Begin(); pi < mesh3d.Points().End(); pi++)
                        meshing.AddPoint(mesh3d[pi], pi);

                    mesh3d.GetIdentifications().GetPairs(0, connectednodes);
                    for (int i = 1; i <= connectednodes.size(); i++)
                        meshing.AddConnectedPair(connectednodes.Get(i));

                    for (int i = 1; i <= mesh3d.GetNOpenElements(); i++) {
                        Element2d hel = mesh3d.OpenElement(i);
                        meshing.AddBoundaryElement(hel);
                    }

                    oldne = mesh3d.GetNE();

                    meshing.GenerateMesh(mesh3d, mpquad);

                    for (int i = oldne + 1; i <= mesh3d.GetNE(); i++)
                        mesh3d.VolumeElement(i).SetIndex(k);

                    LOG_DEBUG("mesh has " << mesh3d.GetNE() << " prism/pyramid elements");

                    mesh3d.FindOpenElements(k);
                }
            }

            if (mesh3d.HasOpenQuads()) {
                LOG_ERROR("mesh has still open quads");
                throw std::runtime_error("Stop meshing since too many attempts");
            }

            if (mp.delaunay && mesh3d.GetNOpenElements()) {
                Meshing3 meshing((const char**) NULL);

                mesh3d.FindOpenElements(k);


                for (PointIndex pi = mesh3d.Points().Begin(); pi < mesh3d.Points().End(); pi++)
                    meshing.AddPoint(mesh3d[pi], pi);

                for (int i = 1; i <= mesh3d.GetNOpenElements(); i++)
                    meshing.AddBoundaryElement(mesh3d.OpenElement(i));

                oldne = mesh3d.GetNE();

                meshing.Delaunay(mesh3d, k, mp);

                for (int i = oldne + 1; i <= mesh3d.GetNE(); i++)
                    mesh3d.VolumeElement(i).SetIndex(k);

                LOG_DEBUG(mesh3d.GetNP() << " points, " << mesh3d.GetNE() << " elements");
            }

            int cntsteps = 0;
            if (mesh3d.GetNOpenElements())
                do {
                    if (multithread.terminate)
                        break;

                    mesh3d.FindOpenElements(k);
                    LOG_DEBUG(mesh3d.GetNOpenElements() << " open faces");
                    cntsteps++;

                    if (cntsteps > mp.maxoutersteps)
                        throw std::runtime_error("Stop meshing since too many attempts");

                    std::string rulefile = "./tetra.rls";
                    LOG_DEBUG("start tetmeshing");

                    //	  Meshing3 meshing(rulefile);
                    Meshing3 meshing(tetrules);

                    Array<int, PointIndex::BASE> glob2loc(mesh3d.GetNP());
                    glob2loc = -1;

                    for (PointIndex pi = mesh3d.Points().Begin(); pi < mesh3d.Points().End(); pi++)
                        if (domain_bbox.IsIn(mesh3d[pi]))
                            glob2loc[pi] =
                                meshing.AddPoint(mesh3d[pi], pi);

                    for (int i = 1; i <= mesh3d.GetNOpenElements(); i++) {
                        Element2d hel = mesh3d.OpenElement(i);
                        for (int j = 0; j < hel.GetNP(); j++)
                            hel[j] = glob2loc[hel[j]];
                        meshing.AddBoundaryElement(hel);
                        // meshing.AddBoundaryElement (mesh3d.OpenElement(i));
                    }

                    oldne = mesh3d.GetNE();

                    mp.giveuptol = 15 + 10 * cntsteps;
                    mp.sloppy = 5;
                    meshing.GenerateMesh(mesh3d, mp);

                    for (ElementIndex ei = oldne; ei < mesh3d.GetNE(); ei++)
                        mesh3d[ei].SetIndex(k);


                    mesh3d.CalcSurfacesOfNode();
                    mesh3d.FindOpenElements(k);

                    // teterrpow = 2;
                    if (mesh3d.GetNOpenElements() != 0) {
                        meshed = 0;
                        LOG_DEBUG(mesh3d.GetNOpenElements() << " open faces found");

                        MeshOptimize3d optmesh(mp);

                        const char * optstr = "mcmstmcmstmcmstmcm";
                        for (size_t j = 1; j <= strlen(optstr); j++) {
                            mesh3d.CalcSurfacesOfNode();
                            mesh3d.FreeOpenElementsEnvironment(2);
                            mesh3d.CalcSurfacesOfNode();

                            switch (optstr[j - 1]) {
                                case 'c': optmesh.CombineImprove(mesh3d, OPT_REST);
                                    break;
                                case 'd': optmesh.SplitImprove(mesh3d, OPT_REST);
                                    break;
                                case 's': optmesh.SwapImprove(mesh3d, OPT_REST);
                                    break;
                                case 't': optmesh.SwapImprove2(mesh3d, OPT_REST);
                                    break;
                                case 'm': mesh3d.ImproveMesh(mp, OPT_REST);
                                    break;
                            }

                        }

                        mesh3d.FindOpenElements(k);
                        LOG_DEBUG("Call remove problem");
                        RemoveProblem(mesh3d, k);
                        mesh3d.FindOpenElements(k);
                    }
                    else {
                        meshed = 1;
                        LOG_INFO("Success !");
                    }
                } while (!meshed);

            LOG_DEBUG(mesh3d.GetNP() << " points, " << mesh3d.GetNE() << " elements");
        }

        mp.maxh = globmaxh;

        MeshQuality3d(mesh3d);

        return MESHING3_OK;
    }

    MESHING3_RESULT OptimizeVolume(MeshingParameters & mp, Mesh & mesh3d)
    {
        int i;

        LOG_DEBUG("Volume Optimization");

        mesh3d.CalcSurfacesOfNode();
        for (i = 1; i <= mp.optsteps3d; i++) {
            if (multithread.terminate)
                break;

            MeshOptimize3d optmesh(mp);

            // teterrpow = mp.opterrpow;
            for (size_t j = 1; j <= strlen(mp.optimize3d); j++) {
                if (multithread.terminate)
                    break;

                switch (mp.optimize3d[j - 1]) {
                    case 'c': optmesh.CombineImprove(mesh3d, OPT_REST);
                        break;
                    case 'd': optmesh.SplitImprove(mesh3d);
                        break;
                    case 's': optmesh.SwapImprove(mesh3d);
                        break;
                    case 't': optmesh.SwapImprove2(mesh3d);
                        break;
                    case 'm': mesh3d.ImproveMesh(mp);
                        break;
                    case 'M': mesh3d.ImproveMesh(mp);
                        break;
                    case 'j': mesh3d.ImproveMeshJacobian(mp);
                        break;
                }
            }
            mesh3d.mglevels = 1;
            MeshQuality3d(mesh3d);
        }

        return MESHING3_OK;
    }

    void RemoveIllegalElements(Mesh & mesh3d)
    {
        int it = 10;
        int nillegal, oldn;

        LOG_DEBUG("Remove Illegal Elements");
        mesh3d.CalcSurfacesOfNode();

        nillegal = mesh3d.MarkIllegalElements();

        MeshingParameters dummymp;
        MeshOptimize3d optmesh(dummymp);
        while (nillegal && (it--) > 0) {
            if (multithread.terminate)
                break;

            LOG_DEBUG(nillegal << " illegal tets");
            optmesh.SplitImprove(mesh3d, OPT_LEGAL);

            mesh3d.MarkIllegalElements(); // test
            optmesh.SwapImprove(mesh3d, OPT_LEGAL);
            mesh3d.MarkIllegalElements(); // test
            optmesh.SwapImprove2(mesh3d, OPT_LEGAL);

            oldn = nillegal;
            nillegal = mesh3d.MarkIllegalElements();

            if (oldn != nillegal)
                it = 10;
        }
        LOG_DEBUG(nillegal << " illegal tets");
    }
}
