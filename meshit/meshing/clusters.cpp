#include <meshit.hpp>

#include "clusters.hpp"

namespace meshit {

    AnisotropicClusters::AnisotropicClusters(const Mesh & amesh)
        : mesh(amesh) { }

    AnisotropicClusters::~AnisotropicClusters() { }

    void AnisotropicClusters::Update()
    {
        const MeshTopology & top = mesh.GetTopology();

        bool hasedges = top.HasEdges();
        bool hasfaces = top.HasFaces();

        if (!hasedges || !hasfaces) return;

        nv = mesh.GetNV();
        ned = top.GetNEdges();
        nfa = top.GetNFaces();
        int nse = mesh.GetNSE();

        cluster_reps.resize(nv + ned + nfa);
        cluster_reps = -1;

        Array<int> llist(nv + ned + nfa);
        llist = 0;

        Array<int> nnums, ednums, fanums;

        for (int i = 1; i <= nse; i++) {
            const Element2d & el = mesh.SurfaceElement(i);
            ELEMENT_TYPE typ = el.GetType();

            top.GetSurfaceElementEdges(i, ednums);
            int fanum = top.GetSurfaceElementFace(i);

            int elnv = top.GetNVertices(typ);
            int elned = ednums.size();

            nnums.resize(elnv + elned + 1);
            for (int j = 1; j <= elnv; j++)
                nnums.Elem(j) = el.PNum(j);
            for (int j = 1; j <= elned; j++)
                nnums.Elem(elnv + j) = nv + ednums.Elem(j);
            nnums.Elem(elnv + elned + 1) = fanum;

            for (int j = 0; j < nnums.size(); j++)
                cluster_reps.Elem(nnums[j]) = nnums[j];
        }
    }
}
