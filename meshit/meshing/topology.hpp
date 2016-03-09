#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

/**************************************************************************/
/* File:   topology.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   27. Apr. 01                                                    */
/**************************************************************************/

#include "../general/array.hpp"
#include "meshclass.hpp"

/*
    Mesh topology
    (Elements, Faces, Edges, Vertices
 */

namespace meshit {

    class MeshTopology
    {
        const Mesh& mesh;
        bool buildedges;
        bool buildfaces;

        std::vector<INDEX_2> edge2vert;
        std::vector<INDEX_4> face2vert;
        TABLE<int, PointIndex::BASE>* vert2surfelement;
        TABLE<int, PointIndex::BASE>* vert2segment;
        int timestamp;

     public:
        MeshTopology(const Mesh& amesh);
        ~MeshTopology();

        void Update();

        int GetNEdges() const
        {
            return edge2vert.size();
        }

        int GetNFaces() const
        {
            return face2vert.size();
        }

        static inline int GetNVertices(ELEMENT_TYPE et);
        static inline int GetNEdges(ELEMENT_TYPE et);

        inline static const ELEMENT_EDGE* GetEdges(ELEMENT_TYPE et);
        inline static const ELEMENT_FACE* GetFaces(ELEMENT_TYPE et);

        void GetSurfaceElementEdges(int elnr, Array<int>& edges) const;

    };

    inline int MeshTopology::GetNVertices(ELEMENT_TYPE et)
    {
        switch (et) {
            case SEGMENT:
            case SEGMENT3:
                return 2;

            case TRIG:
            case TRIG6:
                return 3;

            case QUAD:
            case QUAD6:
            case QUAD8:
                return 4;

            default:
                MESHIT_LOG_ERROR("MeshTopology::GetNVertices, illegal element type " << et);
                return 0;
        }
    }

    inline int MeshTopology::GetNEdges(ELEMENT_TYPE et)
    {
        switch (et) {
            case SEGMENT:
            case SEGMENT3:
                return 1;

            case TRIG:
            case TRIG6:
                return 3;

            case QUAD:
            case QUAD6:
            case QUAD8:
                return 4;

            default:
                MESHIT_LOG_ERROR("MeshTopology::GetNEdges, illegal element type " << et);
                return 0;
        }
    }

    const ELEMENT_EDGE* MeshTopology::GetEdges(ELEMENT_TYPE et)
    {
        static int segm_edges[1][2] = {
                {0, 1}
        };

        static int trig_edges[3][2] = {
                {2, 0},
                {1, 2},
                {0, 1}
        };

        static int quad_edges[4][2] = {
                {0, 1},
                {2, 3},
                {3, 0},
                {1, 2}
        };

        switch (et) {
            case SEGMENT:
            case SEGMENT3:
                return segm_edges;

            case TRIG:
            case TRIG6:
                return trig_edges;

            case QUAD:
            case QUAD6:
            case QUAD8:
                return quad_edges;

            default:
                MESHIT_LOG_ERROR("MeshTopology::GetEdges, illegal element type " << et);
                return nullptr;
        }
    }

    inline const ELEMENT_FACE* MeshTopology::GetFaces(ELEMENT_TYPE et)
    {
        static const int trig_faces[1][4] = {
                {1, 2, 3, 0}
        };
        static const int quad_faces[1][4] = {
                {1, 2, 3, 4}
        };

        switch (et) {
            case TRIG:
            case TRIG6:
                return trig_faces;

            case QUAD:
            case QUAD6:
            case QUAD8:
                return quad_faces;

            default:
                MESHIT_LOG_ERROR("MeshTopology::GetFaces, illegal element type " << et);
                return nullptr;
        }
    }

}

#endif
