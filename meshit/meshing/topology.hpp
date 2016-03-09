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

    struct T_EDGE
    {
        int nr; // 0-based
    };

    struct T_FACE
    {
        int fnr; // 0-based
    };

    class MeshTopology
    {
        const Mesh& mesh;
        bool buildedges;
        bool buildfaces;

        Array<INDEX_2> edge2vert;
        Array<INDEX_4> face2vert;
        Array<T_EDGE[12]> edges;
        Array<T_FACE[6]> faces;
        Array<T_EDGE[4]> surfedges;
        Array<T_EDGE> segedges;
        Array<T_FACE> surffaces;
        Array<INDEX_2> surf2volelement;
        Array<int> face2surfel;
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
        static inline int GetNPoints(ELEMENT_TYPE et);
        static inline int GetNEdges(ELEMENT_TYPE et);
        static inline int GetNFaces(ELEMENT_TYPE et);

        static const Point3d* GetVertices(ELEMENT_TYPE et);
        inline static const ELEMENT_EDGE* GetEdges1(ELEMENT_TYPE et);
        inline static const ELEMENT_EDGE* GetEdges0(ELEMENT_TYPE et);
        inline static const ELEMENT_FACE* GetFaces1(ELEMENT_TYPE et);
        inline static const ELEMENT_FACE* GetFaces0(ELEMENT_TYPE et);

        void GetFaceVertices(int fnr, Array<int>& vertices) const;

        void GetSurfaceElementEdges(int elnr, Array<int>& edges) const;
        int GetSurfaceElementFace(int elnr) const;

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

            case TET:
            case TET10:
                return 4;

            case PYRAMID:
                return 5;

            case PRISM:
            case PRISM12:
                return 6;

            case HEX:
                return 8;

            default:
                MESHIT_LOG_ERROR("MeshTopology::GetNVertices, illegal element type " << et);
                return 0;
        }
    }

    inline int MeshTopology::GetNPoints(ELEMENT_TYPE et)
    {
        switch (et) {
            case SEGMENT:
                return 2;
            case SEGMENT3:
                return 3;

            case TRIG:
                return 3;
            case TRIG6:
                return 6;

            case QUAD:
            case QUAD6:
            case QUAD8:
                return 4;

            case TET:
                return 4;
            case TET10:
                return 10;

            case PYRAMID:
                return 5;

            case PRISM:
            case PRISM12:
                return 6;

            case HEX:
                return 8;

                // default:
                // std::cerr << "Ng_ME_GetNVertices, illegal element type " << et <<std::endl;
        }
        return 0;
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

            case TET:
            case TET10:
                return 6;

            case PYRAMID:
                return 8;

            case PRISM:
            case PRISM12:
                return 9;

            case HEX:
                return 12;

                // default:
                // std::cerr << "Ng_ME_GetNEdges, illegal element type " << et <<std::endl;
        }
        return 0;
    }

    inline int MeshTopology::GetNFaces(ELEMENT_TYPE et)
    {
        switch (et) {
            case SEGMENT:
            case SEGMENT3:
                return 0;

            case TRIG:
            case TRIG6:
                return 1;

            case QUAD:
            case QUAD6:
            case QUAD8:
                return 1;

            case TET:
            case TET10:
                return 4;

            case PYRAMID:
                return 5;

            case PRISM:
            case PRISM12:
                return 5;

            case HEX:
                return 6;

                // default:
                // std::cerr << "Ng_ME_GetNVertices, illegal element type " << et <<std::endl;
        }
        return 0;
    }

    const ELEMENT_EDGE* MeshTopology::GetEdges1(ELEMENT_TYPE et)
    {
        static int segm_edges[1][2] = {
                {1, 2}
        };

        static int trig_edges[3][2] = {
                {3, 1},
                {2, 3},
                {1, 2}
        };

        static int quad_edges[4][2] = {
                {1, 2},
                {3, 4},
                {4, 1},
                {2, 3}
        };


        static int tet_edges[6][2] = {
                {4, 1},
                {4, 2},
                {4, 3},
                {1, 2},
                {1, 3},
                {2, 3}
        };

        static int prism_edges[9][2] = {
                {3, 1},
                {1, 2},
                {3, 2},
                {6, 4},
                {4, 5},
                {6, 5},
                {3, 6},
                {1, 4},
                {2, 5}
        };

        static int pyramid_edges[8][2] = {
                {1, 2},
                {2, 3},
                {1, 4},
                {4, 3},
                {1, 5},
                {2, 5},
                {3, 5},
                {4, 5}
        };

        static int hex_edges[12][2] = {
                {1, 2},
                {3, 4},
                {4, 1},
                {2, 3},
                {5, 6},
                {7, 8},
                {8, 5},
                {6, 7},
                {1, 5},
                {2, 6},
                {3, 7},
                {4, 8},
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

            case TET:
            case TET10:
                return tet_edges;

            case PYRAMID:
                return pyramid_edges;

            case PRISM:
            case PRISM12:
                return prism_edges;

            case HEX:
                return hex_edges;
                // default:
                // std::cerr << "Ng_ME_GetEdges, illegal element type " << et <<std::endl;
        }
        return 0;
    }

    const ELEMENT_EDGE* MeshTopology::GetEdges0(ELEMENT_TYPE et)
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


        static int tet_edges[6][2] = {
                {3, 0},
                {3, 1},
                {3, 2},
                {0, 1},
                {0, 2},
                {1, 2}
        };

        static int prism_edges[9][2] = {
                {2, 0},
                {0, 1},
                {2, 1},
                {5, 3},
                {3, 4},
                {5, 4},
                {2, 5},
                {0, 3},
                {1, 4}
        };

        static int pyramid_edges[8][2] = {
                {0, 1},
                {1, 2},
                {0, 3},
                {3, 2},
                {0, 4},
                {1, 4},
                {2, 4},
                {3, 4}
        };

        static int hex_edges[12][2] = {
                {0, 1},
                {2, 3},
                {3, 0},
                {1, 2},
                {4, 5},
                {6, 7},
                {7, 4},
                {5, 6},
                {0, 4},
                {1, 5},
                {2, 6},
                {3, 7},
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

            case TET:
            case TET10:
                return tet_edges;

            case PYRAMID:
                return pyramid_edges;

            case PRISM:
            case PRISM12:
                return prism_edges;

            case HEX:
                return hex_edges;
                // default:
                // std::cerr << "Ng_ME_GetEdges, illegal element type " << et <<std::endl;
        }
        return 0;
    }

    inline const ELEMENT_FACE* MeshTopology::GetFaces1(ELEMENT_TYPE et)
    {
        static const int trig_faces[1][4] = {
                {1, 2, 3, 0}
        };
        static const int quad_faces[1][4] = {
                {1, 2, 3, 4}
        };

        static const int tet_faces[4][4] = {
                {4, 2, 3, 0},
                {4, 3, 1, 0},
                {4, 1, 2, 0},
                {1, 3, 2, 0}
        };

        static const int prism_faces[5][4] = {
                {1, 3, 2, 0},
                {4, 5, 6, 0},
                {3, 1, 4, 6},
                {1, 2, 5, 4},
                {2, 3, 6, 5}
        };

        static const int pyramid_faces[5][4] = {
                {1, 2, 5, 0},
                {2, 3, 5, 0},
                {3, 4, 5, 0},
                {4, 1, 5, 0},
                {1, 4, 3, 2}
        };

        static const int hex_faces[6][4] = {
                {1, 4, 3, 2},
                {5, 6, 7, 8},
                {1, 2, 6, 5},
                {2, 3, 7, 6},
                {3, 4, 8, 7},
                {4, 1, 5, 8}
        };


        switch (et) {
            case TRIG:
            case TRIG6:
                return trig_faces;

            case QUAD:
            case QUAD6:
            case QUAD8:
                return quad_faces;


            case TET:
            case TET10:
                return tet_faces;

            case PRISM:
            case PRISM12:
                return prism_faces;

            case PYRAMID:
                return pyramid_faces;

            case SEGMENT:
            case SEGMENT3:

            case HEX:
                return hex_faces;

                // default:
                // std::cerr << "Ng_ME_GetVertices, illegal element type " << et <<std::endl;
        }
        return 0;
    }

    inline const ELEMENT_FACE* MeshTopology::GetFaces0(ELEMENT_TYPE et)
    {
        static const int trig_faces[1][4] = {
                {0, 1, 2, -1}
        };
        static const int quad_faces[1][4] = {
                {0, 1, 2, 3}
        };

        static const int tet_faces[4][4] = {
                {3, 1, 2, -1},
                {3, 2, 0, -1},
                {3, 0, 1, -1},
                {0, 2, 1, -1}
        };

        static const int prism_faces[5][4] = {
                {0, 2, 1, -1},
                {3, 4, 5, -1},
                {2, 0, 3, 5},
                {0, 1, 4, 3},
                {1, 2, 5, 4}
        };

        static const int pyramid_faces[5][4] = {
                {0, 1, 4, -1},
                {1, 2, 4, -1},
                {2, 3, 4, -1},
                {3, 0, 4, -1},
                {0, 3, 2, 1}
        };

        static const int hex_faces[6][4] = {
                {0, 3, 2, 1},
                {4, 5, 6, 7},
                {0, 1, 5, 4},
                {1, 2, 6, 5},
                {2, 3, 7, 6},
                {3, 0, 4, 7}
        };


        switch (et) {
            case TRIG:
            case TRIG6:
                return trig_faces;

            case QUAD:
            case QUAD6:
            case QUAD8:
                return quad_faces;


            case TET:
            case TET10:
                return tet_faces;

            case PRISM:
            case PRISM12:
                return prism_faces;

            case PYRAMID:
                return pyramid_faces;

            case SEGMENT:
            case SEGMENT3:

            case HEX:
                return hex_faces;

                // default:
                // std::cerr << "Ng_ME_GetVertices, illegal element type " << et <<std::endl;
        }
        return 0;
    }

}

#endif
