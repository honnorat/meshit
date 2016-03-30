#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

/**************************************************************************/
/* File:   topology.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   27. Apr. 01                                                    */
/**************************************************************************/

#include "meshclass.hpp"

/*
    Mesh topology
    (Elements, Faces, Edges, Vertices
 */

namespace meshit
{
    class Mesh;

    class MeshTopology
    {
        Mesh& mesh;
        bool buildedges;
        bool buildfaces;

        std::vector<INDEX_2> edge2vert;
        std::vector<INDEX_4> face2vert;
        TABLE<int>* vert2surfelement;
        TABLE<size_t>* vert2segment;

     public:
        explicit MeshTopology(Mesh& amesh);
        ~MeshTopology();

        void Update();

        size_t GetNEdges() const
        {
            return edge2vert.size();
        }

        size_t GetNFaces() const
        {
            return face2vert.size();
        }
    };

}  // namespace meshit

#endif
