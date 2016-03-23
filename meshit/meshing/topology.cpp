#include "topology.hpp"

namespace meshit
{
    MeshTopology::MeshTopology(const Mesh& amesh)
        : mesh(amesh)
    {
        buildedges = 1;
        buildfaces = 1;
        vert2surfelement = nullptr;
        vert2segment = nullptr;
    }

    MeshTopology::~MeshTopology()
    {
        delete vert2surfelement;
        delete vert2segment;
    }

    void MeshTopology::Update()
    {
        size_t nse = mesh.GetNSE();
        size_t nseg = mesh.GetNSeg();
        size_t np = mesh.GetNP();
        size_t nv = mesh.GetNV();

        MESHIT_LOG_DEBUG(" UPDATE MESH TOPOLOGY ");
        MESHIT_LOG_DEBUG("nse  = " << nse);
        MESHIT_LOG_DEBUG("nseg = " << nseg);
        MESHIT_LOG_DEBUG("np   = " << np);
        MESHIT_LOG_DEBUG("nv   = " << nv);

        delete vert2surfelement;
        delete vert2segment;

        /*
          generate:
          vertex to element
          vertex to surface element
          vertex to segment
         */
        std::vector<int> cnt(nv, 0);
        for (size_t i = 0; i < nse; i++) {
            const Element2d& el = mesh.SurfaceElement(i);
            for (size_t j = 0; j < 3; j++) {
                cnt[el[j]]++;
            }
        }

        vert2surfelement = new TABLE<int>(cnt);
        for (size_t i = 0; i < nse; i++) {
            const Element2d& el = mesh.SurfaceElement(i);
            for (size_t j = 0; j < 3; j++) {
                vert2surfelement->AddSave(el[j], i + 1);
            }
        }

        cnt.assign(nv, 0);
        for (size_t i = 0; i < nseg; i++) {
            const Segment& seg = mesh.LineSegment(i);
            cnt[seg[0]]++;
            cnt[seg[1]]++;
        }

        vert2segment = new TABLE<size_t>(cnt);
        for (size_t i = 0; i < nseg; i++) {
            const Segment& seg = mesh.LineSegment(i);
            vert2segment->AddSave(seg[0], i);
            vert2segment->AddSave(seg[1], i);
        }

        if (buildedges) {
            // keep existing edges
            cnt.assign(nv, 0);
            for (size_t i = 0; i < edge2vert.size(); i++) {
                cnt[edge2vert[i][0]]++;
            }
            TABLE<size_t> vert2edge(cnt);
            for (size_t i = 0; i < edge2vert.size(); i++) {
                vert2edge.AddSave(edge2vert[i][0], i);
            }

            // ensure all coarse grid and intermediate level edges
            cnt.assign(nv, 0);
            for (size_t i = 0; i < mesh.mlbetweennodes.size(); i++) {
                INDEX_2 parents = Sort(mesh.mlbetweennodes[i]);
                if (parents[0] >= 0)
                    cnt[parents[0]]++;
            }
            TABLE<size_t> vert2vertcoarse(cnt);
            for (size_t i = 0; i < mesh.mlbetweennodes.size(); i++) {
                INDEX_2 parents = Sort(mesh.mlbetweennodes[i]);
                if (parents[0] > 0)
                    vert2vertcoarse.AddSave(parents[0], parents[1]);
            }

            std::vector<int> edgenr(nv);
            std::vector<int> edgeflag(nv, -1);
            std::vector<int> vertex2;

            size_t ned = edge2vert.size();

            for (size_t i = 0; i < nv; i++) {
                vertex2.resize(0);

                for (size_t j = 0; j < vert2edge[i].size(); j++) {
                    int ednr = vert2edge[i][j];
                    int i2 = edge2vert[ednr][1];
                    edgeflag[i2] = i;
                    edgenr[i2] = ednr;
                }

                for (size_t j = 0; j < vert2vertcoarse[i].size(); j++) {
                    int v2 = vert2vertcoarse[i][j];
                    if (edgeflag[v2] < static_cast<int>(i)) {
                        edgeflag[v2] = i;
                        vertex2.push_back(v2);
                    }
                }

                for (size_t j = 0; j < (*vert2surfelement)[i].size(); j++) {
                    int elnr = (*vert2surfelement)[i][j];
                    const Element2d& el = mesh.SurfaceElement(elnr - 1);

                    const ELEMENT_EDGE eledges[3] = {2, 0,
                                                     1, 2,
                                                     0, 1};
                    for (size_t k = 0; k < 3; k++) {
                        INDEX_2 edge(el[eledges[k][0]], el[eledges[k][1]]);
                        edge.Sort();
                        if (edge.I1() != static_cast<int>(i)) continue;

                        if (edgeflag[edge.I2()] < static_cast<int>(i)) {
                            vertex2.push_back(edge.I2());
                            edgeflag[edge.I2()] = static_cast<int>(i);
                        }
                    }
                }

                for (size_t j = 0; j < (*vert2segment)[i].size(); j++) {
                    int elnr = (*vert2segment)[i][j];
                    const Segment& el = mesh.LineSegment(elnr);

                    INDEX_2 edge(el[0], el[1]);
                    edge.Sort();
                    if (edge.I1() != static_cast<int>(i)) continue;

                    if (edgeflag[edge.I2()] < static_cast<int>(i)) {
                        vertex2.push_back(edge.I2());
                        edgeflag[edge.I2()] = static_cast<int>(i);
                    }
                }

                std::sort(vertex2.begin(), vertex2.end());
                for (size_t j = 0; j < vertex2.size(); j++) {
                    edgenr[vertex2[j]] = ++ned;
                    edge2vert.push_back(INDEX_2(i, vertex2[j]));
                }

                for (size_t j = 0; j < (*vert2surfelement)[i].size(); j++) {
                    int elnr = (*vert2surfelement)[i][j];
                    const Element2d& el = mesh.SurfaceElement(elnr - 1);
                    const ELEMENT_EDGE eledges[3] = {2, 0,
                                                     1, 2,
                                                     0, 1};
                    for (size_t k = 0; k < 3; k++) {
                        INDEX_2 edge(el[eledges[k][0]], el[eledges[k][1]]);

                        int edgedir = (edge.I1() > edge.I2());
                        if (edgedir) std::swap(edge.I1(), edge.I2());
                        if (edge.I1() != static_cast<int>(i)) continue;
                    }
                }

                for (size_t j = 0; j < (*vert2segment)[i].size(); j++) {
                    int elnr = (*vert2segment)[i][j];
                    const Segment& el = mesh.LineSegment(elnr);

                    INDEX_2 edge(el[0], el[1]);

                    int edgedir = (edge.I1() > edge.I2());
                    if (edgedir) std::swap(edge.I1(), edge.I2());
                    if (edge.I1() != static_cast<int>(i)) continue;
                }
            }
        }

        // generate faces
        if (buildfaces) {
            size_t oldnfa = face2vert.size();
            size_t nfa = oldnfa;

            cnt.assign(nv, 0);
            for (size_t i = 0; i < face2vert.size(); i++) {
                cnt[face2vert[i][0]]++;
            }
            TABLE<size_t> vert2oldface(cnt);
            for (size_t i = 0; i < face2vert.size(); i++)
                vert2oldface.AddSave(face2vert[i][0], i);

            size_t max_face_on_vertex = 0;
            for (size_t i = 0; i < nv; i++) {
                size_t onv = vert2oldface[i].size() + (*vert2surfelement)[i].size();
                max_face_on_vertex = std::max(onv, max_face_on_vertex);
            }

            for (size_t v = 0; v < nv; v++) {
                size_t first_fa = nfa;

                INDEX_3_CLOSED_HASHTABLE<int> vert2face(2 * max_face_on_vertex + 10);

                for (size_t j = 0; j < vert2oldface[v].size(); j++) {
                    int fnr = vert2oldface[v][j];
                    INDEX_3 face(face2vert[fnr].I1(),
                                 face2vert[fnr].I2(),
                                 face2vert[fnr].I3());
                    vert2face.Set(face, fnr + 1);
                }

                for (size_t pass = 1; pass <= 2; pass++) {
                    if (pass == 2) {
                        for (size_t j = first_fa; j < face2vert.size(); j++) {
                            if (face2vert[j][0] == static_cast<int>(v)) {
                                INDEX_3 face(face2vert[j].I1(),
                                             face2vert[j].I2(),
                                             face2vert[j].I3());
                                vert2face.Set(face, j + 1);
                            } else {
                                break;
                            }
                        }
                    }

                    for (size_t j = 0; j < (*vert2surfelement)[v].size(); j++) {
                        int elnr = (*vert2surfelement)[v][j];
                        const Element2d& el = mesh.SurfaceElement(elnr - 1);

                        INDEX_3 face(el.PNum(1), el.PNum(2), el.PNum(3));
                        if (face.I1() > face.I2()) {
                            std::swap(face.I1(), face.I2());
                        }
                        if (face.I2() > face.I3()) {
                            std::swap(face.I2(), face.I3());
                        }
                        if (face.I1() > face.I2()) {
                            std::swap(face.I1(), face.I2());
                        }

                        if (face.I1() != static_cast<int>(v)) continue;

                        if (!vert2face.Used(face)) {
                            nfa++;
                            vert2face.Set(face, nfa);

                            INDEX_4 hface(face.I1(), face.I2(), face.I3(), 0);
                            face2vert.push_back(hface);
                        }
                    }

                    // sort faces
                    if (pass == 1) {
                        for (size_t i = 0; i < std::max(0UL, nfa - first_fa); i++) {
                            for (size_t j = first_fa + 1; j < nfa - i; j++) {
                                if (face2vert[j] < face2vert[j - 1]) {
                                    std::swap(face2vert[j - 1], face2vert[j]);
                                }
                            }
                        }
                    }
                }
            }
            face2vert.reserve(nfa);
        }
    }
}  // namespace meshit
