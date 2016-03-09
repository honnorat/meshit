#include "topology.hpp"

namespace meshit {

    template<class T>
    void QuickSortRec(FlatArray<T>& data,
                      int left, int right)
    {
        int i = left;
        int j = right;
        T midval = data[(left + right) / 2];

        do {
            while (data[i] < midval) i++;
            while (midval < data[j]) j--;

            if (i <= j) {
                std::swap(data[i], data[j]);
                i++;
                j--;
            }
        } while (i <= j);
        if (left < j) QuickSortRec(data, left, j);
        if (i < right) QuickSortRec(data, i, right);
    }

    template<class T>
    void QuickSort(FlatArray<T>& data)
    {
        if (data.size() > 1)
            QuickSortRec(data, 0, data.size() - 1);
    }

    MeshTopology::MeshTopology(const Mesh& amesh)
            : mesh(amesh)
    {
        buildedges = 1;
        buildfaces = 1;
        vert2surfelement = 0;
        vert2segment = 0;
        timestamp = -1;
    }

    MeshTopology::~MeshTopology()
    {
        delete vert2surfelement;
        delete vert2segment;
    }

    void MeshTopology::Update()
    {
        if (timestamp > mesh.GetTimeStamp()) return;

        int nse = mesh.GetNSE();
        int nseg = mesh.GetNSeg();
        int np = mesh.GetNP();
        int nv = mesh.GetNV();
        int nfa = 0;
        int ned = edge2vert.size();

        std::cerr << " UPDATE MESH TOPOLOGY " << std::endl;
        std::cerr << "nse  = " << nse << std::endl;
        std::cerr << "nseg = " << nseg << std::endl;
        std::cerr << "np   = " << np << std::endl;
        std::cerr << "nv   = " << nv << std::endl;

        delete vert2surfelement;
        delete vert2segment;

        Array<int, PointIndex::BASE> cnt(nv);
        Array<int> vnums;

        /*
          generate:
          vertex to element 
          vertex to surface element
          vertex to segment 
         */
        cnt = 0;
        for (SurfaceElementIndex sei = 0; sei < nse; sei++) {
            const Element2d& el = mesh.SurfaceElement(sei);
            for (int j = 0; j < el.GetNV(); j++)
                cnt[el[j]]++;
        }

        vert2surfelement = new TABLE<int, PointIndex::BASE>(cnt);
        for (SurfaceElementIndex sei = 0; sei < nse; sei++) {
            const Element2d& el = mesh.SurfaceElement(sei);
            for (int j = 0; j < el.GetNV(); j++)
                vert2surfelement->AddSave(el[j], sei + 1);
        }

        cnt = 0;
        for (int i = 1; i <= nseg; i++) {
            const Segment& seg = mesh.LineSegment(i);
            cnt[seg[0]]++;
            cnt[seg[1]]++;
        }

        vert2segment = new TABLE<int, PointIndex::BASE>(cnt);
        for (int i = 1; i <= nseg; i++) {
            const Segment& seg = mesh.LineSegment(i);
            vert2segment->AddSave(seg[0], i);
            vert2segment->AddSave(seg[1], i);
        }

        if (buildedges) {
            segedges.resize(nseg);

            // keep existing edges
            cnt = 0;
            for (int i = 0; i < edge2vert.size(); i++)
                cnt[edge2vert[i][0]]++;
            TABLE<int, PointIndex::BASE> vert2edge(cnt);
            for (int i = 0; i < edge2vert.size(); i++)
                vert2edge.AddSave(edge2vert[i][0], i + 1);

            // ensure all coarse grid and intermediate level edges
            cnt = 0;
            for (int i = mesh.mlbetweennodes.Begin(); i < mesh.mlbetweennodes.End(); i++) {
                INDEX_2 parents = Sort(mesh.mlbetweennodes[i]);
                if (parents[0] >= PointIndex::BASE) cnt[parents[0]]++;
            }
            TABLE<int, PointIndex::BASE> vert2vertcoarse(cnt);
            for (int i = mesh.mlbetweennodes.Begin(); i < mesh.mlbetweennodes.End(); i++) {
                INDEX_2 parents = Sort(mesh.mlbetweennodes[i]);
                if (parents[0] > PointIndex::BASE) vert2vertcoarse.AddSave(parents[0], parents[1]);
            }

            Array<int, PointIndex::BASE> edgenr(nv);
            Array<int, PointIndex::BASE> edgeflag(nv);
            Array<int> vertex2;

            edgeflag = PointIndex::BASE - 1;

            ned = edge2vert.size();

            for (int i = PointIndex::BASE; i < nv + PointIndex::BASE; i++) {
                vertex2.resize(0);

                for (int j = 0; j < vert2edge[i].size(); j++) {
                    int ednr = vert2edge[i][j];
                    int i2 = edge2vert.Get(ednr)[1];
                    edgeflag[i2] = i;
                    edgenr[i2] = ednr;
                }

                for (int j = 0; j < vert2vertcoarse[i].size(); j++) {
                    int v2 = vert2vertcoarse[i][j];
                    if (edgeflag[v2] < i) {
                        edgeflag[v2] = i;
                        vertex2.push_back(v2);
                    }
                }

                for (int j = 0; j < (*vert2surfelement)[i].size(); j++) {
                    int elnr = (*vert2surfelement)[i][j];
                    const Element2d& el = mesh.SurfaceElement(elnr);

                    int neledges = GetNEdges(el.GetType());
                    const ELEMENT_EDGE* eledges = GetEdges(el.GetType());

                    for (int k = 0; k < neledges; k++) {
                        INDEX_2 edge(el[eledges[k][0]], el[eledges[k][1]]);
                        edge.Sort();
                        if (edge.I1() != i) continue;

                        if (edgeflag[edge.I2()] < i) {
                            vertex2.push_back(edge.I2());
                            edgeflag[edge.I2()] = i;
                        }
                    }
                }

                for (int j = 0; j < (*vert2segment)[i].size(); j++) {
                    int elnr = (*vert2segment)[i][j];
                    const Segment& el = mesh.LineSegment(elnr);

                    INDEX_2 edge(el[0], el[1]);
                    edge.Sort();
                    if (edge.I1() != i) continue;

                    if (edgeflag[edge.I2()] < i) {
                        vertex2.push_back(edge.I2());
                        edgeflag[edge.I2()] = i;
                    }
                }

                QuickSort(vertex2);
                for (int j = 0; j < vertex2.size(); j++) {
                    edgenr[vertex2[j]] = ++ned;
                    edge2vert.push_back(INDEX_2(i, vertex2[j]));
                }

                for (int j = 0; j < (*vert2surfelement)[i].size(); j++) {
                    int elnr = (*vert2surfelement)[i][j];
                    const Element2d& el = mesh.SurfaceElement(elnr);

                    int neledges = GetNEdges(el.GetType());
                    const ELEMENT_EDGE* eledges = GetEdges(el.GetType());

                    for (int k = 0; k < neledges; k++) {
                        INDEX_2 edge(el[eledges[k][0]], el[eledges[k][1]]);

                        int edgedir = (edge.I1() > edge.I2());
                        if (edgedir) std::swap(edge.I1(), edge.I2());

                        if (edge.I1() != i) continue;
                    }
                }

                for (int j = 0; j < (*vert2segment)[i].size(); j++) {
                    int elnr = (*vert2segment)[i][j];
                    const Segment& el = mesh.LineSegment(elnr);

                    INDEX_2 edge(el[0], el[1]);

                    int edgedir = (edge.I1() > edge.I2());
                    if (edgedir) std::swap(edge.I1(), edge.I2());

                    if (edge.I1() != i) continue;

                    segedges.Elem(elnr) = edgenr[edge.I2()] - 1;
                }
            }
        }

        // generate faces
        if (buildfaces) {
            surffaces.resize(nse);

            int oldnfa = face2vert.size();

            cnt = 0;
            for (int i = 0; i < face2vert.size(); i++) {
                cnt[face2vert[i][0]]++;
            }
            TABLE<int, PointIndex::BASE> vert2oldface(cnt);
            for (int i = 0; i < face2vert.size(); i++)
                vert2oldface.AddSave(face2vert[i][0], i);

            int max_face_on_vertex = 0;
            for (int i = PointIndex::BASE; i < nv + PointIndex::BASE; i++) {
                int onv = vert2oldface[i].size() + (*vert2surfelement)[i].size();
                max_face_on_vertex = std::max(onv, max_face_on_vertex);
            }

            nfa = oldnfa;
            for (int v = PointIndex::BASE; v < nv + PointIndex::BASE; v++) {
                int first_fa = nfa;

                INDEX_3_CLOSED_HASHTABLE<int> vert2face(2 * max_face_on_vertex + 10);

                for (int j = 0; j < vert2oldface[v].size(); j++) {
                    int fnr = vert2oldface[v][j];
                    INDEX_3 face(face2vert[fnr].I1(),
                                 face2vert[fnr].I2(),
                                 face2vert[fnr].I3());
                    vert2face.Set(face, fnr + 1);
                }

                for (int pass = 1; pass <= 2; pass++) {

                    if (pass == 2) {
                        for (int j = first_fa; j < face2vert.size(); j++) {
                            if (face2vert[j][0] == v) {
                                INDEX_3 face(face2vert[j].I1(),
                                             face2vert[j].I2(),
                                             face2vert[j].I3());
                                vert2face.Set(face, j + 1);
                            } else {
                                break;
                            }
                        }
                    }

                    for (int j = 0; j < (*vert2surfelement)[v].size(); j++) {
                        int elnr = (*vert2surfelement)[v][j];
                        const Element2d& el = mesh.SurfaceElement(elnr);

                        const ELEMENT_FACE* elfaces = GetFaces(el.GetType());

                        if (elfaces[0][3] == 0) { // triangle

                            int facenum;
                            int facedir;

                            INDEX_3 face(el.PNum(elfaces[0][0]),
                                         el.PNum(elfaces[0][1]),
                                         el.PNum(elfaces[0][2]));

                            // std::cout << "face = " << face <<std::endl;

                            facedir = 0;
                            if (face.I1() > face.I2()) {
                                std::swap(face.I1(), face.I2());
                                facedir += 1;
                            }
                            if (face.I2() > face.I3()) {
                                std::swap(face.I2(), face.I3());
                                facedir += 2;
                            }
                            if (face.I1() > face.I2()) {
                                std::swap(face.I1(), face.I2());
                                facedir += 4;
                            }

                            if (face.I1() != v) continue;

                            if (vert2face.Used(face)) {
                                facenum = vert2face.Get(face);
                            } else {
                                nfa++;
                                vert2face.Set(face, nfa);
                                facenum = nfa;

                                INDEX_4 hface(face.I1(), face.I2(), face.I3(), 0);
                                face2vert.push_back(hface);
                            }

                            surffaces.Elem(elnr) = facenum - 1;
                        } else {
                            // quad
                            int facenum;
                            int facedir;

                            INDEX_4Q face4(el.PNum(elfaces[0][0]),
                                           el.PNum(elfaces[0][1]),
                                           el.PNum(elfaces[0][2]),
                                           el.PNum(elfaces[0][3]));

                            facedir = 0;
                            if (std::min(face4.I1(), face4.I2()) > std::min(face4.I4(), face4.I3())) {
                                // z - orientation
                                facedir += 1;
                                std::swap(face4.I1(), face4.I4());
                                std::swap(face4.I2(), face4.I3());
                            }
                            if (std::min(face4.I1(), face4.I4()) > std::min(face4.I2(), face4.I3())) {
                                // x - orientation
                                facedir += 2;
                                std::swap(face4.I1(), face4.I2());
                                std::swap(face4.I3(), face4.I4());
                            }
                            if (face4.I2() > face4.I4()) {
                                facedir += 4;
                                std::swap(face4.I2(), face4.I4());
                            }

                            INDEX_3 face(face4.I1(), face4.I2(), face4.I3());
                            if (face.I1() != v) continue;

                            if (vert2face.Used(face))
                                facenum = vert2face.Get(face);
                            else {
                                nfa++;
                                vert2face.Set(face, nfa);
                                facenum = nfa;

                                INDEX_4 hface(face4.I1(), face4.I2(), face4.I3(), face4.I3());
                                face2vert.push_back(hface);
                            }

                            surffaces.Elem(elnr) = facenum - 1;
                        }
                    }

                    // sort faces
                    if (pass == 1) {
                        for (int i = 0; i < nfa - first_fa; i++) {
                            for (int j = first_fa + 1; j < nfa - i; j++) {
                                if (face2vert[j] < face2vert[j - 1]) {
                                    std::swap(face2vert[j - 1], face2vert[j]);
                                }
                            }
                        }
                    }
                }
            }


            face2vert.reserve(nfa);


            face2surfel.resize(nfa);
            face2surfel = 0;
            for (int i = 1; i <= nse; i++)
                face2surfel.Elem(GetSurfaceElementFace(i)) = i;

            surf2volelement.resize(nse);
            for (int i = 1; i <= nse; i++) {
                surf2volelement.Elem(i)[0] = 0;
                surf2volelement.Elem(i)[1] = 0;
            }
            face2vert.reserve(face2vert.size());

            // face table complete
            Array<short int> face_els(nfa), face_surfels(nfa);
            face_els = 0;
            face_surfels = 0;
            Array<int> hfaces;
            for (int i = 1; i <= nse; i++)
                face_surfels[GetSurfaceElementFace(i) - 1]++;
        }

        timestamp = NextTimeStamp();
    }

    const Point3d* MeshTopology::GetVertices(ELEMENT_TYPE et)
    {
        static Point3d segm_points[] = {Point3d(1, 0, 0),
                                        Point3d(0, 0, 0)};

        static Point3d trig_points[] = {Point3d(1, 0, 0),
                                        Point3d(0, 1, 0),
                                        Point3d(0, 0, 0)};

        static Point3d quad_points[] = {Point3d(0, 0, 0),
                                        Point3d(1, 0, 0),
                                        Point3d(1, 1, 0),
                                        Point3d(0, 1, 0)};

        static Point3d tet_points[] = {Point3d(1, 0, 0),
                                       Point3d(0, 1, 0),
                                       Point3d(0, 0, 1),
                                       Point3d(0, 0, 0)};

        static Point3d pyramid_points[] = {
                Point3d(0, 0, 0),
                Point3d(1, 0, 0),
                Point3d(1, 1, 0),
                Point3d(0, 1, 0),
                Point3d(0, 0, 1 - 1e-7),
        };

        static Point3d prism_points[] = {
                Point3d(1, 0, 0),
                Point3d(0, 1, 0),
                Point3d(0, 0, 0),
                Point3d(1, 0, 1),
                Point3d(0, 1, 1),
                Point3d(0, 0, 1)
        };


        static Point3d hex_points[] = {Point3d(0, 0, 0),
                                       Point3d(1, 0, 0),
                                       Point3d(1, 1, 0),
                                       Point3d(0, 1, 0),
                                       Point3d(0, 0, 1),
                                       Point3d(1, 0, 1),
                                       Point3d(1, 1, 1),
                                       Point3d(0, 1, 1)};


        switch (et) {
            case SEGMENT:
            case SEGMENT3:
                return segm_points;

            case TRIG:
            case TRIG6:
                return trig_points;

            case QUAD:
            case QUAD6:
            case QUAD8:
                return quad_points;

            case TET:
            case TET10:
                return tet_points;

            case PYRAMID:
                return pyramid_points;

            case PRISM:
            case PRISM12:
                return prism_points;

            case HEX:
                return hex_points;
            default:
                std::cerr << "Ng_ME_GetVertices, illegal element type " << et << std::endl;
        }
        return 0;
    }

    int MeshTopology::GetSurfaceElementFace(int elnr) const
    {
        return surffaces.Get(elnr) + 1;
    }

    void MeshTopology::GetFaceVertices(int fnr, Array<int>& vertices) const
    {
        vertices.resize(4);
        for (int i = 0; i < 4; i++)
            vertices[i] = face2vert.Get(fnr)[i];
        if (vertices[3] == 0)
            vertices.resize(3);
    }
}
