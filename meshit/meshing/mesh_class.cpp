#include <cassert>
#include <stdexcept>

#include "mesh_class.hpp"
#include "../geom2d/geometry2d.hpp"
#include "../gprim/geomtest3d.hpp"
#include "../interface/writeuser.hpp"
#include "mesh_generator.hpp"

namespace meshit {

Mesh::Mesh()
    : surfarea(*this)
{
    loc_h_func = nullptr;
    hglob_ = 1e10;
    hmin_ = 0;
}

Mesh::~Mesh()
{
    delete loc_h_func;

    for (size_t i = 0; i < materials.size(); i++) {
        delete[] materials[i];
    }
}

Mesh& Mesh::operator=(const Mesh& mesh2)
{
    points = mesh2.points;
    segments = mesh2.segments;
    elements = mesh2.elements;
    faces = mesh2.faces;

    return *this;
}

void Mesh::BuildFromSplineGeometry(SplineGeometry& geometry, MeshingParameters& mp)
{
    MESHIT_LOG_DEBUG("Generate Mesh from spline geometry");

    // Take grading from MeshingParameters in priority, or from geometry if not user-defined.
    double grading = mp.grading;
    if (grading < 0) {
        grading = geometry.GetGrading();
    }
    if (grading < 0.01) {
        MESHIT_LOG_WARNING("grading is too small: " << grading << ". We reset it to 0.01.");
        grading = 0.01;
    }
    if (grading > 0.9) {
        MESHIT_LOG_WARNING("grading is too large: " << grading << ". We reset it to 0.9.");
        grading = 0.9;
    }
    geometry.SetGrading(grading);

    // Define local meshing sizes
    Box2d bbox = geometry.GetBoundingBox();
    mp.maxh = std::min(mp.maxh, bbox.Diam());

    SetLocalH(bbox, grading);
    hmin_ = mp.minh;
    hglob_ = mp.maxh;

    geometry.PartitionBoundary(*this, mp);

    // Build mesh faces
    size_t nb_faces = 0;
    for (size_t si = 0; si < segments.size(); si++) {
        if ((segments[si].face_left == CONST<PointIndex>::undefined) ||
            (segments[si].face_right == CONST<PointIndex>::undefined)) {
            throw std::runtime_error("Segment " + std::to_string(si) + " is adjacent to an undefined domain !");
        }
        nb_faces = std::max(nb_faces, static_cast<size_t>(segments[si].face_left));
        nb_faces = std::max(nb_faces, static_cast<size_t>(segments[si].face_right));
    }

    faces.resize(nb_faces);
    for (size_t face_index = 0; face_index < nb_faces; face_index++) {
        faces[face_index] = FaceDescriptor(face_index + 1);
    }

    CalcLocalH();

    size_t np_bdy = points.size();  // boundary points
    double h = mp.maxh;

    MeshGenerator mesher(*this, bbox);

    for (DomainIndex dom_nr = 1; dom_nr <= nb_faces; dom_nr++) {
        if (geometry.GetDomainMaxh(dom_nr) > 0) {
            h = geometry.GetDomainMaxh(dom_nr);
        }
        MESHIT_LOG_DEBUG("Meshing domain " << dom_nr << " / " << nb_faces);

        size_t oldnf = elements.size();

        if (dom_nr > 1) {
            mesher.Reset();
        }

        constexpr PointIndex undef_pid = CONST<PointIndex>::undefined;
        std::vector<PointIndex> compress(np_bdy, undef_pid);
        PointIndex cnt = 0;
        for (size_t pi = 0; pi < np_bdy; pi++) {
            mesher.AddPoint(Point2d(points[pi]), pi);
            compress[pi] = cnt++;
        }
        for (size_t si = 0; si < segments.size(); si++) {
            if (segments[si].face_left == dom_nr) {
                mesher.AddBoundaryElement(compress[segments[si][0]], compress[segments[si][1]]);
            }
            if (segments[si].face_right == dom_nr) {
                mesher.AddBoundaryElement(compress[segments[si][1]], compress[segments[si][0]]);
            }
        }

        mesher.GenerateMesh(mp, h, dom_nr);

        for (size_t sei = oldnf; sei < elements.size(); sei++) {
            elements[sei].SetFaceID(dom_nr);
        }

        // astrid
        char* material;
        geometry.GetMaterial(dom_nr, material);
        if (material) {
            SetMaterial(dom_nr, material);
        }
    }

    int hsteps = mp.optsteps2d;

    mp.optimize2d = "smm";
    mp.optsteps2d = hsteps / 2;
    Optimize2d(mp);

    mp.optimize2d = "Smm";
    mp.optsteps2d = (hsteps + 1) / 2;
    Optimize2d(mp);

    mp.optsteps2d = hsteps;

    Compress();
}

size_t Mesh::AddPoint(const Point2d& p, point_type_t type)
{
    points.push_back(MeshPoint(p, type));
    return points.size() - 1;
}

void Mesh::AddSegment(const Segment& s)
{
    PointIndex maxn = std::max(s[0], s[1]);

    if (maxn < static_cast<PointIndex>(points.size())) {
        if (points[s[0]].Type() > EDGE_POINT) points[s[0]].SetType(EDGE_POINT);
        if (points[s[1]].Type() > EDGE_POINT) points[s[1]].SetType(EDGE_POINT);
    }
    segments.push_back(s);
}

void Mesh::AddElement(const Element2d& el)
{
    size_t si = elements.size();
    elements.push_back(el);

    if (el.face_id_ > faces.size()) {
        MESHIT_LOG_ERROR("has no faces: fd.size = " << faces.size() << ", ind = " << el.face_id_);
    }

    Element2d& bref = elements.back();
    FaceDescriptor& faceref = faces[bref.face_id_ - 1];
    bref.next_ = faceref.first_element;
    faceref.first_element = si;

    if (surfarea.Valid()) {
        surfarea.Add(el);
    }
}

void Mesh::Export(const std::string& filetype, const std::string& filename) const
{
    WriteUserFormat(filetype, const_cast<const Mesh&>(*this), filename);
}

void Mesh::Export(const std::string& filetype, std::ostream& os) const
{
    WriteUserFormat(filetype, const_cast<const Mesh&>(*this), os);
}

void Mesh::Export(const std::string& filename) const
{
    WriteGmsh2Format(const_cast<const Mesh&>(*this), filename);
}

void Mesh::Export(std::ostream& os) const
{
    WriteGmsh2Format(const_cast<const Mesh&>(*this), os);
}

void Mesh::Save(const std::string& filename) const
{
    std::ofstream outfile(filename.c_str());
    Save(outfile);
}

void Mesh::Save(std::ostream& outfile) const
{
    outfile.setf(std::ios::fixed, std::ios::floatfield);
    outfile.setf(std::ios::showpoint);
    outfile.precision(15);

    outfile << "mesh2d" << std::endl;
    outfile << "# face_id      np      p1      p2      p3" << std::endl;
    outfile << "surface_elements" << std::endl;
    outfile << elements.size() << std::endl;

    for (size_t sei = 0; sei < elements.size(); sei++) {
        DomainIndex el_index = elements[sei].FaceID();
        outfile << std::setw(9) << faces[el_index - 1].face_id();

        const Element2d& sel = elements[sei];
        outfile << std::setw(8) << 3;
        outfile << std::setw(8) << sel[0];
        outfile << std::setw(8) << sel[1];
        outfile << std::setw(8) << sel[2] << std::endl;
    }

    outfile << "\n\n";
    outfile << "# surfid      p1      p2   dom_l   dom_r\n";
    outfile << "edge_segments" << std::endl << segments.size() << std::endl;

    for (size_t i = 0; i < segments.size(); i++) {
        const Segment& seg = segments[i];
        outfile << std::setw(8) << seg.edge_id;
        outfile << std::setw(8) << seg[0];
        outfile << std::setw(8) << seg[1];
        outfile << std::setw(8) << seg.face_left;
        outfile << std::setw(8) << seg.face_right;
        outfile << std::endl;
    }

    outfile << "\n\n";
    outfile << "#                    X                     Y\n";
    outfile << "points" << std::endl << points.size() << std::endl;
    for (size_t pi = 0; pi < points.size(); pi++) {
        outfile << std::setw(22) << points[pi].X();
        outfile << std::setw(22) << points[pi].Y() << std::endl;
    }

    int cntmat = 0;
    for (size_t i = 0; i < materials.size(); i++) {
        if (materials[i] && strlen(materials[i])) {
            cntmat++;
        }
    }
    if (cntmat) {
        outfile << "materials" << std::endl << cntmat << std::endl;
        for (size_t i = 0; i < materials.size(); i++) {
            if (materials[i] && strlen(materials[i])) outfile << i + 1 << " " << materials[i] << std::endl;
        }
    }
}

void Mesh::Load(const std::string& filename)
{
    std::ifstream infile(filename.c_str());
    if (!infile.good()) throw std::runtime_error("mesh file not found");

    Load(infile);
}

void Mesh::Load(std::istream& infile)
{
    char str[100];
    int n;

    faces.resize(0);

    bool endmesh = false;

    while (infile.good() && !endmesh) {
        infile >> str;
        if (strcmp(str, "surface_elements") == 0) {
            infile >> n;
            MESHIT_LOG_DEBUG(n << " surface elements");

            for (int i = 1; i <= n; i++) {
                DomainIndex surf_id;
                int nep, faceind = 0;

                infile >> surf_id;

                for (size_t j = 0; j < faces.size(); j++) {
                    if (faces[j].face_id() == static_cast<size_t>(surf_id)) {
                        faceind = j + 1;
                    }
                }
                if (!faceind) {
                    faces.push_back(FaceDescriptor(surf_id));
                }

                infile >> nep;

                if (nep != 3) {
                    MESHIT_LOG_FATAL("Mesh::Load: undefined element type. nep = " << nep << ". Aborting");
                    exit(1);
                }

                Element2d tri;
                tri.SetFaceID(surf_id);

                for (size_t j = 0; j < 3; j++) {
                    infile >> tri.PointID(j);
                }
                AddElement(tri);
            }
        }
        if (strcmp(str, "edge_segments") == 0) {
            infile >> n;
            for (int i = 1; i <= n; i++) {
                Segment seg;
                infile >> seg.edge_id >> seg[0] >> seg[1];
                infile >> seg.face_left >> seg.face_right;
                AddSegment(seg);
            }
        }
        if (strcmp(str, "points") == 0) {
            infile >> n;
            MESHIT_LOG_DEBUG(n << " points");
            for (int i = 1; i <= n; i++) {
                Point2d p;
                infile >> p.X() >> p.Y();
                AddPoint(p);
            }
        }
        if (strcmp(str, "materials") == 0) {
            infile >> n;
            MESHIT_LOG_DEBUG(n << " materials");
            for (int i = 1; i <= n; i++) {
                int nr;
                std::string mat;
                infile >> nr >> mat;
                SetMaterial(nr, mat.c_str());
            }
        }
        if (strcmp(str, "endmesh") == 0) {
            endmesh = true;
        }
        str[0] = '\0';
    }

    IndexBoundaryEdges();
}

void Mesh::IndexBoundaryEdges()
{
    segment_ht.reserve(3 * segments.size());

    for (size_t i = 0; i < segments.size(); i++) {
        const Segment& seg = segments[i];
        MeshPoint& mp1 = points[seg[0]];
        MeshPoint& mp2 = points[seg[1]];
        if (mp1.Type() == INNER_POINT) mp1.SetType(EDGE_POINT);
        if (mp2.Type() == INNER_POINT) mp2.SetType(EDGE_POINT);
    }

    for (size_t i = 0; i < segments.size(); i++) {
        const Segment& seg = segments[i];
        IndexPair i2(seg[0], seg[1]);
        segment_ht[i2.Sort()] = i;
    }
}

void Mesh::SetLocalH(const Box2d& bbox, double grading)
{
    delete loc_h_func;
    loc_h_func = new LocalH();
    loc_h_func->Init(bbox, grading);
}

void Mesh::RestrictLocalH(const Point2d& p, double hloc)
{
    assert(loc_h_func);
    loc_h_func->SetH(p, std::max(hloc, hmin_));
}

void Mesh::RestrictLocalHLine(const Point2d& p1, const Point2d& p2, double hloc)
{
    if (hloc < hmin_) {
        hloc = hmin_;
    }
    Vec2d v(p1, p2);
    int steps = static_cast<int>(v.Length() / hloc) + 1;
    v /= static_cast<double>(steps);

    for (int i = 0; i <= steps; i++) {
        RestrictLocalH(p1 + static_cast<double>(i) * v, hloc);
    }
}

double Mesh::GetH(const Point2d& p) const
{
    assert(loc_h_func);
    return std::min(hglob_, loc_h_func->GetH(p));
}

double Mesh::GetMinH(const Point2d& pmin, const Point2d& pmax)
{
    assert(loc_h_func);
    return std::min(hglob_, loc_h_func->GetMinH(pmin, pmax));
}

double Mesh::AverageH(DomainIndex surf_id) const
{
    double maxh = 0, minh = 1e10;
    double hsum = 0;
    size_t n = 0;
    for (size_t i = 0; i < elements.size(); i++) {
        const Element2d& el = Element(i);
        if (surf_id == 0 || el.FaceID() == surf_id) {
            for (size_t j = 0; j < 3; j++) {
                double hi = Dist(points[el.PointID(j)], points[el.PointID((j + 1) % 3)]);
                hsum += hi;
                if (hi > maxh) maxh = hi;
                if (hi < minh) minh = hi;
                n++;
            }
        }
    }
    MESHIT_LOG_DEBUG("minh = " << minh << " avh = " << (hsum / n) << " maxh = " << maxh);
    return (hsum / n);
}

void Mesh::CalcLocalH()
{
    assert(loc_h_func);

    MESHIT_LOG_DEBUG("CalcLocalH: " << points.size() << " points, "
                     << elements.size() << " surface elements.");

    for (size_t i = 0; i < elements.size(); i++) {
        const Element2d& el = elements[i];
        double h_el = 0.0;
        for (size_t j = 0; j < 3; j++) {
            const MeshPoint& p1 = points[el.PointID(j)];
            const MeshPoint& p2 = points[el.PointID((j + 1) % 3)];
            double h_edge = Dist(p1, p2);
            if (h_edge > h_el) h_el = h_edge;
        }
        if (h_el > 0.0) {
            const MeshPoint& p1 = points[el.PointID(0)];
            const MeshPoint& p2 = points[el.PointID(1)];
            const MeshPoint& p3 = points[el.PointID(2)];
            loc_h_func->SetH(Center(p1, p2, p3), h_el);
        }
    }

    for (size_t i = 0; i < segments.size(); i++) {
        const Segment& seg = segments[i];
        const MeshPoint& p1 = points[seg[0]];
        const MeshPoint& p2 = points[seg[1]];
        loc_h_func->SetH(Center(p1, p2), Dist(p1, p2));
    }
}

void Mesh::GetBox(Point2d& pmin, Point2d& pmax) const
{
    if (points.size() == 0) {
        pmin = pmax = Point2d(0, 0);
        return;
    }

    pmin = Point2d(1e10, 1e10);
    pmax = Point2d(-1e10, -1e10);
    for (size_t pi = 0; pi < points.size(); pi++) {
        pmin.SetToMin(Point2d(points[pi]));
        pmax.SetToMax(Point2d(points[pi]));
    }
    if (pmin.X() > 0.5e10) {
        pmin = pmax = Point2d(0, 0);
    }
}

void Mesh::Compress()
{
    std::vector<PointIndex> op2np(points.size());
    std::vector<MeshPoint> hpoints;
    std::vector<bool> pused(points.size(), false);

    for (size_t i = 0; i < elements.size(); i++) {
        if (elements[i].IsDeleted()) {
            elements.erase(elements.begin() + i);
            i--;
        }
    }
    for (size_t i = 0; i < segments.size(); i++) {
        if (segments[i][0] == CONST<PointIndex>::undefined) {
            segments.erase(segments.begin() + i);
            i--;
        }
    }
    for (size_t i = 0; i < elements.size(); i++) {
        const Element2d& el = elements[i];
        for (size_t j = 0; j < 3; j++) {
            pused[el[j]] = true;
        }
    }

    for (size_t i = 0; i < segments.size(); i++) {
        const Segment& seg = segments[i];
        pused[seg[0]] = true;
        pused[seg[1]] = true;
    }

    PointIndex npi = 0;

    for (size_t pi = 0; pi < points.size(); pi++) {
        if (pused[pi]) {
            op2np[pi] = npi++;
            hpoints.push_back(points[pi]);
        } else {
            op2np[pi] = CONST<PointIndex>::undefined;
        }
    }

    points.resize(0);
    for (size_t i = 0; i < hpoints.size(); i++) {
        points.push_back(hpoints[i]);
    }

    for (size_t i = 0; i < elements.size(); i++) {
        Element2d& el = Element(i);
        for (size_t j = 0; j < 3; j++) {
            el[j] = op2np[el[j]];
        }
    }

    for (size_t i = 0; i < segments.size(); i++) {
        Segment& seg = segments[i];
        seg[0] = op2np[seg[0]];
        seg[1] = op2np[seg[1]];
    }
    for (size_t i = 0; i < faces.size(); i++) {
        faces[i].first_element = CONST<ElementIndex>::undefined;
    }

    for (size_t i = 0; i < elements.size(); i++) {
        size_t idx = elements.size() - i - 1;
        DomainIndex ind = elements[idx].FaceID() - 1;
        elements[idx].next_ = faces[ind].first_element;
        faces[ind].first_element = idx;
    }

    IndexBoundaryEdges();
}

int Mesh::CheckOverlappingBoundary()
{
    Point2d pmin, pmax;
    GetBox(pmin, pmax);
    Box3dTree setree(pmin, pmax);
    std::vector<GenericIndex> inters;

    bool overlap = false;

    for (size_t i = 0; i < elements.size(); i++) {
        const Element2d& tri = Element(i);

        Point2d tpmin{points[tri[0]]};
        Point2d tpmax = tpmin;

        for (size_t k = 1; k < 3; k++) {
            tpmin.SetToMin(Point2d(points[tri[k]]));
            tpmax.SetToMax(Point2d(points[tri[k]]));
        }
        Vec2d diag(tpmin, tpmax);

        tpmax = tpmax + 0.1 * diag;
        tpmin = tpmin - 0.1 * diag;

        setree.Insert(tpmin, tpmax, i + 1);
    }

    for (size_t i = 0; i < elements.size(); i++) {
        const Element2d& tri1 = Element(i);

        Point2d tpmin{points[tri1[0]]};
        Point2d tpmax = tpmin;

        for (size_t k = 1; k < 3; k++) {
            tpmin.SetToMin(Point2d(points[tri1[k]]));
            tpmax.SetToMax(Point2d(points[tri1[k]]));
        }

        setree.GetIntersecting(tpmin, tpmax, inters);

        for (size_t j = 0; j < inters.size(); j++) {
            const Element2d& tri2 = Element(inters[j] - 1);

            const Point2d* trip1[3], * trip2[3];
            for (size_t k = 0; k < 3; k++) {
                trip1[k] = &points[tri1.PointID(k)];
                trip2[k] = &points[tri2.PointID(k)];
            }

            if (IntersectTriangleTriangle(&trip1[0], &trip2[0])) {
                overlap = 1;
                MESHIT_LOG_WARNING("Intersecting elements " << i + 1 << " and " << inters[j]);
                MESHIT_LOG_DEBUG(" el1 = " << tri1);
                MESHIT_LOG_DEBUG(" el2 = " << tri2);

                for (size_t k = 0; k < 3; k++) MESHIT_LOG_DEBUG_CONT(tri1.PointID(k) << "  ");
                MESHIT_LOG_DEBUG("");
                for (size_t k = 0; k < 3; k++) MESHIT_LOG_DEBUG_CONT(tri2.PointID(k) << "  ");
                MESHIT_LOG_DEBUG("");

                for (size_t k = 0; k < 3; k++) MESHIT_LOG_DEBUG_CONT(*trip1[k] << "   ");
                MESHIT_LOG_DEBUG("");
                for (size_t k = 0; k < 3; k++) MESHIT_LOG_DEBUG_CONT(*trip2[k] << "   ");
                MESHIT_LOG_DEBUG("");

                MESHIT_LOG_DEBUG("Face1 = " << faces[tri1.FaceID() - 1]);
                MESHIT_LOG_DEBUG("Face2 = " << faces[tri2.FaceID() - 1]);
            }
        }
    }
    return overlap;
}

void Mesh::RebuildSurfaceElementLists()
{
    for (size_t i = 0; i < faces.size(); i++) {
        faces[i].first_element = CONST<ElementIndex>::undefined;
    }
    for (size_t i = 0; i < elements.size(); i++) {
        size_t idx = elements.size() - i - 1;
        DomainIndex ind = elements[idx].FaceID() - 1;
        elements[idx].next_ = faces[ind].first_element;
        faces[ind].first_element = idx;
    }
}

void Mesh::GetElementsOfFace(DomainIndex facenr, std::vector<ElementIndex>& sei) const
{
    sei.clear();

    ElementIndex si = faces[facenr - 1].first_element;
    while (si != CONST<ElementIndex>::undefined) {
        const Element2d& se = Element(si);
        if (se.FaceID() == facenr && se.PointID(0) != CONST<PointIndex>::undefined && !se.IsDeleted()) {
            sei.push_back(si);
        }
        si = se.next_;
    }
}

size_t Mesh::ComputeNVertices()
{
    PointIndex nb_vertices = 0;
    for (size_t i = 0; i < elements.size(); i++) {
        const Element2d& el = Element(i);
        nb_vertices = std::max(nb_vertices, el.PointID(0));
        nb_vertices = std::max(nb_vertices, el.PointID(1));
        nb_vertices = std::max(nb_vertices, el.PointID(2));
    }
    nb_vertices += 1;
    return static_cast<size_t>(nb_vertices);
}

void Mesh::SetMaterial(size_t domnr, const char* mat)
{
    if (domnr > materials.size()) {
        materials.resize(domnr, nullptr);
    }
    materials[domnr - 1] = new char[strlen(mat) + 1];
    strcpy(materials[domnr - 1], mat);
}

void Mesh::PrintMemInfo(std::ostream& ost) const
{
    ost << "Mesh Mem:" << std::endl;

    ost << points.size() << " mesh points, of size " << sizeof(MeshPoint) << " = "
    << points.size() * sizeof(MeshPoint) << std::endl;

    ost << elements.size() << " elements, of size " << sizeof(Element2d) << " = "
    << elements.size() * sizeof(Element2d) << std::endl;
}

}  // namespace meshit
