#ifndef FILE_MESHING2_HPP
#define FILE_MESHING2_HPP

/**************************************************************************/
/* File:   meshing2.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

#include "adfront2.hpp"
#include "meshclass.hpp"
#include "ruler2.hpp"

namespace meshit {
/*
    The basis class for 2D mesh generation.
    Has the method GenerateMesh

    For surface mesh generation, or non-Euklidean meshing,
    derive from Meshing2, and replace transformation.
 */

class Meshing2
{
 public:
    explicit Meshing2(Mesh& amesh, const Box2d& aboundingbox);
    ~Meshing2();

    /// Load rules, either from file, or compiled rules
    void LoadRules(const char* filename);
    void Reset();
    bool GenerateMesh(const MeshingParameters& mp, double gh, int facenr);

    void AddPoint(const Point2d& p, PointIndex globind);
    void AddBoundaryElement(INDEX i1, INDEX i2);
    void SetMaxArea(double amaxarea);

 protected:
    void EndMesh();

    void DefineTransformation(const Point2d& p1, const Point2d& p2);
    void TransformToPlain(const Point2d& locpoint, Point2d& plainpoint, double h);
    /// return 0 .. ok
    /// return >0 .. cannot transform point to true surface
    void TransformFromPlain(Point2d& plainpoint, Point2d& locpoint, double h);

    /** Applies 2D rules.
     Tests all 2D rules */
    int ApplyRules(std::vector<Point2d>& lpoints,
                   std::vector<int>& legalpoints, size_t maxlegalpoint,
                   std::vector<INDEX_2>& llines, size_t maxlegelline,
                   std::vector<Element2d>& elements,
                   std::vector<uint32_t>& dellines,
                   int tolerance, const MeshingParameters& mp);

 protected:
    Mesh& mesh_;
    const Box2d& boundingbox_;
    AdFront2* adfront;  // the current advancing front

    double max_area;

    std::vector<netrule*> rules;  // rules for mesh generation

    // statistics
    std::vector<int> ruleused, canuse, foundmap;

    Vec2d ex, ey;
    Point2d glob_p1;
};

}  // namespace meshit

#endif
