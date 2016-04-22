#ifndef FILE_MESHING2_HPP
#define FILE_MESHING2_HPP

/**************************************************************************/
/* File:   meshing2.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

#include "adfront2.hpp"
#include "mesh_class.hpp"
#include "ruler2.hpp"

namespace meshit {

/*
    The basis class for 2D mesh generation.
    Has the method GenerateMesh
 */

class MeshGenerator
{
 public:
    explicit MeshGenerator(Mesh& amesh, const Box2d& aboundingbox);
    ~MeshGenerator();

    // Load rules, either from file, or compiled rules
    void LoadRules(const char* filename);
    void Reset();
    bool GenerateMesh(const MeshingParameters& mp, double gh, DomainIndex facenr);

    void AddPoint(const Point2d& p, PointIndex globind);
    void AddBoundaryElement(INDEX i1, INDEX i2);
    void SetMaxArea(double amaxarea);

 protected:
    void DefineTransformation(const Point2d& p1, const Point2d& p2);
    void TransformToPlain(const Point2d& locpoint, Point2d& plainpoint, double h);
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
    std::vector<netrule*> rules;  // rules for mesh generation

    // statistics
    std::vector<int> ruleused, canuse, foundmap;
    double max_area;

    Vec2d ex, ey;
    Point2d glob_p1;
};

}  // namespace meshit

#endif
