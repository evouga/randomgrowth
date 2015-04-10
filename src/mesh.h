#ifndef MESH_H
#define MESH_H

#include <Eigen/Core>
#include <vector>

class Mesh
{
public:
    Mesh() {}

    bool loadMesh(const char *filename);
    bool writeMesh(const char *filename);
    bool loadMesh(const Eigen::VectorXd &deformedPositions, const Eigen::Matrix3Xi &faces);

    Eigen::VectorXd::ConstFixedSegmentReturnType<3>::Type vertPos(int vert) const {return deformedPosition_.segment<3>(3*vert);}
    Eigen::Matrix3Xi::ConstColXpr faceVerts(int face) const {assert(face < faces_.cols()); return faces_.col(face);}
    Eigen::Matrix2Xi::ConstColXpr edgeFaces(int edge) const {return edgeFaces_.col(edge);}
    Eigen::Matrix2Xi::ConstColXpr edgeVerts(int edge) const {return edgeVerts_.col(edge);}
    Eigen::VectorXd::ConstFixedSegmentReturnType<4>::Type faceMetric(int face) const {return restMetrics_.segment<4>(4*face);}
    Eigen::VectorXd::FixedSegmentReturnType<3>::Type vertPos(int vert) {return deformedPosition_.segment<3>(3*vert);}
    Eigen::Matrix3Xi::ColXpr faceVerts(int face) {assert(face < faces_.cols()); return faces_.col(face);}
    Eigen::Matrix2Xi::ColXpr edgeFaces(int edge) {return edgeFaces_.col(edge);}
    Eigen::Matrix2Xi::ColXpr edgeVerts(int edge) {return edgeVerts_.col(edge);}
    Eigen::VectorXd::FixedSegmentReturnType<4>::Type faceMetric(int face) {return restMetrics_.segment<4>(4*face);}

    Eigen::Vector4i buildHinge(int edge) const;

    bool isBoundaryVert(int vert) const {return isBoundaryVert_[vert];}
    bool isBoundaryEdge(int edge) const;

    int numVertices() const {return deformedPosition_.size() / 3;}
    int numFaces() const {return faces_.cols();}    
    int numEdges() const {return edgeFaces_.cols();}

    void resetRestMetric();

protected:
    Eigen::VectorXd deformedPosition_;
    Eigen::Matrix3Xi faces_;
    Eigen::Matrix2Xi edgeFaces_;
    Eigen::Matrix2Xi edgeVerts_;
    Eigen::VectorXd restMetrics_;
    std::vector<bool> isBoundaryVert_;
};

#endif // MESH_H

