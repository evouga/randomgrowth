#ifndef MESH_H
#define MESH_H

#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include <Eigen/Core>
#include <Eigen/Sparse>

struct MyTraits : public OpenMesh::DefaultTraits
{
    EdgeTraits
    {
    private:
        double restlen_;
    public:
        EdgeT() : restlen_(0) {}
        double restlen() const {return restlen_;}
        void setRestlen(double l) {restlen_=l;}
    };
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> OMMesh;

class Mesh
{
public:
    Mesh();

    void elasticEnergy(const Eigen::VectorXd &q, const Eigen::VectorXd &g,
                       double &energy,
                       Eigen::VectorXd &gradq,
                       Eigen::VectorXd &gradg,
                       Eigen::SparseMatrix<double> &hessq,
                       Eigen::SparseMatrix<double> &hessg) const;

    void dofsFromGeometry(Eigen::VectorXd &q, Eigen::VectorXd &g) const;
    void dofsToGeometry(const Eigen::VectorXd &q, const Eigen::VectorXd &g);
    void setIntrinsicLengthsToCurrentLengths();
    int numdofs() const;
    int numedges() const;
    double getYoungsModulus() const;
    double getPoissonRatio() const;
    double getThickness() const;
    void setYoungsModulus(double val);
    void setPoissonRatio(double val);
    void setThickness(double val);

    virtual void render(bool showWireframe, bool smoothShade);

    Eigen::Vector3d centroid();
    double radius();
    bool exportOBJ(const char *filename);
    bool importOBJ(const char *filename);

private:
    void edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2);

    OMMesh *mesh_;
    double YoungsModulus_;
    double PoissonRatio_;
    double h_;
};

#endif // MESH_H
