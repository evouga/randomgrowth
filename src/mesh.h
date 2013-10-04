#ifndef MESH_H
#define MESH_H

#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include <Eigen/Core>
#include <Eigen/Sparse>

const double PI = 3.14159265359;

class Controller;

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

struct ProblemParameters
{
    // simulation
    double h;
    double PoissonRatio;
    double YoungsModulus;
    int maxiters;
    int maxlinesearchiters;
    int maxpoweriters;
    double powertol;
    double tol;

    // rendering
    bool showWireframe;
    bool smoothShade;
};

class Mesh
{
public:
    Mesh();

    void elasticEnergy(const Eigen::VectorXd &q, const Eigen::VectorXd &g,
                        double &energyB, double &energyS) const;

    void elasticEnergyG(const Eigen::VectorXd &q, const Eigen::VectorXd &g,
                        double &energyB,
                        double &energyS,
                        Eigen::VectorXd &gradg,
                        Eigen::SparseMatrix<double> &hessg) const;

    void elasticEnergyQ(const Eigen::VectorXd &q, const Eigen::VectorXd &g,
                        double &energyB,
                        double &energyS,
                        Eigen::VectorXd &gradq,
                        Eigen::SparseMatrix<double> &hessq) const;

    enum RelaxationType {RelaxMetric, RelaxEmbedding};

    bool relaxEnergy(Controller &cont, RelaxationType type);

    int numdofs() const;
    int numedges() const;
    const ProblemParameters &getParameters() const;
    void setParameters(ProblemParameters params);

    void render();

    Eigen::Vector3d centroid();
    double radius();
    bool exportOBJ(const char *filename);
    bool importOBJ(const char *filename);

private:
    enum EnergyDerivatives
    {
        NONE = 0,
        Q = 1,
        G = 2
    };

    void dofsFromGeometry(Eigen::VectorXd &q, Eigen::VectorXd &g) const;
    void dofsToGeometry(const Eigen::VectorXd &q, const Eigen::VectorXd &g);
    void setIntrinsicLengthsToCurrentLengths();
    void edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2);
    double triangleInequalityLineSearch(const Eigen::VectorXd &g, const Eigen::VectorXd &dg) const;
    double triangleInequalityLineSearch(double g0, double g1, double g2, double dg0, double dg1, double dg2) const;
    bool largestMagnitudeEigenvalue(const Eigen::SparseMatrix<double> &M, double &eigenvalue);
    bool smallestEigenvalue(const Eigen::SparseMatrix<double> &M, double &eigenvalue);
    double infinityNorm(const Eigen::VectorXd &v) const;

    double strainDensity(int edgeidx) const;
    double vertexStrainDensity(int vertidx) const;

    void elasticEnergy(const Eigen::VectorXd &q, const Eigen::VectorXd &g,
                       double &energyB,
                       double &energyS,
                       Eigen::VectorXd &gradq,
                       Eigen::VectorXd &gradg,
                       Eigen::SparseMatrix<double> &hessq,
                       Eigen::SparseMatrix<double> &hessg,
                       EnergyDerivatives derivs) const;


    Eigen::Vector3d colormap(double val) const;
    Eigen::Vector3d colormap(double val, double max) const;
    Eigen::Vector3d HSLtoRGB(const Eigen::Vector3d &hsl) const;

    OMMesh *mesh_;
    ProblemParameters params_;
};

#endif // MESH_H
