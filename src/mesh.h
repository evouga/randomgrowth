#ifndef MESH_H
#define MESH_H

#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "elasticenergy.h"
#include <QMutex>

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

struct ProblemParameters : public ElasticParameters
{
    // simulation
    int maxiters;
    int maxlinesearchiters;
    double tol;
    double rho;

    double eulerTimestep;
    double dampingCoeff;
    int numEulerIters;

    // rendering
    bool showWireframe;
    bool smoothShade;
};

class Mesh
{
public:
    Mesh();

    enum RelaxationType {RelaxMetric, RelaxEmbedding, FitMetric};

    bool relaxEnergy(Controller &cont, RelaxationType type);
    bool simulate(Controller &cont);

    int numdofs() const;
    int numedges() const;
    const ProblemParameters &getParameters() const;
    void setParameters(ProblemParameters params);

    // Rendering methods. These run concurrently and must all lock the meshLock before reading from the mesh.
    void render();
    Eigen::Vector3d centroid();
    double radius();
    // End rendering methods

    bool exportOBJ(const char *filename);
    bool importOBJ(const char *filename);

    friend class EmbeddingMinimizer;
    friend class MetricFit;
    friend class ImplicitEulerStep;

private:
    void dofsFromGeometry(Eigen::VectorXd &q, Eigen::VectorXd &g) const;
    void dofsToGeometry(const Eigen::VectorXd &q, const Eigen::VectorXd &g);
    void setIntrinsicLengthsToCurrentLengths();
    void edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2);
    double triangleInequalityLineSearch(const Eigen::VectorXd &g, const Eigen::VectorXd &dg) const;
    double triangleInequalityLineSearch(double g0, double g1, double g2, double dg0, double dg1, double dg2) const;
    double infinityNorm(const Eigen::VectorXd &v) const;
    void buildMassMatrix(const Eigen::VectorXd &q, Eigen::SparseMatrix<double> &M) const;
    double barycentricDualArea(const Eigen::VectorXd &q, int vidx) const;
    double faceArea(const Eigen::VectorXd &q, int fidx) const;

    double strainDensity(int edgeidx) const;
    double vertexStrainDensity(int vertidx) const;

    void elasticEnergy(const Eigen::VectorXd &q, const Eigen::VectorXd &g,
                       double &energyB,
                       double &energyS,
                       Eigen::VectorXd &gradq,
                       Eigen::SparseMatrix<double> &hessq, Eigen::SparseMatrix<double> &gradggradq,
                       int derivativesRequested) const;


    Eigen::Vector3d colormap(double val) const;
    Eigen::Vector3d colormap(double val, double max) const;
    Eigen::Vector3d HSLtoRGB(const Eigen::Vector3d &hsl) const;

    OMMesh *mesh_;
    ProblemParameters params_;

    // The rendering thread reads the mesh and its edge data. Any function must lock this before writing to
    // to the mesh. (The rendering thread does not write to the mesh so reads from the worker thread do not
    // need to lock.)
    QMutex meshLock_;
};

#endif // MESH_H
