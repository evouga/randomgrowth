#ifndef MESH_H
#define MESH_H

#include "omtypes.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "elasticenergy.h"
#include <QMutex>

class Controller;

struct ProblemParameters : public ElasticParameters
{
    // simulation
    double rho;

    double eulerTimestep;
    double dampingCoeff;
    int numEulerIters;

    // rendering
    bool showWireframe;
    bool smoothShade;

    // problem
    double growthAmount;
    double maxEdgeStrain;
    double baseGrowthProbability;

    std::string outputDir;

    virtual void dumpParameters(std::ostream &os);
};

class Mesh
{
public:
    Mesh();

    enum RelaxationType {RelaxMetric, RelaxEmbedding, FitMetric};

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

    bool calculateHarmonicModes(const char *filename);

private:
    void dofsFromGeometry(Eigen::VectorXd &q, Eigen::VectorXd &g) const;
    void dofsToGeometry(const Eigen::VectorXd &q, const Eigen::VectorXd &g);
    void setIntrinsicLengthsToCurrentLengths();
    void edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2);
    double triangleInequalityLineSearch(const Eigen::VectorXd &g, const Eigen::VectorXd &dg) const;
    double triangleInequalityLineSearch(double g0, double g1, double g2, double dg0, double dg1, double dg2) const;
    double infinityNorm(const Eigen::VectorXd &v) const;
    void buildMassMatrix(const Eigen::VectorXd &g, Eigen::SparseMatrix<double> &M) const;
    void buildGeometricMassMatrix(const Eigen::VectorXd &g, Eigen::SparseMatrix<double> &M) const;
    void buildInvMassMatrix(const Eigen::VectorXd &g, Eigen::SparseMatrix<double> &M) const;
    double barycentricDualArea(const Eigen::VectorXd &g, int vidx) const;
    double deformedBarycentricDualArea(const Eigen::VectorXd &q, int vidx) const;
    double faceArea(const Eigen::VectorXd &q, int fidx) const;

    double restFaceArea(const Eigen::VectorXd &g, int fidx) const;
    void undeformedDofsFromGeometry(Eigen::VectorXd &undefq, Eigen::VectorXd &undefg) const;
    double vertexAreaRatio(const Eigen::VectorXd &undefq, const Eigen::VectorXd &g, int vidx);

    double intrinsicCotanWeight(int edgeid, const Eigen::VectorXd &g) const;
    double cotanWeight(int edgeid, const Eigen::VectorXd &q) const;
    void buildIntrinsicDirichletLaplacian(const Eigen::VectorXd &g, Eigen::SparseMatrix<double> &L) const;
    void buildExtrinsicDirichletLaplacian(const Eigen::VectorXd &q, Eigen::SparseMatrix<double> &L) const;
    void gaussianCurvature(const Eigen::VectorXd &q, Eigen::VectorXd &K) const;
    void meanCurvature(const Eigen::VectorXd &q, Eigen::VectorXd &Hdensity) const;
    void vertexAreas(const Eigen::VectorXd &q, Eigen::VectorXd &vareas) const;
    Eigen::Vector3d averageNormal(const Eigen::VectorXd &q, int vidx) const;
    Eigen::Vector3d faceNormal(const Eigen::VectorXd &q, int fidx) const;

    void elasticEnergy(const Eigen::VectorXd &q, const Eigen::VectorXd &g,
                       double &energyB,
                       double &energyS,
                       Eigen::VectorXd &gradq,
                       Eigen::SparseMatrix<double> &hessq, Eigen::SparseMatrix<double> &gradggradq,
                       int derivativesRequested) const;

    double vertexStrainEnergy(const Eigen::VectorXd &q, const Eigen::VectorXd &g, int vidx) const;
    double faceStrainEnergy(const Eigen::VectorXd &q, const Eigen::VectorXd &g, int fidx) const;

    void probabilisticallyGrowDisks(const Eigen::VectorXd &q,
                                    const Eigen::VectorXd &g,
                                    double baseprob, double increment, double maxstrain,
                                    Eigen::VectorXd &newg);

    void dumpFrame();


    Eigen::Vector3d colormap(double val) const;
    Eigen::Vector3d colormap(double val, double max) const;
    Eigen::Vector3d HSLtoRGB(const Eigen::Vector3d &hsl) const;

    double randomRange(double min, double max) const;

    OMMesh *mesh_;
    OMMesh *undeformedMesh_;
    int frameno_;
    ProblemParameters params_;

    // The rendering thread reads the mesh and its edge data. Any function must lock this before writing to
    // to the mesh. (The rendering thread does not write to the mesh so reads from the worker thread do not
    // need to lock.)
    QMutex meshLock_;
};

#endif // MESH_H
