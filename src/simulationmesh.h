#ifndef SIMULATIONMESH_H
#define SIMULATIONMESH_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "midedge.h"
#include "mesh.h"
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
    bool constantPressure;
    double constantPressureVal;
    double airLeakCoeff;

    bool constantVelocity;
    double crushTime;
    double crushMass;

    double coneAngle;
    double coneHeight;

    std::string outputDir;

    virtual void dumpParameters(std::ostream &os);
};

class SimulationMesh : public Mesh
{
public:
    SimulationMesh();

    enum RelaxationType {RelaxMetric, RelaxEmbedding, FitMetric};

    bool crush(Controller &cont);

    const ProblemParameters &getParameters() const;
    void setParameters(ProblemParameters params);

    // Rendering methods. These run concurrently and must all lock the meshLock before reading from the mesh.
    void render();
    Eigen::Vector3d centroid();
    double radius();
    // End rendering methods

    bool exportOBJ(const char *filename);
    bool importOBJ(const char *filename);

    void addRandomNoise(double magnitude);
    void printHessianEigenvalues();
    void setConeHeights(double height);
    void setFlatCone(double height);

private:
    void buildMetricInvMassMatrix(Eigen::SparseMatrix<double> &M) const;
    void metricBarycentricDualAreas(Eigen::VectorXd &areas) const;
    void deformedBarycentricDualAreas(Eigen::VectorXd &areas) const;
    double deformedFaceArea(int fidx) const;

    void enforceConstraints(const Eigen::VectorXd &startq,
                            double planeHeight);

    double metricFaceArea(int fidx) const;
    void vertexNormals(Eigen::VectorXd &vertNormals);
    void gaussianCurvature(Eigen::VectorXd &K) const;
    Eigen::Vector3d faceNormal(int fidx) const;

    void elasticEnergy(const Eigen::VectorXd &q, const Eigen::VectorXd &g,
                       double &energyB,
                       double &energyS,
                       Eigen::VectorXd &gradq,
                       Eigen::SparseMatrix<double> &hessq, Eigen::SparseMatrix<double> &gradggradq,
                       int derivativesRequested) const;

    double vertexStrainEnergy(const Eigen::VectorXd &q, const Eigen::VectorXd &g, int vidx) const;
    double faceStrainEnergy(const Eigen::VectorXd &q, const Eigen::VectorXd &g, int fidx) const;

    void dumpFrame();
    void deleteBadFlatConeFaces();

    Eigen::Vector3d colormap(double val) const;
    Eigen::Vector3d colormap(double val, double max) const;
    Eigen::Vector3d HSLtoRGB(const Eigen::Vector3d &hsl) const;

    double randomRange(double min, double max) const;
    double truncatedConeVolume(double startHeight, double curHeight);
    void pressureForce(double pressure, Eigen::VectorXd &F);
    Eigen::Vector3d surfaceAreaNormal(const Eigen::VectorXd &q, int vidx);

    int frameno_;
    ProblemParameters params_;

    // The rendering thread reads the mesh and its edge data. Any function must lock this before writing to
    // to the mesh. (The rendering thread does not write to the mesh so reads from the worker thread do not
    // need to lock.)
    QMutex meshLock_;
};

#endif // SIMULATIONMESH_H
