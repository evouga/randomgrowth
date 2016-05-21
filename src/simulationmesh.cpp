#include "simulationmesh.h"
#include <iomanip>
#include <Eigen/Geometry>
#include <fstream>

using namespace Eigen;
using namespace std;

const double PI = 3.1415926535898;

SimulationMesh::SimulationMesh() : Mesh(), meshLock_(QMutex::Recursive)
{
    params_.scale = 0.085;
    params_.h = .0005; // 0.0027, 0.0021, 0.0009, 0.0007, 0.0004
    params_.YoungsModulus = 2e9;
    params_.PoissonRatio = 0.33;
    params_.rho = 500.0;
    params_.dampingCoeff = 0.001;
    params_.eulerTimestep = 5e-7;
    params_.simTime = 0.1;

    params_.pressureBehavior = params_.PB_CONSTANT;
    params_.constantPressureVal = 1.0; // in atmospheres
    params_.holeRadius = 0.0085;

    params_.constantVelocity = true;
    params_.crushMass = 2.0;
    params_.crushTime = 1.0;
    params_.initialVel = -2.0;

    params_.coneAngle = 0.25;
    params_.coneHeight = 0.330;

    params_.smoothShade = true;
    params_.showWireframe = true;
    params_.outputDir = "output";
    params_.restoreCheckpoint = "";    
}

void ProblemParameters::dumpParameters(ostream &os)
{
    ElasticParameters::dumpParameters(os);
    os << "rho " << rho << endl;
    os << "dampingCoeff " << dampingCoeff << endl;
    os << "eulerTimestep " << eulerTimestep << endl;
    os << "simTime " << simTime << endl;
    os << "pressureBehavior " << pressureBehavior << endl;
    os << "constantPressureVal " << constantPressureVal << endl;
    os << "holeRadius " << holeRadius << endl;
    os << "coneAngle " << coneAngle << endl;
    os << "coneHeight " << coneHeight << endl;
    os << "constantVelocity " << constantVelocity << endl;
    os << "crushTime " << crushTime << endl;
    os << "crushMass " << crushMass << endl;
    os << "initalVel " << initialVel << endl;
}

bool SimulationMesh::exportOBJ(const char *filename)
{
    return writeMesh(filename);
}

bool SimulationMesh::importOBJ(const char *filename)
{
    bool success = true;
    meshLock_.lock();
    {
        success = loadMesh(filename);
    }
    meshLock_.unlock();
    return success;
}

bool SimulationMesh::restoreCheckpoint(DynamicData &dd, const char *filename)
{
    ifstream ifs(filename);
    if(!ifs)
        return false;

    ifs.read(reinterpret_cast<char*>(&dd.iter), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&dd.frameno), sizeof(int));
    int numverts;
    ifs.read(reinterpret_cast<char*>(&numverts), sizeof(int));

    deformedPosition_.resize(numverts);
    dd.v.resize(numverts);

    for(int i=0; i<numverts; i++)
        ifs.read(reinterpret_cast<char*>(&deformedPosition_[i]), sizeof(double));
    for(int i=0; i<numverts; i++)
        ifs.read(reinterpret_cast<char*>(&dd.v[i]), sizeof(double));

    ifs.read(reinterpret_cast<char*>(&dd.initialV), sizeof(double));
    ifs.read(reinterpret_cast<char*>(&dd.planeVel), sizeof(double));
    ifs.read(reinterpret_cast<char*>(&dd.planeZ), sizeof(double));
    return ifs;
}

bool SimulationMesh::saveCheckpoint(const DynamicData &dd, const char *filename)
{
    ofstream ofs(filename);
    if(!ofs)
        return false;

    ofs.write(reinterpret_cast<const char*>(&dd.iter), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&dd.frameno), sizeof(int));

    int numverts = dd.v.size();
    ofs.write(reinterpret_cast<const char*>(&numverts), sizeof(int));

    for(int i=0; i<numverts; i++)
        ofs.write(reinterpret_cast<const char*>(&deformedPosition_[i]), sizeof(double));
    for(int i=0; i<numverts; i++)
        ofs.write(reinterpret_cast<const char*>(&dd.v[i]), sizeof(double));

    ofs.write(reinterpret_cast<const char*>(&dd.initialV), sizeof(double));
    ofs.write(reinterpret_cast<const char*>(&dd.planeVel), sizeof(double));
    ofs.write(reinterpret_cast<const char*>(&dd.planeZ), sizeof(double));
    return ofs;
}

const ProblemParameters &SimulationMesh::getParameters() const
{
    return params_;
}

void SimulationMesh::setParameters(ProblemParameters params)
{
    params_ = params;
}

void SimulationMesh::dumpFrame(const DynamicData &dd)
{
    stringstream ss;
    ss << params_.outputDir;
    ss << "/frame_";
    ss << setfill('0') << setw(8) << dd.frameno << ".obj";
    exportOBJ(ss.str().c_str());

    stringstream ss2;
    ss2 << params_.outputDir;
    ss2 << "/frame_";
    ss2 << setfill('0') << setw(8) << dd.frameno << ".chk";
    saveCheckpoint(dd, ss2.str().c_str());
}

void SimulationMesh::setConeHeights(double height)
{
    int numverts = numVertices();
    for(int i=0; i<numverts; i++)
    {
        Vector3d pos = vertPos(i);
        double newz = height*(1.0 - sqrt(pos[0]*pos[0] + pos[1]*pos[1]));
        vertPos(i)[2] = newz;
    }
    resetRestMetric();
}

Eigen::Vector3d SimulationMesh::faceNormal(int fidx) const
{
    Vector3d qs[3];
    for(int i=0; i<3; i++)
        qs[i] = vertPos(faceVerts(fidx)[i]);
    Vector3d result = (qs[1]-qs[0]).cross(qs[2]-qs[0]);
    return result / result.norm();
}

double SimulationMesh::truncatedConeVolume(double startHeight, double curHeight)
{
    double base = params_.scale*params_.scale*PI;
    double totV = base*params_.scale*startHeight/3.0;
    double topr = params_.scale*(1.0 - curHeight/startHeight);
    double topbase = PI*topr*topr;
    double topV = topbase*params_.scale*(startHeight-curHeight)/3.0;
    return totV-topV;
}

void SimulationMesh::pressureForce(double pressure, VectorXd &F)
{
    F.resize(3*numVertices());

    VectorXd areas(numVertices());
    metricBarycentricDualAreas(areas);

    VectorXd normals;
    vertexNormals(normals);

    for(int v=0; v<numVertices(); v++)
    {
        F.segment<3>(3*v) = pressure*areas[v]*normals.segment<3>(3*v);
    }
}

void SimulationMesh::metricBarycentricDualAreas(VectorXd &areas) const
{
    VectorXd faceareas(numFaces());
    for(int i=0; i<numFaces(); i++)
        faceareas[i] = Midedge::intrinsicArea(*this, i, params_);

    areas.resize(numVertices());
    areas.setZero();
    for(int i=0; i<numFaces(); i++)
    {
        for(int j=0; j<3; j++)
        {
            areas[faceVerts(i)[j]] += faceareas[i];
        }
    }
    areas /= 3.0;
}

void SimulationMesh::vertexNormals(VectorXd &vertNormals)
{
    vertNormals.resize(3*numVertices());
    int nfaces = numFaces();

    VectorXd faceNormals(3*nfaces);
    for(int i=0; i<nfaces; i++)
    {
        Vector3d normal = faceNormal(i);
        normal *= deformedFaceArea(i);
        faceNormals.segment<3>(3*i) = normal;
    }

    vertNormals.setZero();
    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<3; j++)
        {
            vertNormals.segment<3>(3*faceVerts(i)[j]) += faceNormals.segment<3>(3*i);
        }
    }

    for(int i=0; i<numVertices(); i++)
    {
        vertNormals.segment<3>(3*i) /= vertNormals.segment<3>(3*i).norm();
    }
}

double SimulationMesh::deformedFaceArea(int fidx) const
{
    Vector3d q0 = vertPos(faceVerts(fidx)[0]);
    Vector3d q1 = vertPos(faceVerts(fidx)[1]);
    Vector3d q2 = vertPos(faceVerts(fidx)[2]);

    double A = ((q1-q0).cross(q2-q0)).norm();
    return 0.5*params_.scale*params_.scale*A;
}

