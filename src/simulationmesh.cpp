#include "simulationmesh.h"
#include <iomanip>
#include <Eigen/Geometry>
#include <fstream>

using namespace Eigen;
using namespace std;

SimulationMesh::SimulationMesh() : Mesh(), meshLock_(QMutex::Recursive)
{
    params_.scale = 0.085;
    params_.h = .0006; // 0.0027, 0.0021, 0.0009, 0.0007, 0.0004
    params_.YoungsModulus = 2e9;
    params_.PoissonRatio = 0.33;
    params_.rho = 500.0;
    params_.dampingCoeff = 1e-3;
    params_.eulerTimestep = 1e-6;
    params_.numEulerIters = numeric_limits<int>::max();

    params_.pullMag = 10;

    params_.smoothShade = true;
    params_.showWireframe = true;
    params_.outputDir = "output";

    frameno_ = 0;
}

void ProblemParameters::dumpParameters(ostream &os)
{
    ElasticParameters::dumpParameters(os);
    os << "rho " << rho << endl;
    os << "dampingCoeff " << dampingCoeff << endl;
    os << "eulerTimestep " << eulerTimestep << endl;
    os << "numEulerIters " << numEulerIters << endl;
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
        frameno_ = 0;
        energies_.resize(numFaces());
        energies_.setZero();
    }
    meshLock_.unlock();
    return success;
}

const ProblemParameters &SimulationMesh::getParameters() const
{
    return params_;
}

void SimulationMesh::setParameters(ProblemParameters params)
{
    params_ = params;
}

void SimulationMesh::dumpFrame()
{
    stringstream ss;
    ss << params_.outputDir;
    ss << "/frame_";
    ss << setfill('0') << setw(8) << frameno_ << ".obj";
    exportOBJ(ss.str().c_str());

//    stringstream ss2;
//    ss2 << params_.outputDir;
//    ss2 << "/frame_";
//    ss2 << setfill('0') << setw(8) << frameno_ << ".geo";

//    stringstream ss3;
//    ss3 << params_.outputDir;
//    ss3 << "/frame_";
//    ss3 << setfill('0') << setw(8) << frameno_ << ".g";

//    ofstream ofs(ss2.str().c_str());
//    VectorXd q,g;
//    dofsFromGeometry(q, g);
//    VectorXd Hdensity, Kdensity, vareas;
//    meanCurvature(q, Hdensity);
//    gaussianCurvature(q, Kdensity);
//    vertexAreas(q, vareas);

//    for(int i=0; i<(int)mesh_->n_vertices(); i++)
//    {
//        ofs << i << " " << Hdensity[i] << " " << Kdensity[i] << " " << vareas[i] << endl;
//    }

//    ofstream gofs(ss3.str().c_str());
//    gofs << g;

    frameno_++;
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

void SimulationMesh::setFlatCone(double height)
{
//    deleteBadFlatConeFaces();
//    VectorXd q, g;
//    dofsFromGeometry(q, g);
//    int numverts = mesh_->n_vertices();
//    for(int i=0; i<numverts; i++)
//    {
//        Vector3d pos = q.segment<3>(3*i);
//        double r = sqrt(1.0+height*height)*sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
//        double theta = atan2(pos[1], pos[0]) / sqrt(1.0+height*height);
//        q[3*i+0] = r*cos(theta);
//        q[3*i+1] = r*sin(theta);
//        q[3*i+2] = 0;
//    }
//    dofsToGeometry(q);
}

void SimulationMesh::deleteBadFlatConeFaces()
{
//    set<int> badfaces;

//    for(OMMesh::FaceIter fi=mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
//    {
//        for(OMMesh::FaceEdgeIter fei = mesh_->fe_iter(fi.handle()); fei; ++fei)
//        {
//            OMMesh::HalfedgeHandle heh = mesh_->halfedge_handle(fei.handle(),0);
//            OMMesh::VertexHandle topt = mesh_->to_vertex_handle(heh);
//            OMMesh::VertexHandle frompt = mesh_->from_vertex_handle(heh);
//            double p1y = mesh_->point(topt)[1];
//            double p2y = mesh_->point(frompt)[1];
//            double t = -p2y/(p1y-p2y);
//            if(t < 0 || t > 1)
//                continue;
//            double p1x = mesh_->point(topt)[0];
//            double p2x = mesh_->point(frompt)[0];
//            double x = t*p1x + (1-t)*p2x;
//            if(x < 0)
//                badfaces.insert(fi.handle().idx());
//        }
//    }

//    OMMesh *newmesh = new OMMesh;
//    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
//    {
//        OMMesh::Point pt = mesh_->point(vi.handle());
//        newmesh->add_vertex(pt);
//    }

//    for(int i=0; i<(int)mesh_->n_faces(); i++)
//    {
//        if(badfaces.count(i) > 0)
//            continue;

//        OMMesh::FaceHandle fh = mesh_->face_handle(i);
//        vector<OMMesh::VertexHandle> newface;
//        for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fh); fvi; ++fvi)
//        {
//            int idx = fvi.handle().idx();
//            OMMesh::VertexHandle vert = newmesh->vertex_handle(idx);
//            newface.push_back(vert);
//        }
//        newmesh->add_face(newface);
//    }

//    meshLock_.lock();
//    {
//        delete mesh_;
//        mesh_ = newmesh;
//    }
//    meshLock_.unlock();
}

Eigen::Vector3d SimulationMesh::faceNormal(int fidx) const
{
    Vector3d qs[3];
    for(int i=0; i<3; i++)
        qs[i] = vertPos(faceVerts(fidx)[i]);
    Vector3d result = (qs[1]-qs[0]).cross(qs[2]-qs[0]);
    return result / result.norm();
}
