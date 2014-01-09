#include <OpenMesh/Core/IO/MeshIO.hh>
#include "mesh.h"
#include <iomanip>
#include <Eigen/Geometry>
#include <fstream>

using namespace Eigen;
using namespace OpenMesh;
using namespace std;

Mesh::Mesh() : meshLock_(QMutex::Recursive)
{
    params_.scale = 0.1;
    params_.h = .0003;
    params_.YoungsModulus = 4e8;
    params_.PoissonRatio = 0.5;
    params_.rho = 500.0;
    params_.dampingCoeff = 0.0001;
    params_.eulerTimestep = 1e-6;
    params_.numEulerIters = 500000;

    params_.growthAmount = 50;
    params_.maxEdgeStrain = 1e8;
    params_.baseGrowthProbability = 0.5;

    params_.smoothShade = true;
    params_.showWireframe = true;
    params_.outputDir = "output";

    mesh_ = new OMMesh();
    undeformedMesh_ = new OMMesh();
    frameno_ = 0;
}

void ProblemParameters::dumpParameters(ostream &os)
{
    ElasticParameters::dumpParameters(os);
    os << "rho " << rho << endl;
    os << "dampingCoeff " << dampingCoeff << endl;
    os << "eulerTimestep " << eulerTimestep << endl;
    os << "numEulerIters " << numEulerIters << endl;
    os << "growthAmount " << growthAmount << endl;
    os << "baseGrowthProbability " << baseGrowthProbability << endl;
    os << "maxEdgeStrain " << maxEdgeStrain << endl;
}

int Mesh::numdofs() const
{
    return 3*mesh_->n_vertices();
}

int Mesh::numedges() const
{
    return mesh_->n_edges();
}

void Mesh::dofsFromGeometry(Eigen::VectorXd &q, Eigen::VectorXd &g) const
{
    q.resize(numdofs());
    g.resize(numedges());

    for(int i=0; i<(int)mesh_->n_vertices(); i++)
    {
        OMMesh::Point pt = mesh_->point(mesh_->vertex_handle(i));
        for(int j=0; j<3; j++)
            q[3*i+j] = pt[j];
    }

    for(int i=0; i<(int)mesh_->n_edges(); i++)
    {
        g[i] = mesh_->data(mesh_->edge_handle(i)).restlen();
    }
}

void Mesh::undeformedDofsFromGeometry(VectorXd &undefq, VectorXd &undefg) const
{
    undefq.resize(numdofs());

    for(int i=0; i<(int)undeformedMesh_->n_vertices(); i++)
    {
        OMMesh::Point pt = undeformedMesh_->point(undeformedMesh_->vertex_handle(i));
        for(int j=0; j<3; j++)
            undefq[3*i+j] = pt[j];
    }

    for(int i=0; i<(int)undeformedMesh_->n_edges(); i++)
    {
        undefg[i] = undeformedMesh_->data(undeformedMesh_->edge_handle(i)).restlen();
    }
}

void Mesh::dofsToGeometry(const VectorXd &q, const VectorXd &g)
{    
    meshLock_.lock();
    {
        assert(q.size() == numdofs());
        assert(g.size() == numedges());

        for(int i=0; i<(int)mesh_->n_vertices(); i++)
        {
            OMMesh::Point &pt = mesh_->point(mesh_->vertex_handle(i));
            for(int j=0; j<3; j++)
                pt[j] = q[3*i+j];
        }

        for(int i=0; i<(int)mesh_->n_edges(); i++)
        {
            mesh_->data(mesh_->edge_handle(i)).setRestlen(g[i]);
        }
    }
    meshLock_.unlock();
}

void Mesh::edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2)
{
    OMMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle(eh, 0);
    pt1 = mesh_->point(mesh_->from_vertex_handle(heh1));
    pt2 = mesh_->point(mesh_->to_vertex_handle(heh1));
}

bool Mesh::exportOBJ(const char *filename)
{
    OpenMesh::IO::Options opt;
    mesh_->request_face_normals();
    mesh_->request_vertex_normals();
    mesh_->update_normals();
    opt.set(OpenMesh::IO::Options::VertexNormal);
    return OpenMesh::IO::write_mesh(*mesh_, filename, opt);
}

bool Mesh::importOBJ(const char *filename)
{
    bool success = true;
    meshLock_.lock();
    {
        OpenMesh::IO::Options opt;
        mesh_->request_face_normals();
        mesh_->request_vertex_normals();
        opt.set(OpenMesh::IO::Options::VertexNormal);
        success = OpenMesh::IO::read_mesh(*mesh_, filename, opt);
        mesh_->update_normals();

        setIntrinsicLengthsToCurrentLengths();

        delete undeformedMesh_;
        undeformedMesh_ = new OMMesh(*mesh_);

        for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
            mesh_->point(vi.handle())[2] += randomRange(0,0.0001);
    }
    frameno_ = 0;
    meshLock_.unlock();
    return success;
}

const ProblemParameters &Mesh::getParameters() const
{
    return params_;
}

void Mesh::setParameters(ProblemParameters params)
{
    params_ = params;
}

void Mesh::setIntrinsicLengthsToCurrentLengths()
{
    meshLock_.lock();
    {
        for(OMMesh::EdgeIter ei = mesh_->edges_begin(); ei != mesh_->edges_end(); ++ei)
        {
            double length = mesh_->calc_edge_length(ei.handle());
            mesh_->data(ei.handle()).setRestlen(length);
        }
    }       
    meshLock_.unlock();
}

Vector3d Mesh::averageNormal(const VectorXd &q, int vidx) const
{
    OMMesh::VertexHandle vh = mesh_->vertex_handle(vidx);
    Vector3d result;
    int denom = 0;
    result.setZero();
    for(OMMesh::VertexFaceIter vfi = mesh_->vf_iter(vh); vfi; ++vfi)
    {
        result += faceNormal(q, vfi.handle().idx());
        denom++;
    }
    return result/denom;
}

Vector3d Mesh::faceNormal(const VectorXd &q, int fidx) const
{
    int vids[3];
    int idx=0;
    OMMesh::FaceHandle fh = mesh_->face_handle(fidx);
    for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fh); fvi; ++fvi)
    {
        vids[idx++] = fvi.handle().idx();
    }

    Vector3d v0 = q.segment<3>(3*vids[0]);
    Vector3d v1 = q.segment<3>(3*vids[1]);
    Vector3d v2 = q.segment<3>(3*vids[2]);
    Vector3d n = (v0-v1).cross(v2-v1);
    return n/n.norm();
}

double Mesh::infinityNorm(const VectorXd &v) const
{
    double maxval = 0;
    for(int i=0; i<v.size(); i++)
        maxval = std::max(fabs(v[i]), maxval);
    return maxval;
}

double Mesh::randomRange(double min, double max) const
{
    double val = (max-min)*double(rand())/double(RAND_MAX);
    return val+min;
}

void Mesh::dumpFrame()
{
    stringstream ss;
    ss << params_.outputDir;
    ss << "/frame_";
    ss << setfill('0') << setw(8) << frameno_ << ".obj";
    exportOBJ(ss.str().c_str());

    stringstream ss2;
    ss2 << params_.outputDir;
    ss2 << "/frame_";
    ss2 << setfill('0') << setw(8) << frameno_ << ".geo";

    ofstream ofs(ss2.str().c_str());
    VectorXd q,g;
    dofsFromGeometry(q, g);
    VectorXd Hdensity, Kdensity, vareas;
    meanCurvature(q, Hdensity);
    gaussianCurvature(q, Kdensity);
    vertexAreas(q, vareas);

    for(int i=0; i<(int)mesh_->n_vertices(); i++)
    {
        ofs << i << " " << Hdensity[i] << " " << Kdensity[i] << " " << vareas[i] << endl;
    }

    frameno_++;
}
