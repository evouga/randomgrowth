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
    params_.h = .0003; // 0.0027, 0.0021, 0.0009, 0.0007, 0.0004
    params_.YoungsModulus = 4e7;
    params_.PoissonRatio = 0.5;
    params_.rho = 500.0;
    params_.dampingCoeff = 1e-4;
    params_.eulerTimestep = 1e-6;
    params_.numEulerIters = 500000;

    params_.growthAmount = 50;
    params_.maxEdgeStrain = 1e8;
    params_.baseGrowthProbability = 0.5;

    params_.smoothShade = true;
    params_.showWireframe = true;
    params_.outputDir = "output";

    mesh_ = new OMMesh();
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
    if(q.size() != numdofs())
        q.resize(numdofs());
    if(g.size() != numedges())
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

void Mesh::targetMetricFromGeometry(VectorXd &targetg) const
{
    if(targetg.size() != numedges())
        targetg.resize(numedges());

    for(int i=0; i<numedges(); i++)
    {
        targetg[i] = mesh_->data(mesh_->edge_handle(i)).targetlen();
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

void Mesh::targetMetricToGeometry(const VectorXd &targetg)
{
    meshLock_.lock();
    {
        assert(targetg.size() == numedges());

        for(int i=0; i<(int)numedges(); i++)
        {
            mesh_->data(mesh_->edge_handle(i)).setTargetlen(targetg[i]);
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

bool Mesh::importMetric(const char *filename)
{
    bool success = true;
    meshLock_.lock();
    {
        VectorXd q, g;
        dofsFromGeometry(q, g);
        ifstream ifs(filename);
        if(!ifs)
        {
            success = false;
        }
        else
        {
            for(int i=0; i<numedges(); i++)
            {
                double newelen;
                ifs >> newelen;
                if(!ifs)
                {
                    success = false;
                    break;
                }
                g[i] = newelen;
            }
        }

        if(success)
        {
            dofsToGeometry(q, g);
            setNoTargetMetric();
        }

    }
    meshLock_.unlock();
    return success;
}

void Mesh::addRandomNoise(double magnitude)
{
    meshLock_.lock();
    {
        for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
            mesh_->point(vi.handle())[2] += randomRange(-magnitude,magnitude);
    }
    meshLock_.unlock();
}

void Mesh::setNoTargetMetric()
{
    meshLock_.lock();
    {
        VectorXd q, g;
        dofsFromGeometry(q, g);
        targetMetricToGeometry(g);
    }
    meshLock_.unlock();
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
        setNoTargetMetric();
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
    Vector3d n = (v2-v1).cross(v0-v1);
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

    stringstream ss3;
    ss3 << params_.outputDir;
    ss3 << "/frame_";
    ss3 << setfill('0') << setw(8) << frameno_ << ".g";

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

    ofstream gofs(ss3.str().c_str());
    gofs << g;

    frameno_++;
}

void Mesh::setConeHeights(double height)
{
    VectorXd q, g;
    dofsFromGeometry(q, g);
    int numverts = mesh_->n_vertices();
    for(int i=0; i<numverts; i++)
    {
        Vector3d pos = q.segment<3>(3*i);
        double newz = height*(1.0 - sqrt(pos[0]*pos[0] + pos[1]*pos[1]));
        q[3*i+2] = newz;
    }
    dofsToGeometry(q, g);
    setIntrinsicLengthsToCurrentLengths();
    setNoTargetMetric();
}

void Mesh::setFlatCone(double height)
{
    deleteBadFlatConeFaces();
    VectorXd q, g;
    dofsFromGeometry(q, g);
    int numverts = mesh_->n_vertices();
    for(int i=0; i<numverts; i++)
    {
        Vector3d pos = q.segment<3>(3*i);
        double r = sqrt(1.0+height*height)*sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
        double theta = atan2(pos[1], pos[0]) / sqrt(1.0+height*height);
        q[3*i+0] = r*cos(theta);
        q[3*i+1] = r*sin(theta);
        q[3*i+2] = 0;
    }
    dofsToGeometry(q, g);
    setIntrinsicLengthsToCurrentLengths();
    setNoTargetMetric();
}

void Mesh::deleteBadFlatConeFaces()
{
    set<int> badfaces;

    for(OMMesh::FaceIter fi=mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        for(OMMesh::FaceEdgeIter fei = mesh_->fe_iter(fi.handle()); fei; ++fei)
        {
            OMMesh::HalfedgeHandle heh = mesh_->halfedge_handle(fei.handle(),0);
            OMMesh::VertexHandle topt = mesh_->to_vertex_handle(heh);
            OMMesh::VertexHandle frompt = mesh_->from_vertex_handle(heh);
            double p1y = mesh_->point(topt)[1];
            double p2y = mesh_->point(frompt)[1];
            double t = -p2y/(p1y-p2y);
            if(t < 0 || t > 1)
                continue;
            double p1x = mesh_->point(topt)[0];
            double p2x = mesh_->point(frompt)[0];
            double x = t*p1x + (1-t)*p2x;
            if(x < 0)
                badfaces.insert(fi.handle().idx());
        }
    }

    OMMesh *newmesh = new OMMesh;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        OMMesh::Point pt = mesh_->point(vi.handle());
        newmesh->add_vertex(pt);
    }

    for(int i=0; i<(int)mesh_->n_faces(); i++)
    {
        if(badfaces.count(i) > 0)
            continue;

        OMMesh::FaceHandle fh = mesh_->face_handle(i);
        vector<OMMesh::VertexHandle> newface;
        for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fh); fvi; ++fvi)
        {
            int idx = fvi.handle().idx();
            OMMesh::VertexHandle vert = newmesh->vertex_handle(idx);
            newface.push_back(vert);
        }
        newmesh->add_face(newface);
    }

    meshLock_.lock();
    {
        delete mesh_;
        mesh_ = newmesh;
    }
    meshLock_.unlock();
}

Vector3d Mesh::surfaceAreaNormal(const VectorXd &q, int vidx)
{
    OMMesh::VertexHandle vh = mesh_->vertex_handle(vidx);
    Vector3d result;
    result.setZero();

    for(OMMesh::VertexFaceIter vfi = mesh_->vf_iter(vh); vfi; ++vfi)
    {
        double area = faceArea(q, vfi.handle().idx());
        Vector3d N = this->faceNormal(q, vfi.handle().idx());
        result += area*N;
    }
    result /= result.norm();
    return result;
}
