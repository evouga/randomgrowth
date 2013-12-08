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

    params_.growthAmount = 0.05;
    params_.growthRadius = 0.2;
    params_.growthTime = 1000;
    params_.newGrowthRate = 1000;

    params_.smoothShade = true;
    params_.showWireframe = true;

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
    os << "growthRadius " << growthRadius << endl;
    os << "growthTime " << growthTime << endl;
    os << "newGrowthRate " << newGrowthRate << endl;
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

void Mesh::growPlanarDisk(Vector2d center, double radius, double strainincrement, double maxstrain)
{
    assert(mesh_->n_edges() == undeformedMesh_->n_edges());
    meshLock_.lock();
    {
        for(OMMesh::EdgeIter ei = undeformedMesh_->edges_begin(); ei != undeformedMesh_->edges_end(); ++ei)
        {
            OMMesh::HalfedgeHandle heh = undeformedMesh_->halfedge_handle(ei.handle(), 0);
            OMMesh::VertexHandle v1 = undeformedMesh_->from_vertex_handle(heh);
            OMMesh::VertexHandle v2 = undeformedMesh_->to_vertex_handle(heh);
            Vector2d v1p, v2p;
            for(int j=0; j<2; j++)
            {
                v1p[j] = undeformedMesh_->point(v1)[j];
                v2p[j] = undeformedMesh_->point(v2)[j];
            }

            double dist1 = (v1p-center).norm();
            double dist2 = (v2p-center).norm();

            double fracinside = 0;

            if(dist1 > radius && dist2 > radius)
            {
                // both points outside
                continue;
            }
            else if(dist1 <= radius && dist2 <= radius)
            {
                // both points inside
                fracinside = 1;
            }
            else
            {
                if(dist1 > radius)
                    swap(v1p, v2p);

                //(v1-center).(v1-center) + 2t(v1-center).(v2-v1)+t^2(v2-v1).(v2-v1) = radius^2
                double a = (v2p-v1p).dot(v2p-v1p);
                double b = 2.0*(v1p-center).dot(v2p-v1p);
                double c = (v1p-center).dot(v1p-center) - radius*radius;
                double t = (-b + sqrt(b*b-4.0*a*c))/(2.0*a);
                assert(t >= 0 && t <= 1.0);
                Vector2d intpt = v1p + t*(v2p-v1p);
                fracinside = (intpt-v1p).norm()/(v2p-v1p).norm();
            }

            double origlen = (v2p-v1p).norm();
            OMMesh::EdgeHandle defeh = mesh_->edge_handle(ei.handle().idx());
            double curlen = mesh_->data(defeh).restlen();

            double maxlen = maxstrain*origlen;
            double newlen = fracinside*(1.0+strainincrement)*curlen + (1.0-fracinside)*curlen;
            newlen = min(newlen, maxlen);
            mesh_->data(defeh).setRestlen(newlen);
        }
    }
    meshLock_.unlock();
}

double Mesh::strainDensity(int edgeidx) const
{
    OMMesh::EdgeHandle eh = mesh_->edge_handle(edgeidx);
    OMMesh::HalfedgeHandle heh = mesh_->halfedge_handle(eh, 0);
    OMMesh::VertexHandle vh1 = mesh_->to_vertex_handle(heh);
    OMMesh::VertexHandle vh2 = mesh_->from_vertex_handle(heh);
    OMMesh::Point pt1 = mesh_->point(vh1);
    OMMesh::Point pt2 = mesh_->point(vh2);

    Vector3d edgevec(pt1[0]-pt2[0], pt1[1]-pt2[1],pt1[2]-pt2[2]);
    double l = edgevec.norm();
    double L = undeformedMesh_->data(undeformedMesh_->edge_handle(edgeidx)).restlen();//mesh_->data(eh).restlen();
    return (l-L)/L;
}

double Mesh::vertexStrainDensity(int vertidx) const
{
    OMMesh::VertexHandle vh = mesh_->vertex_handle(vertidx);
    double totstrain = 0;
    int numedges = 0;
    for(OMMesh::VertexEdgeIter vei = mesh_->ve_iter(vh); vei; ++vei)
    {
        totstrain += strainDensity(vei.handle().idx());
        numedges ++;
    }
    return totstrain/numedges;
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
    VectorXd Hdensity, Kdensity;
    meanCurvatureDensity(q, Hdensity);
    gaussianCurvatureDensity(q, Kdensity);

    for(int i=0; i<(int)mesh_->n_vertices(); i++)
    {
        ofs << i << " " << Hdensity[i] << " " << Kdensity[i] << endl;
    }

    frameno_++;
}
