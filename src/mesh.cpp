#include <OpenMesh/Core/IO/MeshIO.hh>
#include "mesh.h"

using namespace Eigen;
using namespace OpenMesh;
using namespace std;

Mesh::Mesh()
{
    params_.h = 1;
    params_.YoungsModulus = 1;
    params_.PoissonRatio = 0.5;
    params_.maxiters = 15;
    params_.maxlinesearchiters = 10;
    params_.tol = 1e-6;
    mesh_ = new OMMesh();
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

Vector3d Mesh::centroid()
{
    Vector3d centroid(0,0,0);
    int numpts = mesh_->n_vertices();

    for(int i=0; i<numpts; i++)
    {
        OMMesh::Point pt = mesh_->point(mesh_->vertex_handle(i));
        for(int j=0; j<3; j++)
            centroid[j] += pt[j];
    }
    centroid /= numpts;
    return centroid;
}

double Mesh::radius()
{
    Vector3d cent = centroid();
    int numpts = mesh_->n_vertices();
    double maxradius = 0;
    for(int i=0; i<numpts; i++)
    {
        OMMesh::Point pt = mesh_->point(mesh_->vertex_handle(i));
        Vector3d ept(pt[0],pt[1],pt[2]);
        double radius = (ept-cent).squaredNorm();
        if(radius > maxradius)
            maxradius = radius;
    }
    return sqrt(maxradius);
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
    OpenMesh::IO::Options opt;
    mesh_->request_face_normals();
    mesh_->request_vertex_normals();
    opt.set(OpenMesh::IO::Options::VertexNormal);
    bool success = OpenMesh::IO::read_mesh(*mesh_, filename, opt);
    mesh_->update_normals();

    setIntrinsicLengthsToCurrentLengths();    

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
    for(OMMesh::EdgeIter ei = mesh_->edges_begin(); ei != mesh_->edges_end(); ++ei)
    {
        double length = mesh_->calc_edge_length(ei.handle());
        mesh_->data(ei.handle()).setRestlen(length);
    }
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
    double L = mesh_->data(eh).restlen();
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
