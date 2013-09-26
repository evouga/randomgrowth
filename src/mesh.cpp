#include "mesh.h"
#include "autodifftemplates.h"

using namespace Eigen;
using namespace OpenMesh;

Mesh::Mesh()
{
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
