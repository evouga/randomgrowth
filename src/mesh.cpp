#include <OpenMesh/Core/IO/MeshIO.hh>
#include "mesh.h"
#include <GL/gl.h>


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

void Mesh::render(bool showWireframe, bool smoothShade)
{

    glEnable(GL_LIGHTING);
    glEnable(GL_DITHER);

    glPolygonOffset(1.0, 1.0);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    if(smoothShade)
    {
        glShadeModel(GL_SMOOTH);
    }
    else
    {
        glShadeModel(GL_FLAT);
    }

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);

    static vector<GLfloat> colors;
    static vector<int> indices;
    static vector<GLfloat> pos;
    static vector<GLfloat> normal;

    colors.clear();
    indices.clear();
    pos.clear();
    normal.clear();

    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fi.handle()); fvi; ++fvi)
        {
            Vector3d color(0.0, 186/255., 0.0);
            OMMesh::VertexHandle v = fvi.handle();
            OMMesh::Point pt = mesh_->point(v);
            OMMesh::Point n;
            mesh_->calc_vertex_normal_correct(v, n);
            n.normalize();
            for(int j=0; j<3; j++)
            {
                pos.push_back(pt[j]);
                normal.push_back(n[j]);
                colors.push_back(color[j]);
            }
        }
    }

    glVertexPointer(3, GL_FLOAT, 0, &pos[0]);
    glNormalPointer(GL_FLOAT, 0, &normal[0]);
    glColorPointer(3, GL_FLOAT, 0, &colors[0]);

    int idx=0;
    for (int i=0; i<(int)mesh_->n_faces(); i++)
    {
        for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(mesh_->face_handle(i)); fvi; ++fvi)
        {
            indices.push_back(idx++);
        }
    }
    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisable(GL_POLYGON_OFFSET_FILL);

    if(showWireframe)
    {
        glLineWidth(1.0);
        glBegin(GL_LINES);
        for(OMMesh::ConstEdgeIter ei = mesh_->edges_begin(); ei != mesh_->edges_end(); ++ei)
        {
            glColor3f(0.0, 0.0, 0.0);
            OMMesh::Point pt1, pt2;
            edgeEndpoints(ei.handle(), pt1, pt2);
            glVertex3d(pt1[0], pt1[1], pt1[2]);
            glVertex3d(pt2[0], pt2[1], pt2[2]);
        }
        glEnd();
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
