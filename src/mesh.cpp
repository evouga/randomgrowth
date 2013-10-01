#include <OpenMesh/Core/IO/MeshIO.hh>
#include "mesh.h"
#include <fadbad.h>
#include <fadiff.h>
#include "autodifftemplates.h"
#include <GL/gl.h>


using namespace Eigen;
using namespace OpenMesh;
using namespace std;
using namespace fadbad;

typedef Eigen::Triplet<double> Tr;

Mesh::Mesh() :
    YoungsModulus_(1.0), PoissonRatio_(0.5), h_(1.0)
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

void Mesh::elasticEnergy(const VectorXd &q, const VectorXd &g, double &energy, VectorXd &gradq, VectorXd &gradg, Eigen::SparseMatrix<double> &hessq, Eigen::SparseMatrix<double> &hessg) const
{
    assert(q.size() == numdofs());
    assert(g.size() == numedges());
    energy = 0;
    gradq.resize(numdofs());
    gradq.setZero();
    gradg.resize(numedges());
    gradg.setZero();
    hessq.resize(numdofs(), numdofs());
    hessg.resize(numedges(), numedges());

    vector<Tr> Hqcoeffs;
    vector<Tr> Hgcoeffs;

    // bending energy
    std::cout << "bending energy" << std::endl;

    double bendcoeff = h_*h_*h_*YoungsModulus_/24.0/(1.0+PoissonRatio_);
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        if(mesh_->is_boundary(vi.handle()))
            continue;

        F<F<double> > stencilenergy = 0;
        // bending energy Laplacian term
        double Lcoeff = 1.0/(1.0-PoissonRatio_);
        vector<double> spokelens;
        vector<int> spokeidx;
        vector<double> rightopplens;
        vector<int> rightoppidx;

        vector<OMMesh::Point> nbq;
        vector<int> nbidx;
        for(OMMesh::VertexOHalfedgeIter voh = mesh_->voh_iter(vi.handle()); voh; ++voh)
        {
            OMMesh::HalfedgeHandle heh = voh.handle();
            int eidx = mesh_->edge_handle(heh).idx();
            spokelens.push_back(g[eidx]);
            spokeidx.push_back(eidx);

            OMMesh::VertexOHalfedgeIter nextoh = voh;
            ++nextoh;
            if(!nextoh)
                nextoh = mesh_->voh_iter(vi.handle());

            OMMesh::VertexHandle nextvert = mesh_->to_vertex_handle(nextoh.handle());

            OMMesh::HalfedgeHandle opp = mesh_->next_halfedge_handle(heh);;
            if(mesh_->to_vertex_handle(opp) != nextvert)
            {
                opp = mesh_->prev_halfedge_handle(mesh_->opposite_halfedge_handle(heh));
                assert(mesh_->from_vertex_handle(opp) == nextvert);
            }

            int oidx = mesh_->edge_handle(opp).idx();
            rightopplens.push_back(g[oidx]);
            rightoppidx.push_back(oidx);

            OMMesh::VertexHandle vh = mesh_->to_vertex_handle(heh);
            nbq.push_back(mesh_->point(vh));
            nbidx.push_back(vh.idx());
        }

        int numnbs = (int)spokelens.size();
        int diffqvars = 3*(numnbs+1);
        int diffvars = 3*(numnbs+1) + 2*numnbs;

        F<double> *dspokelens = new F<double>[numnbs];
        F<F<double> > *ddspokelens = new F<F<double> >[numnbs];
        for(int i=0; i<numnbs; i++)
        {
            dspokelens[i] = spokelens[i];
            dspokelens[i].diff(diffqvars+i, diffvars);
            ddspokelens[i] = dspokelens[i];
            ddspokelens[i].diff(diffqvars+i, diffvars);
        }
        F<double> *dopplens = new F<double>[numnbs];
        F<F<double> > *ddopplens = new F<F<double> >[numnbs];
        for(int i=0; i<numnbs; i++)
        {            
            dopplens[i] = rightopplens[i];
            dopplens[i].diff(diffqvars+numnbs+i, diffvars);
            ddopplens[i] = dopplens[i];
            ddopplens[i].diff(diffqvars+numnbs+i, diffvars);
        }

        F<double> dq[3];
        F<F<double> > ddq[3];
        int vidx = vi.handle().idx();
        OMMesh::Point q = mesh_->point(vi.handle());
        for(int i=0; i<3; i++)
        {
            dq[i] = q[i];
            dq[i].diff(i,diffvars);
            ddq[i] = dq[i];
            ddq[i].diff(i, diffvars);
        }

        F<double> *dnbq[3];
        F<F<double> > *ddnbq[3];
        for(int i=0; i<3; i++)
        {
            dnbq[i] = new F<double>[numnbs];
            ddnbq[i] = new F<F<double> >[numnbs];
            for(int j=0; j<numnbs; j++)
            {
                dnbq[i][j] = nbq[j][i];
                dnbq[i][j].diff(3+3*j+i,diffvars);
                ddnbq[i][j] = dnbq[i][j];
                ddnbq[i][j].diff(3+3*j+i,diffvars);
            }
        }


        F<F<double> > Lr[3];
        for(int j=0; j<3; j++)
            Lr[j] = L(ddq[j], numnbs, ddnbq[j], ddspokelens, ddopplens);
        F<F<double> > materialArea = dualarea(numnbs, ddspokelens, ddopplens);

        stencilenergy += bendcoeff*Lcoeff*normSquared(Lr)/materialArea;

        delete[] dspokelens;
        delete[] ddspokelens;
        delete[] dopplens;
        delete[] ddopplens;
        for(int i=0; i<3; i++)
        {
            delete[] dnbq[i];
            delete[] ddnbq[i];
        }

        // bending energy det term
        vector<F<double> *> nbdq;
        vector<F<F<double> > *> nbddq;
        for(int i=0; i<numnbs; i++)
        {
            F<double> *dq = new F<double>[3];
            nbdq.push_back(dq);
            F<F<double> > *ddq = new F<F<double> >[3];
            nbddq.push_back(ddq);
            for(int j=0; j<3; j++)
            {
                dq[j] = nbq[i][j];
                dq[j].diff(3+3*i+j, diffvars);
                ddq[j] = dq[j];
                ddq[j].diff(3+3*i+j, diffvars);
            }
        }

        F<F<double> > Kcurv = K(ddq, nbddq);
        F<F<double> > embeddedArea = dualarea(ddq, nbddq);
        stencilenergy += bendcoeff*Kcurv*embeddedArea*embeddedArea/materialArea;

        for(int i=0; i<numnbs; i++)
        {
            delete[] nbdq[i];
            delete[] nbddq[i];
        }

        energy += stencilenergy.val().val();
        for(int j=0; j<3; j++)
        {
            gradq[3*vidx+j] += stencilenergy.d(j).val();
        }
        for(int i=0; i<numnbs; i++)
        {
            for(int j=0; j<3; j++)
            {
                gradq[3*nbidx[i]+j] += stencilenergy.d(3+3*i+j).val();
            }
            gradg[spokeidx[i]] += stencilenergy.d(diffqvars+i).val();
            gradg[rightoppidx[i]] += stencilenergy.d(diffqvars+numnbs+i).val();
        }

        //TODO Hessian
    }

    std::cout << energy << std::endl;

    std::cout << "stretching energy" << std::endl;
    // Stretching energy
    double stretchcoeff = h_*YoungsModulus_/8.0/(1.0+PoissonRatio_);

    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        F<F<double> > stencilenergy = 0;
        int vidx[3];
        OMMesh::Point pts[3];
        double elens[3];
        int eidx[3];
        int idx=0;
        for(OMMesh::FaceHalfedgeIter fei = mesh_->fh_iter(fi.handle()); fei; ++fei)
        {
            vidx[idx] = mesh_->to_vertex_handle(fei.handle()).idx();
            pts[idx] = mesh_->point(mesh_->to_vertex_handle(fei.handle()));
            int eid = mesh_->edge_handle(fei.handle()).idx();
            eidx[idx] = eid;
            elens[idx] = g[eid];
            idx++;
        }
        assert(idx==3);

        int diffvars = 12;
        int diffqvars = 9;

        F<double> dq[3][3];
        F<F<double> > ddq[3][3];
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                dq[i][j] = pts[i][j];
                dq[i][j].diff(3*i+j, diffvars);
                ddq[i][j] = dq[i][j];
                ddq[i][j].diff(3*i+j, diffvars);
            }
        }

        F<double> delen[3];
        F<F<double> > ddelen[3];
        for(int i=0; i<3; i++)
        {
            delen[i] = elens[i];
            delen[i].diff(diffqvars+i, diffvars);
            ddelen[i] = delen[i];
            ddelen[i].diff(diffqvars+i, diffvars);
        }

        F<F<double> > emblens[3];
        for(int i=0; i<3; i++)
        {
            int previdx = (3+i-1)%3;
            F<F<double> > vec[3];
            diff(ddq[i],ddq[previdx],vec);
            emblens[i] = norm(vec);
        }

        F<F<double> > g[4];
        F<F<double> > a[4];
        metrictri(ddelen[0], ddelen[1], ddelen[2], g);
        metrictri(emblens[0], emblens[1], emblens[2], a);
        F<F<double> > ginva[4];
        invGtimesH(g, a, ginva);
        ginva[0] -= 1.0;
        ginva[3] -= 1.0;
        F<F<double> > matarea = 0.5*sqrt(det(g));
        stencilenergy += matarea*stretchcoeff/(1.0-PoissonRatio_)*tr(ginva);
        stencilenergy += matarea*stretchcoeff*-2.0*det(ginva);

        energy += stencilenergy.val().val();

        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                gradq[3*vidx[i]+j] += stencilenergy.d(3*i+j).val();
            }
            gradg[eidx[i]] += stencilenergy.d(diffqvars + i).val();
        }
        //TODO
    }

    hessq.setFromTriplets(Hqcoeffs.begin(), Hqcoeffs.end());
    hessg.setFromTriplets(Hgcoeffs.begin(), Hgcoeffs.end());
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

    Eigen::VectorXd q(numdofs());
    Eigen::VectorXd g(numedges());
    dofsFromGeometry(q, g);
    Eigen::VectorXd gradq, gradg;
    Eigen::SparseMatrix<double> hessq, hessg;
    double energy;
    std::cout << "elastic energy" << std::endl;
    elasticEnergy(q, g, energy, gradq, gradg, hessq, hessg);
    std::cout << "Energy " << energy << std::endl;
    std::cout << gradg.transpose() << std::endl;
    return success;
}

double Mesh::getPoissonRatio() const
{
    return PoissonRatio_;
}

double Mesh::getYoungsModulus() const
{
    return YoungsModulus_;
}

double Mesh::getThickness() const
{
    return h_;
}

void Mesh::setPoissonRatio(double val)
{
    PoissonRatio_ = val;
}

void Mesh::setYoungsModulus(double val)
{
    YoungsModulus_ = val;
}

void Mesh::setThickness(double val)
{
    h_ = val;
}

void Mesh::setIntrinsicLengthsToCurrentLengths()
{
    for(OMMesh::EdgeIter ei = mesh_->edges_begin(); ei != mesh_->edges_end(); ++ei)
    {
        double length = mesh_->calc_edge_length(ei.handle());
        mesh_->data(ei.handle()).setRestlen(length);
    }
}
