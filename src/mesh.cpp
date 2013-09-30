#include "mesh.h"
#include <fadbad.h>
#include <fadiff.h>
#include "autodifftemplates.h"

using namespace Eigen;
using namespace OpenMesh;
using namespace std;
using namespace fadbad;

typedef Eigen::Triplet<double> Tr;

Mesh::Mesh(double YoungsModulus, double PoissonRatio, double h) :
    YoungsModulus_(YoungsModulus), PoissonRatio_(PoissonRatio), h_(h)
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

    int totvars = numdofs()+numedges();
    F<F<double> > ddenergy;

    // bending energy

    double bendcoeff = h_*h_*h_*YoungsModulus_/24.0/(1.0+PoissonRatio_);
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        if(mesh_->is_boundary(vi.handle()))
            continue;

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
            spokelens.push_back(mesh_->calc_edge_length(heh));
            spokeidx.push_back(mesh_->edge_handle(heh).idx());
            OMMesh::HalfedgeHandle opp = mesh_->next_halfedge_handle(heh);
            rightopplens.push_back(mesh_->calc_edge_length(opp));
            rightoppidx.push_back(mesh_->edge_handle(opp).idx());

            OMMesh::VertexHandle vh = mesh_->to_vertex_handle(heh);
            nbq.push_back(mesh_->point(vh));
            nbidx.push_back(vh.idx());
        }

        int numnbs = (int)spokelens.size();

        F<double> *dspokelens = new F<double>[numnbs];
        F<F<double> > *ddspokelens = new F<F<double> >[numnbs];
        for(int i=0; i<numnbs; i++)
        {
            dspokelens[i] = spokelens[i];
            dspokelens[i].diff(numdofs()+spokeidx[i], totvars);
            ddspokelens[i] = dspokelens[i];
            ddspokelens[i].diff(numdofs()+spokeidx[i], totvars);
        }
        F<double> *dopplens = new F<double>[numnbs];
        F<F<double> > *ddopplens = new F<F<double> >[numnbs];
        for(int i=0; i<numnbs; i++)
        {
            dopplens[i] = rightopplens[i];
            dopplens[i].diff(numdofs()+rightoppidx[i], totvars);
            ddopplens[i] = dopplens[i];
            ddopplens[i].diff(numdofs()+rightoppidx[i], totvars);
        }

        F<double> dq[3];
        F<F<double> > ddq[3];
        int vidx = vi.handle().idx();
        OMMesh::Point q = mesh_->point(vi.handle());
        for(int i=0; i<3; i++)
        {
            dq[i] = q[i];
            dq[i].diff(3*vidx+i,totvars);
            ddq[i] = dq[i];
            ddq[i].diff(3*vidx+i, totvars);
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
                dnbq[i][j].diff(3*nbidx[j]+i,totvars);
                ddnbq[i][j] = dnbq[i][j];
                ddnbq[i][j].diff(3*nbidx[j]+i,totvars);
            }
        }

        F<F<double> > Lr[3];
        for(int j=0; j<3; j++)
            Lr[j] = L(ddq[j], numnbs, ddnbq[j], ddspokelens, ddopplens);
        F<F<double> > materialArea = dualarea(numnbs, ddspokelens, ddopplens);

        ddenergy += bendcoeff*Lcoeff*normSquared(Lr)/materialArea;

        delete[] dspokelens;
        delete[] ddspokelens;
        delete[] dopplens;
        delete[] ddopplens;
        for(int i=0; i<3; i++)
        {
            delete[] dnbq[3];
            delete[] ddnbq[3];
        }

        // bending energy det term
        vector<F<double> *> nbdq;
        vector<F<F<double> > *> nbddq;
        for(int i=0; i<numnbs; i++)
        {
            F<double> *dq = new F<double>[3];
            F<F<double> > *ddq = new F<F<double> >[3];
            for(int j=0; j<3; j++)
            {
                dq[j] = nbq[i][j];
                dq[j].diff(nbidx[i]+j,totvars);
                ddq[j] = dq[j];
                ddq[j].diff(nbidx[i]+j, totvars);
            }
        }

        F<F<double> > Kcurv = K(ddq, nbddq);
        F<F<double> > embeddedArea = dualarea(ddq, nbddq);
        ddenergy += bendcoeff*Kcurv*embeddedArea*embeddedArea/materialArea;

        for(int i=0; i<numnbs; i++)
        {
            delete[] nbdq[i];
            delete[] nbddq[i];
        }
    }

    // Stretching energy
    double stretchcoeff = h_*YoungsModulus_/8.0/(1.0+PoissonRatio_);

    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        int vidx[3];
        OMMesh::Point pts[3];
        double elens[3];
        int eidx[3];
        int idx=0;
        for(OMMesh::FaceHalfedgeIter fei = mesh_->fh_iter(fi.handle()); fei; ++fei)
        {
            vidx[idx] = mesh_->to_vertex_handle(fei.handle()).idx();
            pts[idx] = mesh_->point(mesh_->to_vertex_handle(fei.handle()));
            eidx[idx] = mesh_->edge_handle(fei.handle()).idx();
            elens[idx] = mesh_->calc_edge_length(fei.handle());
            idx++;
        }
        assert(idx==3);

        F<double> dq[3][3];
        F<F<double> > ddq[3][3];
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                dq[i][j] = pts[i][j];
                dq[i][j].diff(3*vidx[i]+j, totvars);
                ddq[i][j] = dq[i][j];
                ddq[i][j].diff(3*vidx[i]+j, totvars);
            }
        }

        F<double> delen[3];
        F<F<double> > ddelen[3];
        for(int i=0; i<3; i++)
        {
            delen[i] = elens[i];
            delen[i].diff(numdofs()+eidx[i], totvars);
            ddelen[i] = delen[i];
            ddelen[i].diff(numdofs()+eidx[i], totvars);
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
        F<F<double> > matarea = sqrt(det(g));
        ddenergy += matarea*stretchcoeff/(1.0-PoissonRatio_)*tr(ginva);
        ddenergy += matarea*stretchcoeff*-2.0*det(ginva);
    }

    energy = ddenergy.val().val();

    // derivatives
    for(int i=0; i<numdofs(); i++)
    {
        gradq[i] = ddenergy.d(i).val();
        for(int j=0; j<numdofs(); j++)
        {
            double secondderiv = ddenergy.d(i).d(j);
            if(secondderiv != 0.0)
                Hqcoeffs.push_back(Tr(i,j,secondderiv));
        }
    }

    for(int i=0; i<numedges(); i++)
    {
        gradg[i] = ddenergy.d(numdofs()+i).val();
        for(int j=0; j<numedges(); j++)
        {
            double secondderiv = ddenergy.d(numdofs()+i).d(numdofs()+j);
            if(secondderiv != 0.0)
                Hgcoeffs.push_back(Tr(i,j,secondderiv));
        }
    }

    hessq.setFromTriplets(Hqcoeffs.begin(), Hqcoeffs.end());
    hessg.setFromTriplets(Hgcoeffs.begin(), Hgcoeffs.end());
}
