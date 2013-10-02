#include "mesh.h"
#include <fadbad.h>
#include <fadiff.h>
#include "autodifftemplates.h"

using namespace std;
using namespace Eigen;
using namespace fadbad;
using namespace OpenMesh;

typedef Eigen::Triplet<double> Tr;


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
        F<F<double> > materialArea = dualbarycentricarea(numnbs, ddspokelens, ddopplens);

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
        F<F<double> > embeddedArea = dualbarycentricarea(ddq, nbddq);
        stencilenergy += -2.0*bendcoeff*Kcurv*embeddedArea*embeddedArea/materialArea;

        for(int i=0; i<numnbs; i++)
        {
            delete[] nbdq[i];
            delete[] nbddq[i];
        }
        if(stencilenergy.val().val() < 0)
            std::cout << "stencil energy " << stencilenergy.val().val() << std::endl;
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

        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                double hess = stencilenergy.d(j).d(k);
                if(hess != 0)
                    Hqcoeffs.push_back(Tr(3*vidx+j,3*vidx+k,hess));
            }
            for(int k=0; k<numnbs; k++)
            {
                for(int l=0; l<3; l++)
                {
                    double hess = stencilenergy.d(j).d(3+3*k+j);
                    if(hess != 0)
                    {
                        Hqcoeffs.push_back(Tr(3*vidx+j,3*nbidx[k]+l,hess));
                        Hqcoeffs.push_back(Tr(3*nbidx[k]+l,3*vidx+j,hess));
                    }
                }
            }
        }
        for(int i=0; i<numnbs; i++)
        {
            for(int j=0; j<numnbs; j++)
            {
                for(int k=0; k<3; k++)
                {
                    for(int l=0; l<3; l++)
                    {
                        double hess = stencilenergy.d(3+3*i+k).d(3+3*j+l);
                        if(hess != 0)
                            Hqcoeffs.push_back(Tr(3*nbidx[i]+k,3*nbidx[j]+l,hess));
                    }
                }
            }
        }

        for(int i=0; i<numnbs; i++)
        {
            for(int j=0; j<numnbs; j++)
            {
                double hess = stencilenergy.d(diffqvars+i).d(diffqvars+j);
                if(hess != 0)
                    Hgcoeffs.push_back(Tr(spokeidx[i],spokeidx[j],hess));
            }
            for(int j=0; j<numnbs; j++)
            {
                double hess = stencilenergy.d(diffqvars+i).d(diffqvars+numnbs+j);
                if(hess != 0)
                {
                    Hgcoeffs.push_back(Tr(spokeidx[i],rightoppidx[j],hess));
                    Hgcoeffs.push_back(Tr(rightoppidx[j],spokeidx[i],hess));
                }
            }
            for(int j=0; j<numnbs; j++)
            {
                double hess = stencilenergy.d(diffqvars+numnbs+i).d(diffqvars+numnbs+j);
                if(hess != 0)
                    Hgcoeffs.push_back(Tr(rightoppidx[i],rightoppidx[j],hess));
            }
        }
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
        stencilenergy += matarea*stretchcoeff/(1.0-PoissonRatio_)*tr(ginva)*tr(ginva);
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

        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    for(int l=0; l<3; l++)
                    {
                        double hess = stencilenergy.d(3*i+k).d(3*j+l);
                        if(hess != 0.0)
                        {
                            Hqcoeffs.push_back(Tr(3*vidx[i]+k,3*vidx[j]+l,hess));
                        }
                    }
                }

                double hess = stencilenergy.d(diffqvars+i).d(diffqvars+j);
                if(hess != 0)
                {
                    Hgcoeffs.push_back(Tr(eidx[i],eidx[j],hess));
                }
            }
        }
    }

    hessq.setFromTriplets(Hqcoeffs.begin(), Hqcoeffs.end());
    hessg.setFromTriplets(Hgcoeffs.begin(), Hgcoeffs.end());
}

bool Mesh::relaxIntrinsicLengths()
{
    VectorXd q(numdofs());
    VectorXd g(numedges());
    dofsFromGeometry(q, g);

    VectorXd dq(numdofs());
    VectorXd dg(numedges());

    SparseMatrix<double> hq,hg;

    double energy;

    elasticEnergy(q, g, energy, dq, dg, hq, hg);
    std::cout << "Initial energy: " << energy << std::endl;

    for(int i=0; i<3; i++)
    {

        SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
        solver.compute(hg);
        VectorXd searchdir = -solver.solve(dg);
        std::cout << "searchdir dot " << searchdir.dot(dg) << endl;
        VectorXd testq(numedges());
        testq.setZero();
        testq[0] += 1.0;

        g += 0.5*searchdir;
        elasticEnergy(q, g, energy, dq, dg, hq, hg);
        std::cout << "new energy " << energy << std::endl;
    }

    return true;
}
