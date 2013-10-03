#include "mesh.h"
#include <fadbad.h>
#include <fadiff.h>
#include <iomanip>
#include "autodifftemplates.h"

using namespace std;
using namespace Eigen;
using namespace fadbad;
using namespace OpenMesh;

typedef Eigen::Triplet<double> Tr;

void Mesh::elasticEnergy(const VectorXd &q, const VectorXd &g, double &energy) const
{
    EnergyDerivatives derivs = NONE;
    VectorXd gradq, gradg;
    SparseMatrix<double> hessq, hessg;
    elasticEnergy(q, g, energy, gradq, gradg, hessq, hessg, derivs);
}

void Mesh::elasticEnergyG(const VectorXd &q, const VectorXd &g, double &energy, VectorXd &gradg, Eigen::SparseMatrix<double> &hessg) const
{
    EnergyDerivatives derivs = G;
    VectorXd gradq;
    SparseMatrix<double> hessq;
    elasticEnergy(q, g, energy, gradq, gradg, hessq, hessg, derivs);
}

void Mesh::elasticEnergyQ(const VectorXd &q, const VectorXd &g, double &energy, VectorXd &gradq, Eigen::SparseMatrix<double> &hessq) const
{
    EnergyDerivatives derivs = Q;
    VectorXd gradg;
    SparseMatrix<double> hessg;
    elasticEnergy(q, g, energy, gradq, gradg, hessq, hessg, derivs);
}

void Mesh::elasticEnergy(const VectorXd &q, const VectorXd &g, double &energy, VectorXd &gradq, VectorXd &gradg, Eigen::SparseMatrix<double> &hessq, Eigen::SparseMatrix<double> &hessg, EnergyDerivatives derivs) const
{
    assert(q.size() == numdofs());
    assert(g.size() == numedges());
    energy = 0;

    if(derivs & Q)
    {
        gradq.resize(numdofs());
        gradq.setZero();
        hessq.resize(numdofs(), numdofs());
    }

    if(derivs & G)
    {
        gradg.resize(numedges());
        gradg.setZero();
        hessg.resize(numedges(), numedges());
    }


    vector<Tr> Hqcoeffs;
    vector<Tr> Hgcoeffs;

    // bending energy
    double bendcoeff = params_.h*params_.h*params_.h*params_.YoungsModulus/24.0/(1.0+params_.PoissonRatio);
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        if(mesh_->is_boundary(vi.handle()))
            continue;

        F<F<double> > stencilenergy = 0;
        // bending energy Laplacian term
        double Lcoeff = 1.0/(1.0-params_.PoissonRatio);
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
        int diffqvars;
        if(derivs & Q)
            diffqvars = 3*(numnbs+1);
        else
            diffqvars = 0;

        int diffvars;
        if(derivs & G)
            diffvars = diffqvars + 2*numnbs;
        else
            diffvars = diffqvars;

        F<double> *dspokelens = new F<double>[numnbs];
        F<F<double> > *ddspokelens = new F<F<double> >[numnbs];
        for(int i=0; i<numnbs; i++)
        {
            dspokelens[i] = spokelens[i];
            if(derivs & G)
                dspokelens[i].diff(diffqvars+i, diffvars);
            ddspokelens[i] = dspokelens[i];
            if(derivs & G)
                ddspokelens[i].diff(diffqvars+i, diffvars);
        }
        F<double> *dopplens = new F<double>[numnbs];
        F<F<double> > *ddopplens = new F<F<double> >[numnbs];
        for(int i=0; i<numnbs; i++)
        {
            dopplens[i] = rightopplens[i];
            if(derivs & G)
                dopplens[i].diff(diffqvars+numnbs+i, diffvars);
            ddopplens[i] = dopplens[i];
            if(derivs & G)
                ddopplens[i].diff(diffqvars+numnbs+i, diffvars);
        }

        F<double> dq[3];
        F<F<double> > ddq[3];
        int vidx = vi.handle().idx();
        OMMesh::Point q = mesh_->point(vi.handle());
        for(int i=0; i<3; i++)
        {
            dq[i] = q[i];
            if(derivs & Q)
                dq[i].diff(i,diffvars);
            ddq[i] = dq[i];
            if(derivs & Q)
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
                if(derivs & Q)
                    dnbq[i][j].diff(3+3*j+i,diffvars);
                ddnbq[i][j] = dnbq[i][j];
                if(derivs & Q)
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
                if(derivs & Q)
                    dq[j].diff(3+3*i+j, diffvars);
                ddq[j] = dq[j];
                if(derivs & Q)
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
        energy += stencilenergy.val().val();

        if(derivs & Q)
        {
            for(int j=0; j<3; j++)
            {
                gradq[3*vidx+j] += stencilenergy.d(j).val();
            }
        }

        for(int i=0; i<numnbs; i++)
        {
            if(derivs & Q)
            {
                for(int j=0; j<3; j++)
                {
                    gradq[3*nbidx[i]+j] += stencilenergy.d(3+3*i+j).val();
                }
            }
            if(derivs & G)
            {
                gradg[spokeidx[i]] += stencilenergy.d(diffqvars+i).val();
                gradg[rightoppidx[i]] += stencilenergy.d(diffqvars+numnbs+i).val();
            }
        }

        if(derivs & Q)
        {
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
        }

        if(derivs & G)
        {
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
    }

    // Stretching energy
    double stretchcoeff = params_.h*params_.YoungsModulus/8.0/(1.0+params_.PoissonRatio);

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

        int diffqvars;
        if(derivs & Q)
            diffqvars = 9;
        else
            diffqvars = 0;

        int diffvars;
        if(derivs & G)
            diffvars = diffqvars + 3;
        else
            diffvars = diffqvars;

        F<double> dq[3][3];
        F<F<double> > ddq[3][3];
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                dq[i][j] = pts[i][j];
                if(derivs & Q)
                    dq[i][j].diff(3*i+j, diffvars);
                ddq[i][j] = dq[i][j];
                if(derivs & Q)
                    ddq[i][j].diff(3*i+j, diffvars);
            }
        }

        F<double> delen[3];
        F<F<double> > ddelen[3];
        for(int i=0; i<3; i++)
        {
            delen[i] = elens[i];
            if(derivs & G)
                delen[i].diff(diffqvars+i, diffvars);
            ddelen[i] = delen[i];
            if(derivs & G)
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
        stencilenergy += matarea*stretchcoeff/(1.0-params_.PoissonRatio)*tr(ginva)*tr(ginva);
        stencilenergy += matarea*stretchcoeff*-2.0*det(ginva);

        energy += stencilenergy.val().val();

        for(int i=0; i<3; i++)
        {
            if(derivs & Q)
            {
                for(int j=0; j<3; j++)
                {
                    gradq[3*vidx[i]+j] += stencilenergy.d(3*i+j).val();
                }
            }
            if(derivs & G)
                gradg[eidx[i]] += stencilenergy.d(diffqvars + i).val();
        }

        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                if(derivs & Q)
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
                }

                if(derivs & G)
                {
                    double hess = stencilenergy.d(diffqvars+i).d(diffqvars+j);
                    if(hess != 0)
                    {
                        Hgcoeffs.push_back(Tr(eidx[i],eidx[j],hess));
                    }
                }
            }
        }
    }

    if(derivs & Q)
        hessq.setFromTriplets(Hqcoeffs.begin(), Hqcoeffs.end());
    if(derivs & G)
        hessg.setFromTriplets(Hgcoeffs.begin(), Hgcoeffs.end());
}

double Mesh::triangleInequalityLineSearch(const VectorXd &g, const VectorXd &dg) const
{
    assert(g.size() == numedges());
    assert(dg.size() == numedges());

    double maxt = std::numeric_limits<double>::infinity();

    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        vector<double> gs;
        vector<double> dgs;
        for(OMMesh::FaceEdgeIter fei = mesh_->fe_iter(fi.handle()); fei; ++fei)
        {
            gs.push_back(g[fei.handle().idx()]);
            dgs.push_back(dg[fei.handle().idx()]);
        }
        assert(gs.size() == 3);
        for(int i=0; i<3; i++)
        {
            int idx[3];
            for(int j=0; j<3; j++)
            {
                idx[j] = (i+j)%3;
            }
            double thismax = triangleInequalityLineSearch(gs[idx[0]], gs[idx[1]], gs[idx[2]],
                    dgs[idx[0]], dgs[idx[1]], dgs[idx[2]]);

            maxt = std::min(maxt, thismax);
        }
    }

    return maxt;
}

double Mesh::triangleInequalityLineSearch(double g0, double g1, double g2, double dg0, double dg1, double dg2) const
{
    double num = g0+g1-g2;
    double denom = dg2-dg0-dg1;
    if(denom == 0)
        return std::numeric_limits<double>::infinity();
    double cand = num/denom;
    if(cand < 0)
        return std::numeric_limits<double>::infinity();
    return cand;
}

bool Mesh::relaxIntrinsicLengths()
{
    VectorXd q(numdofs());
    VectorXd g(numedges());
    dofsFromGeometry(q, g);

    VectorXd dg(numedges());

    SparseMatrix<double> hg;

    double energy;

    elasticEnergyG(q, g, energy, dg, hg);

    for(int i=0; i<params_.maxiters; i++)
    {
        if(dg.norm() < params_.tol)
            break;
        SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
        solver.compute(hg);
        VectorXd searchdir = -solver.solve(dg);
        std::cout << std::fixed << std::setprecision(8) << "Iter " << i+1 << "   E " << energy << "   |dg| " << dg.norm() << "   sd = " << searchdir.dot(dg);

        double stepsize = triangleInequalityLineSearch(g, searchdir);
        stepsize = std::min(1.0, 0.9*stepsize);
        std::cout << std::fixed << std::setprecision(8) << "   h0 = " << stepsize;
        double initialenergy = energy;

        VectorXd newg;
        int lsiters = 0;
        bool abort = false;

        do
        {
            newg = g + stepsize*searchdir;
            elasticEnergyG(q, newg, energy, dg, hg);
            stepsize /= 2.0;
            if(++lsiters > params_.maxlinesearchiters)
            {
                abort = true;
                break;
            }
        }
        while(energy > initialenergy);

        if(abort)
            break;

        g = newg;
        std::cout << std::fixed << std::setprecision(8) << "   h " << stepsize*2.0 << std::endl;
    }
    dofsToGeometry(q, g);
    if(dg.norm() < params_.tol)
    {
        std::cout << "Converged, final E " << energy << std::endl;
        return true;
    }
    std::cout << "Failed to converge" << std::endl;
    return false;
}
