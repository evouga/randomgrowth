#include "mesh.h"
#include <fadbad.h>
#include <fadiff.h>
#include <badiff.h>
#include <iomanip>
#include "autodifftemplates.h"
#include "controller.h"
#include <Eigen/Dense>
#include "newton.h"

using namespace std;
using namespace Eigen;
using namespace fadbad;
using namespace OpenMesh;

void Mesh::elasticEnergy(const VectorXd &q, const VectorXd &g, double &energyB, double &energyS) const
{
    EnergyDerivatives derivs = NONE;
    VectorXd gradq, gradg;
    SparseMatrix<double> hessq, hessg, gradggradq;
    elasticEnergy(q, g, energyB, energyS, gradq, gradg, hessq, hessg, gradggradq, derivs);
}

void Mesh::elasticEnergyG(const VectorXd &q,
                          const VectorXd &g,
                          double &energyB,
                          double &energyS,
                          VectorXd &gradg,
                          Eigen::SparseMatrix<double> &hessg) const
{
    EnergyDerivatives derivs = G;
    VectorXd gradq;
    SparseMatrix<double> hessq, gradggradq;
    elasticEnergy(q, g, energyB, energyS, gradq, gradg, hessq, hessg, gradggradq, derivs);
}

void Mesh::elasticEnergyQ(const VectorXd &q,
                          const VectorXd &g,
                          double &energyB,
                          double &energyS,
                          VectorXd &gradq,
                          Eigen::SparseMatrix<double> &hessq) const
{
    EnergyDerivatives derivs = Q;
    VectorXd gradg;
    SparseMatrix<double> hessg, gradggradq;
    elasticEnergy(q, g, energyB, energyS, gradq, gradg, hessq, hessg, gradggradq, derivs);
}

void Mesh::elasticEnergyGQ(const VectorXd &q, const VectorXd &g, VectorXd &gradq, Eigen::SparseMatrix<double> &gradggradq)
{
    double energyB, energyS;
    VectorXd gradg;
    SparseMatrix<double> hessg, hessq;
    EnergyDerivatives derivs = BOTH;
    elasticEnergy(q, g, energyB, energyS, gradq, gradg, hessq, hessg, gradggradq, derivs);
}

void Mesh::elasticEnergy(const VectorXd &q,
                         const VectorXd &g,
                         double &energyB,
                         double &energyS,
                         VectorXd &gradq,
                         VectorXd &gradg,
                         Eigen::SparseMatrix<double> &hessq,
                         Eigen::SparseMatrix<double> &hessg,
                         Eigen::SparseMatrix<double> &gradggradq,
                         EnergyDerivatives derivs) const
{
    assert(q.size() == numdofs());
    assert(g.size() == numedges());
    energyB = energyS = 0;

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

    if(derivs & G && derivs & Q)
    {
        gradggradq.resize(numedges(), numdofs());
    }


    vector<Tr> Hqcoeffs;
    vector<Tr> Hgcoeffs;
    vector<Tr> dgdqcoeffs;

    // bending energy
    double bendcoeff = params_.h*params_.h*params_.YoungsModulus/24.0/(1.0+params_.PoissonRatio);
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        if(mesh_->is_boundary(vi.handle()))
            continue;

        B<F<double> > stencilenergy = 0;
        // bending energy Laplacian term
        double Lcoeff = 1.0/(1.0-params_.PoissonRatio);
        vector<double> spokelens;
        vector<int> spokeidx;
        vector<double> rightopplens;
        vector<int> rightoppidx;

        vector<Vector3d> nbq;
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
            nbq.push_back(q.segment<3>(3*vh.idx()));
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
        B<F<double> > *ddspokelens = new B<F<double> >[numnbs];
        for(int i=0; i<numnbs; i++)
        {
            dspokelens[i] = spokelens[i];
            if(derivs & G)
                dspokelens[i].diff(diffqvars+i, diffvars);
            ddspokelens[i] = dspokelens[i];
            //if(derivs & G)
            //    ddspokelens[i].diff(diffqvars+i, diffvars);
        }
        F<double> *dopplens = new F<double>[numnbs];
        B<F<double> > *ddopplens = new B<F<double> >[numnbs];
        for(int i=0; i<numnbs; i++)
        {
            dopplens[i] = rightopplens[i];
            if(derivs & G)
                dopplens[i].diff(diffqvars+numnbs+i, diffvars);
            ddopplens[i] = dopplens[i];
            //if(derivs & G)
            //    ddopplens[i].diff(diffqvars+numnbs+i, diffvars);
        }

        F<double> dq[3];
        B<F<double> > ddq[3];
        int vidx = vi.handle().idx();
        Vector3d thisq = q.segment<3>(3*vi.handle().idx());
        for(int i=0; i<3; i++)
        {
            dq[i] = thisq[i];
            if(derivs & Q)
                dq[i].diff(i,diffvars);
            ddq[i] = dq[i];
        }

        F<double> *dnbq[3];
        B<F<double> > *ddnbq[3];
        for(int i=0; i<3; i++)
        {
            dnbq[i] = new F<double>[numnbs];
            ddnbq[i] = new B<F<double> >[numnbs];
            for(int j=0; j<numnbs; j++)
            {
                dnbq[i][j] = nbq[j][i];
                if(derivs & Q)
                    dnbq[i][j].diff(3+3*j+i,diffvars);
                ddnbq[i][j] = dnbq[i][j];
            }
        }

        B<F<double> > test;

        {
            B<F<double> > Lr[3];
            for(int j=0; j<3; j++)
                Lr[j] = L(ddq[j], numnbs, ddnbq[j], ddspokelens, ddopplens);
            B<F<double> > materialArea = dualbarycentricarea(numnbs, ddspokelens, ddopplens);

            //stencilenergy += bendcoeff*Lcoeff*normSquared(Lr)/materialArea;

            // bending energy det term
            vector<B<F<double> > *> nbddq;
            for(int i=0; i<numnbs; i++)
            {
                B<F<double> > *ddq = new B<F<double> >[3];
                nbddq.push_back(ddq);
                for(int j=0; j<3; j++)
                {
                    ddq[j] = ddnbq[j][i];
                }
            }

            B<F<double> > Kcurv = K(ddq, nbddq);
            test = dualbarycentricarea(ddq, nbddq);
            B<F<double> > embeddedArea = dualbarycentricarea(ddq, nbddq);
            stencilenergy += -2.0*bendcoeff*Kcurv*embeddedArea*embeddedArea/materialArea;
            for(int i=0; i<numnbs; i++)
                delete[] nbddq[i];
        }

        VectorXd Dq(numdofs());
        Dq.setZero();
        vector<Tr> hq, dgdq;
        int dr = ElasticEnergy::DR_DQ | ElasticEnergy::DR_HQ | ElasticEnergy::DR_DGDQ;
        cout << "e " << ElasticEnergy::bendTwo(q, g, vidx, nbidx, spokeidx, rightoppidx, Dq, hq, dgdq, params_, dr);
        cout << " " << stencilenergy.val().val() << endl;
        SparseMatrix<double> Hq(numdofs(), numdofs());
        Hq.setFromTriplets(hq.begin(), hq.end());
        SparseMatrix<double> DgDq(numedges(), numdofs());
        DgDq.setFromTriplets(dgdq.begin(), dgdq.end());

//        test.diff(0,1);
//        for(int i=0; i<3; i++)
//            for(int j=0; j<3; j++)
//            {
//                cout << ddq[i].d(0).d(j) << " ";
//            }
//        cout << endl;
//        exit(0);

        stencilenergy.diff(0,1);

        energyB += stencilenergy.val().val();

        if(derivs & Q)
        {
            for(int j=0; j<3; j++)
            {
                gradq[3*vidx+j] += ddq[j].d(0).val();
                cout << "dq " << ddq[j].d(0).val() << " " << Dq[3*vidx+j] << endl;
            }
        }

        for(int i=0; i<numnbs; i++)
        {
            if(derivs & Q)
            {
                for(int j=0; j<3; j++)
                {
                    gradq[3*nbidx[i]+j] += ddnbq[j][i].d(0).val();
                    cout << "dq " << ddnbq[j][i].d(0).val() << " " << Dq[3*nbidx[i]+j] << endl;
                }
            }
            if(derivs & G)
            {
                gradg[spokeidx[i]] += ddspokelens[i].d(0).val();
                gradg[rightoppidx[i]] += ddopplens[i].d(0).val();
            }
        }

        if(derivs & Q)
        {
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    double hess = ddq[j].d(0).d(k);
                    if(hess != 0)
                    {
                        Hqcoeffs.push_back(Tr(3*vidx+j,3*vidx+k,hess));
                        cout << "hcent " << hess << " " << Hq.coeffRef(3*vidx+j, 3*vidx+k) << endl;
                    }
                }
                for(int k=0; k<numnbs; k++)
                {
                    for(int l=0; l<3; l++)
                    {
                        double hess = ddq[j].d(0).d(3+3*k+l);
                        if(hess != 0)
                        {
                            Hqcoeffs.push_back(Tr(3*vidx+j,3*nbidx[k]+l,hess));
                            cout << "hmix " << hess << " " << Hq.coeffRef(3*vidx+j, 3*nbidx[k]+l) << endl;
                            Hqcoeffs.push_back(Tr(3*nbidx[k]+l,3*vidx+j,hess));
                            cout << "hmix " << hess << " " << Hq.coeffRef(3*nbidx[k]+l, 3*vidx+j) << endl;
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
                            double hess = ddnbq[k][i].d(0).d(3+3*j+l);
                            if(hess != 0)
                            {
                                Hqcoeffs.push_back(Tr(3*nbidx[i]+k,3*nbidx[j]+l,hess));
                                if(i==j) cout << "=";
                                cout << "hnb " << hess << " " << Hq.coeffRef(3*nbidx[i]+k, 3*nbidx[j]+l) << endl;
                            }
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
                    double hess = ddspokelens[i].d(0).d(diffqvars+j);
                    if(hess != 0)
                        Hgcoeffs.push_back(Tr(spokeidx[i],spokeidx[j],hess));
                }
                for(int j=0; j<numnbs; j++)
                {
                    double hess = ddspokelens[i].d(0).d(diffqvars+numnbs+j);
                    if(hess != 0)
                    {
                        Hgcoeffs.push_back(Tr(spokeidx[i],rightoppidx[j],hess));
                        Hgcoeffs.push_back(Tr(rightoppidx[j],spokeidx[i],hess));
                    }
                }
                for(int j=0; j<numnbs; j++)
                {
                    double hess = ddopplens[i].d(0).d(diffqvars+numnbs+j);
                    if(hess != 0)
                        Hgcoeffs.push_back(Tr(rightoppidx[i],rightoppidx[j],hess));
                }

                if(derivs & Q)
                {
                    for(int j=0; j<3; j++)
                    {
                        double hess = ddspokelens[i].d(0).d(j);
                        if(hess != 0)
                        {
                            dgdqcoeffs.push_back(Tr(spokeidx[i],3*vidx+j,hess));
                            std::cout << "mix " << hess << " " << DgDq.coeffRef(spokeidx[i], 3*vidx+j) << endl;
                        }
                        hess = ddopplens[i].d(0).d(j);
                        if(hess != 0)
                        {
                            dgdqcoeffs.push_back(Tr(rightoppidx[i],3*vidx+j,hess));
                            std::cout << "mixon " << hess << " " << DgDq.coeffRef(rightoppidx[i],3*vidx+j) << endl;
                        }
                    }

                    for(int j=0; j<numnbs; j++)
                    {
                        for(int k=0; k<3; k++)
                        {
                            double hess = ddspokelens[i].d(0).d(3+3*j+k);
                            if(hess != 0)
                            {
                                dgdqcoeffs.push_back(Tr(spokeidx[i], 3*nbidx[j]+k, hess));
                                std::cout << "mix " << hess << " " << DgDq.coeffRef(spokeidx[i], 3*nbidx[j]+k) << endl;
                            }

                            hess = ddopplens[i].d(0).d(3+3*j+k);
                            if(hess != 0)
                            {
                                dgdqcoeffs.push_back(Tr(rightoppidx[i], 3*nbidx[j]+k, hess));
                                std::cout << "mixoc " << hess << " " << DgDq.coeffRef(rightoppidx[i], 3*nbidx[j]+k) << endl;
                            }
                        }
                    }
                }
            }
        }

        delete[] dspokelens;
        delete[] ddspokelens;
        delete[] dopplens;
        delete[] ddopplens;
        for(int i=0; i<3; i++)
        {
            delete[] dnbq[i];
            delete[] ddnbq[i];
        }

    }

    // Stretching energy

    for(OMMesh::FaceIter it = mesh_->faces_begin(); it != mesh_->faces_end(); ++it)
    {
        int qidx[3];
        int gidx[3];

        int idx=0;
        for(OMMesh::FaceHalfedgeIter fhi = mesh_->fh_iter(it.handle()); fhi; ++fhi)
        {
            assert(idx < 3);
            OMMesh::HalfedgeHandle heh = fhi.handle();
            OMMesh::EdgeHandle eh = mesh_->edge_handle(heh);
            OMMesh::VertexHandle from = mesh_->from_vertex_handle(heh);
            gidx[idx] = eh.idx();
            qidx[(idx+1)%3] = from.idx();
            idx++;
        }
        assert(idx == 3);

        int dr = ElasticEnergy::DR_NONE;
        dr |= ElasticEnergy::DR_DQ;

        if(derivs & Q)
        {
            dr |= ElasticEnergy::DR_HQ;
        }
        if(derivs & Q && derivs & G)
        {
            dr |= ElasticEnergy::DR_DGDQ;
        }
        energyS += ElasticEnergy::stretchingEnergy(q, g, qidx, gidx, gradq, Hqcoeffs, dgdqcoeffs, params_, dr);
    }

    if(derivs & Q)
        hessq.setFromTriplets(Hqcoeffs.begin(), Hqcoeffs.end());
    if(derivs & G)
        hessg.setFromTriplets(Hgcoeffs.begin(), Hgcoeffs.end());
    if(derivs & G && derivs & Q)
        gradggradq.setFromTriplets(dgdqcoeffs.begin(), dgdqcoeffs.end());
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

class EmbeddingMinimizer : public NewtonObjective
{
public:
    EmbeddingMinimizer(Mesh &m, Controller &cont, const VectorXd &g) : m_(m), cont_(cont), g_(g) {}

    virtual double getEnergy(const VectorXd &q) const
    {
        double energyB, energyS;
        VectorXd gradq;
        SparseMatrix<double> hessq;
        m_.elasticEnergyQ(q, g_, energyB, energyS, gradq, hessq);
        return energyB+energyS;
    }

    virtual void getGradient(const VectorXd &q, VectorXd &grad) const
    {
        double energyB, energyS;
        SparseMatrix<double> hessq;
        m_.elasticEnergyQ(q, g_, energyB, energyS, grad, hessq);
    }

    virtual void getHessian(const VectorXd &q, Eigen::SparseMatrix<double> &hess) const
    {
        double energyB, energyS;
        VectorXd grad;
        m_.elasticEnergyQ(q, g_, energyB, energyS, grad, hess);
    }

    virtual void showCurrentIteration(const VectorXd &q) const
    {
        m_.dofsToGeometry(q, g_);
        cont_.updateGL();
    }

private:
    Mesh &m_;
    Controller &cont_;
    const VectorXd &g_;
};

class MetricFit : public NewtonObjective
{
public:
    MetricFit(Mesh &m, Controller &cont, const VectorXd &q) : m_(m), cont_(cont), q_(q) {}

    virtual double getEnergy(const VectorXd &q) const
    {
        VectorXd gradq;
        SparseMatrix<double> dgdq;
        m_.elasticEnergyGQ(q_, q, gradq, dgdq);
        return gradq.norm();
    }

    virtual void getGradient(const VectorXd &q, VectorXd &grad) const
    {
        SparseMatrix<double> dgdq;
        VectorXd gradq;
        m_.elasticEnergyGQ(q_, q, gradq, dgdq);
        grad = dgdq*gradq;
    }

    virtual void getHessian(const VectorXd &q, Eigen::SparseMatrix<double> &hess) const
    {
        SparseMatrix<double> dgdq;
        VectorXd gradq;
        m_.elasticEnergyGQ(q_, q, gradq, dgdq);
        hess = dgdq*dgdq.transpose();
    }

    virtual void showCurrentIteration(const VectorXd &q) const
    {
        m_.dofsToGeometry(q_, q);
        cont_.updateGL();
    }

private:
    Mesh &m_;
    Controller &cont_;
    const VectorXd &q_;
};

bool Mesh::relaxEnergy(Controller &cont, RelaxationType type)
{
    VectorXd q(numdofs());
    VectorXd g(numedges());
    dofsFromGeometry(q, g);

    NewtonObjective *obj = NULL;
    VectorXd guess;

    if(type == RelaxEmbedding)
    {
        obj = new EmbeddingMinimizer(*this, cont, g);
        guess = q;
    }
    else if(type == FitMetric)
    {
        obj = new MetricFit(*this, cont, q);
        guess = g;
    }
    else
        return false;

    Newton n(*obj);
    NewtonParameters params;
    params.tol = params_.tol;
    params.maxiters = params_.maxiters;
    params.lsmaxiters = params_.maxlinesearchiters;
    VectorXd result;
    Newton::SolverStatus ss = n.solve(params, guess, result);
    std::cout << n.solverStatusMessage(ss) << std::endl;

    if(ss != Newton::CONVERGED)
        return false;

    if(type == RelaxEmbedding)
    {
        q = result;
    }
    else if(type == FitMetric)
    {
        g = result;
    }

    dofsToGeometry(q, g);
    cont.updateGL();
    return true;
}

bool Mesh::largestMagnitudeEigenvalue(const Eigen::SparseMatrix<double> &M, double &eigenvalue)
{
    int dim = M.cols();
    VectorXd v(dim);
    v.setRandom();
    v.normalize();
    double curestimate=0;
    int iter=0;

    for(iter=0; iter < params_.maxpoweriters; iter++)
    {
        v = M*v;        
        v.normalize();

        curestimate = v.dot(M*v);
        double err = (M*v-curestimate*v).norm();
        if(err < params_.powertol)
            break;
    }
    eigenvalue = curestimate;
    return iter < params_.maxpoweriters;
}

bool Mesh::smallestEigenvalue(const Eigen::SparseMatrix<double> &M, double &eigenvalue)
{
    double largesteval;
    if(!largestMagnitudeEigenvalue(M, largesteval))
        return false;

    if(largesteval < 0)
    {
        eigenvalue = largesteval;
        return true;
    }

    int dim = M.cols();
    SparseMatrix<double> shift(dim,dim);
    shift.setIdentity();
    shift *= -largesteval;
    SparseMatrix<double> newM = M + shift;
    double newlargest;
    if(!largestMagnitudeEigenvalue(newM, newlargest))
        return false;
    eigenvalue = newlargest+largesteval;
    return true;
}
