#include "mesh.h"
#include <iomanip>
#include "controller.h"
#include <Eigen/Dense>
#include "newton.h"

using namespace std;
using namespace Eigen;
using namespace OpenMesh;

void Mesh::elasticEnergy(const VectorXd &q,
                         const VectorXd &g,
                         double &energyB,
                         double &energyS,
                         VectorXd &gradq,
                         Eigen::SparseMatrix<double> &hessq,
                         Eigen::SparseMatrix<double> &gradggradq,
                         int derivativesRequested) const
{
    assert(q.size() == numdofs());
    assert(g.size() == numedges());
    energyB = energyS = 0;

    if(derivativesRequested & ElasticEnergy::DR_DQ)
    {
        gradq.resize(numdofs());
        gradq.setZero();
        if(derivativesRequested & ElasticEnergy::DR_HQ)
            hessq.resize(numdofs(), numdofs());
    }

    if(derivativesRequested & ElasticEnergy::DR_DGDQ)
    {
        gradggradq.resize(numedges(), numdofs());
    }


    vector<Tr> Hqcoeffs;
    vector<Tr> Hgcoeffs;
    vector<Tr> dgdqcoeffs;

    // bending energy
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        if(mesh_->is_boundary(vi.handle()))
            continue;

        vector<int> spokeidx;
        vector<int> rightoppidx;
        vector<int> nbidx;
        for(OMMesh::VertexOHalfedgeIter voh = mesh_->voh_iter(vi.handle()); voh; ++voh)
        {
            OMMesh::HalfedgeHandle heh = voh.handle();
            int eidx = mesh_->edge_handle(heh).idx();
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
            rightoppidx.push_back(oidx);

            OMMesh::VertexHandle vh = mesh_->to_vertex_handle(heh);
            nbidx.push_back(vh.idx());
        }

        int centidx = vi.handle().idx();

        energyB += ElasticEnergy::bendingEnergy(q, g, centidx, nbidx, spokeidx, rightoppidx, gradq, Hqcoeffs, dgdqcoeffs, params_, derivativesRequested);
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

        energyS += ElasticEnergy::stretchingEnergy(q, g, qidx, gidx, gradq, Hqcoeffs, dgdqcoeffs, params_, derivativesRequested);
    }

    if(derivativesRequested & ElasticEnergy::DR_HQ)
        hessq.setFromTriplets(Hqcoeffs.begin(), Hqcoeffs.end());
    if(derivativesRequested & ElasticEnergy::DR_DGDQ)
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
        SparseMatrix<double> gradggradq;
        int derivs = ElasticEnergy::DR_NONE;
        m_.elasticEnergy(q, g_, energyB, energyS, gradq, hessq, gradggradq, derivs);
        return energyB+energyS;
    }

    virtual double getEnergyAndDerivatives(const VectorXd &q, VectorXd &grad, SparseMatrix<double> &hess) const
    {
        double energyB, energyS;
        SparseMatrix<double> dgdq;
        int derivs = ElasticEnergy::DR_DQ | ElasticEnergy::DR_HQ;
        m_.elasticEnergy(q, g_, energyB, energyS, grad, hess, dgdq, derivs);
        return energyB+energyS;
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

class ImplicitEulerStep : public NewtonObjective
{
public:
    ImplicitEulerStep(Mesh &m, Controller &cont, const VectorXd &g, const VectorXd &qi, const VectorXd &vi, double h) :
        m_(m), cont_(cont), g_(g), qi_(qi), vi_(vi), h_(h) {}

    virtual double getEnergy(const VectorXd &q) const
    {
        double energyB, energyS;
        VectorXd gradq;
        SparseMatrix<double> hessq;
        SparseMatrix<double> gradggradq;
        int derivs = ElasticEnergy::DR_NONE;
        m_.elasticEnergy(qi_+h_*q, g_, energyB, energyS, gradq, hessq, gradggradq, derivs);
        SparseMatrix<double> M;
        m_.buildMassMatrix(qi_, M);

        return 0.5*(q-vi_).dot(M*(q-vi_)) + energyB + energyS;
    }

    virtual double getEnergyAndDerivatives(const VectorXd &q, VectorXd &grad, SparseMatrix<double> &hess) const
    {
        double energyB, energyS;
        VectorXd gradq;
        SparseMatrix<double> hessq;
        SparseMatrix<double> gradggradq;
        int derivs = ElasticEnergy::DR_DQ | ElasticEnergy::DR_HQ;
        m_.elasticEnergy(qi_+h_*q, g_, energyB, energyS, gradq, hessq, gradggradq, derivs);

        SparseMatrix<double> M;
        m_.buildMassMatrix(qi_, M);

        grad = M*(q-vi_) + h_ * gradq;

        hess = M + h_*h_*hessq;

        return 0.5*(q-vi_).dot(M*(q-vi_)) + energyB + energyS;
    }

    virtual void showCurrentIteration(const VectorXd &) const
    {
        m_.dofsToGeometry(qi_, g_);
        cont_.updateGL();
    }

private:
    Mesh &m_;
    Controller &cont_;
    const VectorXd &g_;
    const VectorXd &qi_;
    const VectorXd &vi_;
    double h_;
};

class MetricFit : public NewtonObjective
{
public:
    MetricFit(Mesh &m, Controller &cont, const VectorXd &q) : m_(m), cont_(cont), q_(q) {}

    virtual double getEnergy(const VectorXd &q) const
    {
        VectorXd gradq;
        SparseMatrix<double> dgdq, hessq;
        int derivs = ElasticEnergy::DR_DQ;
        double energyB, energyS;
        m_.elasticEnergy(q_, q, energyB, energyS, gradq, hessq, dgdq, derivs);
        return gradq.norm();
    }

    virtual double getEnergyAndDerivatives(const VectorXd &q, VectorXd &grad, SparseMatrix<double> &hess) const
    {
        SparseMatrix<double> dgdq, hessq;
        VectorXd gradq;
        double energyB, energyS;
        int derivs = ElasticEnergy::DR_DQ | ElasticEnergy::DR_DGDQ;
        m_.elasticEnergy(q_, q, energyB, energyS, gradq, hessq, dgdq, derivs);
        grad = dgdq*gradq;
        hess = dgdq*dgdq.transpose();
        return gradq.norm();
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


bool Mesh::simulate(Controller &cont)
{
    VectorXd q(numdofs());
    VectorXd g(numedges());
    dofsFromGeometry(q, g);

    double h = params_.eulerTimestep;
    VectorXd v(numdofs());
    v.setZero();
    int numsteps = params_.numEulerIters;

    for(int i=0; i<numsteps; i++)
    {
        ImplicitEulerStep obj(*this, cont, g, q, v, h);
        Newton n(obj);
        NewtonParameters params;
        params.tol = params_.tol;
        params.maxiters = params_.maxiters;
        params.lsmaxiters = params_.maxlinesearchiters;
        VectorXd newv(numdofs());
        Newton::SolverStatus ss = n.solve(params, v, newv);
        std::cout << n.solverStatusMessage(ss) << std::endl;

        v = newv*params_.dampingCoeff;
        q += h*v;
    }

    dofsToGeometry(q, g);
    cont.updateGL();
    return true;
}

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

    delete obj;

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

void Mesh::buildMassMatrix(const VectorXd &q, Eigen::SparseMatrix<double> &M) const
{
    M.resize(numdofs(), numdofs());
    vector<Tr> entries;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        int vidx = vi.handle().idx();
        double area = barycentricDualArea(q, vidx);
        for(int i=0; i<3; i++)
            entries.push_back(Tr(3*vidx+i, 3*vidx+i, area));
    }

    M.setFromTriplets(entries.begin(), entries.end());
}

double Mesh::barycentricDualArea(const VectorXd &q, int vidx) const
{
    double result = 0;
    OMMesh::VertexHandle vh = mesh_->vertex_handle(vidx);
    for(OMMesh::VertexFaceIter vfi = mesh_->vf_iter(vh); vfi; ++vfi)
    {
        result += faceArea(q, vfi.handle().idx());
    }
    return result/3.0;
}

double Mesh::faceArea(const VectorXd &q, int fidx) const
{
    FaceHandle fh = mesh_->face_handle(fidx);
    int verts[3];
    int idx=0;
    for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fh); fvi; ++fvi)
    {
        verts[idx++] = fvi.handle().idx();
    }

    Vector3d q0 = q.segment<3>(3*verts[0]);
    Vector3d q1 = q.segment<3>(3*verts[1]);
    Vector3d q2 = q.segment<3>(3*verts[2]);

    double A = ((q1-q0).cross(q2-q0)).norm();
    return 0.5*A;
}

