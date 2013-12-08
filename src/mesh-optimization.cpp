#include "mesh.h"
#include <iomanip>
#include "controller.h"
#include <Eigen/Dense>
#include <fstream>

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

bool Mesh::simulate(Controller &cont)
{
    {
        string name = params_.outputDir + "/parameters";
        ofstream paramfile(name.c_str());
        params_.dumpParameters(paramfile);
    }

    VectorXd q(numdofs());
    VectorXd g(numedges());
    dofsFromGeometry(q, g);

    double h = params_.eulerTimestep;
    VectorXd v(numdofs());
    v.setZero();
    int numsteps = params_.numEulerIters;

    Vector2d o(0,0);

    for(int i=0; i<numsteps; i++)
    {
        if(i%params_.newGrowthRate == 0)
        {
            double x,y;
            do
            {
                x = randomRange(-1.0,1.0);
                y = randomRange(-1.0,1.0);
            }
            while(x*x+y*y > 1);
            o[0] = x;
            o[1] = y;
        }
        q += h*v;
        int derivs = ElasticEnergy::DR_DQ;
        VectorXd gradq;
        double energyB, energyS;
        SparseMatrix<double> hessq, gradggradq;
        elasticEnergy(q, g, energyB, energyS, gradq, hessq, gradggradq, derivs);        

        SparseMatrix<double> Minv;
        buildInvMassMatrix(q, Minv);

        v -= h*Minv*gradq + h*Minv*params_.dampingCoeff*v;
        dofsToGeometry(q, g);
        if(i%10 == 0)
            dumpFrame();

        growPlanarDisk(o, params_.growthRadius, params_.growthAmount/params_.growthTime, 1.0+4.0*params_.growthAmount);
        dofsFromGeometry(q, g);
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
        double mass = area*params_.rho*params_.h*params_.scale*params_.scale*params_.scale;
        if(mesh_->is_boundary(vi.handle()))
            mass = std::numeric_limits<double>::infinity();
        for(int i=0; i<3; i++)
            entries.push_back(Tr(3*vidx+i, 3*vidx+i, mass));
    }

    M.setFromTriplets(entries.begin(), entries.end());
}

void Mesh::buildInvMassMatrix(const VectorXd &q, Eigen::SparseMatrix<double> &Minv) const
{
    Minv.resize(numdofs(), numdofs());
    vector<Tr> entries;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        int vidx = vi.handle().idx();
        double area = barycentricDualArea(q, vidx);
        double invmass = 1.0/area/params_.rho/params_.h/params_.scale/params_.scale/params_.scale;
        if(mesh_->is_boundary(vi.handle()))
            invmass = 0;
        for(int i=0; i<3; i++)
            entries.push_back(Tr(3*vidx+i, 3*vidx+i, invmass));
    }

    Minv.setFromTriplets(entries.begin(), entries.end());
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

void Mesh::gaussianCurvatureDensity(const VectorXd &q, VectorXd &Kdensity) const
{
    int numverts = mesh_->n_vertices();
    Kdensity.resize(numverts);
    Kdensity.setZero();

    for(int i=0; i<numverts; i++)
    {
        OMMesh::VertexHandle vh = mesh_->vertex_handle(i);
        if(mesh_->is_boundary(vh))
            continue;

        Vector3d cent = q.segment<3>(3*i);

        double angsum = 0;

        vector<int> star;

        for(OMMesh::VertexOHalfedgeIter voh = mesh_->voh_iter(vh); voh; ++voh)
        {
            star.push_back(mesh_->to_vertex_handle(voh.handle()).idx());
        }

        int numnbs = star.size();
        for(int j=0; j<numnbs; j++)
        {
            Vector3d v1 = q.segment<3>(3*star[j]);
            Vector3d v2 = q.segment<3>(3*star[(j+1)%numnbs]);
            double sintheta = ((v1-cent).cross(v2-cent)).norm();
            double costheta = (v1-cent).dot(v2-cent);
            angsum += atan2(sintheta, costheta);
        }

        Kdensity[i] = 2*PI-angsum;
        double area = barycentricDualArea(q, i);
        Kdensity[i] /= area;
    }
}

void Mesh::meanCurvatureDensity(const VectorXd &q, VectorXd &Hdensity) const
{
    int numverts = mesh_->n_vertices();
    VectorXd x(numverts), y(numverts), z(numverts);
    for(int i=0; i<numverts; i++)
    {
        x[i] = q[3*i];
        y[i] = q[3*i+1];
        z[i] = q[3*i+2];
    }

    SparseMatrix<double> L;
    buildExtrinsicDirichletLaplacian(q, L);
    VectorXd Hx = L*x;
    VectorXd Hy = L*y;
    VectorXd Hz = L*z;

    Hdensity.resize(numverts);
    for(int i=0; i<numverts; i++)
    {
        double area = barycentricDualArea(q, i);
        Vector3d Hnormal(Hx[i], Hy[i], Hz[i]);
        double sign = 1.0;

        Vector3d avnormal = averageNormal(q, i);
        if(avnormal.dot(Hnormal) < 0)
            sign = -1.0;

        Hdensity[i] = sign * Hnormal.norm()/area;
    }
}

double Mesh::cotanWeight(int edgeid, const VectorXd &q) const
{
    OMMesh::EdgeHandle eh = mesh_->edge_handle(edgeid);
    double weight = 0;
    for(int i=0; i<2; i++)
    {
        OMMesh::HalfedgeHandle heh = mesh_->halfedge_handle(eh,i);

        if(mesh_->is_boundary(heh))
            continue;
        OMMesh::HalfedgeHandle next = mesh_->next_halfedge_handle(heh);

        Vector3d e1, e2;
        OMMesh::VertexHandle oppv = mesh_->to_vertex_handle(next);
        OMMesh::VertexHandle v1 = mesh_->to_vertex_handle(heh);
        OMMesh::VertexHandle v2 = mesh_->from_vertex_handle(heh);
        e1 = q.segment<3>(3*v1.idx())-q.segment<3>(3*oppv.idx());
        e2 = q.segment<3>(3*v2.idx())-q.segment<3>(3*oppv.idx());
        double cosang = e1.dot(e2);
        double sinang = (e1.cross(e2)).norm();
        weight += 0.5*cosang/sinang;
    }
    return weight;
}

void Mesh::buildExtrinsicDirichletLaplacian(const VectorXd &q, Eigen::SparseMatrix<double> &L) const
{
    int numverts = mesh_->n_vertices();
    L.resize(numverts, numverts);
    L.setZero();

    vector<Tr> Lcoeffs;

    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        if(mesh_->is_boundary(vi.handle()))
            continue;

        double totweight = 0;
        OMMesh::VertexHandle cent = vi.handle();
        for(OMMesh::VertexOHalfedgeIter voh = mesh_->voh_iter(cent); voh; ++voh)
        {
            OMMesh::EdgeHandle eh = mesh_->edge_handle(voh.handle());
            OMMesh::VertexHandle to = mesh_->to_vertex_handle(voh.handle());
            double weight = cotanWeight(eh.idx(), q);
            Lcoeffs.push_back(Tr(cent.idx(), to.idx(), weight));
            totweight += weight;
        }
        Lcoeffs.push_back(Tr(cent.idx(), cent.idx(), -totweight));
    }
    L.setFromTriplets(Lcoeffs.begin(), Lcoeffs.end());
}
