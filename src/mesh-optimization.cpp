#include "mesh.h"
#include <iomanip>
#include "controller.h"
#include <Eigen/Dense>
#include <fstream>
#include "midedge.h"
#include "robustleastsquares.h"
#include <Eigen/SPQRSupport>

using namespace std;
using namespace Eigen;
using namespace OpenMesh;

const double PI = 3.1415926535898;

typedef Eigen::Triplet<double> Tr;

void Mesh::elasticEnergy(const VectorXd &q,
                         const VectorXd &g1,
                         const VectorXd &g2,
                         double &energy,
                         VectorXd &gradq,
                         bool derivativesRequested)
{
    assert(q.size() == numdofs());
    assert(g1.size() == 4*mesh_->n_faces());
    assert(g2.size() == 4*mesh_->n_faces());
    VectorXd *derivs = NULL;
    if(derivativesRequested)
    {
        gradq.resize(q.size());
        derivs = &gradq;
    }
    energy = Midedge::elasticEnergy(*mesh_, q, g1, g2, params_, derivs);
}

void Mesh::calculateGrowth(const VectorXd &q,
                     VectorXd &g1,
                     VectorXd &g2)
{
    vector<vector<Tr> > mats;
    for(int i=0; i<4; i++)
    {
        mats.push_back(vector<Tr>());
    }
    Midedge::elasticEnergyFactor(*mesh_, q, g1, g2, params_, mats);

    int numfaces = mesh_->n_faces();

    VectorXd scales(2*numfaces+q.size());
    scales.setZero();


    vector<Tr> Lcoeffs;
    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        int nbs = 0;
        for(OMMesh::FaceFaceIter ffi = mesh_->ff_iter(fi.handle()); ffi; ++ffi)
        {
            Lcoeffs.push_back(Tr(fi.handle().idx(), ffi.handle().idx(), 1.0));
            nbs++;
        }
        Lcoeffs.push_back(Tr(fi.handle().idx(), fi.handle().idx(), -double(nbs)));
    }
    SparseMatrix<double> L(numfaces, numfaces);
    L.setFromTriplets(Lcoeffs.begin(), Lcoeffs.end());

    vector<Tr> Mcoeffs;

    for(int i=0; i<numfaces; i++)
    {
        for(int j=0; j<numfaces; j++)
        {
            if(L.coeff(i, j) != 0)
            {
                Mcoeffs.push_back(Tr(i, j, L.coeff(i,j)));
                Mcoeffs.push_back(Tr(i+numfaces, j+numfaces, L.coeff(i,j)));
            }
        }
    }
    double reg = 1e-5;
    for(int i=0; i<q.size(); i++)
    {
        Mcoeffs.push_back(Tr(2*numfaces+i, 2*numfaces+i, reg));
    }


    vector<Tr> Ncoeffs;

    for(int i=0; i<(int)mats[0].size(); i++)
    {
        {
            Mcoeffs.push_back(Tr(mats[0][i].row() + 2*numfaces, mats[0][i].col(), mats[0][i].value()));
            Mcoeffs.push_back(Tr(mats[0][i].col(), mats[0][i].row() + 2*numfaces, mats[0][i].value()));
            Ncoeffs.push_back(mats[0][i]);
        }
    }
    for(int i=0; i<(int)mats[1].size(); i++)
    {
        {
            Mcoeffs.push_back(Tr(mats[1][i].row() + 2*numfaces, numfaces+mats[1][i].col(), mats[1][i].value()));
            Mcoeffs.push_back(Tr(mats[1][i].col() + numfaces, 2*numfaces+mats[1][i].row(), mats[1][i].value()));
            Ncoeffs.push_back(Tr(mats[1][i].row(), numfaces+mats[1][i].col(), mats[1][i].value()));
        }
    }
    SparseMatrix<double> M(2*numfaces+q.size(), 2*numfaces+q.size());
    M.setFromTriplets(Mcoeffs.begin(), Mcoeffs.end());

    SparseMatrix<double> N(q.size(), 2*numfaces);
    N.setFromTriplets(Ncoeffs.begin(), Ncoeffs.end());

    VectorXd rhs(2*numfaces+q.size());
    VectorXd nrhs(q.size());
    rhs.setZero();
    nrhs.setZero();
    for(int i=2; i<4; i++)
    {
        for(int j=0; j<(int)mats[i].size(); j++)
        {
            {
                rhs[mats[i][j].row()+2*numfaces] -= mats[i][j].value();
                nrhs[mats[i][j].row()] -= mats[i][j].value();
            }
        }
    }

    double oldenergy;
    VectorXd gradq;
    elasticEnergy(q, g1, g2, oldenergy, gradq, true);
    cout << "old residual " << gradq.norm() << endl;

    //VectorXd lsrhs = rhs - M*scales;
    //RobustLeastSquares::solve(M, lsrhs, scales, true);
    RobustLeastSquares::solve(M, rhs, scales, true);

    SPQR<SparseMatrix<double> > ssqr(N);
    cout << ssqr.rank() << " / " << 2*numfaces << endl;

    std::cout << "scales norm " << scales.norm() << std::endl;
    std::cout << "minimum scale ";
    double minval = std::numeric_limits<double>::infinity();
    for(int i=0; i<2*numfaces; i++)
    {
        minval = std::min(minval, scales[i]);
        if(scales[i] < 0.1)
            scales[i] = 0.1;
    }
    cout << minval << endl;

    std::cout << (M*scales-rhs).norm() << std::endl;
    std::cout << (N*scales.segment(0,2*numfaces) - nrhs).norm() << std::endl;

    colors1_.resize(numfaces);
    colors2_.resize(numfaces);
    for(int i=0; i<numfaces; i++)
    {
        colors1_[i] = -1.0+scales[i];
        colors2_[i] = -1.0+scales[i+numfaces];
    }

    for(int i=0; i<numfaces; i++)
        for(int k=0; k<4; k++)
        {
            g1[4*i+k] /= (scales[numfaces+i]);
            g2[4*i+k] /= (scales[i]);
        }

    VectorXd newq;
    dofsFromGeometry(newq);
    double newenergy;
    VectorXd newgradq(newq.size());
    newgradq.setZero();

    newenergy = Midedge::elasticEnergy(*mesh_, newq, g1, g2, params_, &newgradq);
    std::cout << "new residual " << newgradq.norm() << std::endl;
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

bool Mesh::relaxConfiguration(Controller &cont)
{
    VectorXd q(numdofs());
    dofsFromGeometry(q);

    double h = params_.eulerTimestep;
    int numsteps = 50000;
    SparseMatrix<double> Minv;
    buildInvMassMatrix(g1_, g2_, Minv);
    VectorXd v(numdofs());
    v.setZero();

    for(int iter=0; iter<numsteps; iter++)
    {
        std::cout << "iter " << iter;
        VectorXd gradq;
        q += h*v;
        double energy;
        elasticEnergy(q, g1_, g2_, energy, gradq, true);

        std::cout << " energy " << energy << " force magnitude: " << gradq.norm() << " impulse magnitude: " << (h*Minv*gradq).norm() << std::endl;

        v = -h*Minv*gradq;
        dofsToGeometry(q);
    }

    dofsToGeometry(q);

    cont.updateGL();
    return true;
}

bool Mesh::simulateGrowth(Controller &cont)
{
    {
        string name = params_.outputDir + "/parameters";
        ofstream paramfile(name.c_str());
        params_.dumpParameters(paramfile);
    }

    VectorXd q(numdofs());
    dofsFromGeometry(q);

    double h = params_.eulerTimestep;
    VectorXd v(numdofs());
    v.setZero();
    int numsteps = params_.numEulerIters;   

    VectorXd targetg1 = g1_;
    VectorXd targetg2 = g2_;

    setInducedMetric();
    VectorXd origg1 = g1_;
    VectorXd origg2 = g2_;

    int growthsteps = numsteps/10;

    for(int iter=0; iter<numsteps; iter++)
    {
        double t = double(iter);
        if(t > growthsteps)
            t = double(growthsteps);

        t /= double(growthsteps);

        g1_ = (1.0-t)*origg1 + t*targetg1;
        g2_ = (1.0-t)*origg2 + t*targetg2;

        std::cout << "iter " << iter;
        q += h*v;
        VectorXd gradq;
        double energy;
        elasticEnergy(q, g1_, g2_, energy, gradq, true);

        SparseMatrix<double> Minv;
        buildInvMassMatrix(g1_, g2_, Minv);
        std::cout << " energy " << energy << " force magnitude: " << gradq.norm() << " impulse magnitude: " << (h*Minv*gradq).norm() << std::endl;

        v -= h*Minv*gradq + h*Minv*params_.dampingCoeff*v;
        dofsToGeometry(q);
        if(iter%100 == 0)
            dumpFrame();
    }

    dofsToGeometry(q);

    cont.updateGL();
    return true;
}

void Mesh::buildMassMatrix(const VectorXd &g, Eigen::SparseMatrix<double> &M) const
{
    M.resize(numdofs(), numdofs());
    vector<Tr> entries;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        int vidx = vi.handle().idx();
        double area = barycentricDualArea(g, vidx);
        double mass = area*params_.rho*params_.h;
//        if(mesh_->is_boundary(vi.handle()))
//            mass = std::numeric_limits<double>::infinity();
        for(int i=0; i<3; i++)
            entries.push_back(Tr(3*vidx+i, 3*vidx+i, mass));
    }

    M.setFromTriplets(entries.begin(), entries.end());
}

void Mesh::buildInvMassMatrix(const VectorXd &g1, const VectorXd &g2, Eigen::SparseMatrix<double> &Minv) const
{
    Minv.resize(numdofs(), numdofs());
    vector<Tr> entries;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        int vidx = vi.handle().idx();
        double area = barycentricDualArea(g1, vidx);
        area += barycentricDualArea(g2, vidx);
        double invmass = 0.5/area/params_.rho/params_.h/params_.scale/params_.scale;
//        if(mesh_->is_boundary(vi.handle()))
//            invmass = 0;
        for(int i=0; i<3; i++)
            entries.push_back(Tr(3*vidx+i, 3*vidx+i, invmass));
    }

    Minv.setFromTriplets(entries.begin(), entries.end());
}

double Mesh::barycentricDualArea(const VectorXd &g, int vidx) const
{
    double result = 0;
    OMMesh::VertexHandle vh = mesh_->vertex_handle(vidx);
    for(OMMesh::VertexFaceIter vfi = mesh_->vf_iter(vh); vfi; ++vfi)
    {
        result += restFaceArea(g, vfi.handle().idx());
    }
    return result/3.0;
}

double Mesh::deformedBarycentricDualArea(const VectorXd &q, int vidx) const
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

void Mesh::gaussianCurvature(const VectorXd &q, VectorXd &K) const
{
    int numverts = mesh_->n_vertices();
    K.resize(numverts);
    K.setZero();

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

        K[i] = 2*PI-angsum;
        double area =  deformedBarycentricDualArea(q, i);
        K[i] /= area;
    }
}

void Mesh::vertexAreas(const VectorXd &q, VectorXd &vareas) const
{
    int numverts = mesh_->n_vertices();
    vareas.resize(numverts);
    vareas.setZero();
    for(int i=0; i<numverts; i++)
    {
        double area = deformedBarycentricDualArea(q, i);
        vareas[i] = area;
    }
}

void Mesh::meanCurvature(const VectorXd &q, VectorXd &Hdensity) const
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
        Vector3d Hnormal(Hx[i], Hy[i], Hz[i]);
        double sign = 1.0;

        Vector3d avnormal = averageNormal(q, i);
        if(avnormal.dot(Hnormal) < 0)
            sign = -1.0;

        Hdensity[i] = 0.5*sign * Hnormal.norm();
    }
}

double Mesh::intrinsicCotanWeight(int edgeid, const VectorXd &g) const
{
    OMMesh::EdgeHandle eh = mesh_->edge_handle(edgeid);
    double c = g[eh.idx()];
    double weight = 0;
    for(int i=0; i<2; i++)
    {
        OMMesh::HalfedgeHandle heh = mesh_->halfedge_handle(eh,i);

        if(mesh_->is_boundary(heh))
            continue;
        OMMesh::HalfedgeHandle next = mesh_->next_halfedge_handle(heh);
        OMMesh::EdgeHandle edge1 = mesh_->edge_handle(next);
        double a = g[edge1.idx()];

        OMMesh::HalfedgeHandle prev = mesh_->prev_halfedge_handle(heh);
        OMMesh::EdgeHandle edge2 = mesh_->edge_handle(prev);
        double b = g[edge2.idx()];

        double denom = 4.0*a*a*b*b-(a*a+b*b-c*c)*(a*a+b*b-c*c);
        if(denom < 0)
            denom = 0;
        weight += 0.5*(a*a+b*b-c*c)/sqrt(denom);
    }
    return weight;
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

void Mesh::buildIntrinsicDirichletLaplacian(const VectorXd &g, Eigen::SparseMatrix<double> &L) const
{
    int numverts = mesh_->n_vertices();
    L.resize(numverts, numverts);

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
            double weight = intrinsicCotanWeight(eh.idx(), g);
            Lcoeffs.push_back(Tr(cent.idx(), to.idx(), weight));
            totweight += weight;
        }
        Lcoeffs.push_back(Tr(cent.idx(), cent.idx(), -totweight));
    }
    L.setFromTriplets(Lcoeffs.begin(), Lcoeffs.end());
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

double Mesh::restFaceArea(const VectorXd &g, int fidx) const
{
    return Midedge::intrinsicArea(fidx, g, params_);
}

double Mesh::vertexAreaRatio(const VectorXd &undefq, const VectorXd &g, int vidx)
{
    double restarea = 0;
    double curarea = 0;
    OMMesh::VertexHandle vh = mesh_->vertex_handle(vidx);
    for(OMMesh::VertexFaceIter vfi = mesh_->vf_iter(vh); vfi; ++vfi)
    {
        int fidx = vfi.handle().idx();
        restarea += faceArea(undefq, fidx);
        curarea += restFaceArea(g, fidx);
    }

    return curarea/restarea;
}

void Mesh::enforceConstraints(VectorXd &q, const VectorXd &startq, double planeHeight)
{
    int numverts = mesh_->n_vertices();

    for(int i=0; i<numverts; i++)
    {
        double z = q[3*i+2];
        if(z > planeHeight)
            q[3*i+2] = planeHeight;
    }
    // pin boundary
    for(int i=0; i<numverts; i++)
    {
        OMMesh::VertexHandle vh = mesh_->vertex_handle(i);
        if(mesh_->is_boundary(vh))
        {
            Vector3d startpos = startq.segment<3>(3*i);
            double rsq = startpos[0]*startpos[0]+startpos[1]*startpos[1];
            if(rsq > 0.9)
                q[3*i+2] = 0;
        }
    }
}

void Mesh::pressureForce(const VectorXd &q, double pressure, VectorXd &F)
{
    F.resize(q.size());
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        Vector3d normal = this->surfaceAreaNormal(q, vi.handle().idx());
        double area = params_.scale*params_.scale*this->deformedBarycentricDualArea(q, vi.handle().idx());
        for(int i=0; i<3; i++)
            F[3*vi.handle().idx()+i] = pressure*normal[i]*area;
    }
}
