#include "mesh.h"
#include <iomanip>
#include "controller.h"
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;
using namespace OpenMesh;

const double PI = 3.1415926535898;

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

double Mesh::truncatedConeVolume(double startHeight, double curHeight)
{
    double base = params_.scale*params_.scale*PI;
    double totV = base*params_.scale*startHeight/3.0;
    double topr = params_.scale*(1.0 - curHeight/startHeight);
    double topbase = PI*topr*topr;
    double topV = topbase*params_.scale*(startHeight-curHeight)/3.0;
    return totV-topV;
}

bool Mesh::crush(Controller &cont, double coneHeight, double endHeight)
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

    VectorXd startq = q;

    const int crushTime = 100000;

    //double initialV = truncatedConeVolume(coneHeight, coneHeight);
    double airpressure = 30000;//101325.0;

    for(int i=0; i<numsteps; i++)
    {        
        std::cout << "iter " << i << std::endl;

        double t = i;
        if(t > crushTime)
            t = crushTime;
        double planeZ = (t*endHeight + (crushTime-t)*coneHeight)/double(crushTime);

        q += h*v;
        VectorXd gradq;
        std::cout << "computing energy" << std::endl;
        Midedge::elasticEnergy(*mesh_,q,g,params_,&gradq);
        std::cout << "force magnitude: " << gradq.norm() << std::endl;
        SparseMatrix<double> Minv;
        buildInvMassMatrix(g, Minv);
        VectorXd F = -gradq;

        //double curV = truncatedConeVolume(coneHeight, planeZ);
        //double pressure = airpressure*(initialV/curV - 1.0);
        //const double leakconst = 1e-6;
        //initialV = std::max(curV, initialV - leakconst*h*pressure);
        VectorXd pressureF;
        pressureForce(q,airpressure, pressureF);

        //std::cout << "Regular force " << F.norm() << " pressure force " << pressureF.norm() << " pressure " << pressure << " initialV " << initialV << " curV " << curV << std::endl;

        F += pressureF;

        v += h*Minv*F - h*Minv*params_.dampingCoeff*v;
        dofsToGeometry(q);
        if(i%100 == 0)
            dumpFrame();
        // enforce constraints
        enforceConstraints(q, startq, planeZ);
        std::cout << "done" << std::endl;
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
        double mass = area*params_.rho*params_.h*params_.scale*params_.scale*params_.scale;
//        if(mesh_->is_boundary(vi.handle()))
//            mass = std::numeric_limits<double>::infinity();
        for(int i=0; i<3; i++)
            entries.push_back(Tr(3*vidx+i, 3*vidx+i, mass));
    }

    M.setFromTriplets(entries.begin(), entries.end());
}

void Mesh::buildGeometricMassMatrix(const VectorXd &g, Eigen::SparseMatrix<double> &M) const
{
    M.resize(mesh_->n_vertices(), mesh_->n_vertices());
    vector<Tr> entries;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        int vidx = vi.handle().idx();
        double area = barycentricDualArea(g, vidx);
        double mass = area;
//        if(mesh_->is_boundary(vi.handle()))
//            mass = std::numeric_limits<double>::infinity();
        entries.push_back(Tr(vidx, vidx, mass));
    }

    M.setFromTriplets(entries.begin(), entries.end());
}

void Mesh::buildInvMassMatrix(const VectorXd &g, Eigen::SparseMatrix<double> &Minv) const
{
    Minv.resize(numdofs(), numdofs());
    vector<Tr> entries;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        int vidx = vi.handle().idx();
        double area = barycentricDualArea(g, vidx);
        double invmass = 1.0/area/params_.rho/params_.h;
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
        result += Midedge::intrinsicArea(vfi.handle().idx(), g, params_);
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

        Hdensity[i] = sign * Hnormal.norm();
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
    OMMesh::FaceHandle fh = mesh_->face_handle(fidx);
    double gs[3];
    int idx=0;
    for(OMMesh::FaceEdgeIter fei = mesh_->fe_iter(fh); fei; ++fei)
    {
        gs[idx++] = g[fei.handle().idx()];
    }

    double s = 0.5*(gs[0]+gs[1]+gs[2]);
    return sqrt(s*(s-gs[0])*(s-gs[1])*(s-gs[2]));
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

double Mesh::sampleHeight(const Vector2d pos, const VectorXd &q)
{
    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        int verts[3];
        int idx=0;
        for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fi.handle()); fvi; ++fvi)
        {
            verts[idx++] = fvi.handle().idx();
        }
        assert(idx==3);

        Vector2d e1 = q.segment<2>(3*verts[1])-q.segment<2>(3*verts[0]);
        Vector2d e2 = q.segment<2>(3*verts[2])-q.segment<2>(3*verts[0]);

        Matrix2d M;
        M << e1[0], e2[0], e1[1], e2[1];
        Vector2d rhs = pos - q.segment<2>(3*verts[0]);
        Matrix2d Minv = M.inverse();
        Vector2d barycoords = Minv*rhs;
        if(barycoords[0] > 0 && barycoords[1] > 0 && barycoords[0]+barycoords[1] <= 1.0)
        {
            return barycoords[0]*q[3*verts[1]+2] + barycoords[1]*q[3*verts[2]+2] + (1-barycoords[0]-barycoords[1])*q[3*verts[0]+2];
        }
    }

    double closestdist = std::numeric_limits<double>::infinity();
    double result = std::numeric_limits<double>::infinity();
    for(int i=0; i<(int)mesh_->n_vertices(); i++)
    {
        Vector2d vpos = q.segment<2>(3*i);
        double dist = (vpos-pos).norm();
        if(dist < closestdist)
        {
            closestdist = dist;
            result = q[3*i+2];
        }
    }

    return result;
}

void Mesh::symmetrize(int nfold)
{
    VectorXd q, g;
    dofsFromGeometry(q, g);
    VectorXd newq = q;

    double angle = 2*PI/nfold;
    Matrix2d rot;
    rot << cos(angle), sin(angle), -sin(angle), cos(angle);

    for(int i=0; i<(int)mesh_->n_vertices(); i++)
    {
        Vector2d pos = q.segment<2>(3*i);
        double totz = q[3*i+2];
        int numz = 1;
        for(int j=0; j<nfold; j++)
        {
            pos = rot*pos;
            double newz = sampleHeight(pos, q);
            if(newz != numeric_limits<double>::infinity())
            {
                totz += newz;
                numz++;
            }
        }
        newq[3*i+2] = totz/numz;
    }

    dofsToGeometry(newq);
}

void Mesh::enforceConstraints(VectorXd &q, const VectorXd &startq, double planeHeight)
{
    int numverts = mesh_->n_vertices();

    for(int i=0; i<numverts; i++)
    {
        OMMesh::VertexHandle vh = mesh_->vertex_handle(i);
        // pin boundary
        if(mesh_->is_boundary(vh))
        {
            Vector3d startpos = startq.segment<3>(3*i);
            q.segment<3>(3*i) = startpos;
//            double rsq = startpos[0]*startpos[0]+startpos[1]*startpos[1];
//            if(rsq > 0.9)
//                q[3*i+2] = 0;
        }
        // and crush
        double z = q[3*i+2];
        if(z > planeHeight)
        {
            q[3*i+2] = planeHeight;
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
