#include "shellforces.h"
#include "mesh.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

ShellForces::ShellForces(OMMesh &mesh, const VectorXd &undefq, double stretchingStiffness, double bendingStiffness)
    : mesh_(mesh), stretchingStiffness_(stretchingStiffness), bendingStiffness_(bendingStiffness)
{
    undeformedLengths_.resize(mesh_.n_faces());
    precomputedMatrices_.resize(mesh_.n_faces());

    for(OMMesh::FaceIter fi = mesh_.faces_begin(); fi != mesh_.faces_end(); ++fi)
    {
        OMMesh::FaceVertexIter fvi = mesh_.fv_iter(fi.handle());
        int v0 = fvi.handle().idx(); ++fvi;
        int v1 = fvi.handle().idx(); ++fvi;
        int v2 = fvi.handle().idx();

        int faceidx = fi.handle().idx();

        undeformedLengths_[faceidx] = Vector3d( (undefq.segment<3>(3*v1) - undefq.segment<3>(3*v0)).squaredNorm(),
                                          (undefq.segment<3>(3*v2) - undefq.segment<3>(3*v1)).squaredNorm(),
                                          (undefq.segment<3>(3*v0) - undefq.segment<3>(3*v2)).squaredNorm());

        Matrix3d Tm;
        precomputeMatrix(undefq, v0, v1, v2, Tm);
        precomputedMatrices_[faceidx] = Tm;
    }

    f0_.resize(mesh_.n_edges());
    e0_.resize(mesh_.n_edges());
    h0_.resize(mesh_.n_edges());

    for(OMMesh::EdgeIter ei = mesh_.edges_begin(); ei != mesh_.edges_end(); ++ei)
    {
        if(mesh_.is_boundary(ei.handle()))
            continue;
        OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(ei.handle(), 0);
        int p1 = mesh_.from_vertex_handle(heh).idx();
        int p2 = mesh_.to_vertex_handle(heh).idx();
        int q1 = mesh_.opposite_vh(heh).idx();
        int q2 = mesh_.opposite_vh(mesh_.opposite_halfedge_handle(heh)).idx();

        int eidx = ei.handle().idx();

        double f0 = computeAngle(undefq, p1, p2, q1, q2);
        double e0 = (undefq.segment<3>(3*p1) - undefq.segment<3>(3*p2)).norm();
        double A1 = 0.5 * (undefq.segment<3>(3*p2) - undefq.segment<3>(3*p1)).cross(undefq.segment<3>(3*q1) - undefq.segment<3>(3*p1)).norm();
        double A2 = 0.5 * (undefq.segment<3>(3*p2) - undefq.segment<3>(3*p1)).cross(undefq.segment<3>(3*q2) - undefq.segment<3>(3*p1)).norm();
        double h0 = 2.0*(A1+A2)/(3.0 * e0);

        f0_[eidx] = f0;
        e0_[eidx] = e0;
        h0_[eidx] = h0;
    }
}

void ShellForces::precomputeMatrix(const Eigen::VectorXd &undefq, int v0, int v1, int v2, Matrix3d &Tm)
{
    Vector3d uv[3];
    uv[0] = undefq.segment<3>(3*v2) - undefq.segment<3>(3*v1);
    uv[1] = undefq.segment<3>(3*v0) - undefq.segment<3>(3*v2);
    uv[2] = undefq.segment<3>(3*v1) - undefq.segment<3>(3*v0);

    double A = 0.5 * (uv[2].cross(-uv[1])).norm();

    int i, j, k, l, m, n;
    // T_{il1} = (1/8/A^4)*(dot(v_k,v_m)*dot(v_j, v_n)+dot(v_k,v_n)*dot(v_j,v_m))
    for (i = 0; i < 3; i++) {
        j = (i + 1) % 3;
        k = (i + 2) % 3;

        for (l = i; l < 3; l++) {
            m = (l + 1) % 3;
            n = (l + 2) % 3;

            Tm(i,l) = Tm(l,i) = (uv[k].dot(uv[m]) * uv[j].dot(uv[n]) +
                                   uv[k].dot(uv[n]) * uv[j].dot(uv[m]));
        }
    }

    // Scale
    //
    for (i=0; i<3; ++i)
    {
        for (j=0; j<3; ++j)
            Tm(i,j) = 1.0 / (256.0 * (A * A * A)) * Tm(i,j);
    }
}

void ShellForces::getStretchingForce(const VectorXd &q, VectorXd &f)
{
    f.resize(q.size());
    f.setZero();
    for(OMMesh::FaceIter fi = mesh_.faces_begin(); fi != mesh_.faces_end(); ++fi)
    {
        OMMesh::FaceVertexIter fvItr = mesh_.fv_iter(fi.handle());
        int v0 = fvItr.handle().idx(); ++fvItr;
        int v1 = fvItr.handle().idx(); ++fvItr;
        int v2 = fvItr.handle().idx();

        Vector3d force[3];
        int trinum = fi.handle().idx();
        getStencilStretchingForces(q, force, v0, v1, v2, precomputedMatrices_[trinum], trinum);
        f.segment<3>(3*v0) += stretchingStiffness_*force[0];
        f.segment<3>(3*v1) += stretchingStiffness_*force[1];
        f.segment<3>(3*v2) += stretchingStiffness_*force[2];
    }
}

void ShellForces::getStencilStretchingForces(const VectorXd &q, Vector3d *force, int v0, int v1, int v2, Matrix3d &Tm, int trinum)
{
    Vector3d v[3] = { q.segment<3>(3*v2) - q.segment<3>(3*v1),
                      q.segment<3>(3*v0) - q.segment<3>(3*v2),
                      q.segment<3>(3*v1) - q.segment<3>(3*v0) };
    double uv[3];

    uv[0] = undeformedLengths_[trinum][1];
    uv[1] = undeformedLengths_[trinum][2];
    uv[2] = undeformedLengths_[trinum][0];

    // Difference of squared lengths of edges
    //
    double  s[3] = { v[0].squaredNorm() - uv[0],
                   v[1].squaredNorm() - uv[1],
                   v[2].squaredNorm() - uv[2] };

    for(int i=0; i<3; i++)
        force[i].setZero();

    // taking derivative of the energy with respect to degree of freedom pkl
    // (l-th component of vertex k) gives the formula below for the force on
    // degree of freedom pkl
    //
    for (uint k = 0; k < 3; k++) {
        uint k2 = (k + 3 - 2) % 3;
        uint k1 = (k + 3 - 1) % 3;

        for (uint l = 0; l < 3; l++) {
            for (uint i = 0; i < 3; i++) {
                force[k][l] -= 4  * (Tm(i,k2) * s[i] * v[k2][l] -
                                     Tm(i,k1) * s[i] * v[k1][l]);
            }
        }
    }
}



void ShellForces::getBendingForce(const VectorXd &q, VectorXd &f)
{
    f.resize(q.size());
    f.setZero();
    for(OMMesh::EdgeIter ei = mesh_.edges_begin(); ei != mesh_.edges_end(); ++ei)
    {
        if(mesh_.is_boundary(ei.handle()))
            continue;
        OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(ei.handle(), 0);
        int p1 = mesh_.from_vertex_handle(heh).idx();
        int p2 = mesh_.to_vertex_handle(heh).idx();
        int q1 = mesh_.opposite_vh(heh).idx();
        int q2 = mesh_.opposite_vh(mesh_.opposite_halfedge_handle(heh)).idx();

        int eidx = ei.handle().idx();

        double theta = computeAngle(q, p1, p2, q1, q2);
        double theta0 = f0_[eidx];

        Vector3d dp1, dp2, dq1, dq2;
        ComputeDihedralAngleDerivatives(q, p1, p2, q1, q2, dp1, dp2, dq1, dq2);

        Vector3d force[4];
        force[0] = -dp1 * e0_[eidx] / h0_[eidx] * (theta-theta0);
        force[1] = -dp2 * e0_[eidx] / h0_[eidx] * (theta-theta0);
        force[2] = -dq1 * e0_[eidx] / h0_[eidx] * (theta-theta0);
        force[3] = -dq2 * e0_[eidx] / h0_[eidx] * (theta-theta0);

        f.segment<3>(3*p1) += force[0] * bendingStiffness_;
        f.segment<3>(3*p2) += force[1] * bendingStiffness_;
        f.segment<3>(3*q1) += force[2] * bendingStiffness_;
        f.segment<3>(3*q2) += force[3] * bendingStiffness_;
    }    
}

double ShellForces::computeAngle(const VectorXd &q, int p1, int p2, int q1, int q2)
{
    Vector3d ev = q.segment<3>(3*p2) - q.segment<3>(3*p1);

    // compute normals of the two deformed triangles and the angle between them
    //
    Vector3d n1 = ev.cross(q.segment<3>(3*q1) - q.segment<3>(3*p1));
    Vector3d n2 = (q.segment<3>(3*q2) - q.segment<3>(3*p1)).cross(ev);

    double mag = ev.norm();
    ev /= mag;
    //XXXX
    assert(mag >= 1e-8);
    return atan2(n1.cross(n2).dot(ev), n1.dot(n2));
}

void ShellForces::ComputeDihedralAngleDerivatives(const VectorXd &q, int p1, int p2, int q1, int q2,
                                                        Vector3d &del_p1_theta,
                                                        Vector3d &del_p2_theta,
                                                        Vector3d &del_q1_theta,
                                                        Vector3d &del_q2_theta) const
{
    // Edge vectors
    //
    Vector3d v11 = q.segment<3>(3*q1) - q.segment<3>(3*p2);
    Vector3d v12 = q.segment<3>(3*q1) - q.segment<3>(3*p1);

    Vector3d v21 = q.segment<3>(3*q2) - q.segment<3>(3*p2);
    Vector3d v22 = q.segment<3>(3*q2) - q.segment<3>(3*p1);

    Vector3d v = q.segment<3>(3*p2) - q.segment<3>(3*p1);

    double vNorm2 = v.dot(v);
    double vNorm	= sqrt(vNorm2);

    // normal vectors
    Vector3d n1 = v12.cross(v);
    Vector3d n2 = v.cross(v22);

    double n1Norm2 = n1.dot(n1);
    double n2Norm2 = n2.dot(n2);

    if (vNorm < 1e-8 || n1Norm2 < 1e-8 || n2Norm2 < 1e-8)
    {
        del_p1_theta.setZero();
        del_p2_theta.setZero();
        del_q1_theta.setZero();
        del_q2_theta.setZero();
    }
    else
    {
        // gradients of theta
        //
        del_q1_theta = (n1 * vNorm) / n1Norm2;
        del_q2_theta = (n2 * vNorm) / n2Norm2;

        Vector3d F11 = n1 * (v11.dot(v) / (vNorm * n1Norm2));
        Vector3d F12 = n2 * (v21.dot(v) / (vNorm * n2Norm2));

        Vector3d F21 = n1 * (-v12.dot(v) / (vNorm * n1Norm2));
        Vector3d F22 = n2 * (-v22.dot(v) / (vNorm * n2Norm2));

        del_p1_theta = F11 + F12;
        del_p2_theta = F21 + F22;
    }
}

void ShellForces::getForce(const VectorXd &q, VectorXd &f)
{
    VectorXd bf, sf;
    getStretchingForce(q, sf);
    getBendingForce(q, bf);
    f = bf + sf;
}
