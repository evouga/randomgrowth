#include "midedge.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

Midedge::Midedge()
{

}

Eigen::Matrix3d Midedge::crossMatrix(const Vector3d &v)
{
    Matrix3d result;
    result << 0, -v[2], v[1],
            v[2], 0, -v[0],
            -v[1], v[0], 0;
    return result;
}

double Midedge::area(const Vector3d &q0, const Vector3d &q1, const Vector3d &q2)
{
    return 0.5*(q0-q1).cross(q2-q1).norm();
}

void Midedge::Darea(const Vector3d &q0, const Vector3d &q1, const Vector3d &q2, vector<Vector3d> &partials)
{
    Vector3d inner = (q0-q1).cross(q2-q1);
    double denom = 2.0*inner.norm();
    partials.clear();
    partials.push_back(-inner.transpose()*crossMatrix(q2-q1)/denom);
    partials.push_back(inner.transpose()*(-crossMatrix(q0-q1) + crossMatrix(q2-q1))/denom);
    partials.push_back( inner.transpose()*crossMatrix(q0-q1)/denom);
}

Vector3d Midedge::diamondEdgeNormal(const Vector3d &q0, const Vector3d &q1, const Vector3d &q2, const Vector3d &q3)
{
    double a1 = area(q0, q2, q3);
    double a2 = area(q0, q1, q2);
    Vector3d result = a1 * (q1-q0).cross(q2-q0) + a2*(q2-q0).cross(q3-q0);
    return result/result.norm();
}

void Midedge::DdiamondEdgeNormal(const Vector3d &q0, const Vector3d &q1, const Vector3d &q2, const Vector3d &q3, vector<Matrix3d> &partials)
{
    double a1 = area(q0, q2, q3);
    vector<Vector3d> da1;
    Darea(q0, q2, q3, da1);

    double a2 = area(q0, q1, q2);
    vector<Vector3d> da2;
    Darea(q0, q1, q2, da2);

    Vector3d num = a1*(q1-q0).cross(q2-q0) + a2*(q2-q0).cross(q3-q0);
    double denom = num.norm();
    partials.clear();
    Matrix3d dn0 = (q1-q0).cross(q2-q0) * da1[0].transpose()
                        + a1 * crossMatrix(q2-q0)
                        - a1 * crossMatrix(q1-q0)
                        + (q2-q0).cross(q3-q0) * da2[0].transpose()
                        + a2 * crossMatrix(q3-q0)
                        - a2 * crossMatrix(q2-q0);
    partials.push_back(dn0/denom - num*num.transpose()*dn0/denom/denom/denom);
    Matrix3d dn1 = -a1 * crossMatrix(q2-q0)
            + (q2-q0).cross(q3-q0) * da2[1].transpose();
    partials.push_back(dn1/denom - num*num.transpose()*dn1/denom/denom/denom);
    Matrix3d dn2 = (q1-q0).cross(q2-q0) * da1[1].transpose()
            + a1*crossMatrix(q1-q0)
            + (q2-q0).cross(q3-q0) * da2[2].transpose()
            - a2*crossMatrix(q3-q0);
    partials.push_back(dn2/denom - num*num.transpose()*dn2/denom/denom/denom);
    Matrix3d dn3 = (q1-q0).cross(q2-q0) * da1[2].transpose()
            + a2*crossMatrix(q2-q0);
    partials.push_back(dn3/denom - num*num.transpose()*dn3/denom/denom/denom);
}

Vector3d Midedge::g(const Vector3d &q1, const Vector3d &q2, const Vector3d &q3)
{
    Vector3d result;
    result[0] = (q3-q2).dot(q3-q2);
    result[1] = (q3-q2).dot(q1-q3);
    result[2] = (q1-q3).dot(q1-q3);
    return result;
}

void Midedge::Dg(const Vector3d &q1, const Vector3d &q2, const Vector3d &q3, std::vector<Eigen::Matrix3d> &partials)
{
    partials.clear();

    Matrix3d dg1;
    dg1 << 0, 0, 0,
            (q3-q2).transpose(),
            2.0*(q1-q3).transpose();
    partials.push_back(dg1);

    Matrix3d dg2;
    dg2 << -2.0*(q3-q2).transpose(),
            -(q1-q3).transpose(),
            0,0,0;
    partials.push_back(dg2);

    Matrix3d dg3;
    dg3 << 2.0*(q3-q2).transpose(),
            (q1-q3).transpose() - (q3-q2).transpose(),
            -2.0*(q1-q3).transpose();
    partials.push_back(dg3);
}

Vector3d Midedge::b(const Vector3d &q1, const Vector3d &q2, const Vector3d &q3, const Vector3d &n1, const Vector3d &n2, const Vector3d &n3)
{
    Vector3d result;
    result[0] = 2.0*(n3-n2).dot(q3-q2);
    result[1] = 2.0*(n1-n3).dot(q3-q2);
    result[2] = 2.0*(n1-n3).dot(q1-q3);
    return result;
}

void Midedge::Db(const Vector3d &q1, const Vector3d &q2, const Vector3d &q3, const Vector3d &n1, const Vector3d &n2, const Vector3d &n3, std::vector<Matrix3d> &partials)
{
    partials.clear();

    Matrix3d db1;
    db1 << 0, 0, 0,
            0,0,0,
            2.0*(n1-n3).transpose();
    partials.push_back(db1);

    Matrix3d db2;
    db2 << -2.0*(n3-n2).transpose(),
            -2.0*(n1-n3).transpose(),
            0, 0, 0;
    partials.push_back(db2);

    Matrix3d db3;
    db3 << 2.0*(n3-n2).transpose(),
            2.0*(n1-n3).transpose(),
            -2.0*(n1-n3).transpose();
    partials.push_back(db3);

    Matrix3d db4;
    db4 << 0,0,0,
            2.0*(q3-q2).transpose(),
            2.0*(q1-q3).transpose();
    partials.push_back(db4);

    Matrix3d db5;
    db5 << -2.0*(q3-q2).transpose(),
            0, 0, 0,
            0, 0, 0;
    partials.push_back(db5);

    Matrix3d db6;
    db6 << 2.0*(q3-q2).transpose(),
            -2.0*(q3-q2).transpose(),
            -2.0*(q1-q3).transpose();
    partials.push_back(db6);
}

Vector3d Midedge::faceNormal(const Vector3d &q0, const Vector3d &q1, const Vector3d &q2)
{
    Vector3d result = (q1-q0).cross(q2-q0);
    return result/result.norm();
}

void Midedge::DfaceNormal(const Vector3d &q0, const Vector3d &q1, const Vector3d &q2, std::vector<Matrix3d> &partials)
{
    Vector3d n = (q1-q0).cross(q2-q0);
    double denom = n.norm();
    partials.clear();

    Matrix3d dn0 = crossMatrix(q2-q0) - crossMatrix(q1-q0);
    partials.push_back( dn0/denom - n*n.transpose()*dn0/denom/denom/denom );

    Matrix3d dn1 = -crossMatrix(q2-q0);
    partials.push_back( dn1/denom - n*n.transpose()*dn1/denom/denom/denom );

    Matrix3d dn2 = crossMatrix(q1-q0);
    partials.push_back( dn2/denom - n*n.transpose()*dn2/denom/denom/denom );
}

Vector3d Midedge::edgeNormal(const OMMesh &mesh, int edgeid, const VectorXd &q)
{
    OMMesh::EdgeHandle eh = mesh.edge_handle(edgeid);
    OMMesh::HalfedgeHandle heh0 = mesh.halfedge_handle(eh,0);
    if(mesh.is_boundary(eh))
    {
        if(mesh.is_boundary(heh0))
            heh0 = mesh.opposite_halfedge_handle(heh0);
        Vector3d qs[3];
        qs[0] = q.segment<3>(3*mesh.from_vertex_handle(heh0).idx());
        qs[1] = q.segment<3>(3*mesh.to_vertex_handle(heh0).idx());
        qs[2] = q.segment<3>(3*mesh.to_vertex_handle(mesh.next_halfedge_handle(heh0)).idx());
        return faceNormal(qs[0], qs[1], qs[2]);
    }
    else
    {
        Vector3d qs[4];
        qs[0] = q.segment<3>(3*mesh.to_vertex_handle(heh0).idx());
        qs[1] = q.segment<3>(3*mesh.to_vertex_handle(mesh.next_halfedge_handle(heh0)).idx());
        qs[2] = q.segment<3>(3*mesh.from_vertex_handle(heh0).idx());
        qs[3] = q.segment<3>(3*mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(heh0))).idx());
        return diamondEdgeNormal(qs[0], qs[1], qs[2], qs[3]);
    }
}

void Midedge::DedgeNormal(const OMMesh &mesh, int edgeid, const VectorXd &q, const Vector3d &prefactor, VectorXd &partials)
{
    OMMesh::EdgeHandle eh = mesh.edge_handle(edgeid);
    OMMesh::HalfedgeHandle heh0 = mesh.halfedge_handle(eh,0);
    partials.resize(q.size());
    partials.setZero();
    if(mesh.is_boundary(eh))
    {
        if(mesh.is_boundary(heh0))
            heh0 = mesh.opposite_halfedge_handle(heh0);
        Vector3d qs[3];
        int idx[3];
        idx[0] = mesh.from_vertex_handle(heh0).idx();
        qs[0] = q.segment<3>(3*idx[0]);
        idx[1] = mesh.to_vertex_handle(heh0).idx();
        qs[1] = q.segment<3>(3*idx[1]);
        idx[2] = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh0)).idx();
        qs[2] = q.segment<3>(3*idx[2]);
        vector<Matrix3d> partials3d;
        DfaceNormal(qs[0], qs[1], qs[2], partials3d);
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    partials[3*idx[i]+k] += partials3d[i](j,k)*prefactor[j];
                }
            }
        }
    }
    else
    {
        Vector3d qs[4];
        int idx[4];
        idx[0] = mesh.to_vertex_handle(heh0).idx();
        qs[0] = q.segment<3>(3*idx[0]);
        idx[1] = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh0)).idx();
        qs[1] = q.segment<3>(3*idx[1]);
        idx[2] = mesh.from_vertex_handle(heh0).idx();
        qs[2] = q.segment<3>(3*idx[2]);
        idx[3] = mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(heh0))).idx();
        qs[3] = q.segment<3>(3*idx[3]);
        vector<Matrix3d> partials3d;
        DdiamondEdgeNormal(qs[0], qs[1], qs[2], qs[3], partials3d);
        for(int i=0; i<4; i++)
        {
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    partials[3*idx[i]+k] += partials3d[i](j,k)*prefactor[j];
                }
            }
        }
    }
}

Vector4d Midedge::g(const OMMesh &mesh, int faceid, const VectorXd &q)
{
    OMMesh::FaceHandle fh = mesh.face_handle(faceid);
    Vector3d qs[3];
    int idx[3];
    int curidx=0;
    for(OMMesh::ConstFaceHalfedgeIter fhi = mesh.cfh_iter(fh); fhi; ++fhi)
    {
        idx[curidx] = mesh.from_vertex_handle(fhi.handle()).idx();
        qs[curidx] = q.segment<3>(3*idx[curidx]);
        curidx++;
    }
    Vector3d g3 = g(qs[0], qs[1], qs[2]);
    Vector4d result;
    result[0] = g3[0];
    result[1] = result[2] = g3[1];
    result[3] = g3[2];
    return result;
}

void Midedge::Dg(const OMMesh &mesh, int faceid, const VectorXd &q, const Vector4d &prefactor, VectorXd &partials)
{
    OMMesh::FaceHandle fh = mesh.face_handle(faceid);
    Vector3d qs[3];
    int idx[3];
    int curidx=0;
    for(OMMesh::ConstFaceHalfedgeIter fhi = mesh.cfh_iter(fh); fhi; ++fhi)
    {
        idx[curidx] = mesh.from_vertex_handle(fhi.handle()).idx();
        qs[curidx] = q.segment<3>(3*idx[curidx]);
        curidx++;
    }
    vector<Matrix3d> partials3;
    Dg(qs[0], qs[1], qs[2], partials3);
    partials.resize(q.size());
    partials.setZero();

    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            partials[3*idx[i]+j] += partials3[i](0,j)*prefactor[0];
            partials[3*idx[i]+j] += partials3[i](1,j)*prefactor[1];
            partials[3*idx[i]+j] += partials3[i](1,j)*prefactor[2];
            partials[3*idx[i]+j] += partials3[i](2,j)*prefactor[3];
        }
    }
}

Vector4d Midedge::b(const OMMesh &mesh, int faceid, const VectorXd &q)
{
    Vector3d ns[3];
    int idx[3];
    Vector3d qs[3];
    int curidx = 0;
    OMMesh::FaceHandle fh = mesh.face_handle(faceid);
    for(OMMesh::ConstFaceHalfedgeIter fhi = mesh.cfh_iter(fh); fhi; ++fhi)
    {
        idx[curidx] = mesh.from_vertex_handle(fhi.handle()).idx();
        qs[curidx] = q.segment<3>(3*idx[curidx]);
        int edgeid = mesh.edge_handle(fhi.handle()).idx();
        ns[(2+curidx)%3] = edgeNormal(mesh, edgeid, q);
        curidx++;
    }
    Vector3d b3 = b(qs[0], qs[1], qs[2], ns[0], ns[1], ns[2]);
    Vector4d result(b3[0], b3[1], b3[1], b3[2]);
    return result;
}

void Midedge::Db(const OMMesh &mesh, int faceid, const VectorXd &q, const Vector4d &prefactor, VectorXd &partials)
{
    Vector3d ns[3];
    int idx[3];
    Vector3d qs[3];
    int curidx = 0;
    int edgeids[3];
    OMMesh::FaceHandle fh = mesh.face_handle(faceid);
    for(OMMesh::ConstFaceHalfedgeIter fhi = mesh.cfh_iter(fh); fhi; ++fhi)
    {
        idx[curidx] = mesh.from_vertex_handle(fhi.handle()).idx();
        qs[curidx] = q.segment<3>(3*idx[curidx]);
        edgeids[curidx] = mesh.edge_handle(fhi.handle()).idx();
        ns[(2+curidx)%3] = edgeNormal(mesh, edgeids[curidx], q);
        curidx++;
    }
    vector<Matrix3d> partials3;
    Db(qs[0], qs[1], qs[2], ns[0], ns[1], ns[2], partials3);
    Matrix<double, 4, 3> conv;
    conv << 1,0,0,
            0,1,0,
            0,1,0,
            0,0,1;
    partials.resize(q.size());
    partials.setZero();
    for(int i=0; i<3; i++)
    {
        Vector3d prod = (prefactor.transpose()*conv)*partials3[i];
        for(int j=0; j<3; j++)
            partials[3*idx[i]+j] += prod[j];
    }
    for(int i=0; i<3; i++)
    {
        Vector3d prod = (prefactor.transpose()*conv)*partials3[3+i];
        VectorXd npart;
        DedgeNormal(mesh, edgeids[(1+i)%3], q, prod, npart);
        partials += npart;
    }
}

Vector4d Midedge::matMult(const Vector4d &m1, const Vector4d &m2)
{
    Vector4d result;
    result << m1[0]*m2[0]+m1[1]*m2[2],
            m1[0]*m2[1]+m1[1]*m2[3],
            m1[2]*m2[0]+m1[3]*m2[2],
            m1[2]*m2[1]+m1[3]*m2[3];
    return result;
}

void Midedge::DmatMult(const Vector4d &m1, const Vector4d &m2, Matrix4d &m1partials, Matrix4d &m2partials)
{
    m1partials << m2[0], m2[2], 0, 0,
            m2[1], m2[3], 0, 0,
            0, 0, m2[0], m2[2],
            0, 0, m2[1], m2[3];

    m2partials << m1[0], 0, m1[1], 0,
            0, m1[0], 0, m1[1],
            m1[2], 0, m1[3], 0,
            0, m1[2], 0, m1[3];

    return;
}

double Midedge::trace(const Vector4d &m1)
{
    return m1[0] + m1[3];
}

void Midedge::Dtrace(Vector4d &partials)
{
    partials << 1, 0, 0, 1;
}

Vector4d Midedge::gbar(const OMMesh &mesh, int faceid, const VectorXd &g)
{
    double gi[3];
    int idx=0;
    OMMesh::FaceHandle fh = mesh.face_handle(faceid);
    for(OMMesh::ConstFaceHalfedgeIter fhi = mesh.cfh_iter(fh); fhi; ++fhi)
    {
        gi[idx++] = g[mesh.edge_handle(fhi.handle()).idx()];
    }
    Vector4d result;
    result << gi[1]*gi[1],
            -0.5*(gi[1]*gi[1]+gi[2]*gi[2]-gi[0]*gi[0]),
            -0.5*(gi[1]*gi[1]+gi[2]*gi[2]-gi[0]*gi[0]),
            gi[2]*gi[2];
    return result;
}

double Midedge::det(const Vector4d &m)
{
    return m[0]*m[3]-m[1]*m[2];
}

void Midedge::Ddet(const Vector4d &m, Vector4d &partials)
{
    partials << m[3], -m[2], -m[1], m[0];
}

Vector4d Midedge::matInv(const Vector4d &m)
{
    double detval = det(m);
    Vector4d result(m[3], -m[1], -m[2], m[0]);
    result /= detval;
    return result;
}

void Midedge::DmatInv(const Vector4d &m, Matrix4d &partials)
{
    Vector4d adj(m[3], -m[1], -m[2], m[0]);
    Vector4d ddet;
    Ddet(m, ddet);
    double detval = det(m);
    partials << 0, 0, 0, 1,
            0, -1, 0, 0,
            0, 0, -1, 0,
            1, 0, 0, 0;
    partials /= detval;
    partials -= adj*ddet.transpose()/detval/detval;

}

double Midedge::H(const OMMesh &mesh, int faceid, const VectorXd &q)
{
    return 0.5*trace(matMult(matInv(g(mesh, faceid, q)), b(mesh, faceid, q)));
}

void Midedge::DH(const OMMesh &mesh, int faceid, const VectorXd &q, VectorXd &partials)
{
    Vector4d tracepartials;
    Dtrace(tracepartials);
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, faceid, q);
    Matrix4d invpartials;
    DmatInv(gmat, invpartials);
    Matrix4d m1partials, m2partials;
    DmatMult(matInv(gmat), bmat, m1partials, m2partials);
    Vector4d prefactor = ((tracepartials.transpose()*m1partials)*invpartials);
    Dg(mesh, faceid, q, prefactor, partials);
    Vector4d prefactor2 = (tracepartials.transpose()*m2partials);
    VectorXd bpartials;
    Db(mesh, faceid, q, prefactor2, bpartials);
    partials += bpartials;
    partials *= 0.5;
}

double Midedge::K(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q)
{
    return det(matMult(matInv(g(mesh, faceid, q)), b(mesh, faceid, q)));
}

void Midedge::DK(const OMMesh &mesh, int faceid, const VectorXd &q, VectorXd &partials)
{
    Vector4d detpartials;
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, faceid, q);
    Vector4d ginv = matInv(gmat);
    Ddet(matMult(ginv, bmat), detpartials);
    Matrix4d invpartials;
    DmatInv(gmat, invpartials);
    Matrix4d m1partials, m2partials;
    DmatMult(ginv, bmat, m1partials, m2partials);
    Vector4d prefactor1 = ((detpartials.transpose()*m1partials)*invpartials);
    Dg(mesh, faceid, q, prefactor1, partials);
    Vector4d prefactor2 = detpartials.transpose()*m2partials;
    VectorXd bpartials;
    Db(mesh, faceid, q, prefactor2, bpartials);
    partials += bpartials;
}

Vector4d Midedge::c(const OMMesh &mesh, int faceid, const VectorXd &q)
{
    return 2.0*H(mesh, faceid, q)*b(mesh, faceid, q) - K(mesh, faceid, q)*g(mesh, faceid, q);
}

void Midedge::Dc(const OMMesh &mesh, int faceid, const VectorXd &q, const Vector4d &prefactor, VectorXd &partials)
{
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, faceid, q);
    VectorXd db;
    Db(mesh, faceid, q, prefactor, db);
    VectorXd dg;
    Dg(mesh, faceid, q, prefactor, dg);
    partials = 2.0*H(mesh, faceid, q)*db - K(mesh, faceid, q)*dg;
    VectorXd dh;
    DH(mesh, faceid, q, dh);
    VectorXd dk;
    DK(mesh, faceid, q, dk);
    double coeff1 = 2.0*prefactor.dot(bmat);
    partials += coeff1*dh;
    double coeff2 = prefactor.dot(gmat);
    partials -= coeff2*dk;
}


double Midedge::elasticEnergyOne(const OMMesh &mesh, int faceid, const VectorXd &q, const VectorXd &gbar1, const VectorXd &gbar2, const ElasticParameters &params)
{
    Vector4d I(1.0, 0, 0, 1.0);
    Vector4d gbarinv1 = matInv(gbar(mesh, faceid, gbar1));
    Vector4d gbarinv2 = matInv(gbar(mesh, faceid, gbar2));
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, faceid, q);
    Vector4d cmat = c(mesh, faceid, q);

    double A1 = intrinsicArea(mesh, faceid, gbar1, params);
    double A2 = intrinsicArea(mesh, faceid, gbar2, params);

    //std::cout << trace(matMult(gbarinv1,gmat)-I) << std::endl;

    double result = params.h/2.0 * (
                A2*trace(matMult(gbarinv2, gmat) - I)*trace(matMult(gbarinv2, gmat) - I)
                + A1*trace(matMult(gbarinv1, gmat) - I)*trace(matMult(gbarinv1, gmat) - I)
                );
    result += params.h*params.h/4.0/params.scale * (
                A2*trace(matMult(gbarinv2, gmat)-I)*trace(matMult(gbarinv2, bmat))
                - A1*trace(matMult(gbarinv1, gmat)-I)*trace(matMult(gbarinv1, bmat))
                );
    result += params.h*params.h*params.h/24.0/params.scale/params.scale * (
                A2*2.0*trace(matMult(gbarinv2, gmat)-I)*trace(matMult(gbarinv2, cmat))
                + A2*trace(matMult(gbarinv2, bmat))*trace(matMult(gbarinv2, bmat))
                + A1*2.0*trace(matMult(gbarinv1, gmat)-I)*trace(matMult(gbarinv1, cmat))
                + A1*trace(matMult(gbarinv1, bmat))*trace(matMult(gbarinv1, bmat))
                );
    return result;
}

void Midedge::DelasticEnergyOne(const OMMesh &mesh, int faceid, const VectorXd &q, const VectorXd &gbar1, const VectorXd &gbar2, const ElasticParameters &params, Eigen::VectorXd &result)
{
    Vector4d I(1.0, 0, 0, 1.0);
    Vector4d gbarinv1 = matInv(gbar(mesh, faceid, gbar1));
    Vector4d gbarinv2 = matInv(gbar(mesh, faceid, gbar2));
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, faceid, q);
    Vector4d cmat = c(mesh, faceid, q);
    Vector4d dtrace;
    Dtrace(dtrace);

    Matrix4d dgbar2m;
    Matrix4d dgbar1m;
    Matrix4d dummy;
    DmatMult(gbarinv2, gmat, dummy, dgbar2m);
    DmatMult(gbarinv1, gmat, dummy, dgbar1m);

    double A1 = intrinsicArea(mesh, faceid, gbar1, params);
    double A2 = intrinsicArea(mesh, faceid, gbar2, params);

    Vector4d gbar2prefactor = dtrace.transpose()*dgbar2m;
    VectorXd gbar2dg;
    Dg(mesh, faceid, q, gbar2prefactor, gbar2dg);
    Vector4d gbar1prefactor = dtrace.transpose()*dgbar1m;
    VectorXd gbar1dg;
    Dg(mesh, faceid, q, gbar1prefactor, gbar1dg);

    VectorXd gbar2db;
    Db(mesh, faceid, q, gbar2prefactor, gbar2db);
    VectorXd gbar1db;
    Db(mesh, faceid, q, gbar1prefactor, gbar1db);

    VectorXd gbar2dc;
    Dc(mesh, faceid, q, gbar2prefactor, gbar2dc);
    VectorXd gbar1dc;
    Dc(mesh, faceid, q, gbar1prefactor, gbar1dc);

    result = params.h/2.0 * (
                A2*2.0*trace(matMult(gbarinv2, gmat) - I)*gbar2dg
                + A1*2.0*trace(matMult(gbarinv1, gmat) - I)*gbar1dg
                );
    result += params.h*params.h/4.0/params.scale * (
                A2*trace(matMult(gbarinv2, gmat)-I)*gbar2db
                + A2*trace(matMult(gbarinv2, bmat))*gbar2dg
                - A1*trace(matMult(gbarinv1, gmat)-I)*gbar1db
                - A1*trace(matMult(gbarinv1, bmat))*gbar1dg
                );
    result += params.h*params.h*params.h/24.0/params.scale/params.scale * (
                A2*2.0*trace(matMult(gbarinv2, gmat)-I)*gbar2dc
                + A2*2.0*trace(matMult(gbarinv2, cmat))*gbar2dg
                + A2*2.0*trace(matMult(gbarinv2, bmat))*gbar2db
                + A1*2.0*trace(matMult(gbarinv1, gmat)-I)*gbar1dc
                + A1*2.0*trace(matMult(gbarinv1, cmat))*gbar1dg
                + A1*2.0*trace(matMult(gbarinv1, bmat))*gbar1db
                );
}

double Midedge::elasticEnergyTwo(const OMMesh &mesh, int faceid, const VectorXd &q, const VectorXd &gbar1, const VectorXd &gbar2, const ElasticParameters &params)
{
    Vector4d I(1.0, 0, 0, 1.0);
    Vector4d gbarinv1 = matInv(gbar(mesh, faceid, gbar1));
    Vector4d gbarinv2 = matInv(gbar(mesh, faceid, gbar2));
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, faceid, q);
    Vector4d cmat = c(mesh, faceid, q);

    double A1 = intrinsicArea(mesh, faceid, gbar1, params);
    double A2 = intrinsicArea(mesh, faceid, gbar2, params);


    double result = params.h/2.0 * (
                A2*trace(matMult( matMult(gbarinv2, gmat) - I, matMult(gbarinv2, gmat) - I))
                + A1*trace(matMult( matMult(gbarinv1, gmat) - I, matMult(gbarinv1, gmat) - I))
                );
    result += params.h*params.h/4.0/params.scale * (
                A2*trace(matMult( matMult(gbarinv2, gmat)-I, matMult(gbarinv2, bmat) ))
                - A1*trace(matMult( matMult(gbarinv1, gmat)-I, matMult(gbarinv1, bmat) ))
                );
    result += params.h*params.h*params.h/24.0/params.scale/params.scale * (
                A2*2.0*trace(matMult( matMult(gbarinv2, gmat)-I, matMult(gbarinv2, cmat) ))
                + A2*trace(matMult( matMult(gbarinv2, bmat), matMult(gbarinv2, bmat) ))
                + A1*2.0*trace(matMult( matMult(gbarinv1, gmat)-I, matMult(gbarinv1, cmat) ))
                + A1*trace(matMult( matMult(gbarinv1, bmat), matMult(gbarinv1, bmat) ))
                );

    return result;
}

void Midedge::DelasticEnergyTwo(const OMMesh &mesh, int faceid, const VectorXd &q, const VectorXd &gbar1, const VectorXd &gbar2, const ElasticParameters &params, Eigen::VectorXd &result)
{
    Vector4d I(1.0, 0, 0, 1.0);
    Vector4d gbarinv1 = matInv(gbar(mesh, faceid, gbar1));
    Vector4d gbarinv2 = matInv(gbar(mesh, faceid, gbar2));
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, faceid, q);
    Vector4d cmat = c(mesh, faceid, q);
    Vector4d dtrace;
    Dtrace(dtrace);

    double A1 = intrinsicArea(mesh, faceid, gbar1, params);
    double A2 = intrinsicArea(mesh, faceid, gbar2, params);

    Matrix4d dgbar2m;
    Matrix4d dgbar1m;
    Matrix4d dummy;
    DmatMult(gbarinv2, gmat, dummy, dgbar2m);
    DmatMult(gbarinv1, gmat, dummy, dgbar1m);

    Matrix4d term11, term12;
    DmatMult(matMult(gbarinv2, gmat)-I, matMult(gbarinv2, gmat)-I, term12, dummy);
    DmatMult(matMult(gbarinv1, gmat)-I, matMult(gbarinv1, gmat)-I, term11, dummy);

    Vector4d prefactor12 = (dtrace.transpose()*term12)*dgbar2m;
    VectorXd d12;
    Dg(mesh, faceid, q, prefactor12, d12);

    Vector4d prefactor11 = (dtrace.transpose()*term11)*dgbar1m;
    VectorXd d11;
    Dg(mesh, faceid, q, prefactor11, d11);

    result = params.h/2.0 * (
                A2*2.0*d12
                + A1*2.0*d11
                );

    Matrix4d term22a, term22b, term21a, term21b;
    DmatMult(matMult(gbarinv2, bmat),matMult(gbarinv2, gmat)-I,term22a, term22b);
    DmatMult(matMult(gbarinv1, bmat),matMult(gbarinv1, gmat)-I,term21a, term21b);

    Vector4d prefactor22a = (dtrace.transpose()*term22a)*dgbar2m;
    VectorXd d22a;
    Db(mesh, faceid, q, prefactor22a, d22a);

    Vector4d prefactor22b = (dtrace.transpose()*term22b)*dgbar2m;
    VectorXd d22b;
    Dg(mesh, faceid, q, prefactor22b, d22b);

    Vector4d prefactor21a = (dtrace.transpose()*term21a)*dgbar1m;
    VectorXd d21a;
    Db(mesh, faceid, q, prefactor21a, d21a);

    Vector4d prefactor21b = (dtrace.transpose()*term21b)*dgbar1m;
    VectorXd d21b;
    Dg(mesh, faceid, q, prefactor21b, d21b);

    result += params.h*params.h/4.0/params.scale * (
                A2*d22a
                + A2*d22b
                - A1*d21a
                - A1*d21b
                );
    Matrix4d term32a, term32b, term32c, term31a, term31b, term31c;
    DmatMult(matMult(gbarinv2, cmat), matMult(gbarinv2, gmat)-I, term32a, term32b);
    DmatMult(matMult(gbarinv2, bmat), matMult(gbarinv2, bmat), term32c, dummy);
    DmatMult(matMult(gbarinv1, cmat), matMult(gbarinv1, gmat)-I, term31a, term31b);
    DmatMult(matMult(gbarinv1, bmat), matMult(gbarinv1, bmat), term31c, dummy);

    Vector4d prefactor32a = (dtrace.transpose()*term32a)*dgbar2m;
    VectorXd d32a;
    Dc(mesh, faceid, q, prefactor32a, d32a);

    Vector4d prefactor32b = (dtrace.transpose()*term32b)*dgbar2m;
    VectorXd d32b;
    Dg(mesh, faceid, q, prefactor32b, d32b);

    Vector4d prefactor32c = (dtrace.transpose()*term32c)*dgbar2m;
    VectorXd d32c;
    Db(mesh, faceid, q, prefactor32c, d32c);


    Vector4d prefactor31a = (dtrace.transpose()*term31a)*dgbar1m;
    VectorXd d31a;
    Dc(mesh, faceid, q, prefactor31a, d31a);

    Vector4d prefactor31b = (dtrace.transpose()*term31b)*dgbar1m;
    VectorXd d31b;
    Dg(mesh, faceid, q, prefactor31b, d31b);

    Vector4d prefactor31c = (dtrace.transpose()*term31c)*dgbar1m;
    VectorXd d31c;
    Db(mesh, faceid, q, prefactor31c, d31c);

    result += params.h*params.h*params.h/24.0/params.scale/params.scale * (
                A2*2.0*d32a
                + A2*2.0*d32b
                + A2*2.0*d32c
                + A1*2.0*d31a
                + A1*2.0*d31b
                + A1*2.0*d31c
                );
}

double Midedge::elasticEnergy(const OMMesh &mesh, const VectorXd &q, const VectorXd &gbar1, const VectorXd &gbar2, const ElasticParameters &params)
{
    double result = 0;
    for(OMMesh::FaceIter fi = mesh.faces_begin(); fi != mesh.faces_end(); ++fi)
    {
        result += params.YoungsModulus*params.PoissonRatio/(1.0-params.PoissonRatio*params.PoissonRatio)*elasticEnergyOne(mesh, fi.handle().idx(), q, gbar1, gbar2, params);
        result += params.YoungsModulus/(1.0+params.PoissonRatio)*elasticEnergyTwo(mesh, fi.handle().idx(), q, gbar1, gbar2, params);
    }
    return result;
}

void Midedge::DelasticEnergy(const OMMesh &mesh, const VectorXd &q, const VectorXd &gbar1, const VectorXd &gbar2, const ElasticParameters &params, VectorXd &result)
{
    result.resize(q.size());
    result.setZero();
    for(OMMesh::FaceIter fi = mesh.faces_begin(); fi != mesh.faces_end(); ++fi)
    {
        VectorXd tmp;
        DelasticEnergyOne(mesh, fi.handle().idx(), q, gbar1, gbar2, params, tmp);
        result += params.YoungsModulus*params.PoissonRatio/(1.0-params.PoissonRatio*params.PoissonRatio)*tmp;
        DelasticEnergyTwo(mesh, fi.handle().idx(), q, gbar1, gbar2, params, tmp);
        result += params.YoungsModulus/(1.0+params.PoissonRatio)*tmp;
    }
}

double Midedge::intrinsicArea(const OMMesh &mesh, int faceid, const VectorXd &gbar, const ElasticParameters &params)
{
    double gs[3];
    int idx=0;
    OMMesh::FaceHandle fh = mesh.face_handle(faceid);
    for(OMMesh::ConstFaceEdgeIter fei = mesh.cfe_iter(fh); fei; ++fei)
    {
        gs[idx++] = gbar[fei.handle().idx()];
    }
    double s = 0.5*(gs[0]+gs[1]+gs[2]);
    return params.scale*params.scale*0.5*sqrt(s*(s-gs[0])*(s-gs[1])*(s-gs[2]));
}
