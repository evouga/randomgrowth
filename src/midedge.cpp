#include "midedge.h"
#include <Eigen/Dense>
#include <GL/gl.h>

using namespace Eigen;
using namespace std;

Midedge::Midedge()
{

}

Eigen::Matrix3d Midedge::crossMatrix(const Vector3d &v)
{
    Matrix3d result;
    result(0,0)=result(1,1)=result(2,2) = 0;
    result(0,1) = -v[2];
    result(0,2) = v[1];
    result(1,0) = v[2];
    result(2,0) = -v[1];
    result(1,2) = -v[0];
    result(2,1) = v[0];
    /*result << 0, -v[2], v[1],
            v[2], 0, -v[0],
            -v[1], v[0], 0;*/
    return result;
}

double Midedge::area(const Vector3d &q0, const Vector3d &q1, const Vector3d &q2)
{
    return 0.5*(q0-q1).cross(q2-q1).norm();
}

void Midedge::Darea(const Vector3d &q0, const Vector3d &q1, const Vector3d &q2, Vector3d *partials)
{
    Vector3d inner = (q0-q1).cross(q2-q1);
    double denom = 2.0*inner.norm();
    partials[0] = -inner.transpose()*crossMatrix(q2-q1)/denom;
    partials[1] = inner.transpose()*(-crossMatrix(q0-q1) + crossMatrix(q2-q1))/denom;
    partials[2] = inner.transpose()*crossMatrix(q0-q1)/denom;
}

Vector3d Midedge::diamondEdgeNormal(const Vector3d &q0, const Vector3d &q1, const Vector3d &q2, const Vector3d &q3)
{
    double a1 = area(q0, q2, q3);
    double a2 = area(q0, q1, q2);
    assert(a1 > 0);
    assert(a2 > 0);
    Vector3d result = a1 * (q1-q0).cross(q2-q0) + a2*(q2-q0).cross(q3-q0);
    assert(result.norm() > 0);
    return result/result.norm();
}

void Midedge::DdiamondEdgeNormal(const Vector3d &q0, const Vector3d &q1, const Vector3d &q2, const Vector3d &q3, std::vector<Matrix3d> &partials)
{
    double a1 = area(q0, q2, q3);
    assert(a1 > 0);
    assert(!isnan(a1));
    Vector3d da1[3];
    Darea(q0, q2, q3, da1);

    double a2 = area(q0, q1, q2);
    assert(a2 > 0);
    assert(!isnan(a2));
    Vector3d da2[3];
    Darea(q0, q1, q2, da2);

    Vector3d num = a1*(q1-q0).cross(q2-q0) + a2*(q2-q0).cross(q3-q0);
    double denom = num.norm();    
    assert(denom > 0);
    Matrix3d dn0 = (q1-q0).cross(q2-q0) * da1[0].transpose()
                        + a1 * crossMatrix(q2-q1)
                        + (q2-q0).cross(q3-q0) * da2[0].transpose()
                        + a2 * crossMatrix(q3-q2);
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

void Midedge::Db(const Vector3d &q1, const Vector3d &q2, const Vector3d &q3, const Vector3d &n1, const Vector3d &n2, const Vector3d &n3, Matrix3d *partials)
{
    for(int j=0; j<3; j++)
    {
        partials[0](0,j) = 0;
        partials[0](1,j) = 0;
        partials[0](2,j) = 2.0*(n1-n3)[j];

        partials[1](0,j) = -2.0*(n3-n2)[j];
        partials[1](1,j) = -2.0*(n1-n3)[j];
        partials[1](2,j) = 0;

        partials[2](0,j) = 2.0*(n3-n2)[j];
        partials[2](1,j) = 2.0*(n1-n3)[j];
        partials[2](2,j) = -2.0*(n1-n3)[j];

        partials[3](0,j) = 0;
        partials[3](1,j) = 2.0*(q3-q2)[j];
        partials[3](2,j) = 2.0*(q1-q3)[j];

        partials[4](0,j) = -2.0*(q3-q2)[j];
        partials[4](1,j) = 0;
        partials[4](2,j) = 0;

        partials[5](0,j) = 2.0*(q3-q2)[j];
        partials[5](1,j) = -2.0*(q3-q2)[j];
        partials[5](2,j) = -2.0*(q1-q3)[j];
    }
}

Vector3d Midedge::faceNormal(const Vector3d &q0, const Vector3d &q1, const Vector3d &q2)
{
    Vector3d result = (q1-q0).cross(q2-q0);
    for(int i=0; i<3; i++)
        assert(!isnan(result[i]));
    return result/result.norm();
}

void Midedge::DfaceNormal(const Vector3d &q0, const Vector3d &q1, const Vector3d &q2, std::vector<Matrix3d> &partials)
{
    Vector3d n = (q1-q0).cross(q2-q0);
    double denom = n.norm();
    assert(denom > 0);

    Matrix3d dn0 = crossMatrix(q2-q0) - crossMatrix(q1-q0);
    partials.push_back( dn0/denom - n*n.transpose()*dn0/denom/denom/denom );

    Matrix3d dn1 = -crossMatrix(q2-q0);
    partials.push_back( dn1/denom - n*n.transpose()*dn1/denom/denom/denom );

    Matrix3d dn2 = crossMatrix(q1-q0);
    partials.push_back( dn2/denom - n*n.transpose()*dn2/denom/denom/denom );
}

Vector3d Midedge::edgeNormal(int edgeid, const EdgeNormalDerivatives &nderivs)
{
    return nderivs.normals[edgeid];
}

void Midedge::DedgeNormal(int edgeid, const EdgeNormalDerivatives &dnormals, const Vector3d &prefactor, VectorXd &partials)
{
    for(int i=dnormals.startidx[edgeid]; i<dnormals.startidx[edgeid+1]; i++)
    {
        partials.segment<3>(3*dnormals.vertidx[i]) += prefactor.transpose()*dnormals.partials[i];
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

Vector4d Midedge::b(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const VectorXd &q)
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
        ns[(2+curidx)%3] = edgeNormal(edgeid, nderivs);
        curidx++;
    }
    Vector3d b3 = b(qs[0], qs[1], qs[2], ns[0], ns[1], ns[2]);
    Vector4d result(b3[0], b3[1], b3[1], b3[2]);

    return result;
}

void Midedge::Db(const OMMesh &mesh, const EdgeNormalDerivatives &dnormals, int faceid, const VectorXd &q, const Vector4d &prefactor, VectorXd &partials)
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
        ns[(2+curidx)%3] = edgeNormal(edgeids[curidx], dnormals);
        curidx++;
    }
    Matrix3d partials3[6];
    Db(qs[0], qs[1], qs[2], ns[0], ns[1], ns[2], partials3);

    Vector3d condensedPrefactor(prefactor[0], prefactor[1]+prefactor[2], prefactor[3]);
    for(int i=0; i<3; i++)
    {
        Vector3d prod = condensedPrefactor.transpose()*partials3[i];
        for(int j=0; j<3; j++)
            partials[3*idx[i]+j] += prod[j];
    }
    for(int i=0; i<3; i++)
    {
        Vector3d prod = condensedPrefactor.transpose()*partials3[3+i];
        DedgeNormal(edgeids[(1+i)%3], dnormals, prod, partials);
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

const Vector4d Midedge::Dtrace(1.0, 0.0, 0.0, 1.0);

Vector4d Midedge::inducedG(const OMMesh &mesh, int faceid, const VectorXd &q)
{
    double gi[3];
    int idx=0;
    OMMesh::FaceHandle fh = mesh.face_handle(faceid);
    for(OMMesh::ConstFaceHalfedgeIter fhi = mesh.cfh_iter(fh); fhi; ++fhi)
    {
        Vector3d topt = q.segment<3>(3*mesh.to_vertex_handle(fhi.handle()).idx());
        Vector3d frompt = q.segment<3>(3*mesh.from_vertex_handle(fhi.handle()).idx());
        gi[idx++] = (topt-frompt).norm();
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

double Midedge::H(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const VectorXd &q)
{
    return 0.5*trace(matMult(matInv(g(mesh, faceid, q)), b(mesh, nderivs, faceid, q)));
}

void Midedge::DH(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const VectorXd &q, double prefactor, VectorXd &partials)
{
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, nderivs, faceid, q);
    Matrix4d invpartials;
    DmatInv(gmat, invpartials);
    Matrix4d m1partials, m2partials;
    DmatMult(matInv(gmat), bmat, m1partials, m2partials);
    Vector4d prefactor1 = 0.5*prefactor*((Dtrace.transpose()*m1partials)*invpartials);
    Dg(mesh, faceid, q, prefactor1, partials);
    Vector4d prefactor2 = 0.5*prefactor*(Dtrace.transpose()*m2partials);
    Db(mesh, nderivs, faceid, q, prefactor2, partials);
}

double Midedge::K(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q)
{
    Vector4d bmat = b(mesh, nderivs, faceid, q);
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d ginv = matInv(gmat);
    Vector4d ginvb = matMult(ginv, bmat);
    return det(ginvb);
}

void Midedge::DK(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const VectorXd &q, double prefactor, VectorXd &partials)
{
    Vector4d detpartials;
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, nderivs, faceid, q);
    Vector4d ginv = matInv(gmat);
    Ddet(matMult(ginv, bmat), detpartials);
    Matrix4d invpartials;
    DmatInv(gmat, invpartials);
    Matrix4d m1partials, m2partials;
    DmatMult(ginv, bmat, m1partials, m2partials);
    Vector4d prefactor1 = prefactor*((detpartials.transpose()*m1partials)*invpartials);
    Dg(mesh, faceid, q, prefactor1, partials);
    Vector4d prefactor2 = prefactor*detpartials.transpose()*m2partials;
    Db(mesh, nderivs, faceid, q, prefactor2, partials);
}

Vector4d Midedge::c(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const VectorXd &q)
{
    return 2.0*H(mesh, nderivs, faceid, q)*b(mesh, nderivs, faceid, q) - K(mesh, nderivs, faceid, q)*g(mesh, faceid, q);
}

void Midedge::Dc(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const VectorXd &q, const Vector4d &prefactor, VectorXd &partials)
{
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, nderivs, faceid, q);
    Db(mesh, nderivs, faceid, q, 2.0*H(mesh, nderivs, faceid, q)*prefactor, partials);
    Dg(mesh, faceid, q, - K(mesh, nderivs, faceid, q)*prefactor, partials);
    DH(mesh, nderivs, faceid, q, 2.0*prefactor.dot(bmat), partials);
    DK(mesh, nderivs, faceid, q,  -prefactor.dot(gmat), partials);
}


double Midedge::elasticEnergyOne(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const VectorXd &q, const VectorXd &gbar, const ElasticParameters &params)
{
    Vector4d I(1.0, 0, 0, 1.0);
    Vector4d gbarinv = matInv(gbar.segment<4>(4*faceid));
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, nderivs, faceid, q);
    Vector4d cmat = c(mesh, nderivs, faceid, q);

    double A = intrinsicArea(faceid, gbar, params);
    for(int i=0; i<4; i++)
    {
        assert(!isnan(gmat[i]));
        assert(!isnan(bmat[i]));
        assert(!isnan(cmat[i]));
    }
    assert(!isnan(A));

    double result = params.h*A*trace(matMult(gbarinv, gmat) - I)*trace(matMult(gbarinv, gmat) - I);

    result += params.h*params.h*params.h/24.0/params.scale/params.scale * (
                A*4.0*trace(matMult(gbarinv, gmat)-I)*trace(matMult(gbarinv, cmat))
                + A*2.0*trace(matMult(gbarinv, bmat))*trace(matMult(gbarinv, bmat))
                );

    assert(!isnan(result));
    return result;
}

void Midedge::DelasticEnergyOne(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const VectorXd &q, const VectorXd &gbar, const ElasticParameters &params, double prefactor, Eigen::VectorXd &result)
{
    Vector4d I(1.0, 0, 0, 1.0);
    Vector4d gbarinv = matInv(gbar.segment<4>(4*faceid));
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, nderivs, faceid, q);
    Vector4d cmat = c(mesh, nderivs, faceid, q);

    Matrix4d dgbarm;
    Matrix4d dummy;
    DmatMult(gbarinv, gmat, dummy, dgbarm);

    double A = intrinsicArea(faceid, gbar, params);

    double coeff1 = params.h;
    double coeff3 = params.h*params.h*params.h/24.0/params.scale/params.scale;

    Vector4d gbarprefactor = Dtrace.transpose()*dgbarm;

    Vector4d totalgprefactor = (2.0*coeff1*A*trace(matMult(gbarinv, gmat) - I)
                                    +2.0*coeff3*A*2.0*trace(matMult(gbarinv, cmat))
                                    )*gbarprefactor;

    Dg(mesh, faceid, q, prefactor*totalgprefactor, result);

    Vector4d totalbprefactor = (coeff3*A*4.0*trace(matMult(gbarinv, bmat))
                                     )*gbarprefactor;

    Db(mesh, nderivs, faceid, q, prefactor*totalbprefactor, result);

    Vector4d totalcprefactor = coeff3*A*4.0*trace(matMult(gbarinv, gmat)-I)*gbarprefactor;
    Dc(mesh, nderivs, faceid, q, prefactor*totalcprefactor, result);
}

double Midedge::elasticEnergyTwo(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const VectorXd &q, const VectorXd &gbar, const ElasticParameters &params)
{
    Vector4d I(1.0, 0, 0, 1.0);
    Vector4d gbarinv = matInv(gbar.segment<4>(4*faceid));
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, nderivs, faceid, q);
    Vector4d cmat = c(mesh, nderivs, faceid, q);

    double A = intrinsicArea(faceid, gbar, params);

    double result = params.h *
                A*trace(matMult( matMult(gbarinv, gmat) - I, matMult(gbarinv, gmat) - I));
    result += params.h*params.h*params.h/24.0/params.scale/params.scale * (
                A*4.0*trace(matMult( matMult(gbarinv, gmat)-I, matMult(gbarinv, cmat) ))
                + A*trace(matMult( matMult(gbarinv, bmat), matMult(gbarinv, bmat) ))
                );

    assert(!isnan(result));
    return result;
}

void Midedge::DelasticEnergyTwo(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const VectorXd &q, const VectorXd &gbar, const ElasticParameters &params, double prefactor, Eigen::VectorXd &result)
{
    Vector4d I(1.0, 0, 0, 1.0);
    Vector4d gbarinv = matInv(gbar.segment<4>(4*faceid));
    Vector4d gmat = g(mesh, faceid, q);
    Vector4d bmat = b(mesh, nderivs, faceid, q);
    Vector4d cmat = c(mesh, nderivs, faceid, q);

    double A = intrinsicArea(faceid, gbar, params);
    double coeff1 = params.h;
    double coeff3 = params.h*params.h*params.h/24.0/params.scale/params.scale;

    Matrix4d dgbarm;
    Matrix4d dummy;
    DmatMult(gbarinv, gmat, dummy, dgbarm);

    Matrix4d term1;
    DmatMult(matMult(gbarinv, gmat)-I, matMult(gbarinv, gmat)-I, term1, dummy);

    Vector4d prefactor1 = (Dtrace.transpose()*term1)*dgbarm;

    Matrix4d term2a, term2b;
    DmatMult(matMult(gbarinv, bmat),matMult(gbarinv, gmat)-I,term2a, term2b);

    Matrix4d term3a, term3b, term3c;
    DmatMult(matMult(gbarinv, cmat), matMult(gbarinv, gmat)-I, term3a, term3b);
    DmatMult(matMult(gbarinv, bmat), matMult(gbarinv, bmat), term3c, dummy);

    Vector4d prefactor3a = (Dtrace.transpose()*term3a)*dgbarm;

    Vector4d prefactor3b = (Dtrace.transpose()*term3b)*dgbarm;

    Vector4d prefactor3c = (Dtrace.transpose()*term3c)*dgbarm;

    Vector4d totalgprefactor = coeff1*2.0*A*prefactor1
            + coeff3*A*4.0*prefactor3b
            ;

    Dg(mesh, faceid, q, prefactor*totalgprefactor, result);

    Vector4d totalbprefactor = coeff3*A*4.0*prefactor3c
            ;

    Db(mesh, nderivs, faceid, q, prefactor*totalbprefactor, result);

    Vector4d totalcprefactor = coeff3*A*4.0*prefactor3a;

    Dc(mesh, nderivs, faceid, q, prefactor*totalcprefactor, result);
}

void Midedge::visualizeNormals(const OMMesh &mesh, const VectorXd &q, double scale)
{
    EdgeNormalDerivatives nderivs;
    gatherEdgeNormalDerivatives(mesh, q, nderivs);
    glBegin(GL_LINES);
    glColor3f(1.0,0.0,1.0);
    for(OMMesh::EdgeIter eh = mesh.edges_begin(); eh != mesh.edges_end(); ++eh)
    {
        OMMesh::HalfedgeHandle heh = mesh.halfedge_handle(eh.handle(),0);
        int src = mesh.from_vertex_handle(heh).idx();
        int dst = mesh.to_vertex_handle(heh).idx();

        Vector3d pos = 0.5*(q.segment<3>(3*src)+q.segment<3>(3*dst));
        Vector3d n = nderivs.normals[eh.handle().idx()];
        n *= scale;
        glVertex3d(pos[0],pos[1],pos[2]);
        glVertex3d(pos[0]+n[0],pos[1]+n[1],pos[2]+n[2]);
    }
    glEnd();
}

double Midedge::elasticEnergy(const OMMesh &mesh, const VectorXd &q, const VectorXd &gbar, const ElasticParameters &params, VectorXd *derivs)
{
    double result = 0;
    if(derivs)
    {
        if(derivs->size() != q.size())
            derivs->resize(q.size());
        derivs->setZero();
    }
    EdgeNormalDerivatives nderivs;
    gatherEdgeNormalDerivatives(mesh, q, nderivs);


    for(OMMesh::FaceIter fi = mesh.faces_begin(); fi != mesh.faces_end(); ++fi)
    {
        result += params.YoungsModulus/8.0*params.PoissonRatio/(1.0-params.PoissonRatio*params.PoissonRatio)*elasticEnergyOne(mesh, nderivs, fi.handle().idx(), q, gbar, params);
        result += params.YoungsModulus/8.0/(1.0+params.PoissonRatio)*elasticEnergyTwo(mesh, nderivs, fi.handle().idx(), q, gbar, params);

        if(derivs)
        {
            DelasticEnergyOne(mesh, nderivs, fi.handle().idx(), q, gbar, params, params.YoungsModulus/8.0*params.PoissonRatio/(1.0-params.PoissonRatio*params.PoissonRatio), *derivs);
            DelasticEnergyTwo(mesh, nderivs, fi.handle().idx(), q, gbar, params, params.YoungsModulus/8.0/(1.0+params.PoissonRatio), *derivs);
        }
    }

    return result;
}

double Midedge::intrinsicArea(int faceid, const VectorXd &gbar, const ElasticParameters &params)
{
    double gdet = det(gbar.segment<4>(4*faceid));
    assert(gdet >= 0);
    return params.scale*params.scale*0.5*sqrt(det(gbar.segment<4>(4*faceid)));
}

void Midedge::gatherEdgeNormalDerivatives(const OMMesh &mesh, const VectorXd &q, EdgeNormalDerivatives &dnormals)
{
    dnormals.partials.clear();
    dnormals.vertidx.clear();
    dnormals.startidx.clear();

    for(OMMesh::EdgeIter eh = mesh.edges_begin(); eh != mesh.edges_end(); ++eh)
    {
        dnormals.startidx.push_back(dnormals.partials.size());
        if(mesh.is_boundary(eh.handle()))
        {
            OMMesh::HalfedgeHandle heh = mesh.halfedge_handle(eh.handle(),0);
            if(mesh.is_boundary(heh))
                heh = mesh.halfedge_handle(eh.handle(),1);
            assert(!mesh.is_boundary(heh));

            Vector3d qs[3];
            int idx[3];
            idx[0] = mesh.from_vertex_handle(heh).idx();
            qs[0] = q.segment<3>(3*idx[0]);
            idx[1] = mesh.to_vertex_handle(heh).idx();
            qs[1] = q.segment<3>(3*idx[1]);
            idx[2] = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)).idx();
            qs[2] = q.segment<3>(3*idx[2]);
            DfaceNormal(qs[0], qs[1], qs[2], dnormals.partials);
            dnormals.vertidx.push_back(idx[0]);
            dnormals.vertidx.push_back(idx[1]);
            dnormals.vertidx.push_back(idx[2]);
            dnormals.normals.push_back(faceNormal(qs[0], qs[1], qs[2]));
        }
        else
        {
            OMMesh::HalfedgeHandle heh = mesh.halfedge_handle(eh.handle(),0);
            Vector3d qs[4];
            int idx[4];
            idx[0] = mesh.to_vertex_handle(heh).idx();
            qs[0] = q.segment<3>(3*idx[0]);
            idx[1] = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)).idx();
            qs[1] = q.segment<3>(3*idx[1]);
            idx[2] = mesh.from_vertex_handle(heh).idx();
            qs[2] = q.segment<3>(3*idx[2]);
            idx[3] = mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(heh))).idx();
            qs[3] = q.segment<3>(3*idx[3]);
            DdiamondEdgeNormal(qs[0], qs[1], qs[2], qs[3], dnormals.partials);
            dnormals.vertidx.push_back(idx[0]);
            dnormals.vertidx.push_back(idx[1]);
            dnormals.vertidx.push_back(idx[2]);
            dnormals.vertidx.push_back(idx[3]);
            dnormals.normals.push_back(diamondEdgeNormal(qs[0], qs[1], qs[2], qs[3]));
        }
    }
    dnormals.startidx.push_back(dnormals.partials.size());
}

void Midedge::gaussianCurvature(const OMMesh &mesh, const VectorXd &q, VectorXd &Ks)
{
    int nfaces = (int)mesh.n_faces();
    Ks.resize(nfaces);
    EdgeNormalDerivatives nderivs;
    gatherEdgeNormalDerivatives(mesh, q, nderivs);
    for(int i=0; i<(int)mesh.n_faces(); i++)
    {
        Ks[i] = K(mesh, nderivs, i, q);
    }
}

void Midedge::meanCurvature(const OMMesh &mesh, const VectorXd &q, VectorXd &Hs)
{
    int nfaces = (int)mesh.n_faces();
    Hs.resize(nfaces);
    EdgeNormalDerivatives nderivs;
    gatherEdgeNormalDerivatives(mesh, q, nderivs);
    for(int i=0; i<(int)mesh.n_faces(); i++)
    {
        Hs[i] = H(mesh, nderivs, i, q);
    }
}
