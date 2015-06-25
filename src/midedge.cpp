#include "midedge.h"
#include "mesh.h"
#include <Eigen/Dense>
#include <GL/gl.h>
#include <iostream>

using namespace Eigen;
using namespace std;

const Vector3d &EdgeNormalDerivatives::edgeNormal(int face, int vertno) const
{
    return normals[3*face+vertno];
}

void EdgeNormalDerivatives::DedgeNormal(int face, int vertno, const Vector3d &prefactor, VectorXd &partials) const
{
    int numpartials = (int)edgeNormalsPartials[3*face+vertno].size();
    for(int i=0; i<numpartials; i++)
    {
        for(int j=0; j<3; j++)
        {
            #pragma omp atomic
            partials[3*edgeNormalsPartialIndices[3*face+vertno][i] + j] += (prefactor.transpose()*edgeNormalsPartials[3*face+vertno][i])[j];
        }
    }
}



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
    assert(!std::isnan(a1));
    Vector3d da1[3];
    Darea(q0, q2, q3, da1);

    double a2 = area(q0, q1, q2);
    assert(a2 > 0);
    assert(!std::isnan(a2));
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
        assert(!std::isnan(result[i]));
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

Vector4d Midedge::g(const Mesh &mesh, const VectorXd &q, int faceid)
{
    Vector3d qs[3];
    for(int i=0; i<3; i++)
    {
        qs[i] = q.segment<3>(3*mesh.faceVerts(faceid)[i]);
    }

    Vector3d g3 = g(qs[0], qs[1], qs[2]);
    Vector4d result;
    result[0] = g3[0];
    result[1] = result[2] = g3[1];
    result[3] = g3[2];
    return result;
}

void Midedge::Dg(const Mesh &mesh, const VectorXd &q, int faceid, const Vector4d &prefactor, VectorXd &partials)
{

    Vector3d qs[3];
    int idx[3];
    for(int i=0; i<3; i++)
    {
        idx[i] = mesh.faceVerts(faceid)[i];
        qs[i] = q.segment<3>(3*idx[i]);
    }
    vector<Matrix3d> partials3;
    Dg(qs[0], qs[1], qs[2], partials3);

    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            #pragma omp atomic
            partials[3*idx[i]+j] += partials3[i](0,j)*prefactor[0];
            #pragma omp atomic
            partials[3*idx[i]+j] += partials3[i](1,j)*prefactor[1];
            #pragma omp atomic
            partials[3*idx[i]+j] += partials3[i](1,j)*prefactor[2];
            #pragma omp atomic
            partials[3*idx[i]+j] += partials3[i](2,j)*prefactor[3];
        }
    }
}

Vector4d Midedge::b(const Mesh &mesh, const VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid)
{
    Vector3d ns[3];
    Vector3d qs[3];
    for(int i=0; i<3; i++)
    {
        qs[i] = q.segment<3>(3*mesh.faceVerts(faceid)[i]);
        ns[i] = nderivs.edgeNormal(faceid, i);
    }
    Vector3d b3 = b(qs[0], qs[1], qs[2], ns[0], ns[1], ns[2]);
    Vector4d result(b3[0], b3[1], b3[1], b3[2]);

    return result;
}

void Midedge::Db(const Mesh &mesh, const VectorXd &q, const EdgeNormalDerivatives &dnormals, int faceid, const Vector4d &prefactor, VectorXd &partials)
{
    Vector3d ns[3];
    int idx[3];
    Vector3d qs[3];
    for(int i=0; i<3; i++)
    {
        idx[i] = mesh.faceVerts(faceid)[i];
        qs[i] = q.segment<3>(3*idx[i]);
        ns[i] = dnormals.edgeNormal(faceid, i);
    }

    Matrix3d partials3[6];
    Db(qs[0], qs[1], qs[2], ns[0], ns[1], ns[2], partials3);

    Vector3d condensedPrefactor(prefactor[0], prefactor[1]+prefactor[2], prefactor[3]);
    for(int i=0; i<3; i++)
    {
        Vector3d prod = condensedPrefactor.transpose()*partials3[i];
        for(int j=0; j<3; j++)
        {
            #pragma omp atomic
            partials[3*idx[i]+j] += prod[j];
        }
    }
    for(int i=0; i<3; i++)
    {
        Vector3d prod = condensedPrefactor.transpose()*partials3[3+i];
        dnormals.DedgeNormal(faceid, i, prod, partials);
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

//Vector4d Midedge::inducedG(const Mesh &mesh, int faceid)
//{
//    double gi[3];
//    int idx=0;
//    OMMesh::FaceHandle fh = mesh.face_handle(faceid);
//    for(OMMesh::ConstFaceHalfedgeIter fhi = mesh.cfh_iter(fh); fhi; ++fhi)
//    {
//        Vector3d topt = q.segment<3>(3*mesh.to_vertex_handle(fhi.handle()).idx());
//        Vector3d frompt = q.segment<3>(3*mesh.from_vertex_handle(fhi.handle()).idx());
//        gi[idx++] = (topt-frompt).norm();
//    }
//    Vector4d result;
//    result << gi[1]*gi[1],
//            -0.5*(gi[1]*gi[1]+gi[2]*gi[2]-gi[0]*gi[0]),
//            -0.5*(gi[1]*gi[1]+gi[2]*gi[2]-gi[0]*gi[0]),
//            gi[2]*gi[2];
//    return result;
//}

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

double Midedge::H(const Mesh &mesh, const VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid)
{
    return 0.5*trace(matMult(matInv(g(mesh, q, faceid)), b(mesh, q, nderivs, faceid)));
}

void Midedge::DH(const Mesh &mesh, const VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid, double prefactor, VectorXd &partials)
{
    Vector4d gmat = g(mesh, q, faceid);
    Vector4d bmat = b(mesh, q, nderivs, faceid);
    Matrix4d invpartials;
    DmatInv(gmat, invpartials);
    Matrix4d m1partials, m2partials;
    DmatMult(matInv(gmat), bmat, m1partials, m2partials);
    Vector4d prefactor1 = 0.5*prefactor*((Dtrace.transpose()*m1partials)*invpartials);
    Dg(mesh, q, faceid, prefactor1, partials);
    Vector4d prefactor2 = 0.5*prefactor*(Dtrace.transpose()*m2partials);
    Db(mesh, q, nderivs, faceid, prefactor2, partials);
}

double Midedge::K(const Mesh &mesh, const VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid)
{
    Vector4d bmat = b(mesh, q, nderivs, faceid);
    Vector4d gmat = g(mesh, q, faceid);
    Vector4d ginv = matInv(gmat);
    Vector4d ginvb = matMult(ginv, bmat);
    return det(ginvb);
}

void Midedge::DK(const Mesh &mesh, const VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid, double prefactor, VectorXd &partials)
{
    Vector4d detpartials;
    Vector4d gmat = g(mesh, q, faceid);
    Vector4d bmat = b(mesh, q, nderivs, faceid);
    Vector4d ginv = matInv(gmat);
    Ddet(matMult(ginv, bmat), detpartials);
    Matrix4d invpartials;
    DmatInv(gmat, invpartials);
    Matrix4d m1partials, m2partials;
    DmatMult(ginv, bmat, m1partials, m2partials);
    Vector4d prefactor1 = prefactor*((detpartials.transpose()*m1partials)*invpartials);
    Dg(mesh, q, faceid, prefactor1, partials);
    Vector4d prefactor2 = prefactor*detpartials.transpose()*m2partials;
    Db(mesh, q, nderivs, faceid, prefactor2, partials);
}

Vector4d Midedge::c(const Mesh &mesh, const VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid)
{
    return 2.0*H(mesh, q, nderivs, faceid)*b(mesh, q, nderivs, faceid) - K(mesh, q, nderivs, faceid)*g(mesh, q, faceid);
}

void Midedge::Dc(const Mesh &mesh, const VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid, const Vector4d &prefactor, VectorXd &partials)
{
    Vector4d gmat = g(mesh, q, faceid);
    Vector4d bmat = b(mesh, q, nderivs, faceid);
    Db(mesh, q, nderivs, faceid, 2.0*H(mesh, q, nderivs, faceid)*prefactor, partials);
    Dg(mesh, q ,faceid, - K(mesh, q, nderivs, faceid)*prefactor, partials);
    DH(mesh, q, nderivs, faceid, 2.0*prefactor.dot(bmat), partials);
    DK(mesh, q, nderivs, faceid, -prefactor.dot(gmat), partials);
}


double Midedge::elasticEnergyOne(const Mesh &mesh, const VectorXd &q, const PrecomputedFaceQuantities &data, int faceid, const ElasticParameters &params)
{
    double A = data.intrinsicArea;
    double trgminusI = data.trgminusI;
    double trb = data.trb;
    double trc = data.trc;

    assert(!std::isnan(A));

    double result = params.h*A*trgminusI*trgminusI;

    result += params.h*params.h*params.h/24.0/params.scale/params.scale * (
                A*4.0*trgminusI*trc
                + A*2.0*trb*trb
                );

    assert(!std::isnan(result));
    return result;
}

void Midedge::DelasticEnergyOne(const Mesh &mesh, const VectorXd &q, const EdgeNormalDerivatives &nderivs, const PrecomputedFaceQuantities &data, int faceid, const ElasticParameters &params, double prefactor, Eigen::VectorXd &result)
{
    Vector4d I(1.0, 0, 0, 1.0);
    const Vector4d &gbarinv = data.gbarinv;
    const Vector4d &gmat = data.g;

    Matrix4d dgbarm;
    Matrix4d dummy;
    DmatMult(gbarinv, gmat, dummy, dgbarm);

    double A = data.intrinsicArea;
    double trgminusI = data.trgminusI;
    double trb = data.trb;
    double trc = data.trc;

    double coeff1 = params.h;
    double coeff3 = params.h*params.h*params.h/24.0/params.scale/params.scale;

    Vector4d gbarprefactor = Dtrace.transpose()*dgbarm;

    Vector4d totalgprefactor = (2.0*coeff1*A*trgminusI
                                    +2.0*coeff3*A*2.0*trc
                                    )*gbarprefactor;

    Dg(mesh, q, faceid, prefactor*totalgprefactor, result);

    Vector4d totalbprefactor = (coeff3*A*4.0*trb
                                     )*gbarprefactor;

    Db(mesh, q, nderivs, faceid, prefactor*totalbprefactor, result);

    Vector4d totalcprefactor = coeff3*A*4.0*trgminusI*gbarprefactor;
    Dc(mesh, q, nderivs, faceid, prefactor*totalcprefactor, result);
}

double Midedge::elasticEnergyTwo(const Mesh &mesh, const VectorXd &q, const PrecomputedFaceQuantities &data, int faceid, const ElasticParameters &params)
{
    const Vector4d I(1.0, 0, 0, 1.0);
    const Vector4d &gbarinv = data.gbarinv;
    const Vector4d &gmat = data.g;
    const Vector4d &bmat = data.b;
    const Vector4d &cmat = data.c;

    double A = data.intrinsicArea;

    double result = params.h *
                A*trace(matMult( matMult(gbarinv, gmat) - I, matMult(gbarinv, gmat) - I));
    result += params.h*params.h*params.h/24.0/params.scale/params.scale * (
                A*4.0*trace(matMult( matMult(gbarinv, gmat)-I, matMult(gbarinv, cmat) ))
                + A*trace(matMult( matMult(gbarinv, bmat), matMult(gbarinv, bmat) ))
                );

    assert(!std::isnan(result));
    return result;
}

void Midedge::DelasticEnergyTwo(const Mesh &mesh, const VectorXd &q, const EdgeNormalDerivatives &nderivs, const PrecomputedFaceQuantities &data, int faceid, const ElasticParameters &params, double prefactor, Eigen::VectorXd &result)
{
    const Vector4d I(1.0, 0, 0, 1.0);
    const Vector4d &gbarinv = data.gbarinv;
    const Vector4d &gmat = data.g;
    const Vector4d &bmat = data.b;
    const Vector4d &cmat = data.c;

    double A = data.intrinsicArea;
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

    Dg(mesh, q, faceid, prefactor*totalgprefactor, result);

    Vector4d totalbprefactor = coeff3*A*4.0*prefactor3c
            ;

    Db(mesh, q, nderivs, faceid, prefactor*totalbprefactor, result);

    Vector4d totalcprefactor = coeff3*A*4.0*prefactor3a;

    Dc(mesh, q, nderivs, faceid, prefactor*totalcprefactor, result);
}

double Midedge::elasticEnergy(const Mesh &mesh, const VectorXd &q, const ElasticParameters &params, VectorXd *derivs, VectorXd *energies)
{
    double total = 0;
    if(derivs)
    {
        if(derivs->size() != 3*mesh.numVertices())
            derivs->resize(3*mesh.numVertices());
        derivs->setZero();
    }
    if(energies)
        energies->resize(mesh.numFaces());
    EdgeNormalDerivatives nderivs;
    precomputeEdgeNormalDerivatives(mesh, q, nderivs);

    const Vector4d I(1.0, 0.0, 0.0, 1.0);

    int nfaces = mesh.numFaces();
    #pragma omp parallel for
    for(int i=0; i<nfaces; i++)
    {        
        PrecomputedFaceQuantities data;

        data.gbarinv = matInv(mesh.faceMetric(i));
        data.b = b(mesh, q, nderivs, i);
        data.c = c(mesh, q, nderivs, i);
        data.g = g(mesh, q, i);
        data.intrinsicArea = intrinsicArea(mesh, i, params);
        data.trgminusI = trace(matMult(data.gbarinv, data.g) - I);
        data.trb = trace(matMult(data.gbarinv, data.b));
        data.trc = trace(matMult(data.gbarinv, data.c));

        double result = params.YoungsModulus/8.0*params.PoissonRatio/(1.0-params.PoissonRatio*params.PoissonRatio)*elasticEnergyOne(mesh, q, data, i, params);
        result += params.YoungsModulus/8.0/(1.0+params.PoissonRatio)*elasticEnergyTwo(mesh, q, data, i, params);

        if(derivs)
        {
            DelasticEnergyOne(mesh, q, nderivs, data, i, params, params.YoungsModulus/8.0*params.PoissonRatio/(1.0-params.PoissonRatio*params.PoissonRatio), *derivs);
            DelasticEnergyTwo(mesh, q, nderivs, data, i, params, params.YoungsModulus/8.0/(1.0+params.PoissonRatio), *derivs);
        }

        #pragma omp atomic
        total += result;

        if(energies)
            (*energies)[i] = result;
    }

    return total;
}

double Midedge::intrinsicArea(const Mesh &mesh, int faceid, const ElasticParameters &params)
{
    double gdet = det(mesh.faceMetric(faceid));
    assert(gdet >= 0);
    return 0.5*params.scale*params.scale*sqrt(gdet);
}

void Midedge::precomputeEdgeNormalDerivatives(const Mesh &mesh, const VectorXd &q, EdgeNormalDerivatives &dnormals)
{
    int nentries = 3*mesh.numFaces();
    dnormals.edgeNormalsPartials.resize(nentries);
    dnormals.edgeNormalsPartialIndices.resize(nentries);
    dnormals.normals.resize(nentries);

    int nedges = mesh.numEdges();
    #pragma omp parallel for
    for(int i=0; i<nedges; i++)
    {
        Vector4i hinge = mesh.buildHinge(i);
        if(mesh.isBoundaryEdge(i))
        {
            int face = mesh.edgeFaces(i)[0];
            int offset = -1;
            for(int j=0; j<3; j++)
            {
                if(mesh.edgeVerts(i)[1] == mesh.faceVerts(face)[(j+1)%3])
                    offset = j;
            }
            assert(offset != -1);
            Vector3d qs[3];
            for(int j=0; j<3; j++)
                qs[j] = q.segment<3>(3*hinge[j]);
            DfaceNormal(qs[0], qs[1], qs[2], dnormals.edgeNormalsPartials[3*face+offset]);
            dnormals.edgeNormalsPartialIndices[3*face+offset].push_back(hinge[0]);
            dnormals.edgeNormalsPartialIndices[3*face+offset].push_back(hinge[1]);
            dnormals.edgeNormalsPartialIndices[3*face+offset].push_back(hinge[2]);
            dnormals.normals[3*face+offset] = faceNormal(qs[0], qs[1], qs[2]);
        }
        else
        {
            for(int side=0; side<2; side++)
            {
                int face = mesh.edgeFaces(i)[side];
                int offset = -1;
                for(int j=0; j<3; j++)
                {
                    if(mesh.edgeVerts(i)[1-side] == mesh.faceVerts(face)[(j+1)%3])
                        offset = j;
                }
                assert(offset != -1);
                Vector3d qs[4];
                for(int j=0; j<4; j++)
                    qs[j] = q.segment<3>(3*hinge[j]);
                DdiamondEdgeNormal(qs[0], qs[1], qs[2], qs[3], dnormals.edgeNormalsPartials[3*face+offset]);
                for(int j=0; j<4; j++)
                    dnormals.edgeNormalsPartialIndices[3*face+offset].push_back(hinge[j]);
                dnormals.normals[3*face+offset] =  diamondEdgeNormal(qs[0], qs[1], qs[2], qs[3]);
            }
        }
    }    
}

void Midedge::gaussianCurvature(const Mesh &mesh, const VectorXd &q, VectorXd &Ks)
{
    int nfaces = mesh.numFaces();
    Ks.resize(nfaces);
    EdgeNormalDerivatives nderivs;
    precomputeEdgeNormalDerivatives(mesh, q, nderivs);
    for(int i=0; i<nfaces; i++)
    {
        Ks[i] = K(mesh, q, nderivs, i);
    }
}

void Midedge::meanCurvature(const Mesh &mesh, const VectorXd &q, VectorXd &Hs)
{
    int nfaces = mesh.numFaces();
    Hs.resize(nfaces);
    EdgeNormalDerivatives nderivs;
    precomputeEdgeNormalDerivatives(mesh, q, nderivs);
    for(int i=0; i<nfaces; i++)
    {
        Hs[i] = H(mesh, q, nderivs, i);
    }
}
