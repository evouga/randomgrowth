#include "elasticenergy.h"
#include <iostream>

using namespace Eigen;

double ElasticEnergy::stretchingEnergy(const VectorXd &qs,
                                       const VectorXd &gs,
                                       int *qidx,
                                       int *gidx,
                                       VectorXd &dq,
                                       std::vector<Tr> &hq,
                                       std::vector<Tr> &dgdq,
                                       const ElasticParameters &params,
                                       int derivsRequested)
{
    return stretchOne(qs, gs, qidx, gidx, dq, hq, dgdq, params, derivsRequested)
            + stretchTwo(qs, gs, qidx, gidx, dq, hq, dgdq, params, derivsRequested);
}

double ElasticEnergy::stretchOne(const VectorXd &qs,
                                 const VectorXd &gs,
                                 int *qidx,
                                 int *gidx,
                                 VectorXd &dq,
                                 std::vector<Tr> &hq,
                                 std::vector<Tr> &dgdq,
                                 const ElasticParameters &params,
                                 int derivsRequested)
{
    // Y/(8(1+v)(1-v)) sqrt(det g) tr(g^-1 a - I)^2

    double coeff = params.YoungsModulus/8.0/(1.0-params.PoissonRatio*params.PoissonRatio);

    double g[3];
    Vector3d q[3];
    for(int i=0; i<3; i++)
    {
        g[i] = gs[gidx[i]];
        q[i] = qs.segment<3>(3*qidx[i]);
    }

    double detg = g[0]*g[0]*g[1]*g[1] - (g[0]*g[0]+g[1]*g[1]-g[2]*g[2])*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2])/4.0;

    double A = (q[2]-q[1]).dot(2.0*g[1]*g[1]*(q[2]-q[1])+(g[2]*g[2]-g[0]*g[0]-g[1]*g[1])*(q[2]-q[0]))
            +(q[2]-q[0]).dot(2.0*g[0]*g[0]*(q[2]-q[0])+(g[2]*g[2]-g[0]*g[0]-g[1]*g[1])*(q[2]-q[1]))
            -4.0*detg;

    if(derivsRequested != DR_NONE)
    {
        Vector3d dAdq[3];
        dAdq[0] = 2.0*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2])*(q[2]-q[1])-4.0*g[0]*g[0]*(q[2]-q[0]);
        dAdq[1] = 2.0*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2])*(q[2]-q[0])-4.0*g[1]*g[1]*(q[2]-q[1]);
        dAdq[2] = 2.0*(g[1]*g[1]+g[2]*g[2]-g[0]*g[0])*(q[2]-q[1])+2.0*(g[0]*g[0]+g[2]*g[2]-g[1]*g[1])*(q[2]-q[0]);

        if(derivsRequested & DR_DQ)
        {
            double dqcoeff = coeff*A/(4.0*detg*sqrt(detg));

            for(int i=0; i<3; i++)
                dq.segment<3>(3*qidx[i]) += dqcoeff*dAdq[i];
        }

        if(derivsRequested & DR_HQ)
        {
            for(int i=0; i<3; i++)
            {
                Matrix3d dqidqi;
                dqidqi.setIdentity();
                dqidqi *= A*g[i]*g[i];
                dqidqi += 0.25*(dAdq[i]*dAdq[i].transpose());
                dqidqi *= coeff/(detg*sqrt(detg));

                for(int j=0; j<3; j++)
                    for(int k=0; k<3; k++)
                    {
                        hq.push_back(Tr(3*qidx[i]+j, 3*qidx[i]+k, dqidqi(k,j)));
                    }
            }

            Matrix3d dq0dq1;
            dq0dq1.setIdentity();
            dq0dq1 *= 0.5*(g[2]*g[2]-g[0]*g[0]-g[1]*g[1])*A;
            dq0dq1 += 0.25*dAdq[1]*dAdq[0].transpose();
            dq0dq1 *= coeff/(detg*sqrt(detg));
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    hq.push_back(Tr(3*qidx[0]+j, 3*qidx[1]+k, dq0dq1(k,j)));
                    hq.push_back(Tr(3*qidx[1]+k, 3*qidx[0]+j, dq0dq1(k,j)));
                }
            }

            for(int i=0; i<2; i++)
            {
                Matrix3d dq2dqi;
                dq2dqi.setIdentity();
                dq2dqi *= 0.5*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2]-2.0*g[i]*g[i])*A;
                dq2dqi += 0.25*dAdq[i]*dAdq[2].transpose();
                dq2dqi *= coeff/(detg*sqrt(detg));
                for(int j=0; j<3; j++)
                    for(int k=0; k<3; k++)
                    {
                        hq.push_back(Tr(3*qidx[2]+j, 3*qidx[i]+k, dq2dqi(k,j)));
                        hq.push_back(Tr(3*qidx[i]+k, 3*qidx[2]+j, dq2dqi(k,j)));
                    }
            }
        }

        if(derivsRequested & DR_DGDQ)
        {
            Vector3d dgdetg(g[0]*(g[1]*g[1]+g[2]*g[2]-g[0]*g[0]),
                    g[1]*(g[0]*g[0]+g[2]*g[2]-g[1]*g[1]),
                    g[2]*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2]));

            Vector3d dgA = -4.0*dgdetg;
            dgA[0] += 4.0*(g[0]*(q[0]-q[2]).dot(q[0]-q[1]));
            dgA[1] += 4.0*(g[1]*(q[1]-q[2]).dot(q[1]-q[0]));
            dgA[2] += 4.0*(g[2]*(q[2]-q[1]).dot(q[2]-q[0]));

            Matrix3d dgdqA[3];
            dgdqA[0].col(0) = 4.0*g[0]*(q[0]-q[2])+4.0*g[0]*(q[0]-q[1]);
            dgdqA[0].col(1) = 4.0*g[1]*(q[2]-q[1]);
            dgdqA[0].col(2) = 4.0*g[2]*(q[1]-q[2]);

            dgdqA[1].col(0) = 4.0*g[0]*(q[2]-q[0]);
            dgdqA[1].col(1) = 4.0*g[1]*(q[1]-q[0])+4.0*g[1]*(q[1]-q[2]);
            dgdqA[1].col(2) = 4.0*g[2]*(q[0]-q[2]);

            dgdqA[2].col(0) = 4.0*g[0]*(q[1]-q[0]);
            dgdqA[2].col(1) = 4.0*g[1]*(q[0]-q[1]);
            dgdqA[2].col(2) = 4.0*g[2]*(q[2]-q[0]) + 4.0*g[2]*(q[2]-q[1]);

            for(int i=0; i<3; i++)
            {
                Matrix3d dgdqi = A/(4.0*detg*sqrt(detg)) * dgdqA[i];
                dgdqi += 1.0/(4.0*detg*sqrt(detg)) * dAdq[i]*dgA.transpose();
                dgdqi -= 3.0*A/(8.0 * detg * detg * sqrt(detg)) * dAdq[i]*dgdetg.transpose();
                dgdqi *= coeff;

                for(int j=0; j<3; j++)
                    for(int k=0; k<3; k++)
                    {
                        dgdq.push_back(Tr(gidx[j], 3*qidx[i]+k, dgdqi(k,j)));
                    }
            }
        }
    }

    return coeff*0.5*A*A/(4.0*detg*sqrt(detg));
}

double ElasticEnergy::stretchTwo(const VectorXd &qs,
                                 const VectorXd &gs,
                                 int *qidx,
                                 int *gidx,
                                 VectorXd &dq,
                                 std::vector<Tr> &hq,
                                 std::vector<Tr> &dgdq,
                                 const ElasticParameters &params,
                                 int derivsRequested)
{
    // -Y/(8(1+v)) sqrt(det g) det(g^-1 a - I)

    double coeff = -params.YoungsModulus/(8.0*(1.0+params.PoissonRatio));

    double g[3];
    Vector3d q[3];
    for(int i=0; i<3; i++)
    {
        g[i] = gs[gidx[i]];
        q[i] = qs.segment<3>(3*qidx[i]);
    }

    double detg = g[0]*g[0]*g[1]*g[1] - (g[0]*g[0]+g[1]*g[1]-g[2]*g[2])*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2])/4.0;

    double A = (q[2]-q[1]).dot(q[2]-q[1])*(q[2]-q[0]).dot(q[2]-q[0])-(q[2]-q[1]).dot(q[2]-q[0])*(q[2]-q[1]).dot(q[2]-q[0]);
    double B = (g[0]*g[0]+g[1]*g[1]-g[2]*g[2])*(q[2]-q[1]).dot(q[2]-q[0]) - g[0]*g[0]*(q[2]-q[0]).dot(q[2]-q[0]) - g[1]*g[1]*(q[2]-q[1]).dot(q[2]-q[1]);

    if(derivsRequested != DR_NONE)
    {
        Vector3d dqA[3];
        dqA[0] = 2.0*(q[2]-q[1]).dot(q[2]-q[0])*(q[2]-q[1]) - 2.0*(q[2]-q[1]).dot(q[2]-q[1])*(q[2]-q[0]);
        dqA[1] = 2.0*(q[2]-q[1]).dot(q[2]-q[0])*(q[2]-q[0]) - 2.0*(q[2]-q[0]).dot(q[2]-q[0])*(q[2]-q[1]);
        dqA[2] = 2.0*(q[2]-q[0]).dot(q[1]-q[0])*(q[2]-q[1]) - 2.0*(q[2]-q[1]).dot(q[1]-q[0])*(q[2]-q[0]);

        Vector3d dqB[3];
        dqB[0] = (g[2]*g[2]-g[0]*g[0]-g[1]*g[1])*(q[2]-q[1]) + 2.0*g[0]*g[0]*(q[2]-q[0]);
        dqB[1] = (g[2]*g[2]-g[0]*g[0]-g[1]*g[1])*(q[2]-q[0]) + 2.0*g[1]*g[1]*(q[2]-q[1]);
        dqB[2] = (g[1]*g[1]-g[2]*g[2]-g[0]*g[0])*(q[2]-q[0]) + (g[0]*g[0]-g[1]*g[1]-g[2]*g[2])*(q[2]-q[1]);

        if(derivsRequested & DR_DQ)
        {
            for(int i=0; i<3; i++)
                dq.segment<3>(3*qidx[i]) += coeff*(dqA[i]+dqB[i])/(sqrt(detg));
        }

        if(derivsRequested & DR_HQ)
        {

            Matrix3d dqdq0A[3];
            dqdq0A[0].setIdentity();
            dqdq0A[0] *= 2.0*(q[2]-q[1]).dot(q[2]-q[1]);
            dqdq0A[0] -= 2.0*(q[2]-q[1])*(q[2]-q[1]).transpose();

            dqdq0A[1].setIdentity();
            dqdq0A[1] *= -2.0*(q[2]-q[1]).dot(q[2]-q[0]);
            dqdq0A[1] += 4.0*(q[2]-q[0])*(q[2]-q[1]).transpose()-2.0*(q[2]-q[1])*(q[2]-q[0]).transpose();

            dqdq0A[2].setIdentity();
            dqdq0A[2] *= 2.0*(q[2]-q[1]).dot(q[1]-q[0]);
            dqdq0A[2] += 2.0*(q[2]-q[1])*(q[2]-q[1]).transpose();
            dqdq0A[2] += 2.0*(q[2]-q[1])*(q[2]-q[0]).transpose();
            dqdq0A[2] -= 4.0*(q[2]-q[0])*(q[2]-q[1]).transpose();

            Matrix3d dqdq1A[3];
            dqdq1A[0].setIdentity();
            dqdq1A[0] *= -2.0*(q[2]-q[1]).dot(q[2]-q[0]);
            dqdq1A[0] += 4.0*(q[2]-q[1])*(q[2]-q[0]).transpose();
            dqdq1A[0] += -2.0*(q[2]-q[0])*(q[2]-q[1]).transpose();

            dqdq1A[1].setIdentity();
            dqdq1A[1] *= 2.0*(q[2]-q[0]).dot(q[2]-q[0]);
            dqdq1A[1] += -2.0*(q[2]-q[0])*(q[2]-q[0]).transpose();

            dqdq1A[2].setIdentity();
            dqdq1A[2] *= 2.0*(q[2]-q[0]).dot(q[0]-q[1]);
            dqdq1A[2] += 2.0*(q[2]-q[0])*(q[2]-q[0]).transpose();
            dqdq1A[2] += 2.0*(q[2]-q[0])*(q[2]-q[1]).transpose();
            dqdq1A[2] += -4.0*(q[2]-q[1])*(q[2]-q[0]).transpose();

            Matrix3d dqdq2A[3];
            dqdq2A[0].setIdentity();
            dqdq2A[0] *= 2.0*(q[2]-q[1]).dot(q[1]-q[0]);
            dqdq2A[0] += -2.0*(q[2]-q[1])*(q[2]-q[0]).transpose();
            dqdq2A[0] += -2.0*(q[2]-q[1])*(q[1]-q[0]).transpose();
            dqdq2A[0] += 2.0*(q[2]-q[0])*(q[2]-q[1]).transpose();

            dqdq2A[1].setIdentity();
            dqdq2A[1] *= -2.0*(q[2]-q[0]).dot(q[1]-q[0]);
            dqdq2A[1] += 2.0*(q[2]-q[1])*(q[2]-q[0]).transpose();
            dqdq2A[1] += -2.0*(q[2]-q[0])*(q[2]-q[1]).transpose();
            dqdq2A[1] += 2.0*(q[2]-q[0])*(q[1]-q[0]).transpose();

            dqdq2A[2].setIdentity();
            dqdq2A[2] *= 2.0*(q[1]-q[0]).dot(q[1]-q[0]);
            dqdq2A[2] += -2.0*(q[1]-q[0])*(q[1]-q[0]).transpose();

            Matrix3d dqdq0B[3];
            dqdq0B[0].setIdentity();
            dqdq0B[0] *= -2.0*g[0]*g[0];

            dqdq0B[1].setIdentity();
            dqdq0B[1] *= (g[0]*g[0]+g[1]*g[1]-g[2]*g[2]);

            dqdq0B[2].setIdentity();
            dqdq0B[2] *= (g[0]*g[0]+g[2]*g[2]-g[1]*g[1]);

            Matrix3d dqdq1B[3];
            dqdq1B[0].setIdentity();
            dqdq1B[0] *= (g[0]*g[0]+g[1]*g[1]-g[2]*g[2]);

            dqdq1B[1].setIdentity();
            dqdq1B[1] *= -2.0*g[1]*g[1];

            dqdq1B[2].setIdentity();
            dqdq1B[2] *= (g[1]*g[1]+g[2]*g[2]-g[0]*g[0]);

            Matrix3d dqdq2B[3];
            dqdq2B[0].setIdentity();
            dqdq2B[0] *= (g[0]*g[0]+g[2]*g[2]-g[1]*g[1]);

            dqdq2B[1].setIdentity();
            dqdq2B[1] *= (g[1]*g[1]+g[2]*g[2]-g[0]*g[0]);

            dqdq2B[2].setIdentity();
            dqdq2B[2] *= -2.0*g[2]*g[2];

            for(int i=0; i<3; i++)
            {
                Matrix3d dqdq0E = coeff*(dqdq0A[i]+dqdq0B[i])/(sqrt(detg));
                Matrix3d dqdq1E = coeff*(dqdq1A[i]+dqdq1B[i])/(sqrt(detg));
                Matrix3d dqdq2E = coeff*(dqdq2A[i]+dqdq2B[i])/(sqrt(detg));
                for(int j=0; j<3; j++)
                {
                    for(int k=0; k<3; k++)
                    {
                        hq.push_back(Tr(3*qidx[i]+j, 3*qidx[0]+k, dqdq0E(k,j)));
                        hq.push_back(Tr(3*qidx[i]+j, 3*qidx[1]+k, dqdq1E(k,j)));
                        hq.push_back(Tr(3*qidx[i]+j, 3*qidx[2]+k, dqdq2E(k,j)));
                    }
                }
            }
        }

        if(derivsRequested & DR_DGDQ)
        {
            Matrix3d dgdqB[3];

            dgdqB[0].col(0) = 2.0*g[0]*(q[2]-q[0]) + 2.0*g[0]*(q[1]-q[0]);
            dgdqB[0].col(1) = -2.0*g[1]*(q[2]-q[1]);
            dgdqB[0].col(2) = 2.0*g[2]*(q[2]-q[1]);

            dgdqB[1].col(0) = -2.0*g[0]*(q[2]-q[0]);
            dgdqB[1].col(1) = 2.0*g[1]*(q[2]-q[1]) + 2.0*g[1]*(q[0]-q[1]);
            dgdqB[1].col(2) = 2.0*g[2]*(q[2]-q[0]);

            dgdqB[2].col(0) = -2.0*g[0]*(q[1]-q[0]);
            dgdqB[2].col(1) = 2.0*g[1]*(q[1]-q[0]);
            dgdqB[2].col(2) = 2.0*g[2]*(q[0]-q[2]) + 2.0*g[2]*(q[1]-q[2]);

            Vector3d dgdetg(g[0]*(g[1]*g[1]+g[2]*g[2]-g[0]*g[0]),
                    g[1]*(g[0]*g[0]+g[2]*g[2]-g[1]*g[1]),
                    g[2]*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2]));


            for(int i=0; i<3; i++)
            {
                Matrix3d dgdqi = -1.0/(2.0*detg*sqrt(detg))*(dqA[i]+dqB[i])*dgdetg.transpose();
                dgdqi += 1.0/(sqrt(detg)) * dgdqB[i];
                dgdqi *= coeff;
                for(int j=0; j<3; j++)
                {
                    for(int k=0; k<3; k++)
                    {
                        dgdq.push_back(Tr(gidx[j], 3*qidx[i]+k, dgdqi(k,j)));
                    }
                }
            }
        }
    }

    return coeff*sqrt(detg) + coeff*(A+B)/sqrt(detg);
}

double ElasticEnergy::bendOne(const VectorXd &qs, const VectorXd &gs, int centqidx, const std::vector<int> &nbqidx, const std::vector<int> &spokegidx, const std::vector<int> &oppgidx, VectorXd &dq, std::vector<Tr> &hq, std::vector<Tr> &dgdq, const ElasticParameters &params, int derivsRequested)
{
    int numnbs = nbqidx.size();
    assert((int)spokegidx.size() == numnbs);
    assert((int)oppgidx.size() == numnbs);

    double coeff = params.YoungsModulus*params.h*params.h/(24.0*(1.0-params.PoissonRatio*params.PoissonRatio));

    std::vector<double> A;
    std::vector<double> thetas;
    std::vector<double> psis;
    for(int i=0; i<numnbs; i++)
    {
        int nextid = (i+1)%numnbs;
        int previd = (i+numnbs-1)%numnbs;
        double gi = gs[spokegidx[i]];
        double gn = gs[spokegidx[nextid]];
        double gin = gs[oppgidx[i]];
        double gp = gs[spokegidx[previd]];
        double gpi = gs[oppgidx[previd]];
        double area = 0.5*sqrt(gi*gi*gn*gn - (gi*gi+gn*gn-gin*gin)*(gi*gi+gn*gn-gin*gin)/4.0);
        A.push_back(area);

        double thetanum = gp*gp+gpi*gpi-gi*gi;
        double thetadenom = sqrt(4.0*gp*gp*gpi*gpi - (gp*gp+gpi*gpi-gi*gi)*(gp*gp+gpi*gpi-gi*gi));

        thetas.push_back(thetanum/thetadenom);

        double psinum = gn*gn+gin*gin-gi*gi;
        double psidenom = sqrt(4.0*gn*gn*gin*gin - (gn*gn+gin*gin-gi*gi)*(gn*gn+gin*gin-gi*gi));

        psis.push_back(psinum/psidenom);
    }

    double rtdetg = 0;
    for(int i=0; i<numnbs; i++)
        rtdetg += A[i];
    rtdetg /= 3.0;

    Vector3d B(0,0,0);
    for(int i=0; i<numnbs; i++)
    {
        B += 0.5*(thetas[i]+psis[i])*(qs.segment<3>(3*nbqidx[i])-qs.segment<3>(3*centqidx));
    }

    if(derivsRequested != DR_NONE)
    {
        double cotsums = 0;
        for(int i=0; i<numnbs; i++)
        {
            if(derivsRequested & DR_DQ)
                dq.segment<3>(3*nbqidx[i]) += coeff*(thetas[i]+psis[i])/rtdetg * B;
            cotsums += thetas[i]+psis[i];
        }

        if(derivsRequested & DR_DQ)
            dq.segment<3>(3*centqidx) += -coeff*cotsums/rtdetg*B;

        if(derivsRequested & DR_HQ)
        {
            for(int i=0; i<numnbs; i++)
            {
                for(int j=0; j<numnbs; j++)
                {
                    double hess = coeff*(thetas[i]+psis[i])*(thetas[j]+psis[j])/(2.0*rtdetg);

                    for(int k=0; k<3; k++)
                    {
                        hq.push_back(Tr(3*nbqidx[i]+k, 3*nbqidx[j]+k, hess));
                    }
                }

                double hess = -coeff/(2.0*rtdetg)*(thetas[i]+psis[i])*cotsums;
                for(int k=0; k<3; k++)
                {
                    hq.push_back(Tr(3*nbqidx[i]+k, 3*centqidx+k, hess));
                    hq.push_back(Tr(3*centqidx+k, 3*nbqidx[i]+k, hess));
                }
            }

            double hess = coeff*cotsums*cotsums/(2.0*rtdetg);

            for(int k=0; k<3; k++)
            {
                hq.push_back(Tr(3*centqidx+k, 3*centqidx+k, hess));
            }
        }

        if(derivsRequested & DR_DGDQ)
        {
            VectorXd dgspokedrtdeg(numnbs);
            VectorXd dgoppdrtdeg(numnbs);
            dgspokedrtdeg.setZero();
            dgoppdrtdeg.setZero();

            MatrixXd dgdcottheta(numnbs, numnbs);
            MatrixXd dgdcotpsi(numnbs, numnbs);
            dgdcottheta.setZero();
            dgdcotpsi.setZero();

            MatrixXd dgoppdcottheta(numnbs, numnbs);
            MatrixXd dgoppdcotpsi(numnbs, numnbs);
            dgoppdcottheta.setZero();
            dgoppdcotpsi.setZero();

            for(int i=0; i<numnbs; i++)
            {
                int nextid = (i+1)%numnbs;
                int previd = (i+numnbs-1)%numnbs;
                double gi = gs[spokegidx[i]];
                double gn = gs[spokegidx[nextid]];
                double gp = gs[spokegidx[previd]];
                double gin = gs[oppgidx[i]];
                double gpn = gs[oppgidx[previd]];

                double spoke = 1.0/24.0*(gi*(gn*gn+gin*gin-gi*gi)/A[i] + gi*(gp*gp+gpn*gpn-gi*gi)/A[previd]);
                double opp = 1.0/24.0*(gin*(gi*gi+gn*gn-gin*gin)/A[i]);
                dgspokedrtdeg[i] = spoke;
                dgoppdrtdeg[i] = opp;

                double thetadenom = 4.0*gp*gp*gpn*gpn-(gp*gp+gpn*gpn-gi*gi)*(gp*gp+gpn*gpn-gi*gi);
                double psidenom = 4.0*gn*gn*gin*gin-(gn*gn+gin*gin-gi*gi)*(gn*gn+gin*gin-gi*gi);

                dgdcottheta(previd, i) = 2.0*gp*(1.0/sqrt(thetadenom) + (gp*gp+gpn*gpn-gi*gi)*(gp*gp-gpn*gpn-gi*gi)/(thetadenom*sqrt(thetadenom)));
                dgdcottheta(i, i) = -2.0*gi*(1.0/sqrt(thetadenom) + (gp*gp+gpn*gpn-gi*gi)*(gp*gp+gpn*gpn-gi*gi)/(thetadenom*sqrt(thetadenom)));
                dgoppdcottheta(previd, i) = 2.0*gpn*(1.0/sqrt(thetadenom) + (gp*gp+gpn*gpn-gi*gi)*(gpn*gpn-gp*gp-gi*gi)/(thetadenom*sqrt(thetadenom)));

                dgdcotpsi(nextid, i) = 2.0*gn*(1.0/sqrt(psidenom) + (gn*gn+gin*gin-gi*gi)*(gn*gn-gin*gin-gi*gi)/(psidenom*sqrt(psidenom)));
                dgdcotpsi(i, i) = -2.0*gi*(1.0/sqrt(psidenom) + (gn*gn+gin*gin-gi*gi)*(gn*gn+gin*gin-gi*gi)/(psidenom*sqrt(psidenom)));
                dgoppdcotpsi(i, i) = 2.0*gin*(1.0/sqrt(psidenom) + (gn*gn+gin*gin-gi*gi)*(gin*gin-gn*gn-gi*gi)/(psidenom*sqrt(psidenom)));
            }

            MatrixXd dgdB(3,numnbs);
            dgdB.setZero();

            MatrixXd dgoppdB(3,numnbs);
            dgoppdB.setZero();

            for(int j=0; j<numnbs; j++)
            {
                dgdB += 0.5*(qs.segment<3>(3*nbqidx[j])-qs.segment<3>(3*centqidx))
                        * (dgdcottheta.col(j)+dgdcotpsi.col(j)).transpose();
                dgoppdB += 0.5*(qs.segment<3>(3*nbqidx[j])-qs.segment<3>(3*centqidx))
                        * (dgoppdcottheta.col(j)+dgoppdcotpsi.col(j)).transpose();
            }

            for(int j=0; j<numnbs; j++)
            {
                MatrixXd dgdqj = B*dgspokedrtdeg.transpose();
                dgdqj *= -1.0/(rtdetg);
                dgdqj += dgdB;
                dgdqj *= (thetas[j]+psis[j]);
                dgdqj += B*(dgdcottheta.col(j)+dgdcotpsi.col(j)).transpose();
                dgdqj *= coeff/rtdetg;

                MatrixXd dgoppdqj = B*dgoppdrtdeg.transpose();
                dgoppdqj *= -1.0/(rtdetg);
                dgoppdqj += dgoppdB;
                dgoppdqj *= (thetas[j]+psis[j]);
                dgoppdqj += B*(dgoppdcottheta.col(j)+dgoppdcotpsi.col(j)).transpose();
                dgoppdqj *= coeff/rtdetg;

                for(int i=0; i<3; i++)
                    for(int k=0; k<numnbs; k++)
                    {
                        dgdq.push_back(Tr(spokegidx[k],3*nbqidx[j]+i,dgdqj(i,k)));
                        dgdq.push_back(Tr(oppgidx[k], 3*nbqidx[j]+i,dgoppdqj(i,k)));
                    }
            }

            MatrixXd dgdqc = B*dgspokedrtdeg.transpose();
            dgdqc *= -1.0/rtdetg;
            dgdqc += dgdB;
            dgdqc *= cotsums;

            MatrixXd dgoppdqc = B*dgoppdrtdeg.transpose();
            dgoppdqc *= -1.0/rtdetg;
            dgoppdqc += dgoppdB;
            dgoppdqc *= cotsums;
            for(int j=0; j<numnbs; j++)
            {
                dgdqc += B*(dgdcottheta.col(j)+dgdcotpsi.col(j)).transpose();
                dgoppdqc += B*(dgoppdcottheta.col(j)+dgoppdcotpsi.col(j)).transpose();
            }
            dgdqc *= -1.0*coeff/rtdetg;
            dgoppdqc *= -1.0*coeff/rtdetg;
            for(int i=0; i<3; i++)
            {
                for(int k=0; k<numnbs; k++)
                {
                    dgdq.push_back(Tr(spokegidx[k], 3*centqidx+i, dgdqc(i,k)));
                    dgdq.push_back(Tr(oppgidx[k], 3*centqidx+i, dgoppdqc(i,k)));
                }
            }
        }
    }

    return coeff*B.dot(B)/rtdetg;
}
