#include "elasticenergy.h"

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
