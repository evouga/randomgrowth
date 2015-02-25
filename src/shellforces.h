#ifndef SHELLENERGY_H
#define SHELLENERGY_H

#include <vector>
#include <Eigen/Core>
#include "mesh.h"

class ShellForces
{
public:
    ShellForces(OMMesh &mesh, const Eigen::VectorXd &undefq, double stretchingStiffness, double bendingStiffness);

    void getForce(const Eigen::VectorXd &q, Eigen::VectorXd &f);

private:
    void getStencilStretchingForces(const Eigen::VectorXd &q,
                                    Eigen::Vector3d *force,
                                    int v0, int v1, int v2,
                                    Eigen::Matrix3d &Tm, int trinum);

    void precomputeMatrix(const Eigen::VectorXd &undefq,
                          int v0, int v1, int v2,
                          Eigen::Matrix3d &Tm);

    void getStretchingForce(const Eigen::VectorXd &q, Eigen::VectorXd &f);
    void getBendingForce(const Eigen::VectorXd &q, Eigen::VectorXd &f);
    double computeAngle(const Eigen::VectorXd &q, int p1, int p2, int q1, int q2);
    void ComputeDihedralAngleDerivatives(const Eigen::VectorXd &q, int p1, int p2, int q1, int q2,
                                                            Eigen::Vector3d &del_p1_theta,
                                                            Eigen::Vector3d &del_p2_theta,
                                                            Eigen::Vector3d &del_q1_theta,
                                                            Eigen::Vector3d &del_q2_theta) const;

    std::vector<Eigen::Matrix3d> precomputedMatrices_;
    std::vector<Eigen::Vector3d> undeformedLengths_;
    OMMesh &mesh_;
    std::vector<double> f0_;
    std::vector<double> e0_;
    std::vector<double> h0_;

    double stretchingStiffness_;
    double bendingStiffness_;
};

#endif // SHELLENERGY_H
