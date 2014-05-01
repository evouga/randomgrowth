#ifndef MIDEDGE_H
#define MIDEDGE_H

#include <Eigen/Core>
#include <vector>
#include "omtypes.h"

struct ElasticParameters
{
    double h;
    double PoissonRatio;
    double YoungsModulus;
    double scale;

    virtual void dumpParameters(std::ostream &os)
    {
        os << "scale " << scale << std::endl;
        os << "h " << h << std::endl;
        os << "YoungsModulus " << YoungsModulus << std::endl;
        os << "PoissonRatio " << PoissonRatio << std::endl;
    }

    virtual ~ElasticParameters() {}
};

class Midedge
{
public:
    Midedge();

    static double elasticEnergy(const OMMesh &mesh, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params);
    static void DelasticEnergy(const OMMesh &mesh, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params, Eigen::VectorXd &result);

private:
    static Eigen::Matrix3d crossMatrix(const Eigen::Vector3d &v);

    static Eigen::Vector4d gbar(const OMMesh &mesh, int faceid, const Eigen::VectorXd &g);

    static Eigen::Vector4d matMult(const Eigen::Vector4d &m1, const Eigen::Vector4d &m2);
    static void DmatMult(const Eigen::Vector4d &m1, const Eigen::Vector4d &m2, Eigen::Matrix4d &m1partials, Eigen::Matrix4d &m2partials);

    static double H(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q);
    static void DH(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q, Eigen::VectorXd &partials);

    static double K(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q);
    static void DK(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q, Eigen::VectorXd &partials);

    static double trace(const Eigen::Vector4d &m1);
    static void Dtrace(Eigen::Vector4d &partials);

    static double det(const Eigen::Vector4d &m);
    static void Ddet(const Eigen::Vector4d &m, Eigen::Vector4d &partials);

    static Eigen::Vector4d matInv(const Eigen::Vector4d &m);
    static void DmatInv(const Eigen::Vector4d &m, Eigen::Matrix4d &partials);

    static double area(const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2);
    static void Darea(const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, std::vector<Eigen::Vector3d> &partials);

    static Eigen::Vector3d diamondEdgeNormal(const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3);
    static void DdiamondEdgeNormal(const Eigen::Vector3d &q0,const Eigen::Vector3d &q1,const Eigen::Vector3d &q2,const Eigen::Vector3d &q3, std::vector<Eigen::Matrix3d> &partials);

    static Eigen::Vector3d faceNormal(const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2);
    static void DfaceNormal(const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, std::vector<Eigen::Matrix3d> &partials);

    static Eigen::Vector3d g(const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3);
    static void Dg(const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3, std::vector<Eigen::Matrix3d> &partials);

    static Eigen::Vector3d b(const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3, const Eigen::Vector3d &n1, const Eigen::Vector3d &n2,const Eigen::Vector3d &n3);    
    static void Db(const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3, const Eigen::Vector3d &n1, const Eigen::Vector3d &n2,const Eigen::Vector3d &n3, std::vector<Eigen::Matrix3d> &partials);

    static Eigen::Vector3d edgeNormal(const OMMesh &mesh, int edgeid, const Eigen::VectorXd &q);
    static void DedgeNormal(const OMMesh &mesh, int edgeid, const Eigen::VectorXd &q, const Eigen::Vector3d &prefactor, Eigen::VectorXd &partials);

    static Eigen::Vector4d g(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q);
    static void Dg(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q, const Eigen::Vector4d &prefactor, Eigen::VectorXd &partials);

    static Eigen::Vector4d b(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q);
    static void Db(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q, const Eigen::Vector4d &prefactor, Eigen::VectorXd &partials);

    static Eigen::Vector4d c(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q);
    static void Dc(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q, const Eigen::Vector4d &prefactor, Eigen::VectorXd &partials);

    static double elasticEnergyOne(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params);
    static void DelasticEnergyOne(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params, Eigen::VectorXd &result);

    static double elasticEnergyTwo(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params);
    static void DelasticEnergyTwo(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params, Eigen::VectorXd &result);

    static double intrinsicArea(const OMMesh &mesh, int faceid, const Eigen::VectorXd &gbar, const ElasticParameters &params);
};

#endif // MIDEDGE_H
