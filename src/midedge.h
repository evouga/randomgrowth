#ifndef MIDEDGE_H
#define MIDEDGE_H

#include <Eigen/Core>
#include <vector>
#include "omtypes.h"
#include <Eigen/Sparse>

typedef Eigen::Triplet<double> Tr;

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

struct EdgeNormalDerivatives
{
    std::vector<int> vertidx;
    std::vector<Eigen::Matrix3d> partials;
    std::vector<int> startidx;
    std::vector<Eigen::Vector3d> normals;
};

class Midedge
{
public:
    Midedge();

    static void elasticEnergies(const OMMesh &mesh, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params, Eigen::VectorXd &energies);
    static double elasticEnergy(const OMMesh &mesh, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params, Eigen::VectorXd *derivs);
    static void elasticEnergyFactor(const OMMesh &mesh, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params, std::vector<std::vector<Tr> > &mats);
    static Eigen::Vector4d inducedG(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q);
    static void visualizeNormals(const OMMesh &mesh, const Eigen::VectorXd &q, double scale);
    static void gaussianCurvature(const OMMesh &mesh, const Eigen::VectorXd &q, Eigen::VectorXd &Ks);
    static void meanCurvature(const OMMesh &mesh, const Eigen::VectorXd &q, Eigen::VectorXd &Hs);
    static double intrinsicArea(int faceid, const Eigen::VectorXd &gbar, const ElasticParameters &params);

private:
    static Eigen::Matrix3d crossMatrix(const Eigen::Vector3d &v);

    static Eigen::Vector4d matMult(const Eigen::Vector4d &m1, const Eigen::Vector4d &m2);
    static void DmatMult(const Eigen::Vector4d &m1, const Eigen::Vector4d &m2, Eigen::Matrix4d &m1partials, Eigen::Matrix4d &m2partials);

    static double H(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q);
    static void DH(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q, double prefactor, Eigen::VectorXd &partials);

    static double K(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q);
    static void DK(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q, double prefactor, Eigen::VectorXd &partials);

    static double trace(const Eigen::Vector4d &m1);
    const static Eigen::Vector4d Dtrace;

    static double det(const Eigen::Vector4d &m);
    static void Ddet(const Eigen::Vector4d &m, Eigen::Vector4d &partials);

    static Eigen::Vector4d matInv(const Eigen::Vector4d &m);
    static void DmatInv(const Eigen::Vector4d &m, Eigen::Matrix4d &partials);

    static double area(const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2);
    static void Darea(const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, Eigen::Vector3d *partials);

    static Eigen::Vector3d diamondEdgeNormal(const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3);
    static void DdiamondEdgeNormal(const Eigen::Vector3d &q0,const Eigen::Vector3d &q1,const Eigen::Vector3d &q2,const Eigen::Vector3d &q3, std::vector<Eigen::Matrix3d> &partials);

    static Eigen::Vector3d faceNormal(const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2);
    static void DfaceNormal(const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, std::vector<Eigen::Matrix3d> &partials);

    static Eigen::Vector3d g(const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3);
    static void Dg(const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3, std::vector<Eigen::Matrix3d> &partials);

    static Eigen::Vector3d b(const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3, const Eigen::Vector3d &n1, const Eigen::Vector3d &n2,const Eigen::Vector3d &n3);    
    static void Db(const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3, const Eigen::Vector3d &n1, const Eigen::Vector3d &n2,const Eigen::Vector3d &n3, Eigen::Matrix3d *partials);

    static Eigen::Vector3d edgeNormal(int edgeid, const EdgeNormalDerivatives &nderivs);
    static void DedgeNormal(int edgeid, const EdgeNormalDerivatives &dnormas, const Eigen::Vector3d &prefactor, Eigen::VectorXd &partials);

    static Eigen::Vector4d g(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q);
    static void Dg(const OMMesh &mesh, int faceid, const Eigen::VectorXd &q, const Eigen::Vector4d &prefactor, Eigen::VectorXd &partials);

    static Eigen::Vector4d b(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q);
    static void Db(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q, const Eigen::Vector4d &prefactor, Eigen::VectorXd &partials);

    static Eigen::Vector4d c(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q);
    static void Dc(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q, const Eigen::Vector4d &prefactor, Eigen::VectorXd &partials);

    static double elasticEnergyOne(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params);
    static void DelasticEnergyOne(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params, double prefactor, Eigen::VectorXd &result);

    static double elasticEnergyTwo(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params);
    static void DelasticEnergyTwo(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params, double prefactor, Eigen::VectorXd &result);    

    static void DelasticEnergyTwoFactor(const OMMesh &mesh, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::VectorXd &q, const Eigen::VectorXd &gbar1, const Eigen::VectorXd &gbar2, const ElasticParameters &params, double prefactor, std::vector<std::vector<Tr> > &mats);

    static void gatherEdgeNormalDerivatives(const OMMesh &mesh, const Eigen::VectorXd &q, EdgeNormalDerivatives &dnormals);
};

#endif // MIDEDGE_H
