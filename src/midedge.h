#ifndef MIDEDGE_H
#define MIDEDGE_H

#include <Eigen/Core>
#include <vector>
#include <Eigen/Sparse>

typedef Eigen::Triplet<double> Tr;
class Mesh;

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

struct PrecomputedFaceQuantities
{
    double intrinsicArea;
    double trgminusI, trb, trc;
    Eigen::Vector4d gbarinv, g, b, c;
};

struct EdgeNormalDerivatives
{
    std::vector<std::vector<Eigen::Matrix3d> > edgeNormalsPartials;
    std::vector<std::vector<int> > edgeNormalsPartialIndices;
    std::vector<Eigen::Vector3d> normals;

    const Eigen::Vector3d &edgeNormal(int face, int vertno) const;
    void DedgeNormal(int face, int vertno, const Eigen::Vector3d &prefactor, Eigen::VectorXd &partials) const;
};

class Midedge
{
public:
    Midedge();

    static double elasticEnergy(const Mesh &mesh, const Eigen::VectorXd &q, const ElasticParameters &params, Eigen::VectorXd *derivs, Eigen::VectorXd *energies);
    static void gaussianCurvature(const Mesh &mesh, const Eigen::VectorXd &q, Eigen::VectorXd &Ks);
    static void meanCurvature(const Mesh &mesh, const Eigen::VectorXd &q, Eigen::VectorXd &Hs);
    static double intrinsicArea(const Mesh &mesh, int faceid, const ElasticParameters &params);
    static Eigen::Vector4d g(const Mesh &mesh, const Eigen::VectorXd &q, int faceid);

private:
    static Eigen::Matrix3d crossMatrix(const Eigen::Vector3d &v);

    static Eigen::Vector4d matMult(const Eigen::Vector4d &m1, const Eigen::Vector4d &m2);
    static void DmatMult(const Eigen::Vector4d &m1, const Eigen::Vector4d &m2, Eigen::Matrix4d &m1partials, Eigen::Matrix4d &m2partials);

    static double H(const Mesh &mesh, const Eigen::VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid);
    static void DH(const Mesh &mesh, const Eigen::VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid, double prefactor, Eigen::VectorXd &partials);

    static double K(const Mesh &mesh, const Eigen::VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid);
    static void DK(const Mesh &mesh, const Eigen::VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid, double prefactor, Eigen::VectorXd &partials);

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
    static void Dg(const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3, Eigen::Matrix3d &out1, Eigen::Matrix3d &out2, Eigen::Matrix3d &out3);

    static Eigen::Vector3d b(const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3, const Eigen::Vector3d &n1, const Eigen::Vector3d &n2,const Eigen::Vector3d &n3);    
    static void Db(const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3, const Eigen::Vector3d &n1, const Eigen::Vector3d &n2,const Eigen::Vector3d &n3, Eigen::Matrix3d *partials);

    static void Dg(const Mesh &mesh, const Eigen::VectorXd &q, int faceid, const Eigen::Vector4d &prefactor, Eigen::VectorXd &partials);

    static Eigen::Vector4d b(const Mesh &mesh, const Eigen::VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid);
    static void Db(const Mesh &mesh, const Eigen::VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::Vector4d &prefactor, Eigen::VectorXd &partials);

    static Eigen::Vector4d c(const Mesh &mesh, const Eigen::VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid);
    static void Dc(const Mesh &mesh, const Eigen::VectorXd &q, const EdgeNormalDerivatives &nderivs, int faceid, const Eigen::Vector4d &prefactor, Eigen::VectorXd &partials);

    static double elasticEnergyOne(const Mesh &mesh, const Eigen::VectorXd &q, const PrecomputedFaceQuantities &data, int faceid, const ElasticParameters &params);
    static void DelasticEnergyOne(const Mesh &mesh, const Eigen::VectorXd &q, const EdgeNormalDerivatives &nderivs, const PrecomputedFaceQuantities &data, int faceid, const ElasticParameters &params, double prefactor, Eigen::VectorXd &result);

    static double elasticEnergyTwo(const Mesh &mesh, const Eigen::VectorXd &q, const PrecomputedFaceQuantities &data, int faceid, const ElasticParameters &params);
    static void DelasticEnergyTwo(const Mesh &mesh, const Eigen::VectorXd &q, const EdgeNormalDerivatives &nderivs, const PrecomputedFaceQuantities &data, int faceid, const ElasticParameters &params, double prefactor, Eigen::VectorXd &result);

    static void precomputeEdgeNormalDerivatives(const Mesh &mesh, const Eigen::VectorXd &q, EdgeNormalDerivatives &data);
};

#endif // MIDEDGE_H
