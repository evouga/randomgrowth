#ifndef ELASTICENERGY_H
#define ELASTICENERGY_H

#include <Eigen/Core>
#include <Eigen/Sparse>

const double PI = 3.14159265359;

typedef Eigen::Triplet<double> Tr;



struct ElasticParameters
{
    double h;
    double PoissonRatio;
    double YoungsModulus;
};


class ElasticEnergy
{
public:
    enum DerivativesRequested {DR_NONE  = 0,
                               DR_DQ    = 1,
                               DR_HQ    = 2,
                               DR_DGDQ  = 4};

    static double stretchingEnergy(const Eigen::VectorXd &qs, const Eigen::VectorXd &gs,
                                   int *qidx, int *gidx,
                                   Eigen::VectorXd &dq,
                                   std::vector<Tr> &hq,
                                   std::vector<Tr> &dgdq,
                                   const ElasticParameters &params,
                                   int derivsRequested);

    static double bendingEnergy(const Eigen::VectorXd &qs, const Eigen::VectorXd &gs,
                          int centqidx, const std::vector<int> &nbqidx,
                          const std::vector<int> &spokegidx, const std::vector<int> &oppgidx,
                          Eigen::VectorXd &dq,
                          std::vector<Tr> &hq,
                          std::vector<Tr> &dgdq,
                          const ElasticParameters &params,
                          int derivsRequested);

private:

    static double stretchOne(const Eigen::VectorXd &qs, const Eigen::VectorXd &gs,
                      int *qidx, int *gidx,
                      Eigen::VectorXd &dq,
                      std::vector<Tr> &hq,
                      std::vector<Tr> &dgdq,
                      const ElasticParameters &params,
                      int derivsRequested);

    static double stretchTwo(const Eigen::VectorXd &qs, const Eigen::VectorXd &gs,
                      int *qidx, int *gidx,
                      Eigen::VectorXd &dq,
                      std::vector<Tr> &hq,
                      std::vector<Tr> &dgdq,
                      const ElasticParameters &params,
                      int derivsRequested);

    static double bendOne(const Eigen::VectorXd &qs, const Eigen::VectorXd &gs,
                          int centqidx, const std::vector<int> &nbqidx,
                          const std::vector<int> &spokegidx, const std::vector<int> &oppgidx,
                          Eigen::VectorXd &dq,
                          std::vector<Tr> &hq,
                          std::vector<Tr> &dgdq,
                          const ElasticParameters &params,
                          int derivsRequested);

    static double bendTwo(const Eigen::VectorXd &qs, const Eigen::VectorXd &gs,
                          int centqidx, const std::vector<int> &nbqidx,
                          const std::vector<int> &spokegidx, const std::vector<int> &oppgidx,
                          Eigen::VectorXd &dq,
                          std::vector<Tr> &hq,
                          std::vector<Tr> &dgdq,
                          const ElasticParameters &params,
                          int derivsRequested);
};

#endif // ELASTICENERGY_H
