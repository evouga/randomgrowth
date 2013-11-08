#include "newton.h"
#include <iostream>
#include <iomanip>

using namespace Eigen;

Newton::Newton(const NewtonObjective &no) : no_(no)
{
}

Newton::SolverStatus Newton::solve(NewtonParameters params, const VectorXd &guess, VectorXd &result)
{
    VectorXd q = guess;

    result.resize(q.size());
    result.setZero();

    int iter=0;

    for(; iter<params.maxiters; iter++)
    {
        SparseMatrix<double> hessian;
        VectorXd gradient;
        double energy = no_.getEnergyAndDerivatives(q, gradient, hessian);

        if(isnan(energy))
            return BAD_INPUT;

        if(gradient.norm() < params.tol)
            break;

        SimplicialLDLT<SparseMatrix<double> > solver;
        solver.compute(hessian);
        VectorXd searchdir = -solver.solve(gradient);


        std::cout << std::fixed << std::setprecision(8) << "Iter " << iter+1
                  << "   E " << energy
                  << "   |dq| " << gradient.norm()
                  << "   sd = " << searchdir.dot(gradient);

        double stepsize = 1.0;
        double initialenergy = energy;

        double eps = searchdir.dot(gradient) >= 0 ? params.lmfactor : 0;

        while(searchdir.dot(gradient) >= 0)
        {
            SparseMatrix<double> shift(hessian.rows(), hessian.cols());
            shift.setIdentity();
            shift *= eps;
            SparseMatrix<double> newh = hessian + shift;
            SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
            solver.compute(newh);
            searchdir = -solver.solve(gradient);
            eps *= 2;
        }
        std::cout << std::fixed << std::setprecision(8) << "   lm " << eps/2.0;

        int lsiters = 0;

        VectorXd newq;

        do
        {
            newq = q + stepsize*searchdir;

            if(++lsiters > params.lsmaxiters)
            {
                std::cout << std::endl;
                return LSITERS_EXCEEDED;
            }

            energy = no_.getEnergy(newq);
            stepsize /= 2.0;
        }
        while(energy > initialenergy || isnan(energy));

        q = newq;
        no_.showCurrentIteration(q);

        std::cout << std::fixed << std::setprecision(8) << "   h " << stepsize*2.0 << std::endl;
    }

    if(iter == params.maxiters)
        return MAXITERS_EXCEEDED;

    result = q;

    return CONVERGED;
}

std::string Newton::solverStatusMessage(SolverStatus status) const
{
    switch(status)
    {
    case CONVERGED:
        return std::string("Converged");
    case MAXITERS_EXCEEDED:
        return std::string("Maximum outer iterations exceeded");
    case LSITERS_EXCEEDED:
        return std::string("Maximum line search iterations exceeded");
    case BAD_INPUT:
        return std::string("Bad function data supplied by objective functor");
    case NONE:
    default:
        return std::string("Invalid Solver Status");
    }
}
