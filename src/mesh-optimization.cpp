#include "mesh.h"
#include <iomanip>
#include "controller.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>

using namespace std;
using namespace Eigen;

const double PI = 3.1415926535898;

void SimulationMesh::rescaleConeFromAngle()
{
    params_.scale = params_.coneHeight * tan(params_.coneAngle);
}

bool SimulationMesh::crush(Controller &cont)
{
    {
        string name = params_.outputDir + "/parameters";
        ofstream paramfile(name.c_str());
        params_.dumpParameters(paramfile);
    }

    rescaleConeFromAngle();

    double coneHeight = params_.coneHeight / params_.scale;
    double endHeight = 0.5;

    setConeHeights(coneHeight);

    double h = params_.eulerTimestep;
    DynamicData dd;
    dd.v.resize(3*numVertices());
    dd.v.setZero();
    int numsteps = params_.simTime / params_.eulerTimestep;

    VectorXd startq = deformedPosition_;

    const int crushTime = params_.crushTime*numsteps;

    const double atm = 101325.0;

    dd.initialV = truncatedConeVolume(coneHeight, coneHeight);
    dd.planeZ = coneHeight;
    dd.planeVel = params_.initialVel;

    SparseMatrix<double> Minv;
    buildMetricInvMassMatrix(Minv);

    dd.iter = 0;
    dd.frameno = 0;

    if(params_.restoreCheckpoint.length())
    {
        restoreCheckpoint(dd,params_.restoreCheckpoint.c_str());
    }

    ofstream ofs("numbers.txt");

    for(; dd.iter < numsteps; dd.iter++)
    {
        if(dd.iter%(numsteps/100) == 0)
        {
            dd.frameno++;
            dumpFrame(dd);
        }

        std::cout << "iter " << dd.iter << std::endl;

        if(params_.constantVelocity)
        {
            double t = dd.iter;
            if(t > crushTime)
                t = crushTime;
            dd.planeZ = (t*endHeight + (crushTime-t)*coneHeight)/double(crushTime);
        }
        else
        {
            dd.planeZ += h*dd.planeVel/params_.scale;
        }

        //double crushAmt = (coneHeight-endHeight)*sqrt( double(i) / double(crushTime) );
        //double planeZ = coneHeight - crushAmt;

        deformedPosition_ += h*dd.v/params_.scale;

        dd.v += - h*Minv*params_.dampingCoeff*dd.v;

        VectorXd gradq;
        std::cout << "computing energy" << std::endl;
        double elastice = Midedge::elasticEnergy(*this,params_,&gradq);
        std::cout << "force magnitude: " << gradq.norm() << " plane height " << dd.planeZ << std::endl;
        VectorXd F = -gradq;

        double airpressure = 0;
        double curV = truncatedConeVolume(coneHeight, dd.planeZ);
        if(params_.pressureBehavior == params_.PB_CONSTANT)
            airpressure = (params_.constantPressureVal-1.0)*atm;
        else if(params_.pressureBehavior == params_.PB_SCALING)
        {
            airpressure = (dd.initialV/curV - 1.0)*atm;
        }
        else
        {
            double airrho = 1.225;
            double origradius = params_.scale;
            double ucone = dd.planeVel;
            double HH0 = (coneHeight-dd.planeZ)/coneHeight;
            airpressure = airrho * ucone * ucone * origradius * origradius * origradius * origradius * HH0 * HH0 * HH0 * HH0 / 16.0 / params_.holeRadius / params_.holeRadius / params_.holeRadius / params_.holeRadius;
        }
        VectorXd pressureF;
        pressureForce(airpressure, pressureF);

        if(!params_.constantVelocity)
        {
            double curradius = params_.scale * (coneHeight-dd.planeZ)/coneHeight;
            double origradius = params_.scale;
            double surfaceArea = PI*(curradius*curradius + origradius*origradius);
            double massF = surfaceArea*airpressure;
            massF -= params_.crushMass * 9.8;
            dd.planeVel += h*massF/params_.crushMass;
            cout << "Plane force " << massF << " " << surfaceArea*airpressure << " " << dd.planeVel << endl;
        }

        //std::cout << "Regular force " << F.norm() << " pressure force " << pressureF.norm() << " pressure " << pressure << " initialV " << initialV << " curV " << curV << std::endl;

        F += pressureF;

        dd.v += h*Minv*F;//
        // enforce constraints
        enforceConstraints(startq, dd.planeZ);
        std::cout << "done" << std::endl;

        ofs << dd.iter << " " << elastice << " " << pressureF.dot(dd.v)*h << " " << params_.crushMass*dd.planeZ*params_.scale*9.8 << endl;
    }

    cont.updateGL();
    return true;
}

void SimulationMesh::buildMetricInvMassMatrix(Eigen::SparseMatrix<double> &Minv) const
{
    VectorXd areas;
    metricBarycentricDualAreas(areas);
    int numdofs = 3*numVertices();
    Minv.resize(numdofs, numdofs);
    vector<Tr> entries;
    for(int vidx=0; vidx<numVertices(); vidx++)
    {
        double area = areas[vidx];
        double invmass = 1.0/area/params_.rho/params_.h;
//        if(mesh_->is_boundary(vi.handle()))
//            invmass = 0;
        for(int i=0; i<3; i++)
            entries.push_back(Tr(3*vidx+i, 3*vidx+i, invmass));
    }

    Minv.setFromTriplets(entries.begin(), entries.end());
}

void SimulationMesh::buildMetricMassMatrix(Eigen::SparseMatrix<double> &M) const
{
    VectorXd areas;
    metricBarycentricDualAreas(areas);
    int numdofs = 3*numVertices();
    M.resize(numdofs, numdofs);
    vector<Tr> entries;
    for(int vidx=0; vidx<numVertices(); vidx++)
    {
        double area = areas[vidx];
        double mass = area*params_.rho*params_.h;
        for(int i=0; i<3; i++)
            entries.push_back(Tr(3*vidx+i, 3*vidx+i, mass));
    }

    M.setFromTriplets(entries.begin(), entries.end());
}

void SimulationMesh::deformedBarycentricDualAreas(VectorXd &areas) const
{
    VectorXd faceareas(numFaces());
    for(int i=0; i<numFaces(); i++)
        faceareas[i] = deformedFaceArea(i);

    areas.resize(numVertices());
    areas.setZero();
    for(int i=0; i<numFaces(); i++)
    {
        for(int j=0; j<3; j++)
        {
            areas[faceVerts(i)[j]] += faceareas[i];
        }
    }
    areas /= 3.0;
}

void SimulationMesh::enforceConstraints(const VectorXd &startq, double planeHeight)
{
    for(int i=0; i<numVertices(); i++)
    {
        // pin boundary
        if(isBoundaryVert(i))
        {
            Vector3d startpos = startq.segment<3>(3*i);
            vertPos(i) = startpos;
        }
        // and crush
        double z = vertPos(i)[2];
        if(z > planeHeight)
        {
            vertPos(i)[2] = planeHeight;
        }
    }
}

