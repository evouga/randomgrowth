#include "midedge.h"
#include "mesh.h"
#include "simulationmesh.h"

#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 3)
        return -1;

    int iter=1;
    while(true)
    {
        stringstream ss;
        ss << argv[2];
        ss << std::setfill('0') << std::setw(8) << iter;
        ss << ".chk";
        iter++;

        DynamicData dd;
        SimulationMesh sm;
        sm.loadMesh(argv[1]);

        ProblemParameters params;
        params.scale = 0.085;
        params.h = 0.0002;
        params.YoungsModulus = 2e9;
        params.PoissonRatio = 0.33;
        params.rho = 500;
        params.dampingCoeff = 0.001;
        params.eulerTimestep = 2e-7;
        params.simTime = 0.1;
        params.pressureBehavior = params.PB_CONSTANT;
        params.constantPressureVal = 1;
        params.holeRadius = 0.0085;
        params.coneAngle = 0.25;
        params.coneHeight = 0.33;
        params.constantVelocity = 0;
        params.crushTime = 1;
        params.crushMass = 2;
        params.initialVel = -2;
        params.scale = params.coneHeight * tan(params.coneAngle);
        double coneHeight = params.coneHeight / params.scale;
        sm.setConeHeights(coneHeight);

        if(!sm.restoreCheckpoint(dd, ss.str().c_str()))
            return 0;

        double elasticenergy = Midedge::elasticEnergy(sm, params, NULL);

        double curV = sm.truncatedConeVolume(coneHeight, dd.planeZ);
        double initialV = sm.truncatedConeVolume(coneHeight,coneHeight);
        const double atm = 101325.0;
        double airpressure = (initialV/curV - 1.0)*atm;
        {
        double airrho = 1.225;
        double origradius = params.scale;
        double ucone = dd.planeVel;
        double HH0 = (coneHeight-dd.planeZ)/coneHeight;
        airpressure = airrho * ucone * ucone * origradius * origradius * origradius * origradius * HH0 * HH0 * HH0 * HH0 / 16.0 / params.holeRadius / params.holeRadius / params.holeRadius / params.holeRadius;
        }

        double curradius = params.scale * (coneHeight-dd.planeZ)/coneHeight;
        double origradius = params.scale;
        const double PI = 3.1415926535;
        double surfaceArea = PI*(curradius*curradius + origradius*origradius);
        double massF = surfaceArea*airpressure;

        double pressurework = massF*dd.planeVel*params.simTime/100.0;

        double massenergy = 0.5*params.crushMass*dd.planeVel*dd.planeVel + 9.8*params.crushMass*dd.planeZ;

        cout << elasticenergy << " " << pressurework << " " << massenergy << " " << dd.planeZ/coneHeight << endl;
    }
}
