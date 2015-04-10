#include "mesh.h"
#include <iomanip>
#include "controller.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>

using namespace std;
using namespace Eigen;

const double PI = 3.1415926535898;

double SimulationMesh::truncatedConeVolume(double startHeight, double curHeight)
{
    double base = params_.scale*params_.scale*PI;
    double totV = base*params_.scale*startHeight/3.0;
    double topr = params_.scale*(1.0 - curHeight/startHeight);
    double topbase = PI*topr*topr;
    double topV = topbase*params_.scale*(startHeight-curHeight)/3.0;
    return totV-topV;
}

bool SimulationMesh::pull(Controller &cont)
{
    {
        string name = params_.outputDir + "/parameters";
        ofstream paramfile(name.c_str());
        params_.dumpParameters(paramfile);
    }

    double h = params_.eulerTimestep;
    VectorXd v(3*numVertices());
    v.setZero();
    int numsteps = params_.numEulerIters;
    double pullMag = params_.pullMag;

    for(int i=0; i<numVertices(); i++)
        deformedPosition_[3*i+2] += randomRange(-1e-10, 1e-10);

    for(int i=0; i<numsteps; i++)
    {        
        std::cout << "iter " << i << std::endl;

        deformedPosition_ += h*v;
        VectorXd gradq;
        std::cout << "computing energy" << std::endl;
        Midedge::elasticEnergy(*this,params_,&gradq,&energies_);

        SparseMatrix<double> Minv;
        buildMetricInvMassMatrix(Minv);
        VectorXd F = -gradq;
        F[6] += pullMag;
        F[15] -= pullMag;
        std::cout << "force magnitude: " << F.norm() << std::endl;
        //double curV = truncatedConeVolume(coneHeight, planeZ);
        //double pressure = airpressure*(initialV/curV - 1.0);
        //const double leakconst = 1e-6;
        //initialV = std::max(curV, initialV - leakconst*h*pressure);

        //std::cout << "Regular force " << F.norm() << " pressure force " << pressureF.norm() << " pressure " << pressure << " initialV " << initialV << " curV " << curV << std::endl;

        VectorXd dampv = h*Minv*params_.dampingCoeff*v;
        for(int i=0; i<v.size(); i++)
        {
            if(fabs(dampv[i]) > fabs(v[i]))
                dampv[i] = v[i];
        }

        v -= dampv;

        v += h*Minv*F;
        if(i%100 == 0)
            dumpFrame();

        std::cout << "done" << std::endl;
    }

    cont.updateGL();
    return true;
}

void SimulationMesh::metricBarycentricDualAreas(VectorXd &areas) const
{
    VectorXd faceareas(numFaces());
    for(int i=0; i<numFaces(); i++)
        faceareas[i] = Midedge::intrinsicArea(*this, i, params_);

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

double SimulationMesh::deformedFaceArea(int fidx) const
{
    Vector3d q0 = vertPos(faceVerts(fidx)[0]);
    Vector3d q1 = vertPos(faceVerts(fidx)[1]);
    Vector3d q2 = vertPos(faceVerts(fidx)[2]);

    double A = ((q1-q0).cross(q2-q0)).norm();
    return 0.5*params_.scale*params_.scale*A;
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

void SimulationMesh::pressureForce(double pressure, VectorXd &F)
{
    F.resize(3*numVertices());

    VectorXd areas(numVertices());
    deformedBarycentricDualAreas(areas);

    VectorXd normals;
    vertexNormals(normals);

    for(int v=0; v<numVertices(); v++)
    {
        F.segment<3>(3*v) = pressure*areas[v]*normals.segment<3>(3*v);
    }
}

double SimulationMesh::randomRange(double min, double max) const
{
    double val = (max-min)*double(rand())/double(RAND_MAX);
    return val+min;
}
