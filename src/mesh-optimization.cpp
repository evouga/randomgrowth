#include "mesh.h"
#include <iomanip>
#include "controller.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <lbfgs.h>
#include <math.h>

using namespace std;
using namespace Eigen;

const double PI = 3.1415926535898;
ofstream distances;
ofstream gnorms;

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
{
    SimulationMesh* simPtr = (SimulationMesh*)instance;
    int i;
    lbfgsfloatval_t fx = 0.0;
    VectorXd defPos(n);
    for(i = 0; i < n; i++) {
        defPos[i] = x[i];
    }
    VectorXd gradq(n);
    fx = Midedge::elasticEnergy(*simPtr, defPos, simPtr->getParameters(), &gradq, &simPtr->energies_);
    double pullMag = simPtr->getParameters().pullMag;
    gradq[6] -= pullMag;
    gradq[15] += pullMag;
    fx -= pullMag * x[6];
    fx += pullMag * x[15];
    for(i = 0; i < n; i++) {
        g[i] = gradq[i];
    }
    return fx;
}
static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
    for(int i=0; i < n; i++) {
        SimulationMesh* simPtr = (SimulationMesh*)instance;
        simPtr->deformedPosition_[i] = x[i];
    }

    if(k % 1000 == 0) {
        SimulationMesh* simPtr = (SimulationMesh*)instance;
        printf("PULL MAG = %f\n", simPtr->getParameters().pullMag);
        printf("Iteration %d:\n", k);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);

        cout << "energy: " << fx << endl;
        printf("\n");
    }
    return 0;
}


double SimulationMesh::truncatedConeVolume(double startHeight, double curHeight)
{
    double base = params_.scale*params_.scale*PI;
    double totV = base*params_.scale*startHeight/3.0;
    double topr = params_.scale*(1.0 - curHeight/startHeight);
    double topbase = PI*topr*topr;
    double topV = topbase*params_.scale*(startHeight-curHeight)/3.0;
    return totV-topV;
}

void SimulationMesh::opt()
{
    int i, ret = 0;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(3*this->numVertices());
    lbfgs_parameter_t param;
    if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
    }
    /* Initialize the variables. */
    for (i = 0;i < 3*this->numVertices();i ++) {
        x[i] = deformedPosition_[i];
    }
    /* Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&param);
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;

    /*
        Start the L-BFGS optimization; this will invoke the callback functions
        evaluate() and progress() when necessary.
     */
    ret = lbfgs(3*this->numVertices(), x, &fx, evaluate, progress, this, &param);

    cout << "L-BFGS optimization terminated with status code = " << ret << endl;
    /* Report the result. */

    for (i = 0;i < 3*this->numVertices();i ++) {
        deformedPosition_[i] = x[i];
    }

    lbfgs_free(x);
}

// currently unused code
//void SimulationMesh::sim()
//{
//    double h = params_.eulerTimestep;
//    VectorXd v(3*numVertices());
//    v.setZero();
//    int numsteps = params_.numEulerIters;
//    int numsteps = 500;
//    double pullMag = params_.pullMag;

//    for(int i=0; i<numsteps; i++)
//    {
//        std::cout << "iter " << i << std::endl;

//        deformedPosition_ += h*v;
//        VectorXd gradq;
//        std::cout << "computing energy" << std::endl;
//        double fx = Midedge::elasticEnergy(*this, deformedPosition_, params_,&gradq,&energies_);

//        SparseMatrix<double> Minv;
//        buildMetricInvMassMatrix(Minv);
//        VectorXd F = -gradq;
//        F[6] += pullMag;
//        F[15] -= pullMag;
//        std::cout << "force magnitude: " << F.norm() << std::endl;

//        fx -= pullMag * deformedPosition_[6];
//        fx += pullMag * deformedPosition_[15];

//        cout << "energy: " << fx << endl;
//        outfile << fx << "\n";

//        VectorXd dampv = h*Minv*params_.dampingCoeff*v;
//        for(int j=0; j<v.size(); j++)
//        {
//            if(fabs(dampv[j]) > fabs(v[j]))
//                dampv[j] = v[j];
//        }

//        v -= dampv;

//        v += h*Minv*F;

//       std::cout << "done" << std::endl;
//    }
//    opt();
//}

bool SimulationMesh::pull(Controller &cont)
{
    {
        string name = params_.outputDir + "/parameters";
        ofstream paramfile(name.c_str());
        params_.dumpParameters(paramfile);
    }

    cout << "we're running!" << endl;

    distances.open("../output/"+cont.meshName+"/distances.txt");
    gnorms.open("../output/"+cont.meshName+"/gnorms.txt");

    srand(1234);
    for(int i=0; i<numVertices(); i++) {
        deformedPosition_[3*i+2] += randomRange(-1e-10, 1e-10);
        deformedPosition_[3*i+2] += 0.05 * sin(PI * deformedPosition_[3*i]);
//        deformedPosition_[3*i+2] += 0.05 * (deformedPosition_[3*i] - 1) * (deformedPosition_[3*i] - 1);
    }
    resetRestMetric();

    do {

    opt();

    VectorXd gradq(3*numVertices());
    Midedge::elasticEnergy(*this, deformedPosition_, params_,&gradq,&energies_);
    double pullMag = params_.pullMag;
    gradq[6] -= pullMag;
    gradq[15] += pullMag;

    cout << "final gnorm: " << gradq.norm() << endl;

    double finalDist = sqrt(pow(deformedPosition_[6] - deformedPosition_[15], 2) +
                    pow(deformedPosition_[7] - deformedPosition_[16], 2) +
                    pow(deformedPosition_[8] - deformedPosition_[17], 2));

    cout << "final distance: " << finalDist << endl;
    distances << finalDist << endl;
    gnorms << gradq.norm() << endl;
    cont.exportOBJ("../output/"+cont.meshName+"/objs/"+cont.meshName+"_pulled_"+std::to_string(params_.pullMag)+".obj");

//    cont.importOBJ("../meshes/rectangle.obj");
    if(params_.pullMag < 2)
        params_.pullMag += 0.25;
    else
        params_.pullMag += 0.5;

    } while(params_.pullMag <= 10);

    distances.close();
    gnorms.close();

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
