#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include "../../../src/omtypes.h"
#include <iostream>
#include <Eigen/Core>
#include <fstream>
#include <Eigen/Sparse>
#include <sstream>
#include <iomanip>

typedef Eigen::Triplet<double> Tr;
using namespace Eigen;
using namespace std;

string craftOuputFilename(const char *folder, int num)
{
    stringstream ss;
    ss << folder << "/mode_" << setfill('0') << setw(8) << num << ".obj";
    return ss.str();
}


bool readSpectrumData(const char *filename, const OMMesh &mesh, VectorXd &eigenvalues, MatrixXd &eigenmodes)
{
    ifstream ifs(filename);
    if(!ifs)
        return false;

    int numverts = mesh.n_vertices();
    int nummodes;
    ifs >> nummodes;
    eigenvalues.resize(nummodes);
    eigenmodes.resize(3*numverts, nummodes);
    for(int i=0; i<nummodes; i++)
    {
        ifs >> eigenvalues[i];
    }

    for(int i=0; i<nummodes; i++)
    {
        for(int j=0; j<3*numverts; j++)
            ifs >> eigenmodes.coeffRef(j, i);
    }
    return ifs;
}

int main(int argc, char *argv[])
{
    if(argc != 6)
    {
        cerr << "Usage: rendermodes (undeformed .obj) (spectrum data file) (num modes) (max amplitude) (output folder)" << endl;
        return -1;
    }

    OMMesh sourcemesh;
    OpenMesh::IO::Options opt;
    if(!OpenMesh::IO::read_mesh(sourcemesh, argv[1], opt))
    {
        cerr << "Couldn't' read initial .obj file: " << argv[1] << endl;
        return -1;
    }

    cout << "Reading spectrum data" << endl;
    VectorXd eigenvalues;
    MatrixXd eigenmodes;
    if(!readSpectrumData(argv[2], sourcemesh, eigenvalues, eigenmodes))
    {
        cerr << "Couldn't read spectrum data: " << argv[2] << endl;
        return -1;
    }

    int nummodes = (int)strtol(argv[3], NULL, 10);
    nummodes = min(nummodes, (int)eigenmodes.cols());

    double ampl = strtod(argv[4], NULL);
    int numverts = sourcemesh.n_vertices();

    for(int mode = 0; mode < nummodes; mode++)
    {
        OMMesh thismesh = sourcemesh;
        VectorXd evec = eigenmodes.col(mode);
        double maxval = 0;
        for(int j=0; j<(int)3*numverts; j++)
        {
            maxval = max(maxval, fabs(evec[j]));
        }
        evec *= ampl/maxval;
        for(OMMesh::VertexIter vh = thismesh.vertices_begin(); vh != thismesh.vertices_end(); ++vh)
        {
            for(int j=0; j<3; j++)
                thismesh.point(vh.handle())[j] += evec[3*vh.handle().idx()+j];
        }

        string outname = craftOuputFilename(argv[5], mode);
        if(OpenMesh::IO::write_mesh(thismesh, outname.c_str(), opt))
        {
            cout << "Wrote mode " << mode << endl;
        }
        else
        {
            cerr << "Couldn't write mode " << mode << " to file " << outname << endl;
            return -1;
        }
    }
    return 0;
}
