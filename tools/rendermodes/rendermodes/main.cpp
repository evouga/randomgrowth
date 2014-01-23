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


bool readSpectrumData(const char *filename, const OMMesh &mesh, VectorXd &eigenvalues, MatrixXd &eigenmodes, MatrixXd &harmonicfuncs, SparseMatrix<double> &M)
{
    ifstream ifs(filename);
    if(!ifs)
        return false;

    int numverts = mesh.n_vertices();
    int nummodes;
    ifs >> nummodes;
    eigenvalues.resize(nummodes);
    eigenmodes.resize(numverts, nummodes);
    for(int i=0; i<nummodes; i++)
    {
        ifs >> eigenvalues[i];
        for(int j=0; j<numverts; j++)
            ifs >> eigenmodes.coeffRef(j, i);
    }
    int numbdry;
    ifs >> numbdry;
    harmonicfuncs.resize(numverts, numbdry);
    for(int i=0; i<numbdry; i++)
    {
        int bdid;
        ifs >> bdid;
        OMMesh::VertexHandle bdvert = mesh.vertex_handle(bdid);
        if(!mesh.is_boundary(bdvert))
            return false;
        for(int j=0; j<numverts; j++)
            ifs >> harmonicfuncs.coeffRef(j, i);
    }

    int nummassentries;
    ifs >> nummassentries;
    if(nummassentries != numverts)
        return false;

    vector<Tr> Mentries;
    for(int i=0; i<numverts; i++)
    {
        double m;
        ifs >> m;
        Mentries.push_back(Tr(i,i,m));
    }
    M.resize(numverts, numverts);
    M.setFromTriplets(Mentries.begin(), Mentries.end());

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
    MatrixXd harmonics;
    SparseMatrix<double> M;
    if(!readSpectrumData(argv[2], sourcemesh, eigenvalues, eigenmodes, harmonics, M))
    {
        cerr << "Couldn't read spectrum data: " << argv[2] << endl;
        return -1;
    }

    int nummodes = (int)strtol(argv[3], NULL, 10);
    nummodes = min(nummodes, (int)eigenmodes.cols());

    double ampl = strtod(argv[4], NULL);
    for(int mode = 0; mode < nummodes; mode++)
    {
        VectorXd evec = eigenmodes.col(mode);
        double maxval = 0;
        for(int j=0; j<(int)sourcemesh.n_vertices(); j++)
        {
            maxval = max(maxval, fabs(evec[j]));
        }
        evec *= ampl/maxval;
        for(OMMesh::VertexIter vh = sourcemesh.vertices_begin(); vh != sourcemesh.vertices_end(); ++vh)
        {
            double z = evec[vh.handle().idx()];
            sourcemesh.point(vh.handle())[2] = z;
        }

        string outname = craftOuputFilename(argv[5], mode);
        if(OpenMesh::IO::write_mesh(sourcemesh, outname.c_str(), opt))
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
