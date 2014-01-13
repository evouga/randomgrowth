#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include "../../../src/omtypes.h"
#include <iostream>
#include <Eigen/Core>
#include <fstream>
#include <Eigen/Sparse>
#include <sstream>
#include <iomanip>
#include <Eigen/Geometry>

typedef Eigen::Triplet<double> Tr;

using namespace std;
using namespace Eigen;

double computeArea(const OMMesh &mesh)
{
    double totarea = 0;
    for(OMMesh::ConstFaceIter cfi = mesh.faces_begin(); cfi != mesh.faces_end(); ++cfi)
    {
        OMMesh::ConstFaceVertexIter cfvi = mesh.cfv_iter(cfi.handle());
        OMMesh::Point p0 = mesh.point(cfvi.handle());
        ++cfvi;
        OMMesh::Point p1 = mesh.point(cfvi.handle());
        ++cfvi;
        OMMesh::Point p2 = mesh.point(cfvi.handle());
        Vector3d e1, e2;
        for(int j=0; j<3; j++)
        {
            e1[j] = p1[j]-p0[j];
            e2[j] = p2[j]-p0[j];
        }
        totarea += 0.5*(e1.cross(e2)).norm();
    }
    return totarea;
}

string craftInputFilename(const char *folder, int num)
{
    stringstream ss;
    ss << folder << "/frame_" << setfill('0') << setw(8) << num << ".obj";
    return ss.str();
}

string craftOuputFilename(const char *folder, int num)
{
    stringstream ss;
    ss << folder << "/modecoeffs_" << setfill('0') << setw(8) << num << ".dat";
    return ss.str();
}

void decomposeIntoEigenbasis(const MatrixXd &eigenmodes, const SparseMatrix<double> &M, const VectorXd &q, VectorXd &modecoeffs)
{
    modecoeffs = eigenmodes.transpose()*(M*q);
}

void removeHarmonicComponents(const OMMesh &mesh, const MatrixXd &harmonics, VectorXd &q)
{
    int col = 0;
    for(OMMesh::VertexIter vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
    {
        if(!mesh.is_boundary(vi.handle()))
            continue;
        VectorXd mode = harmonics.col(col);
        int vidx = vi.handle().idx();
        assert(fabs(mode[vidx]-1.0) < 1e-8);
        q -= q[vidx]*mode;
        col++;
    }
    assert(col == harmonics.cols());
}

bool checkCompatible(OMMesh &mesh1, OMMesh &mesh2)
{
    int numv = mesh1.n_vertices();
    int nume = mesh1.n_edges();
    int numf = mesh1.n_faces();

    if(numv != (int)mesh2.n_vertices())
        return false;
    if(nume != (int)mesh2.n_edges())
        return false;
    if(numf != (int)mesh2.n_faces())
        return false;

    return true;
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
    if(argc != 5)
    {
        cerr << "Usage: decompose (.obj folder) (undeformed .obj) (spectrum data file) (output folder)" << endl;
        return -1;
    }

    OMMesh sourcemesh;
    OpenMesh::IO::Options opt;
    if(!OpenMesh::IO::read_mesh(sourcemesh, argv[2], opt))
    {
        cerr << "Couldn't' read initial .obj file: " << argv[2] << endl;
        return -1;
    }

    cout << "Reading spectrum data" << endl;
    VectorXd eigenvalues;
    MatrixXd eigenmodes;
    MatrixXd harmonics;
    SparseMatrix<double> M;
    if(!readSpectrumData(argv[3], sourcemesh, eigenvalues, eigenmodes, harmonics, M))
    {
        cerr << "Couldn't read spectrum data: " << argv[3] << endl;
        return -1;
    }

    string areafilename = string(argv[4]) + "/area.dat";
    ofstream areafile(areafilename.c_str());
    if(!areafile)
    {
        cerr << "Couldn't open area output file " << areafilename << endl;
        return -1;
    }

    int curframe = 0;
    while(true)
    {
        cout << "Processing frame " << curframe << ":";
        string filename = craftInputFilename(argv[1], curframe);
        {
            ifstream teststream(filename.c_str());
            if(!teststream)
            {
                cout << "  Doesn't exist, terminating" << endl;
                return 0;
            }
        }

        OMMesh displacedmesh;
        if(!OpenMesh::IO::read_mesh(displacedmesh, filename.c_str(), opt))
        {
            cerr << "  Couldn't read mesh " << filename << endl;
            return -1;
        }
        if(!checkCompatible(sourcemesh, displacedmesh))
        {
            cerr << "  Initial and displaced mesh do not have compatible topology!" << endl;
            return -1;
        }

        cout << "  Computing vertical displacements" << endl;
        int numverts = sourcemesh.n_vertices();
        VectorXd q(numverts);
        for(int i=0; i<numverts; i++)
        {
            OMMesh::VertexHandle srcvert = sourcemesh.vertex_handle(i);
            OMMesh::VertexHandle dstvert = displacedmesh.vertex_handle(i);
            OMMesh::Point srcpt = sourcemesh.point(srcvert);
            OMMesh::Point dstpt = displacedmesh.point(dstvert);
            q[i] = dstpt[2]-srcpt[2];
        }

        cout << "  Removing harmonic components" << endl;
        removeHarmonicComponents(sourcemesh, harmonics, q);

        cout << "  Decomposing into eigenbasis" << endl;
        int nummodes = eigenmodes.cols();
        VectorXd modecoeffs(nummodes);
        decomposeIntoEigenbasis(eigenmodes, M, q, modecoeffs);
        string outname = craftOuputFilename(argv[4], curframe);
        ofstream ofs(outname.c_str());
        if(!ofs)
        {
            cerr << "  Couldn't open output file " << outname << endl;
            return -1;
        }
        ofs << modecoeffs.transpose() << endl;

        cout << "  Computing area: ";
        double area = computeArea(displacedmesh);
        cout << area << endl;
        areafile << area << endl;

        curframe++;
    }
    return 0;
}
