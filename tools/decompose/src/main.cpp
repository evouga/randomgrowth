#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include "../../../src/omtypes.h"
#include <iostream>
#include <Eigen/Core>
#include <fstream>
#include <Eigen/Sparse>

typedef Eigen::Triplet<double> Tr;

using namespace std;
using namespace Eigen;

void decomposeIntoEigenbasis(const MatrixXd &eigenmodes, const SparseMatrix<double> &M, const VectorXd &q, VectorXd &modecoeffs)
{
    int modes = eigenmodes.cols();
    modecoeffs.resize(modes);
    for(int i=0; i<modes; i++)
    {
        VectorXd mode = eigenmodes.col(i);
        modecoeffs[i] = q.dot(M*mode);
    }
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
    if(argc != 4)
    {
        cerr << "Usage: decompose (initial .obj) (displaced .obj) (spectrum data file)" << endl;
        return -1;
    }

    OMMesh sourcemesh;
    OMMesh displacedmesh;
    OpenMesh::IO::Options opt;
    sourcemesh.request_face_normals();
    sourcemesh.request_vertex_normals();
    displacedmesh.request_face_normals();
    displacedmesh.request_vertex_normals();
    opt.set(OpenMesh::IO::Options::VertexNormal);
    if(!OpenMesh::IO::read_mesh(sourcemesh, argv[1], opt))
    {
        cerr << "Couldn't' read initial .obj file: " << argv[1] << endl;
        return -1;
    }
    if(!OpenMesh::IO::read_mesh(displacedmesh, argv[2], opt))
    {
        cerr << "Couldn't read displaced .obj file: " << argv[2] << endl;
        return -1;
    }
    sourcemesh.update_normals();
    displacedmesh.update_normals();

    if(!checkCompatible(sourcemesh, displacedmesh))
    {
        cerr << "Initial and displaced mesh do not have compatible topology!" << endl;
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

    cout << "Computing vertical displacements" << endl;
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
    for(int i=0; i<numverts; i++)
    {
        OMMesh::Point &pt = sourcemesh.point(sourcemesh.vertex_handle(i));
        pt[2] = q[i];
    }
    OpenMesh::IO::write_mesh(sourcemesh, "vertdispl.obj", opt);

    cout << "Removing harmonic components" << endl;
    removeHarmonicComponents(sourcemesh, harmonics, q);
    for(int i=0; i<numverts; i++)
    {
        OMMesh::Point &pt = sourcemesh.point(sourcemesh.vertex_handle(i));
        pt[2] = q[i];
    }
    OpenMesh::IO::write_mesh(sourcemesh, "nonharmonic.obj", opt);
    cout << "Decomposing into eigenbasis" << endl;
    int nummodes = eigenmodes.cols();
    VectorXd modecoeffs(nummodes);
    decomposeIntoEigenbasis(eigenmodes, M, q, modecoeffs);
    ofstream ofs("modecoeffs.dat");
    ofs << modecoeffs.transpose() << endl;
    cout << "Done" << endl;
    return 0;
}
