#include "mesh.h"
#include <vector>
#include <fstream>
#include <map>
#include <iostream>
#include "midedge.h"

using namespace Eigen;
using namespace std;

bool Mesh::isBoundaryEdge(int edge) const
{
    return edgeFaces(edge)[1] == -1;
}

Vector4i Mesh::buildHinge(int edge) const
{
    Vector4i result;
    result[0] = edgeVerts(edge)[0];
    result[2] = edgeVerts(edge)[1];
    for(int i=0; i<3; i++)
    {
        int candidate = faceVerts(edgeFaces(edge)[0])[i];
        if(candidate != result[0] && candidate != result[2])
        {
            result[1] = candidate;
            break;
        }
    }

    if(edgeFaces(edge)[1] == -1)
    {
        result[3] = -1;
    }
    else
    {
        for(int i=0; i<3; i++)
        {
            int candidate = faceVerts(edgeFaces(edge)[1])[i];
            if(candidate != result[0] && candidate != result[2])
            {
                result[3] = candidate;
                break;
            }
        }
    }

    return result;
}

bool Mesh::loadMesh(const char *filename)
{
    ifstream ifs(filename);
    if(!ifs)
        return false;

    vector<double> pos;
    vector<int> faceidx;

    while(true)
    {
        char c;
        ifs >> c;
        if(!ifs)
            break;
        if(c == '#')
        {
            ifs.ignore(std::numeric_limits<int>::max(), '\n');
        }
        if(c == 'v')
        {
            for(int i=0; i<3; i++)
            {
                double p;
                ifs >> p;
                pos.push_back(p);
            }
        }
        if(c == 'f')
        {
            for(int i=0; i<3; i++)
            {
                int vert;
                ifs >> vert;
                faceidx.push_back(vert);
            }
        }
    }

    int nverts = pos.size()/3;
    int nfaces = faceidx.size()/3;

    VectorXd verts(3*nverts);
    for(int i=0; i<3*nverts; i++)
        verts[i] = pos[i];

    Matrix3Xi faces(3, nfaces);
    for(int i=0; i<nfaces; i++)
        for(int j=0; j<3; j++)
        {
            int id = faceidx[3*i+j]-1;
            if(id < 0 || id >= nverts)
                return false;
            faces.coeffRef(j, i) = faceidx[3*i+j]-1;
        }

    return loadMesh(verts, faces);
}

bool Mesh::loadMesh(const VectorXd &deformedPositions, const Matrix3Xi &faces)
{
    deformedPosition_ = deformedPositions;
    faces_ = faces;

    std::vector<std::pair<int, int> > edgef;
    std::map<std::pair<int, int>, int> edgemap;

    for(int i=0; i<numFaces(); i++)
    {
        for(int j=0; j<3; j++)
        {
            int v1 = faceVerts(i)[j];
            int v2 = faceVerts(i)[(j+1)%3];
            pair<int, int> check(v1, v2);
            std::map<pair<int, int>, int>::iterator it = edgemap.find(check);
            if(it == edgemap.end())
            {
                pair<int, int> verts(v2,v1);
                edgemap[verts] = (int)edgef.size();
                pair<int, int> faces(i, -1);
                edgef.push_back(faces);
            }
            else
            {
                int edgeid = it->second;
                edgef[edgeid].second = i;
            }
        }
    }
    isBoundaryVert_.resize(numVertices());
    for(int i=0; i<numVertices(); i++)
        isBoundaryVert_[i] = false;

    int nedges = edgef.size();
    edgeFaces_.resize(2, nedges);
    edgeVerts_.resize(2, nedges);

    for(map<pair<int, int>, int>::iterator it = edgemap.begin(); it != edgemap.end(); ++it)
    {
        int id = it->second;
        edgeVerts_.coeffRef(0, id) = it->first.first;
        edgeVerts_.coeffRef(1, id) = it->first.second;
        edgeFaces_.coeffRef(0, id) = edgef[id].first;
        edgeFaces_.coeffRef(1, id) = edgef[id].second;
        if(edgef[id].second == -1)
        {
            isBoundaryVert_[it->first.first] = true;
            isBoundaryVert_[it->first.second] = true;
        }
    }


    resetRestMetric();
    return true;
}

void Mesh::resetRestMetric()
{
    restMetrics_.resize(4*numFaces());
    for(int i=0; i<numFaces(); i++)
        restMetrics_.segment<4>(4*i) = Midedge::g(*this, this->deformedPosition_, i);
}

bool Mesh::writeMesh(const char *filename)
{
    ofstream ofs(filename);
    if(!ofs)
        return false;

    for(int i=0; i<deformedPosition_.size()/3; i++)
    {
        ofs << "v ";
        for(int j=0; j<3; j++)
            ofs << deformedPosition_[3*i+j] << " ";
        ofs << endl;
    }

    for(int i=0; i<faces_.cols(); i++)
    {
        ofs << "f ";
        for(int j=0; j<3; j++)
            ofs << faces_.coeff(j, i)+1 << " ";
        ofs << endl;
    }
    return ofs;
}
