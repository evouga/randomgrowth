#include "simulationmesh.h"
#include <GL/gl.h>
#include <iostream>

using namespace std;
using namespace Eigen;

const double PI = 3.1415926535;

void SimulationMesh::vertexNormals(VectorXd &vertNormals)
{
    vertNormals.resize(3*numVertices());
    int nfaces = numFaces();

    VectorXd faceNormals(3*nfaces);
    for(int i=0; i<nfaces; i++)
    {
        Vector3d normal = faceNormal(i);
        normal *= deformedFaceArea(i);
        faceNormals.segment<3>(3*i) = normal;
    }

    vertNormals.setZero();
    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<3; j++)
        {
            vertNormals.segment<3>(3*faceVerts(i)[j]) += faceNormals.segment<3>(3*i);
        }
    }

    for(int i=0; i<numVertices(); i++)
    {
        vertNormals.segment<3>(3*i) /= vertNormals.segment<3>(3*i).norm();
    }
}

void SimulationMesh::render()
{
    meshLock_.lock();
    {
        Eigen::VectorXd normals;
        vertexNormals(normals);
        glEnable(GL_LIGHTING);
        glEnable(GL_DITHER);

        glPolygonOffset(1.0, 1.0);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
        glEnable ( GL_COLOR_MATERIAL );

        if(params_.smoothShade)
        {
            glShadeModel(GL_SMOOTH);
        }
        else
        {
            glShadeModel(GL_FLAT);
        }

        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);

        static vector<GLfloat> colors;
        static vector<int> indices;
        static vector<GLfloat> pos;
        static vector<GLfloat> normal;

        colors.clear();
        indices.clear();
        pos.clear();
        normal.clear();

        VectorXd H;
        Midedge::meanCurvature(*this, H);
        //gaussianCurvature(q, H);
        //cout << H << std::endl;

        for(int face=0; face<faces_.cols(); face++)
        {
            for(int vert=0; vert<3; vert++)
            {
                int vertid = faceVerts(face)[vert];
                Vector3d color = colormap(H[face], 10.0);
                Vector3d pt = vertPos(vertid);
                Vector3d n = normals.segment<3>(3*vertid);
                for(int j=0; j<3; j++)
                {
                    pos.push_back(pt[j]);
                    normal.push_back(n[j]);
                    colors.push_back(color[j]);
                }
            }
        }

        glVertexPointer(3, GL_FLOAT, 0, &pos[0]);
        glNormalPointer(GL_FLOAT, 0, &normal[0]);
        glColorPointer(3, GL_FLOAT, 0, &colors[0]);

        int idx=0;

        for (int i=0; i<3*numFaces(); i++)
        {
            indices.push_back(idx++);
        }
        glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);

        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisable(GL_POLYGON_OFFSET_FILL);
        glDisable(GL_LIGHTING);

        if(params_.showWireframe)
        {
            glLineWidth(1.0);
            glBegin(GL_LINES);
            for(int i=0; i<numFaces(); i++)
            {
                for(int j=0; j<3; j++)
                {
                    glColor3f(0,0,0);
                    Vector3d pt1 = vertPos(faceVerts(i)[j]);
                    Vector3d pt2 = vertPos(faceVerts(i)[(j+1)%3]);
                    glVertex3d(pt1[0], pt1[1], pt1[2]);
                    glVertex3d(pt2[0], pt2[1], pt2[2]);
                }
            }
            glEnd();
        }
    }
    meshLock_.unlock();
}

Vector3d SimulationMesh::colormap(double val, double max) const
{
    double mapped = (max+val)/(2.0*max);
    return colormap(mapped);
}

Vector3d SimulationMesh::colormap(double val) const
{
    Vector3d hsl;
    hsl[0] = (1.0-val)*4.0*PI/3.0;
    hsl[1] = 1.0;
    hsl[2] = 0.5*val;
    return HSLtoRGB(hsl);
}

Vector3d SimulationMesh::HSLtoRGB(const Vector3d &hsl) const
{
    double chroma = (1.0 - fabs(2.0*hsl[2]-1.0))*hsl[1];
    double Hprime = hsl[0]*6.0/(2.0*PI);
    double hmod = fmod(Hprime, 2.0);
    double x = chroma*(1.0 - fabs(hmod-1.0));

    Vector3d rgb;

    if(Hprime < 1.0)
    {
        rgb[0] = chroma;
        rgb[1] = x;
        rgb[2] = 0;
    }
    else if(Hprime < 2.0)
    {
        rgb[0] = x;
        rgb[1] = chroma;
        rgb[2] = 0;
    }
    else if(Hprime < 3.0)
    {
        rgb[0] = 0;
        rgb[1] = chroma;
        rgb[2] = x;
    }
    else if(Hprime < 4.0)
    {
        rgb[0] = 0;
        rgb[1] = x;
        rgb[2] = chroma;
    }
    else if(Hprime < 5.0)
    {
        rgb[0] = x;
        rgb[1] = 0;
        rgb[2] = chroma;
    }
    else
    {
        rgb[0] = chroma;
        rgb[1] = 0;
        rgb[2] = x;
    }

    double m = hsl[2]-0.5*chroma;
    for(int i=0; i<3; i++)
        rgb[i] += m;
    return rgb;
}

Vector3d SimulationMesh::centroid()
{
    Vector3d centroid(0,0,0);

    meshLock_.lock();
    {
        int numpts = numVertices();

        for(int i=0; i<numpts; i++)
        {
            Vector3d pt = vertPos(i);
            for(int j=0; j<3; j++)
                centroid[j] += pt[j];
        }
        centroid /= numpts;
    }
    meshLock_.unlock();
    return centroid;
}

double SimulationMesh::radius()
{
    double maxradius = 0;

    meshLock_.lock();
    {
        Vector3d cent = centroid();
        int numpts = numVertices();
        for(int i=0; i<numpts; i++)
        {
            Vector3d pt = vertPos(i);
            Vector3d ept(pt[0],pt[1],pt[2]);
            double radius = (ept-cent).squaredNorm();
            if(radius > maxradius)
                maxradius = radius;
        }
    }
    meshLock_.unlock();
    return sqrt(maxradius);
}
