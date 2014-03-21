#ifndef OMTYPES_H
#define OMTYPES_H

#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"

struct MyTraits : public OpenMesh::DefaultTraits
{
    typedef OpenMesh::Vec3d Point; // use double-values points
    typedef OpenMesh::Vec3d Normal; // use double-values points

    EdgeTraits
    {
    private:
        double restlen_;
        double targetlen_;
    public:
        EdgeT() : restlen_(0) {}
        double restlen() const {return restlen_;}
        void setRestlen(double l) {restlen_=l;}

        double targetlen() const {return targetlen_;}
        void setTargetlen(double l) {targetlen_=l;}
    };
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> OMMesh;


#endif // OMTYPES_H
