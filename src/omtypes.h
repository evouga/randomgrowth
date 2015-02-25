#ifndef OMTYPES_H
#define OMTYPES_H

#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"

struct MyTraits : public OpenMesh::DefaultTraits
{
    typedef OpenMesh::Vec3d Point; // use double-values points
    typedef OpenMesh::Vec3d Normal; // use double-values points
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> OMMesh;


#endif // OMTYPES_H
