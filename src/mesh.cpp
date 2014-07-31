#include <OpenMesh/Core/IO/MeshIO.hh>
#include "mesh.h"
#include <iomanip>
#include <Eigen/Geometry>
#include <fstream>

using namespace Eigen;
using namespace OpenMesh;
using namespace std;

Mesh::Mesh() : meshLock_(QMutex::Recursive)
{
    params_.scale = 0.1;
    params_.h = .003;
    params_.YoungsModulus = 4e7;
    params_.PoissonRatio = 0.0;
    params_.rho = 500.0;
    params_.dampingCoeff = 1e-4;
    params_.eulerTimestep = 1e-6;
    params_.numEulerIters = 500000;

    params_.growthAmount = 50;
    params_.maxEdgeStrain = 1e8;
    params_.baseGrowthProbability = 0.5;

    params_.smoothShade = true;
    params_.showWireframe = true;
    params_.outputDir = "output";
    params_.colorCutoff = 1.0;
    params_.colorMode = ProblemParameters::CRM_STRAIN;

    mesh_ = new OMMesh();
    frameno_ = 0;
}

void ProblemParameters::dumpParameters(ostream &os)
{
    ElasticParameters::dumpParameters(os);
    os << "rho " << rho << endl;
    os << "dampingCoeff " << dampingCoeff << endl;
    os << "eulerTimestep " << eulerTimestep << endl;
    os << "numEulerIters " << numEulerIters << endl;
    os << "growthAmount " << growthAmount << endl;
    os << "baseGrowthProbability " << baseGrowthProbability << endl;
    os << "maxEdgeStrain " << maxEdgeStrain << endl;
}

int Mesh::numdofs() const
{
    return 3*mesh_->n_vertices();
}

int Mesh::numedges() const
{
    return mesh_->n_edges();
}

void Mesh::dofsFromGeometry(Eigen::VectorXd &q) const
{
    if(q.size() != numdofs())
        q.resize(numdofs());

    for(int i=0; i<(int)mesh_->n_vertices(); i++)
    {
        OMMesh::Point pt = mesh_->point(mesh_->vertex_handle(i));
        for(int j=0; j<3; j++)
            q[3*i+j] = pt[j];
    }
}

void Mesh::dofsToGeometry(const VectorXd &q)
{    
    meshLock_.lock();
    {
        assert(q.size() == numdofs());

        for(int i=0; i<(int)mesh_->n_vertices(); i++)
        {
            OMMesh::Point &pt = mesh_->point(mesh_->vertex_handle(i));
            for(int j=0; j<3; j++)
                pt[j] = q[3*i+j];
        }
    }
    meshLock_.unlock();
}

void Mesh::edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2)
{
    OMMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle(eh, 0);
    pt1 = mesh_->point(mesh_->from_vertex_handle(heh1));
    pt2 = mesh_->point(mesh_->to_vertex_handle(heh1));
}

bool Mesh::exportOBJ(const char *filename)
{
    OpenMesh::IO::Options opt;
    mesh_->request_face_normals();
    mesh_->request_vertex_normals();
    mesh_->update_normals();
    opt.set(OpenMesh::IO::Options::VertexNormal);
    return OpenMesh::IO::write_mesh(*mesh_, filename, opt);
}

void Mesh::addRandomNoise(double magnitude)
{
    meshLock_.lock();
    {
        for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
            mesh_->point(vi.handle())[2] += randomRange(-magnitude,magnitude);
    }
    meshLock_.unlock();
}

bool Mesh::importOBJ(const char *filename)
{
    bool success = true;
    meshLock_.lock();
    {
        OpenMesh::IO::Options opt;
        mesh_->request_face_normals();
        mesh_->request_vertex_normals();
        mesh_->request_vertex_texcoords2D();
        opt.set(OpenMesh::IO::Options::VertexNormal);
        opt.set(OpenMesh::IO::Options::VertexTexCoord);
        success = OpenMesh::IO::read_mesh(*mesh_, filename, opt);
        mesh_->update_normals();

        colors1_.resize(mesh_->n_faces());
        colors2_.resize(mesh_->n_faces());
        colors1_.setZero();
        colors2_.setZero();
        VectorXd q;
        dofsFromGeometry(q);
        metricFromEdgeLengths(q, g1_);
        metricFromEdgeLengths(q, g2_);
    }
    frameno_ = 0;
    meshLock_.unlock();
    return success;
}

const ProblemParameters &Mesh::getParameters() const
{
    return params_;
}

void Mesh::setParameters(ProblemParameters params)
{
    params_ = params;
}

Vector3d Mesh::averageNormal(const VectorXd &q, int vidx) const
{
    OMMesh::VertexHandle vh = mesh_->vertex_handle(vidx);
    Vector3d result;
    int denom = 0;
    result.setZero();
    for(OMMesh::VertexFaceIter vfi = mesh_->vf_iter(vh); vfi; ++vfi)
    {
        result += faceNormal(q, vfi.handle().idx());
        denom++;
    }
    return result/denom;
}

Vector3d Mesh::faceNormal(const VectorXd &q, int fidx) const
{
    int vids[3];
    int idx=0;
    OMMesh::FaceHandle fh = mesh_->face_handle(fidx);
    for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fh); fvi; ++fvi)
    {
        vids[idx++] = fvi.handle().idx();
    }

    Vector3d v0 = q.segment<3>(3*vids[0]);
    Vector3d v1 = q.segment<3>(3*vids[1]);
    Vector3d v2 = q.segment<3>(3*vids[2]);
    Vector3d n = (v2-v1).cross(v0-v1);
    return n/n.norm();
}

double Mesh::infinityNorm(const VectorXd &v) const
{
    double maxval = 0;
    for(int i=0; i<v.size(); i++)
        maxval = std::max(fabs(v[i]), maxval);
    return maxval;
}

double Mesh::randomRange(double min, double max) const
{
    double val = (max-min)*double(rand())/double(RAND_MAX);
    return val+min;
}

void Mesh::dumpFrame()
{
    stringstream ss;
    ss << params_.outputDir;
    ss << "/frame_";
    ss << setfill('0') << setw(8) << frameno_ << ".obj";
    exportOBJ(ss.str().c_str());   
    frameno_++;
}

void Mesh::setConeHeights(double height)
{
    VectorXd q;
    dofsFromGeometry(q);
    int numverts = mesh_->n_vertices();
    for(int i=0; i<numverts; i++)
    {
        Vector3d pos = q.segment<3>(3*i);
        double newz = height*(1.0 - sqrt(pos[0]*pos[0] + pos[1]*pos[1]));
        q[3*i+2] = newz;
    }
    dofsToGeometry(q);
}

void Mesh::setCylinder()
{
    VectorXd q;
    dofsFromGeometry(q);
    int numverts = mesh_->n_vertices();
    for(int i=0; i<numverts; i++)
    {
        Vector3d pos = q.segment<3>(3*i);
        double newz = sqrt(1.1-pos[1]*pos[1]);
        q[3*i+2] = newz;
    }
    dofsToGeometry(q);
}

void Mesh::flatten()
{
    VectorXd q;
    dofsFromGeometry(q);
    int numverts = mesh_->n_vertices();
    for(int i=0; i<numverts; i++)
    {
        q[3*i+2] *= 1e-8;
    }
    dofsToGeometry(q);
}

void Mesh::flattenFromUV(double scale)
{
    VectorXd q;
    dofsFromGeometry(q);
    int numverts = mesh_->n_vertices();
    vector<double> oldareas;
    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        oldareas.push_back(faceArea(q, fi.handle().idx()));
    }

    for(int i=0; i<numverts; i++)
    {
        OMMesh::TexCoord2D texcoord = mesh_->texcoord2D(mesh_->vertex_handle(i));
        q[3*i] = scale*texcoord[0];
        q[3*i+1] = scale*texcoord[1];
        q[3*i+2] = 0;
    }

    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        int idx=0;
        int verts[3];
        for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fi.handle()); fvi; ++fvi)
        {
            verts[idx++] = fvi.handle().idx();
        }

        Vector3d n = (q.segment<3>(3*verts[2])-q.segment<3>(3*verts[1])).cross(q.segment<3>(3*verts[0])-q.segment<3>(3*verts[1]));
        n /= n.norm();

        assert(fabs(fabs(n[2]) - 1.0) < 1e-14);

        if(n[2] < 0)
            cout << "inverted face " << fi.handle().idx() << endl;
    }

    dofsToGeometry(q);

    setInducedMetric();
    if(colors1_.size() == mesh_->n_faces())
    {
        for(int i=0; i<(int)mesh_->n_faces(); i++)
        {
            double newarea = faceArea(q, i);
            double ratio = oldareas[i]/newarea;
            for(int k=0; k<4; k++)
            {
                g1_[4*i+k] *= ratio/(1.0+colors2_[i]);
                g2_[4*i+k] *= ratio/(1.0+colors1_[i]);
            }
            colors2_[i] = 1.0/ratio*(1.0+colors2_[i])-1.0;
            colors1_[i] = 1.0/ratio*(1.0+colors1_[i])-1.0;
        }
    }

}


void Mesh::setFlatCone(double height)
{
    deleteBadFlatConeFaces();
    VectorXd q;
    dofsFromGeometry(q);
    int numverts = mesh_->n_vertices();
    for(int i=0; i<numverts; i++)
    {
        Vector3d pos = q.segment<3>(3*i);
        double r = sqrt(1.0+height*height)*sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
        double theta = atan2(pos[1], pos[0]) / sqrt(1.0+height*height);
        q[3*i+0] = r*cos(theta);
        q[3*i+1] = r*sin(theta);
        q[3*i+2] = 0;
    }
    dofsToGeometry(q);
}

void Mesh::deleteBadFlatConeFaces()
{
    set<int> badfaces;

    for(OMMesh::FaceIter fi=mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        for(OMMesh::FaceEdgeIter fei = mesh_->fe_iter(fi.handle()); fei; ++fei)
        {
            OMMesh::HalfedgeHandle heh = mesh_->halfedge_handle(fei.handle(),0);
            OMMesh::VertexHandle topt = mesh_->to_vertex_handle(heh);
            OMMesh::VertexHandle frompt = mesh_->from_vertex_handle(heh);
            double p1y = mesh_->point(topt)[1];
            double p2y = mesh_->point(frompt)[1];
            double t = -p2y/(p1y-p2y);
            if(t < 0 || t > 1)
                continue;
            double p1x = mesh_->point(topt)[0];
            double p2x = mesh_->point(frompt)[0];
            double x = t*p1x + (1-t)*p2x;
            if(x < 0)
                badfaces.insert(fi.handle().idx());
        }
    }

    OMMesh *newmesh = new OMMesh;
    newmesh->request_vertex_texcoords2D();
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        OMMesh::Point pt = mesh_->point(vi.handle());
        OMMesh::VertexHandle newvert = newmesh->add_vertex(pt);
        newmesh->set_texcoord2D(newvert, mesh_->texcoord2D(vi.handle()));
    }

    for(int i=0; i<(int)mesh_->n_faces(); i++)
    {
        if(badfaces.count(i) > 0)
            continue;

        OMMesh::FaceHandle fh = mesh_->face_handle(i);
        vector<OMMesh::VertexHandle> newface;
        for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fh); fvi; ++fvi)
        {
            int idx = fvi.handle().idx();
            OMMesh::VertexHandle vert = newmesh->vertex_handle(idx);
            newface.push_back(vert);
        }
        newmesh->add_face(newface);
    }

    meshLock_.lock();
    {
        delete mesh_;
        mesh_ = newmesh;
    }
    meshLock_.unlock();
}

Vector3d Mesh::surfaceAreaNormal(const VectorXd &q, int vidx)
{
    OMMesh::VertexHandle vh = mesh_->vertex_handle(vidx);
    Vector3d result;
    result.setZero();

    for(OMMesh::VertexFaceIter vfi = mesh_->vf_iter(vh); vfi; ++vfi)
    {
        double area = faceArea(q, vfi.handle().idx());
        Vector3d N = this->faceNormal(q, vfi.handle().idx());
        result += area*N;
    }
    result /= result.norm();
    return result;
}

void Mesh::metricFromEdgeLengths(const VectorXd &q, VectorXd &g)
{
    g.resize(4*mesh_->n_faces());
    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        Vector4d gpart = Midedge::inducedG(*mesh_, fi.handle().idx(), q);
        for(int j=0; j<4; j++)
            g[4*fi.handle().idx()+j] = gpart[j];
    }
}

void Mesh::setInducedMetric()
{
    VectorXd q;
    dofsFromGeometry(q);
    metricFromEdgeLengths(q, g1_);
    metricFromEdgeLengths(q, g2_);
}

void Mesh::setEquilibriumMetric()
{
    VectorXd q;
    dofsFromGeometry(q);
    calculateGrowth(q, g1_, g2_);
}

void Mesh::swapYandZ()
{
    VectorXd q;
    dofsFromGeometry(q);
    for(int i=0; i<(int)mesh_->n_vertices(); i++)
        swap(q[3*i+1],q[3*i+2]);
    dofsToGeometry(q);
}

void Mesh::swapXandZ()
{
    VectorXd q;
    dofsFromGeometry(q);
    for(int i=0; i<(int)mesh_->n_vertices(); i++)
        swap(q[3*i],q[3*i+2]);
    dofsToGeometry(q);
}

void Mesh::reflectY()
{
    VectorXd q;
    dofsFromGeometry(q);
    for(int i=0; i<(int)mesh_->n_vertices(); i++)
        q[3*i+1] = -q[3*i+1];
    dofsToGeometry(q);
}

void Mesh::subdivideLinear()
{
    OMMesh *newmesh = new OMMesh();
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        OMMesh::Point pt;
        pt = mesh_->point(vi.handle());
        newmesh->add_vertex(pt);
    }
    for(OMMesh::EdgeIter ei = mesh_->edges_begin(); ei != mesh_->edges_end(); ++ei)
    {
        OMMesh::HalfedgeHandle heh = mesh_->halfedge_handle(ei.handle(),0);
        OMMesh::Point pt1 = mesh_->point(mesh_->from_vertex_handle(heh));
        OMMesh::Point pt2 = mesh_->point(mesh_->to_vertex_handle(heh));
        OMMesh::Point midpt;
        midpt = 0.5*(pt1+pt2);
        newmesh->add_vertex(midpt);
    }

    int oldverts = mesh_->n_vertices();
    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        int vertidx[3];
        int edgeidx[3];
        int idx=0;
        for(OMMesh::FaceHalfedgeIter fhi = mesh_->fh_iter(fi.handle()); fhi; ++fhi)
        {
            vertidx[idx] = mesh_->from_vertex_handle(fhi.handle()).idx();
            edgeidx[idx] = mesh_->edge_handle(fhi.handle()).idx();
            idx++;
        }
        vector<OMMesh::VertexHandle> toadd;
        toadd.push_back(newmesh->vertex_handle(vertidx[0]));
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[0]));
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[2]));
        newmesh->add_face(toadd);

        toadd.clear();
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[0]));
        toadd.push_back(newmesh->vertex_handle(vertidx[1]));
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[1]));
        newmesh->add_face(toadd);

        toadd.clear();
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[2]));
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[1]));
        toadd.push_back(newmesh->vertex_handle(vertidx[2]));
        newmesh->add_face(toadd);

        toadd.clear();
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[1]));
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[2]));
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[0]));
        newmesh->add_face(toadd);
    }

    meshLock_.lock();
    {
        int numoldfaces = mesh_->n_faces();
        int numnewfaces = newmesh->n_faces();
        delete mesh_;
        mesh_ = newmesh;
        setInducedMetric();
        if(colors1_.size() == numoldfaces)
        {
            VectorXd newcolors1(numnewfaces);
            VectorXd newcolors2(numnewfaces);
            for(int i=0; i<numoldfaces; i++)
            {
                for(int k=0; k<4; k++)
                {
                    newcolors1[4*i+k] = colors1_[i];
                    newcolors2[4*i+k] = colors2_[i];
                }
            }
            colors1_.resize(numnewfaces);
            colors2_.resize(numnewfaces);

            colors1_ = newcolors1;
            colors2_ = newcolors2;

            for(int i=0; i<numnewfaces; i++)
                for(int k=0; k<4; k++)
                {
                    g1_[4*i+k] /= (1.0+colors2_[i]);
                    g2_[4*i+k] /= (1.0+colors1_[i]);
                }
        }

    }
    meshLock_.unlock();
}

void Mesh::subdivideLoop()
{
    OMMesh *newmesh = new OMMesh();
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        OMMesh::Point pt;
        if(mesh_->is_boundary(vi.handle()))
        {
            pt = 3.0/4.0*mesh_->point(vi.handle());
            vector<OMMesh::Point> bdpts;
            for(OMMesh::VertexOHalfedgeIter voh = mesh_->voh_iter(vi.handle()); voh; ++voh)
            {
                if(mesh_->is_boundary(voh.handle()) || mesh_->is_boundary(mesh_->opposite_halfedge_handle(voh.handle())))
                    bdpts.push_back(mesh_->point(mesh_->to_vertex_handle(voh.handle())));
            }
            assert(bdpts.size() == 2);

            for(int i=0; i<(int)bdpts.size(); i++)
                pt += 1.0/8.0 * bdpts[i];
        }
        else
        {
            vector<OMMesh::Point> nbpts;
            for(OMMesh::VertexVertexIter vvi = mesh_->vv_iter(vi.handle()); vvi; ++vvi)
            {
                nbpts.push_back(mesh_->point(vvi.handle()));
            }
            double beta;
            if(nbpts.size() == 3)
                beta = 3.0/16.0;
            else
                beta = 3.0/(8.0*nbpts.size());

            pt = (1.0 - nbpts.size()*beta) * mesh_->point(vi.handle());
            for(int i=0; i<(int)nbpts.size(); i++)
                pt += beta*nbpts[i];
        }
        newmesh->add_vertex(pt);
    }
    for(OMMesh::EdgeIter ei = mesh_->edges_begin(); ei != mesh_->edges_end(); ++ei)
    {
        OMMesh::HalfedgeHandle heh = mesh_->halfedge_handle(ei.handle(),0);
        OMMesh::Point pt1 = mesh_->point(mesh_->from_vertex_handle(heh));
        OMMesh::Point pt2 = mesh_->point(mesh_->to_vertex_handle(heh));
        OMMesh::Point midpt;
        if(mesh_->is_boundary(ei.handle()))
            midpt = 0.5*(pt1+pt2);
        else
        {
            OMMesh::Point opppt1 = mesh_->point(mesh_->to_vertex_handle(mesh_->next_halfedge_handle(heh)));
            OMMesh::Point opppt2 = mesh_->point(mesh_->to_vertex_handle(mesh_->next_halfedge_handle(mesh_->opposite_halfedge_handle(heh))));
            midpt = 3.0/8.0*pt1 + 3.0/8.0*pt2 + 1.0/8.0*opppt1 + 1.0/8.0*opppt2;
        }
        newmesh->add_vertex(midpt);
    }

    int oldverts = mesh_->n_vertices();
    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        int vertidx[3];
        int edgeidx[3];
        int idx=0;
        for(OMMesh::FaceHalfedgeIter fhi = mesh_->fh_iter(fi.handle()); fhi; ++fhi)
        {
            vertidx[idx] = mesh_->from_vertex_handle(fhi.handle()).idx();
            edgeidx[idx] = mesh_->edge_handle(fhi.handle()).idx();
            idx++;
        }
        vector<OMMesh::VertexHandle> toadd;
        toadd.push_back(newmesh->vertex_handle(vertidx[0]));
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[0]));
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[2]));
        newmesh->add_face(toadd);

        toadd.clear();
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[0]));
        toadd.push_back(newmesh->vertex_handle(vertidx[1]));
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[1]));
        newmesh->add_face(toadd);

        toadd.clear();
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[2]));
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[1]));
        toadd.push_back(newmesh->vertex_handle(vertidx[2]));
        newmesh->add_face(toadd);

        toadd.clear();
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[1]));
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[2]));
        toadd.push_back(newmesh->vertex_handle(oldverts+edgeidx[0]));
        newmesh->add_face(toadd);
    }

    meshLock_.lock();
    {
        int numoldfaces = mesh_->n_faces();
        int numnewfaces = newmesh->n_faces();
        delete mesh_;
        mesh_ = newmesh;
        setInducedMetric();
        if(colors1_.size() == numoldfaces)
        {
            VectorXd newcolors1(numnewfaces);
            VectorXd newcolors2(numnewfaces);
            for(int i=0; i<numoldfaces; i++)
            {
                for(int k=0; k<4; k++)
                {
                    newcolors1[4*i+k] = colors1_[i];
                    newcolors2[4*i+k] = colors2_[i];
                }
            }
            colors1_.resize(numnewfaces);
            colors2_.resize(numnewfaces);
            for(int i=0; i<numnewfaces;i++)
            {
                double av1 = newcolors1[i];
                double av2 = newcolors2[i];
                int nbs = 0;
                OMMesh::FaceHandle fh = mesh_->face_handle(i);
                for(OMMesh::FaceFaceIter ffi = mesh_->ff_iter(fh); ffi; ++ffi)
                {
                    av1 += newcolors1[ffi.handle().idx()];
                    av2 += newcolors2[ffi.handle().idx()];
                    nbs++;
                }
                av1 /= (nbs+1.0);
                av2 /= (nbs+1.0);
                colors1_[i] = av1;
                colors2_[i] = av2;
            }

            for(int i=0; i<numnewfaces; i++)
                for(int k=0; k<4; k++)
                {
                    g1_[4*i+k] /= (1.0+colors2_[i]);
                    g2_[4*i+k] /= (1.0+colors1_[i]);
                }
        }

    }
    meshLock_.unlock();
}

void Mesh::deleteSmallUVFaces(double minarea)
{
    meshLock_.lock();
    {
        VectorXd q;
        dofsFromGeometry(q);
        int numfaces = int(mesh_->n_faces());
        vector<int> newfacemap;
        set<int> usedverts;
        for(int i=0; i<numfaces; i++)
        {
            OMMesh::FaceHandle fh = mesh_->face_handle(i);
            int idx=0;
            Vector2d uv[3];
            for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fh); fvi; ++fvi)
            {
                OMMesh::TexCoord2D tex = mesh_->texcoord2D(fvi.handle());
                uv[idx][0] = tex[0];
                uv[idx][1] = tex[1];
                idx++;
            }
            Vector2d e1 = uv[1]-uv[0];
            Vector2d e2 = uv[2]-uv[0];

            double area = fabs(e1[1]*e2[0]-e1[0]*e2[1]);
            if(area > minarea)
            {
                newfacemap.push_back(i);
                for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fh); fvi; ++fvi)
                    usedverts.insert(fvi.handle().idx());
            }
        }

        int numnewfaces = newfacemap.size();
        vector<int> newvertmap;
        vector<int> invvertmap;
        for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
        {
            if(usedverts.count(vi.handle().idx()) > 0)
            {
                invvertmap.push_back(newvertmap.size());
                newvertmap.push_back(vi.handle().idx());
            }
            else
            {
                invvertmap.push_back(-1);
            }
        }

        if(colors1_.size() == numfaces)
        {
            VectorXd newcolors1(numnewfaces);
            VectorXd newcolors2(numnewfaces);
            for(int i=0; i<numnewfaces; i++)
            {
                newcolors1[i] = colors1_[newfacemap[i]];
                newcolors2[i] = colors2_[newfacemap[i]];
            }
            colors1_ = newcolors1;
            colors2_ = newcolors2;
        }

        if(g1_.size() == 4*numfaces)
        {
            VectorXd newg1(4*numnewfaces);
            VectorXd newg2(4*numnewfaces);
            for(int i=0; i<numnewfaces; i++)
            {
                for(int j=0; j<4; j++)
                {
                    newg1[4*i+j] = g1_[4*newfacemap[i]+j];
                    newg2[4*i+j] = g2_[4*newfacemap[i]+j];
                }
            }
            g1_ = newg1;
            g2_ = newg2;
        }

        OMMesh *newmesh = new OMMesh();
        newmesh->request_vertex_texcoords2D();

        int numnewverts = newvertmap.size();

        for(int i=0; i<numnewverts; i++)
        {
            OMMesh::VertexHandle vh = mesh_->vertex_handle(newvertmap[i]);
            OMMesh::Point pt = mesh_->point(vh);
            OMMesh::VertexHandle newvert = newmesh->add_vertex(pt);
            newmesh->set_texcoord2D(newvert, mesh_->texcoord2D(vh));
        }
        for(int i=0; i<numnewfaces; i++)
        {
            vector<OMMesh::VertexHandle> newface;
            for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(mesh_->face_handle(newfacemap[i])); fvi; ++fvi)
            {
                newface.push_back(newmesh->vertex_handle(invvertmap[fvi.handle().idx()]));
            }
            newmesh->add_face(newface);
        }
        delete mesh_;
        mesh_ = newmesh;
    }
    meshLock_.unlock();
}
