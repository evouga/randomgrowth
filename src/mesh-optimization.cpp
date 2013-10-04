#include "mesh.h"
#include <fadbad.h>
#include <fadiff.h>
#include <badiff.h>
#include <iomanip>
#include "autodifftemplates.h"
#include "controller.h"

using namespace std;
using namespace Eigen;
using namespace fadbad;
using namespace OpenMesh;

typedef Eigen::Triplet<double> Tr;

void Mesh::elasticEnergy(const VectorXd &q, const VectorXd &g, double &energyB, double &energyS) const
{
    EnergyDerivatives derivs = NONE;
    VectorXd gradq, gradg;
    SparseMatrix<double> hessq, hessg;
    elasticEnergy(q, g, energyB, energyS, gradq, gradg, hessq, hessg, derivs);
}

void Mesh::elasticEnergyG(const VectorXd &q,
                          const VectorXd &g,
                          double &energyB,
                          double &energyS,
                          VectorXd &gradg,
                          Eigen::SparseMatrix<double> &hessg) const
{
    EnergyDerivatives derivs = G;
    VectorXd gradq;
    SparseMatrix<double> hessq;
    elasticEnergy(q, g, energyB, energyS, gradq, gradg, hessq, hessg, derivs);
}

void Mesh::elasticEnergyQ(const VectorXd &q,
                          const VectorXd &g,
                          double &energyB,
                          double &energyS,
                          VectorXd &gradq,
                          Eigen::SparseMatrix<double> &hessq) const
{
    EnergyDerivatives derivs = Q;
    VectorXd gradg;
    SparseMatrix<double> hessg;
    elasticEnergy(q, g, energyB, energyS, gradq, gradg, hessq, hessg, derivs);
}

void Mesh::elasticEnergy(const VectorXd &q,
                         const VectorXd &g,
                         double &energyB,
                         double &energyS,
                         VectorXd &gradq,
                         VectorXd &gradg,
                         Eigen::SparseMatrix<double> &hessq,
                         Eigen::SparseMatrix<double> &hessg,
                         EnergyDerivatives derivs) const
{
    assert(q.size() == numdofs());
    assert(g.size() == numedges());
    energyB = energyS = 0;

    if(derivs & Q)
    {
        gradq.resize(numdofs());
        gradq.setZero();
        hessq.resize(numdofs(), numdofs());
    }

    if(derivs & G)
    {
        gradg.resize(numedges());
        gradg.setZero();
        hessg.resize(numedges(), numedges());
    }


    vector<Tr> Hqcoeffs;
    vector<Tr> Hgcoeffs;

    // bending energy
    double bendcoeff = params_.h*params_.h*params_.YoungsModulus/24.0/(1.0+params_.PoissonRatio);
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        if(mesh_->is_boundary(vi.handle()))
            continue;

        B<F<double> > stencilenergy = 0;
        // bending energy Laplacian term
        double Lcoeff = 1.0/(1.0-params_.PoissonRatio);
        vector<double> spokelens;
        vector<int> spokeidx;
        vector<double> rightopplens;
        vector<int> rightoppidx;

        vector<Vector3d> nbq;
        vector<int> nbidx;
        for(OMMesh::VertexOHalfedgeIter voh = mesh_->voh_iter(vi.handle()); voh; ++voh)
        {
            OMMesh::HalfedgeHandle heh = voh.handle();
            int eidx = mesh_->edge_handle(heh).idx();
            spokelens.push_back(g[eidx]);
            spokeidx.push_back(eidx);

            OMMesh::VertexOHalfedgeIter nextoh = voh;
            ++nextoh;
            if(!nextoh)
                nextoh = mesh_->voh_iter(vi.handle());

            OMMesh::VertexHandle nextvert = mesh_->to_vertex_handle(nextoh.handle());

            OMMesh::HalfedgeHandle opp = mesh_->next_halfedge_handle(heh);;
            if(mesh_->to_vertex_handle(opp) != nextvert)
            {
                opp = mesh_->prev_halfedge_handle(mesh_->opposite_halfedge_handle(heh));
                assert(mesh_->from_vertex_handle(opp) == nextvert);
            }

            int oidx = mesh_->edge_handle(opp).idx();
            rightopplens.push_back(g[oidx]);
            rightoppidx.push_back(oidx);

            OMMesh::VertexHandle vh = mesh_->to_vertex_handle(heh);
            nbq.push_back(q.segment<3>(3*vh.idx()));
            nbidx.push_back(vh.idx());
        }

        int numnbs = (int)spokelens.size();
        int diffqvars;
        if(derivs & Q)
            diffqvars = 3*(numnbs+1);
        else
            diffqvars = 0;

        int diffvars;
        if(derivs & G)
            diffvars = diffqvars + 2*numnbs;
        else
            diffvars = diffqvars;

        F<double> *dspokelens = new F<double>[numnbs];
        B<F<double> > *ddspokelens = new B<F<double> >[numnbs];
        for(int i=0; i<numnbs; i++)
        {
            dspokelens[i] = spokelens[i];
            if(derivs & G)
                dspokelens[i].diff(diffqvars+i, diffvars);
            ddspokelens[i] = dspokelens[i];
            //if(derivs & G)
            //    ddspokelens[i].diff(diffqvars+i, diffvars);
        }
        F<double> *dopplens = new F<double>[numnbs];
        B<F<double> > *ddopplens = new B<F<double> >[numnbs];
        for(int i=0; i<numnbs; i++)
        {
            dopplens[i] = rightopplens[i];
            if(derivs & G)
                dopplens[i].diff(diffqvars+numnbs+i, diffvars);
            ddopplens[i] = dopplens[i];
            //if(derivs & G)
            //    ddopplens[i].diff(diffqvars+numnbs+i, diffvars);
        }

        F<double> dq[3];
        B<F<double> > ddq[3];
        int vidx = vi.handle().idx();
        Vector3d thisq = q.segment<3>(3*vi.handle().idx());
        for(int i=0; i<3; i++)
        {
            dq[i] = thisq[i];
            if(derivs & Q)
                dq[i].diff(i,diffvars);
            ddq[i] = dq[i];
            //if(derivs & Q)
            //    ddq[i].diff(i, diffvars);
        }

        F<double> *dnbq[3];
        B<F<double> > *ddnbq[3];
        for(int i=0; i<3; i++)
        {
            dnbq[i] = new F<double>[numnbs];
            ddnbq[i] = new B<F<double> >[numnbs];
            for(int j=0; j<numnbs; j++)
            {
                dnbq[i][j] = nbq[j][i];
                if(derivs & Q)
                    dnbq[i][j].diff(3+3*j+i,diffvars);
                ddnbq[i][j] = dnbq[i][j];
                //if(derivs & Q)
                //    ddnbq[i][j].diff(3+3*j+i,diffvars);
            }
        }

        {
            B<F<double> > Lr[3];
            for(int j=0; j<3; j++)
                Lr[j] = L(ddq[j], numnbs, ddnbq[j], ddspokelens, ddopplens);
            B<F<double> > materialArea = dualbarycentricarea(numnbs, ddspokelens, ddopplens);

            stencilenergy += bendcoeff*Lcoeff*normSquared(Lr)/materialArea;

            // bending energy det term
            vector<B<F<double> > *> nbddq;
            for(int i=0; i<numnbs; i++)
            {
                B<F<double> > *ddq = new B<F<double> >[3];
                nbddq.push_back(ddq);
                for(int j=0; j<3; j++)
                {
                    ddq[j] = ddnbq[j][i];
                }
            }

            B<F<double> > Kcurv = K(ddq, nbddq);
            B<F<double> > embeddedArea = dualbarycentricarea(ddq, nbddq);
            stencilenergy += -2.0*bendcoeff*Kcurv*embeddedArea*embeddedArea/materialArea;
            for(int i=0; i<numnbs; i++)
                delete[] nbddq[i];
        }

        stencilenergy.diff(0,1);

        energyB += stencilenergy.val().val();

        if(derivs & Q)
        {
            for(int j=0; j<3; j++)
            {
                gradq[3*vidx+j] += ddq[j].d(0).val();//->stencilenergy.d(j).val();
            }
        }

        for(int i=0; i<numnbs; i++)
        {
            if(derivs & Q)
            {
                for(int j=0; j<3; j++)
                {
                    gradq[3*nbidx[i]+j] += ddnbq[j][i].d(0).val();//stencilenergy.d(3+3*i+j).val();
                }
            }
            if(derivs & G)
            {
                gradg[spokeidx[i]] += ddspokelens[i].d(0).val();// stencilenergy.d(diffqvars+i).val();
                gradg[rightoppidx[i]] += ddopplens[i].d(0).val();// stencilenergy.d(diffqvars+numnbs+i).val();
            }
        }

        if(derivs & Q)
        {
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    double hess = ddq[j].d(0).d(k);//stencilenergy.d(j).d(k);
                    if(hess != 0)
                        Hqcoeffs.push_back(Tr(3*vidx+j,3*vidx+k,hess));
                }
                for(int k=0; k<numnbs; k++)
                {
                    for(int l=0; l<3; l++)
                    {
                        double hess = ddq[j].d(0).d(3+3*k+l);//stencilenergy.d(j).d(3+3*k+l);
                        if(hess != 0)
                        {
                            Hqcoeffs.push_back(Tr(3*vidx+j,3*nbidx[k]+l,hess));
                            Hqcoeffs.push_back(Tr(3*nbidx[k]+l,3*vidx+j,hess));
                        }
                    }
                }
            }

            for(int i=0; i<numnbs; i++)
            {
                for(int j=0; j<numnbs; j++)
                {
                    for(int k=0; k<3; k++)
                    {
                        for(int l=0; l<3; l++)
                        {
                            double hess = ddnbq[k][i].d(0).d(3+3*j+l);//stencilenergy.d(3+3*i+k).d(3+3*j+l);
                            if(hess != 0)
                                Hqcoeffs.push_back(Tr(3*nbidx[i]+k,3*nbidx[j]+l,hess));
                        }
                    }
                }
            }
        }

        if(derivs & G)
        {
            for(int i=0; i<numnbs; i++)
            {
                for(int j=0; j<numnbs; j++)
                {
                    double hess = ddspokelens[i].d(0).d(diffqvars+j);// stencilenergy.d(diffqvars+i).d(diffqvars+j);
                    if(hess != 0)
                        Hgcoeffs.push_back(Tr(spokeidx[i],spokeidx[j],hess));
                }
                for(int j=0; j<numnbs; j++)
                {
                    double hess = ddspokelens[i].d(0).d(diffqvars+numnbs+j);//stencilenergy.d(diffqvars+i).d(diffqvars+numnbs+j);
                    if(hess != 0)
                    {
                        Hgcoeffs.push_back(Tr(spokeidx[i],rightoppidx[j],hess));
                        Hgcoeffs.push_back(Tr(rightoppidx[j],spokeidx[i],hess));
                    }
                }
                for(int j=0; j<numnbs; j++)
                {
                    double hess = ddopplens[i].d(0).d(diffqvars+numnbs+j);// stencilenergy.d(diffqvars+numnbs+i).d(diffqvars+numnbs+j);
                    if(hess != 0)
                        Hgcoeffs.push_back(Tr(rightoppidx[i],rightoppidx[j],hess));
                }
            }
        }

        delete[] dspokelens;
        delete[] ddspokelens;
        delete[] dopplens;
        delete[] ddopplens;
        for(int i=0; i<3; i++)
        {
            delete[] dnbq[i];
            delete[] ddnbq[i];
        }

    }

    // Stretching energy
    double stretchcoeff = params_.YoungsModulus/8.0/(1.0+params_.PoissonRatio);

    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        B<F<double> > stencilenergy = 0;
        int vidx[3];
        Vector3d pts[3];
        double elens[3];
        int eidx[3];
        int idx=0;
        for(OMMesh::FaceHalfedgeIter fei = mesh_->fh_iter(fi.handle()); fei; ++fei)
        {
            vidx[idx] = mesh_->to_vertex_handle(fei.handle()).idx();
            pts[idx] = q.segment<3>(3*mesh_->to_vertex_handle(fei.handle()).idx());
            int eid = mesh_->edge_handle(fei.handle()).idx();
            eidx[idx] = eid;
            elens[idx] = g[eid];
            idx++;
        }
        assert(idx==3);

        int diffqvars;
        if(derivs & Q)
            diffqvars = 9;
        else
            diffqvars = 0;

        int diffvars;
        if(derivs & G)
            diffvars = diffqvars + 3;
        else
            diffvars = diffqvars;

        F<double> dq[3][3];
        B<F<double> > ddq[3][3];
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                dq[i][j] = pts[i][j];
                if(derivs & Q)
                    dq[i][j].diff(3*i+j, diffvars);
                ddq[i][j] = dq[i][j];
                //if(derivs & Q)
                //    ddq[i][j].diff(3*i+j, diffvars);
            }
        }

        F<double> delen[3];
        B<F<double> > ddelen[3];
        for(int i=0; i<3; i++)
        {
            delen[i] = elens[i];
            if(derivs & G)
                delen[i].diff(diffqvars+i, diffvars);
            ddelen[i] = delen[i];
            //if(derivs & G)
            //    ddelen[i].diff(diffqvars+i, diffvars);
        }

        {
            B<F<double> > emblens[3];
            for(int i=0; i<3; i++)
            {
                int previdx = (3+i-1)%3;
                B<F<double> > vec[3];
                diff(ddq[i],ddq[previdx],vec);
                emblens[i] = norm(vec);
            }

            B<F<double> > g[4];
            B<F<double> > a[4];
            metrictri(ddelen[0], ddelen[1], ddelen[2], g);
            metrictri(emblens[0], emblens[1], emblens[2], a);
            B<F<double> > ginva[4];
            invGtimesH(g, a, ginva);
            ginva[0] -= 1.0;
            ginva[3] -= 1.0;
            B<F<double> > matarea = 0.5*sqrt(det(g));
            stencilenergy += matarea*stretchcoeff/(1.0-params_.PoissonRatio)*tr(ginva)*tr(ginva);
            stencilenergy += matarea*stretchcoeff*-2.0*det(ginva);
        }

        stencilenergy.diff(0,1);

        energyS += stencilenergy.val().val();

        for(int i=0; i<3; i++)
        {
            if(derivs & Q)
            {
                for(int j=0; j<3; j++)
                {
                    gradq[3*vidx[i]+j] += ddq[i][j].d(0).val();//stencilenergy.d(3*i+j).val();
                }
            }
            if(derivs & G)
                gradg[eidx[i]] += ddelen[i].d(0).val();//stencilenergy.d(diffqvars + i).val();
        }

        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                if(derivs & Q)
                {
                    for(int k=0; k<3; k++)
                    {
                        for(int l=0; l<3; l++)
                        {
                            double hess = ddq[i][k].d(0).d(3*j+l);//stencilenergy.d(3*i+k).d(3*j+l);
                            if(hess != 0.0)
                            {
                                Hqcoeffs.push_back(Tr(3*vidx[i]+k,3*vidx[j]+l,hess));
                            }
                        }
                    }
                }

                if(derivs & G)
                {
                    double hess = ddelen[i].d(0).d(diffqvars+j);//stencilenergy.d(diffqvars+i).d(diffqvars+j);
                    if(hess != 0)
                    {
                        Hgcoeffs.push_back(Tr(eidx[i],eidx[j],hess));
                    }
                }
            }
        }
    }

    if(derivs & Q)
        hessq.setFromTriplets(Hqcoeffs.begin(), Hqcoeffs.end());
    if(derivs & G)
        hessg.setFromTriplets(Hgcoeffs.begin(), Hgcoeffs.end());
}

double Mesh::triangleInequalityLineSearch(const VectorXd &g, const VectorXd &dg) const
{
    assert(g.size() == numedges());
    assert(dg.size() == numedges());

    double maxt = std::numeric_limits<double>::infinity();

    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        vector<double> gs;
        vector<double> dgs;
        for(OMMesh::FaceEdgeIter fei = mesh_->fe_iter(fi.handle()); fei; ++fei)
        {
            gs.push_back(g[fei.handle().idx()]);
            dgs.push_back(dg[fei.handle().idx()]);
        }
        assert(gs.size() == 3);
        for(int i=0; i<3; i++)
        {
            int idx[3];
            for(int j=0; j<3; j++)
            {
                idx[j] = (i+j)%3;
            }
            double thismax = triangleInequalityLineSearch(gs[idx[0]], gs[idx[1]], gs[idx[2]],
                    dgs[idx[0]], dgs[idx[1]], dgs[idx[2]]);

            maxt = std::min(maxt, thismax);
        }
    }

    return maxt;
}

double Mesh::triangleInequalityLineSearch(double g0, double g1, double g2, double dg0, double dg1, double dg2) const
{
    double num = g0+g1-g2;
    double denom = dg2-dg0-dg1;
    if(denom == 0)
        return std::numeric_limits<double>::infinity();
    double cand = num/denom;
    if(cand < 0)
        return std::numeric_limits<double>::infinity();
    return cand;
}

bool Mesh::relaxEnergy(Controller &cont, RelaxationType type)
{
    VectorXd q(numdofs());
    VectorXd g(numedges());
    dofsFromGeometry(q, g);

    VectorXd gradient;

    SparseMatrix<double> hessian;

    double energyB, energyS;

    if(type == RelaxEmbedding)
        elasticEnergyQ(q, g, energyB, energyS, gradient, hessian);
    else if(type == RelaxMetric)
        elasticEnergyG(q, g, energyB, energyS, gradient, hessian);

    for(int i=0; i<params_.maxiters; i++)
    {
        if(gradient.norm() < params_.tol)
            break;
        SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
        solver.compute(hessian);
        VectorXd searchdir = -solver.solve(gradient);

        std::cout << std::fixed << std::setprecision(8) << "Iter " << i+1 << "   Eb " << energyB << "   Es " << energyS << "   |dq| " << gradient.norm() << "   sd = " << searchdir.dot(gradient);


        double stepsize = 1.0;

        if(type == RelaxMetric)
        {
            stepsize = triangleInequalityLineSearch(g, searchdir);
            stepsize = std::min(1.0, 0.9*stepsize);
        }

        double initialenergy = energyB+energyS;

        VectorXd newdofs;
        int lsiters = 0;
        if(searchdir.dot(gradient) >= 0)
            searchdir = -gradient;

        bool abort = searchdir.dot(gradient) > 0;

        if(!abort)
        {
            do
            {
                if(type == RelaxEmbedding)
                {
                    newdofs = q + stepsize*searchdir;
                    elasticEnergyQ(newdofs, g, energyB, energyS, gradient, hessian);
                }
                else if(type == RelaxMetric)
                {
                    newdofs = g + stepsize*searchdir;
                    elasticEnergyG(q, newdofs, energyB, energyS, gradient, hessian);
                }

                stepsize /= 2.0;
                if(++lsiters > params_.maxlinesearchiters)
                {
                    abort = true;
                    break;
                }
            }
            while(energyB+energyS > initialenergy);
        }

        if(abort)
        {
            std::cout << endl;
            break;
        }

        if(type == RelaxEmbedding)
            q = newdofs;
        else if(type == RelaxMetric)
            g = newdofs;

        std::cout << std::fixed << std::setprecision(8) << "   h " << stepsize*2.0 << std::endl;

        dofsToGeometry(q, g);
        cont.updateGL();
    }

    if(gradient.norm() < params_.tol)
    {
        std::cout << "Converged, final energies " << energyB << ", " << energyS << std::endl;

        double eval;
        std::cout << smallestEigenvalue(hessian, eval) << " ";
        cout << eval << endl;
        return true;
    }
    std::cout << "Failed to converge" << std::endl;
    return false;
}

bool Mesh::largestMagnitudeEigenvalue(const Eigen::SparseMatrix<double> &M, double &eigenvalue)
{
    int dim = M.cols();
    VectorXd v(dim);
    v.setRandom();
    v.normalize();
    double oldestimate = std::numeric_limits<double>::infinity();
    double curestimate=0;

    for(int iter=0; iter < params_.maxpoweriters; iter++)
    {
        v = M*v;
        v.normalize();
        curestimate = v.dot(M*v);
        if(fabs(curestimate-oldestimate) < params_.powertol)
            break;
        oldestimate = curestimate;
    }
    eigenvalue = curestimate;
    return fabs(curestimate-oldestimate) < params_.powertol;
}

bool Mesh::smallestEigenvalue(const Eigen::SparseMatrix<double> &M, double &eigenvalue)
{
    double largesteval;
    if(!largestMagnitudeEigenvalue(M, largesteval))
        return false;

    if(largesteval < 0)
    {
        eigenvalue = largesteval;
        return true;
    }

    int dim = M.cols();
    SparseMatrix<double> shift(dim,dim);
    shift.setIdentity();
    shift *= -largesteval;
    SparseMatrix<double> newM = M + shift;
    double newlargest;
    if(!largestMagnitudeEigenvalue(newM, newlargest))
        return false;
    eigenvalue = newlargest+largesteval;
    return true;
}
