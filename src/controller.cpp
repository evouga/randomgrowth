#include "controller.h"
#include "mainwindow.h"
#include "mesh.h"
#include <string>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <QMessageBox>
#include <fstream>

using namespace std;
using namespace Eigen;

Controller::Controller(MainWindow &mw) : mw_(mw), m_()
{
    ProblemParameters params = m_.getParameters();
    mw_.setParameters(params);
}

void Controller::quit()
{
    QMetaObject::invokeMethod(&mw_, "close");
}

void Controller::renderMesh()
{
    m_.render();
}

void Controller::exportOBJ(string filename)
{
    if(!m_.exportOBJ(filename.c_str()))
    {
        QString msg = "Couldn't write file " + QString::fromStdString(filename) + ". Save failed.";
        QMessageBox::warning(&mw_, "Couldn't Write File", msg, QMessageBox::Ok);
        return;
    }
}

void Controller::importOBJ(string filename)
{
    if(!m_.importOBJ(filename.c_str()))
    {
        string err = string("Couldn't load OBJ: ") + filename;
        QMetaObject::invokeMethod(&mw_, "showError", Q_ARG(std::string, err));
    }
    else
    {
        centerCamera();
    }
}

void Controller::updateParameters(ProblemParameters params)
{    
    m_.setParameters(params);
    updateGL();
}

void Controller::crush()
{
    m_.crush(*this, 3.885, 0.5);
}

void Controller::centerCamera()
{
    Vector3d centroid = m_.centroid();
    double radius = m_.radius();
    QMetaObject::invokeMethod(&mw_, "centerCamera", Q_ARG(Eigen::Vector3d, centroid), Q_ARG(double, radius));
}

void Controller::updateGL()
{
    QMetaObject::invokeMethod(&mw_, "repaintMesh");
}

void Controller::addNoise()
{
    m_.addRandomNoise(1e-3);
    updateGL();
}

void Controller::symmetrize()
{
    m_.symmetrize(3);
}

void Controller::makeCone()
{
    m_.setConeHeights(3.885);
}

void Controller::makeFlatCone()
{
    m_.setFlatCone(3.885);
}
