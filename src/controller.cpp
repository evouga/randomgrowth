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
    double Youngs = m_.getYoungsModulus();
    double Poisson = m_.getPoissonRatio();
    double h = m_.getThickness();
    mw_.setParameters(Youngs, Poisson, h);
}

void Controller::quit()
{
    mw_.close();
}

void Controller::renderMesh()
{
    bool showWireframe = mw_.showWireframe();
    bool smoothShade = mw_.smoothShade();
    m_.render(showWireframe, smoothShade);
}

void Controller::getSceneBounds(Eigen::Vector3d &center, double &radius)
{
    center = m_.centroid();
    radius = m_.radius();
}

void Controller::exportOBJ(const char *filename)
{
    if(!m_.exportOBJ(filename))
    {
        QString msg = "Couldn't write file " + QString(filename) + ". Save failed.";
        QMessageBox::warning(&mw_, "Couldn't Write File", msg, QMessageBox::Ok);
        return;
    }
}

void Controller::importOBJ(const char *filename)
{
    if(!m_.importOBJ(filename))
    {
        string err = string("Couldn't load OBJ: ") + filename;
        mw_.showError(err);
    }
    else
    {
        mw_.centerCamera();
    }
}

void Controller::updateParameters()
{
    double Youngs, Poisson, h;
    mw_.getParameters(Youngs, Poisson, h);
    m_.setYoungsModulus(Youngs);
    m_.setPoissonRatio(Poisson);
    m_.setThickness(h);
}
