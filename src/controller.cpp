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
    ProblemParameters params = mw_.getParameters();
    m_.setParameters(params);
}

void Controller::findMetric()
{
    m_.relaxIntrinsicLengths();
}

void Controller::relaxEmbedding()
{
    m_.relaxEmbedding();
}
