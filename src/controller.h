#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "simulationmesh.h"
#include <Eigen/Core>
#include <QObject>
#include <string>

class MainWindow;

class Controller : public QObject
{
    Q_OBJECT

public:
    Controller(MainWindow &mw);

    void renderMesh();

    std::string meshName = "single_slit";

public slots:
    void exportOBJ(std::string filename);
    void importOBJ(std::string filename);
    void updateParameters(ProblemParameters params);
    void quit();
    void centerCamera();
    void updateGL();
    void makeCone();
    void makeFlatCone();
    void pull();

private:    
    MainWindow &mw_;
    SimulationMesh m_;
};

#endif // CONTROLLER_H
