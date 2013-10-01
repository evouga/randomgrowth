#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "mesh.h"
#include <Eigen/Core>

class MainWindow;

class Controller
{
public:
    Controller(MainWindow &mw);

    void renderMesh();
    void quit();
    void getSceneBounds(Eigen::Vector3d &center, double &radius);
    void exportOBJ(const char *filename);
    void importOBJ(const char *filename);
    void updateParameters();

private:
    MainWindow &mw_;
    Mesh m_;
};

#endif // CONTROLLER_H
