#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "mesh.h"
#include <Eigen/Core>
#include <QObject>

class MainWindow;

class Controller : public QObject
{
    Q_OBJECT

public:
    Controller(MainWindow &mw);

    void renderMesh();

public slots:
    void exportOBJ(std::string filename);
    void importOBJ(std::string filename);
    void importMetric(std::string filename);
    void updateParameters(ProblemParameters params);
    void relaxEmbedding();
    void extremizeWithNewton();
    void quit();
    void centerCamera();
    void updateGL();
    void addNoise();
    void setNoTargetMetric();
    void setNegativeCurvatureTargetMetric();
    void symmetrize();
    void printHessianEigenvalues();
    void makeCone();
    void makeFlatCone();
    void crush();

private:    
    MainWindow &mw_;
    Mesh m_;
};

#endif // CONTROLLER_H
