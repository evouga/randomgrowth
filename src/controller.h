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
    void updateParameters(ProblemParameters params);
    void relaxEmbedding();
    void quit();
    void centerCamera();
    void updateGL();
    void addNoise();
    void makeCone();
    void makeCylinder();
    void makeFlatCone();
    void flattenMesh();
    void setInducedMetric();
    void setEquilibriumMetric();
    void swapYandZ();
    void swapXandZ();
    void reflectY();
    void subdivideLinear();
    void subdivideLoop();
    void flattenFromUV();
    void deleteSmallUVFaces();
    void relaxConfiguration();

private:    
    MainWindow &mw_;
    Mesh m_;
};

#endif // CONTROLLER_H
