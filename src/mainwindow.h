#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <string>
#include "simulationmesh.h"

class Controller;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void setController(Controller &cont);

public slots:
    void setParameters(ProblemParameters params);
    void showError(std::string error);
    void centerCamera(Eigen::Vector3d centroid, double radius);
    void repaintMesh();
    void tick();

private slots:
    void on_actionExit_triggered();

    void on_actionReset_Camera_triggered();

    void on_actionTake_Screenshot_triggered();

    void on_wireframeCheckBox_clicked();

    void on_smoothShadeCheckBox_clicked();

    void on_actionExport_OBJ_triggered();

    void on_actionImport_OBJ_triggered();

    void on_poissonRatioEdit_textEdited(const QString &arg1);

    void on_youngsModulusEdit_textEdited(const QString &arg1);

    void on_thicknessEdit_textEdited(const QString &arg1);

    void on_densityEdit_textEdited(const QString &arg1);

    void on_dampingCoeffEdit_textEdited(const QString &arg1);

    void on_eulerTimestepEdit_textEdited(const QString &arg1);

    void on_scaleEdit_textEdited(const QString &arg1);

    void on_outputEdit_textEdited(const QString &arg1);

    void on_actionImport_Metric_triggered();

    void on_actionMake_Cone_triggered();

    void on_crushButton_clicked();

    void on_constantPressureRadio_toggled(bool checked);

    void on_scalingPressureRadio_toggled(bool checked);

    void on_constantPressureEdit_textEdited(const QString &arg1);

    void on_airLeakEdit_textEdited(const QString &arg1);

    void on_coneAngleEdit_textEdited(const QString &arg1);

    void on_coneHeightEdit_textEdited(const QString &arg1);

    void on_crushTimeEdit_textEdited(const QString &arg1);

    void on_constantVelocityRadio_toggled(bool checked);

    void on_simulatedCrushingRadio_toggled(bool checked);

    void on_crushMassEdit_textEdited(const QString &arg1);

    void on_simTimeEdit_textChanged(const QString &arg1);

    void on_initialVelEdit_textEdited(const QString &arg1);

    void on_restoreEdit_textEdited(const QString &arg1);

    void on_holeRadiusRadio_toggled(bool checked);

    void on_holeRadiusEdit_textEdited(const QString &arg1);

private:
    void updateGL();
    void saveScreenshot();
    void saveScreenshot(const std::string &filename);
    ProblemParameters getParameters();
    std::string launchImportOBJDialog();
    std::string launchImportMetricDialog();

    Ui::MainWindow *ui;
    Controller *cont_;
    QTimer *repainttimer_;
};

#endif // MAINWINDOW_H
