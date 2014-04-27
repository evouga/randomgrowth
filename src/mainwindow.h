#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <string>
#include "mesh.h"

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

    void on_relaxEmbeddingButton_clicked();

    void on_densityEdit_textEdited(const QString &arg1);

    void on_dampingCoeffEdit_textEdited(const QString &arg1);

    void on_eulerTimestepEdit_textEdited(const QString &arg1);

    void on_eulerItersEdit_textEdited(const QString &arg1);

    void on_growthAmountEdit_textEdited(const QString &arg1);

    void on_scaleEdit_textEdited(const QString &arg1);

    void on_outputEdit_textEdited(const QString &arg1);

    void on_maxStrainEdit_textEdited(const QString &arg1);

    void on_baseProbabilityEdit_textEdited(const QString &arg1);

    void on_actionImport_Metric_triggered();

    void on_actionAdd_Noise_triggered();

    void on_actionSet_No_Target_Metric_triggered();

    void on_actionSet_Negative_K_Target_Metric_triggered();

    void on_actionMake_Cone_triggered();

    void on_actionMake_Flat_Cone_triggered();

    void on_actionSet_Current_Lengths_as_Intrinsic_triggered();

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
