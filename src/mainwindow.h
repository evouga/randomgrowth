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

    void on_actionAdd_Noise_triggered();

    void on_actionMake_Cone_triggered();

    void on_actionMake_Flat_Cone_triggered();

    void on_actionMake_Cylinder_triggered();

    void on_actionSet_Induced_Metric_triggered();

    void on_actionSet_Equilibrium_Metric_triggered();

    void on_actionFlatten_triggered();

    void on_actionSwap_Y_and_Z_triggered();

    void on_colorCutoffEdit_textEdited(const QString &arg1);

    void on_actionSwap_X_and_Z_triggered();

    void on_colorMenuComboBox_currentIndexChanged(int index);

    void on_actionSet_from_UV_triggered();

    void on_actionDelete_Small_Faces_triggered();

    void on_actionRelax_Configuration_triggered();

    void on_actionReflect_Y_triggered();

    void on_actionLinear_triggered();

    void on_actionLoop_triggered();

private:
    void updateGL();
    void saveScreenshot();
    void saveScreenshot(const std::string &filename);
    ProblemParameters getParameters();
    std::string launchImportOBJDialog();

    Ui::MainWindow *ui;
    Controller *cont_;
    QTimer *repainttimer_;
    const static QString colorMenuOptionsNames_[ProblemParameters::CRM_SIZE];
    const static ProblemParameters::colorRenderMode colorMenuOptionsVals_[ProblemParameters::CRM_SIZE];
};

#endif // MAINWINDOW_H
