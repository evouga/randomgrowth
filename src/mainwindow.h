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

    void saveScreenshot();
    void saveScreenshot(const std::string &filename);
    void showError(const std::string &error);
    void centerCamera();
    bool showWireframe() const;
    bool smoothShade() const;
    void repaintMesh();
    std::string launchImportOBJDialog();
    void setParameters(const ProblemParameters &params);
    ProblemParameters getParameters();

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

    void on_findMetricButton_clicked();

private:
    void updateGL();
    Ui::MainWindow *ui;
    Controller *cont_;
};

#endif // MAINWINDOW_H
