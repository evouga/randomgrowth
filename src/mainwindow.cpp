#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "controller.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QDir>
#include <QDateTime>
#include <QPixmap>

using namespace Eigen;
using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    cont_(NULL)
{
    ui->setupUi(this);
    repainttimer_ = new QTimer(this);
    repainttimer_->start(300);
    connect(repainttimer_, SIGNAL(timeout()), this, SLOT(tick()));
}

MainWindow::~MainWindow()
{
    delete ui;
    delete repainttimer_;
}

void MainWindow::setController(Controller &cont)
{
    cont_ = &cont;
    ui->GLwidget->setController(cont);
}

void MainWindow::tick()
{
    repaintMesh();
}

void MainWindow::showError(string error)
{
    QMessageBox::warning(this, tr("VViewer"),
                                    QString(error.c_str()),
                                    QMessageBox::Ok, QMessageBox::NoButton);
}

void MainWindow::centerCamera(Vector3d centroid, double radius)
{
    ui->GLwidget->centerCamera(centroid, radius);
    updateGL();
}


void MainWindow::saveScreenshot()
{
    QDateTime dateTime = QDateTime::currentDateTime();
    QString dateTimeString = dateTime.toString("dd_MM_yy_hh_mm_ss_zzz");
    string filename = "screen_" + string(dateTimeString.toStdString()) + ".png";
    saveScreenshot(filename);
}

void MainWindow::saveScreenshot(const string &filename)
{
    updateGL();
    QImage img = ui->GLwidget->grabFrameBuffer(true);

    //QPixmap p = this->grab();
    QString curdir = QDir::currentPath();
    string fullname = string(curdir.toStdString()) + "/output/" + filename;
    //p.save(QString::fromUtf8(fullname.c_str()));
    img.save(QString::fromUtf8(fullname.c_str()));
}

string MainWindow::launchImportOBJDialog()
{
    string filename = QFileDialog::getOpenFileName(this,
                                                   tr("Open Mesh"),
                                                   "",
                                                   tr("Mesh Files (*.obj *.ply)")).toStdString();

    return filename;
}

string MainWindow::launchImportMetricDialog()
{
    string filename = QFileDialog::getOpenFileName(this,
                                                   tr("Open Metric"),
                                                   "",
                                                   tr("Metric Files (*.g)")).toStdString();

    return filename;
}

void MainWindow::repaintMesh()
{
    updateGL();
}

void MainWindow::updateGL()
{
    ui->GLwidget->updateGL();
}

void MainWindow::setParameters(ProblemParameters params)
{
    ui->youngsModulusEdit->setText(QString::number(params.YoungsModulus));
    ui->poissonRatioEdit->setText(QString::number(params.PoissonRatio));
    ui->thicknessEdit->setText(QString::number(params.h));
    ui->wireframeCheckBox->setChecked(params.showWireframe);
    ui->smoothShadeCheckBox->setChecked(params.smoothShade);
    ui->densityEdit->setText(QString::number(params.rho));
    ui->dampingCoeffEdit->setText(QString::number(params.dampingCoeff));
    ui->simTimeEdit->setText(QString::number(params.simTime));
    ui->eulerTimestepEdit->setText(QString::number(params.eulerTimestep));
    ui->scaleEdit->setText(QString::number(params.scale));
    ui->outputEdit->setText(QString::fromStdString(params.outputDir));
    ui->constantPressureRadio->setChecked(params.pressureBehavior == params.PB_CONSTANT);
    ui->scalingPressureRadio->setChecked(params.pressureBehavior == params.PB_SCALING);
    ui->holeRadiusRadio->setChecked(params.pressureBehavior == params.PB_LEAKY);
    ui->constantPressureEdit->setText(QString::number(params.constantPressureVal));
    ui->holeRadiusEdit->setText(QString::number(params.holeRadius));
    ui->coneAngleEdit->setText(QString::number(params.coneAngle));
    ui->coneHeightEdit->setText(QString::number(params.coneHeight));
    ui->crushTimeEdit->setText(QString::number(params.crushTime));
    ui->constantVelocityRadio->setChecked(params.constantVelocity);
    ui->simulatedCrushingRadio->setChecked(!params.constantVelocity);
    ui->crushMassEdit->setText(QString::number(params.crushMass));
    ui->initialVelEdit->setText(QString::number(params.initialVel));
    ui->restoreEdit->setText(QString::fromStdString(params.restoreCheckpoint));
}

ProblemParameters MainWindow::getParameters()
{
    ProblemParameters result;
    result.YoungsModulus = ui->youngsModulusEdit->text().toDouble();
    result.PoissonRatio  = ui->poissonRatioEdit->text().toDouble();
    result.h = ui->thicknessEdit->text().toDouble();
    result.showWireframe = ui->wireframeCheckBox->isChecked();
    result.smoothShade = ui->smoothShadeCheckBox->isChecked();
    result.rho = ui->densityEdit->text().toDouble();
    result.dampingCoeff = ui->dampingCoeffEdit->text().toDouble();
    result.eulerTimestep = ui->eulerTimestepEdit->text().toDouble();
    result.simTime = ui->simTimeEdit->text().toDouble();
    result.scale = ui->scaleEdit->text().toDouble();
    result.outputDir = ui->outputEdit->text().toStdString();
    if(ui->constantPressureRadio->isChecked())
    {
        result.pressureBehavior = result.PB_CONSTANT;
    }
    else if(ui->scalingPressureRadio->isChecked())
    {
        result.pressureBehavior = result.PB_SCALING;
    }
    else
        result.pressureBehavior = result.PB_LEAKY;

    result.constantPressureVal = ui->constantPressureEdit->text().toDouble();
    result.holeRadius = ui->holeRadiusEdit->text().toDouble();
    result.coneAngle = ui->coneAngleEdit->text().toDouble();
    result.coneHeight = ui->coneHeightEdit->text().toDouble();
    result.crushTime = ui->crushTimeEdit->text().toDouble();
    result.constantVelocity = ui->constantPressureRadio->isChecked();
    result.crushMass = ui->crushMassEdit->text().toDouble();
    result.initialVel = ui->initialVelEdit->text().toDouble();
    result.restoreCheckpoint = ui->restoreEdit->text().toStdString();
    return result;
}

void MainWindow::on_actionExit_triggered()
{
    assert(cont_);
    QMetaObject::invokeMethod(cont_, "quit");
}

void MainWindow::on_actionReset_Camera_triggered()
{
    QMetaObject::invokeMethod(cont_, "centerCamera");
    updateGL();
}

void MainWindow::on_actionTake_Screenshot_triggered()
{
    saveScreenshot();
    updateGL();
}

void MainWindow::on_wireframeCheckBox_clicked()
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_smoothShadeCheckBox_clicked()
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_actionExport_OBJ_triggered()
{
    QFileDialog savedialog(this, "Export 3D Geometry", ".", "Mesh Files (*.obj)");
    savedialog.setFileMode(QFileDialog::AnyFile);
    savedialog.setDefaultSuffix("obj");
    savedialog.setViewMode(QFileDialog::List);
    savedialog.setAcceptMode(QFileDialog::AcceptSave);
    if(savedialog.exec())
    {
        QStringList filenames = savedialog.selectedFiles();
        if(filenames.size() > 0)
        {
            QString filename = filenames[0];
            QMetaObject::invokeMethod(cont_, "exportOBJ", Q_ARG(std::string, filename.toStdString()));
        }
    }
}

void MainWindow::on_actionImport_OBJ_triggered()
{
    string filename = launchImportOBJDialog();
    QMetaObject::invokeMethod(cont_, "importOBJ", Q_ARG(std::string, filename));
}

void MainWindow::on_poissonRatioEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_youngsModulusEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_thicknessEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_densityEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_dampingCoeffEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_eulerTimestepEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_scaleEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_outputEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_actionImport_Metric_triggered()
{
    string filename = launchImportMetricDialog();
    QMetaObject::invokeMethod(cont_, "importMetric", Q_ARG(std::string, filename));
}

void MainWindow::on_actionMake_Cone_triggered()
{
    QMetaObject::invokeMethod(cont_, "makeCone");
}

void MainWindow::on_crushButton_clicked()
{
    QMetaObject::invokeMethod(cont_, "crush");
}

void MainWindow::on_constantPressureRadio_toggled(bool )
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_scalingPressureRadio_toggled(bool )
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_constantPressureEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_airLeakEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_coneAngleEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_coneHeightEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_crushTimeEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_constantVelocityRadio_toggled(bool )
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_simulatedCrushingRadio_toggled(bool )
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_crushMassEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_simTimeEdit_textChanged(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_initialVelEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_restoreEdit_textEdited(const QString &)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_holeRadiusRadio_toggled(bool )
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}

void MainWindow::on_holeRadiusEdit_textEdited(const QString &arg1)
{
    QMetaObject::invokeMethod(cont_, "updateParameters", Q_ARG(ProblemParameters, getParameters()));
}
