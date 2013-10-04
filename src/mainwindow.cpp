#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "controller.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QDir>
#include <QDateTime>
#include <QPixmap>

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    cont_(NULL)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setController(Controller &cont)
{
    cont_ = &cont;
    ui->GLwidget->setController(cont);
}

void MainWindow::showError(const string &error)
{
    QMessageBox::warning(this, tr("VViewer"),
                                    QString(error.c_str()),
                                    QMessageBox::Ok, QMessageBox::NoButton);
}

void MainWindow::centerCamera()
{
    ui->GLwidget->centerCamera();
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

bool MainWindow::showWireframe() const
{
    return ui->wireframeCheckBox->isChecked();
}

bool MainWindow::smoothShade() const
{
    return ui->smoothShadeCheckBox->isChecked();
}

void MainWindow::repaintMesh()
{
    updateGL();
}

void MainWindow::updateGL()
{
    ui->GLwidget->updateGL();
}

void MainWindow::setParameters(const ProblemParameters &params)
{
    ui->youngsModulusEdit->setText(QString::number(params.YoungsModulus));
    ui->poissonRatioEdit->setText(QString::number(params.PoissonRatio));
    ui->thicknessEdit->setText(QString::number(params.h));
    ui->maxitersEdit->setText(QString::number(params.maxiters));
    ui->maxlsitersEdit->setText(QString::number(params.maxlinesearchiters));
    ui->tolEdit->setText(QString::number(params.tol));
    ui->maxpoweritersEdit->setText(QString::number(params.maxpoweriters));
    ui->poweritertolEdit->setText(QString::number(params.powertol));
}

ProblemParameters MainWindow::getParameters()
{
    ProblemParameters result;
    result.YoungsModulus = ui->youngsModulusEdit->text().toDouble();
    result.PoissonRatio  = ui->poissonRatioEdit->text().toDouble();
    result.h = ui->thicknessEdit->text().toDouble();
    result.maxiters = ui->maxitersEdit->text().toInt();
    result.maxlinesearchiters = ui->maxlsitersEdit->text().toInt();
    result.tol = ui->tolEdit->text().toDouble();
    result.powertol = ui->poweritertolEdit->text().toDouble();
    result.maxpoweriters = ui->maxpoweritersEdit->text().toInt();
    return result;
}

void MainWindow::on_actionExit_triggered()
{
    assert(cont_);
    cont_->quit();
}

void MainWindow::on_actionReset_Camera_triggered()
{
    centerCamera();
    updateGL();
}

void MainWindow::on_actionTake_Screenshot_triggered()
{
    saveScreenshot();
    updateGL();
}

void MainWindow::on_wireframeCheckBox_clicked()
{
    updateGL();
}

void MainWindow::on_smoothShadeCheckBox_clicked()
{
    updateGL();
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
            cont_->exportOBJ(filename.toStdString().c_str());
        }
    }
}

void MainWindow::on_actionImport_OBJ_triggered()
{
    string filename = launchImportOBJDialog();
    cont_->importOBJ(filename.c_str());
    updateGL();
}

void MainWindow::on_poissonRatioEdit_textEdited(const QString &)
{
    cont_->updateParameters();
}

void MainWindow::on_youngsModulusEdit_textEdited(const QString &)
{
    cont_->updateParameters();
}

void MainWindow::on_thicknessEdit_textEdited(const QString &)
{
    cont_->updateParameters();
}

void MainWindow::on_findMetricButton_clicked()
{
    cont_->findMetric();
    updateGL();
}

void MainWindow::on_maxitersEdit_textEdited(const QString &)
{
    cont_->updateParameters();
}

void MainWindow::on_maxlsitersEdit_textEdited(const QString &)
{
    cont_->updateParameters();
}

void MainWindow::on_tolEdit_textEdited(const QString &)
{
    cont_->updateParameters();
}

void MainWindow::on_relaxEmbeddingButton_clicked()
{
    cont_->relaxEmbedding();
}

void MainWindow::on_maxpoweritersEdit_textEdited(const QString &)
{
    cont_->updateParameters();
}

void MainWindow::on_poweritertolEdit_textEdited(const QString &)
{
    cont_->updateParameters();
}
