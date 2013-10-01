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

void MainWindow::setParameters(double Youngs, double Poisson, double h)
{
    ui->youngsModulusEdit->setText(QString::number(Youngs));
    ui->poissonRatioEdit->setText(QString::number(Poisson));
    ui->thicknessEdit->setText(QString::number(h));
}

void MainWindow::getParameters(double &Youngs, double &Poisson, double &h)
{
    Youngs = ui->youngsModulusEdit->text().toDouble();
    Poisson = ui->poissonRatioEdit->text().toDouble();
    h = ui->thicknessEdit->text().toDouble();
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
