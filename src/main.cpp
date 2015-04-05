#include "mainwindow.h"
#include "controller.h"
#include <QApplication>
#include <QDesktopWidget>
#include <QThread>
#include <QMetaClassInfo>
#include <string>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    qRegisterMetaType<std::string>("std::string");
    qRegisterMetaType<Eigen::Vector3d>("Eigen::Vector3d");
    qRegisterMetaType<ProblemParameters>("ProblemParameters");

    QThread workthread;
    workthread.start();
    MainWindow window;
    Controller cont(window);
    cont.moveToThread(&workthread);

    window.setController(cont);

    int desktopArea = QApplication::desktop()->width() *
                     QApplication::desktop()->height();
    int widgetArea = window.width() * window.height();
    if (((float)widgetArea / (float)desktopArea) < 0.75f)
        window.show();
    else
        window.showMaximized();

    int ret = app.exec();
    workthread.quit();
    workthread.wait();
    return ret;
}
