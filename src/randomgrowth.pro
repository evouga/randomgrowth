#-------------------------------------------------
#
# Project created by QtCreator 2013-09-26T16:48:11
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = randomgrowth
TEMPLATE = app

INCLUDEPATH += ../ext/eigen
QMAKE_CXX = g++-4.8
QMAKE_CXXFLAGS += -g -fopenmp -std=c++11
LIBS += -lpng -lGL -lGLU -fopenmp


SOURCES += main.cpp\
        mainwindow.cpp \
    glwidget.cpp \
    zoomer.cpp \
    translator.cpp \
    rotator.cpp \
    camera.cpp \
    yimage.cpp \
    controller.cpp \
    mesh-rendering.cpp \
    mesh-optimization.cpp \
    midedge.cpp \
    simulationmesh.cpp \
    mesh.cpp

HEADERS  += mainwindow.h \
    glwidget.h \
    zoomer.h \
    translator.h \
    rotator.h \
    camera.h \
    yimage.h \
    controller.h \
    midedge.h \
    simulationmesh.h \
    mesh.h

FORMS    += mainwindow.ui
