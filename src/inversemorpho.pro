#-------------------------------------------------
#
# Project created by QtCreator 2013-09-26T16:48:11
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = inversemorpho
TEMPLATE = app

INCLUDEPATH += ../ext/openmesh ../ext/eigen
QMAKE_LIBDIR += ../ext/openmesh/build/Build/lib/OpenMesh
LIBS += -lOpenMeshCore


SOURCES += main.cpp\
        mainwindow.cpp \
    mesh.cpp

HEADERS  += mainwindow.h \
    mesh.h \
    autodifftemplates.h

FORMS    += mainwindow.ui
