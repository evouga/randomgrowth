#-------------------------------------------------
#
# Project created by QtCreator 2013-09-26T16:48:11
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = midedge
TEMPLATE = app

INCLUDEPATH += ../ext/openmesh ../ext/eigen /usr/include/suitesparse
QMAKE_LIBDIR += ../ext/openmesh/build/Build/lib/OpenMesh
QMAKE_CXXFLAGS += -g -fno-omit-frame-pointer
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3 #-DNDEBUG
LIBS += -lOpenMeshCore -lpng -lGL -lGLU -lspqr -lcholmod


SOURCES += main.cpp\
        mainwindow.cpp \
    mesh.cpp \
    glwidget.cpp \
    zoomer.cpp \
    translator.cpp \
    rotator.cpp \
    camera.cpp \
    yimage.cpp \
    controller.cpp \
    mesh-rendering.cpp \
    mesh-optimization.cpp \
    elasticenergy.cpp \
    midedge.cpp \
    robustleastsquares.cpp

HEADERS  += mainwindow.h \
    mesh.h \
    glwidget.h \
    zoomer.h \
    translator.h \
    rotator.h \
    camera.h \
    yimage.h \
    controller.h \
    elasticenergy.h \
    omtypes.h \
    midedge.h \
    robustleastsquares.h

FORMS    += mainwindow.ui
