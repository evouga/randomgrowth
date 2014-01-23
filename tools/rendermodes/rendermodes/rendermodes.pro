#-------------------------------------------------
#
# Project created by QtCreator 2014-01-13T16:45:58
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = rendermodes
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

INCLUDEPATH += ../../../ext/openmesh ../../../ext/eigen
QMAKE_LIBDIR += ../../../ext/openmesh/build/Build/lib/OpenMesh
QMAKE_CXXFLAGS += -g
LIBS += -Wl,-rpath=../../../ext/openmesh/build/Build/lib/OpenMesh -lOpenMeshCore


SOURCES += main.cpp
