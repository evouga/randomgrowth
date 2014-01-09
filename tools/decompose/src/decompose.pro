#-------------------------------------------------
#
# Project created by QtCreator 2014-01-09T15:08:21
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = decompose
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

INCLUDEPATH += ../../../ext/openmesh ../../../ext/eigen
QMAKE_LIBDIR += ../../../ext/openmesh/build/Build/lib/OpenMesh
QMAKE_CXXFLAGS += -g
LIBS += -Wl,-rpath=../../../ext/openmesh/build/Build/lib/OpenMesh -lOpenMeshCore


SOURCES += main.cpp
