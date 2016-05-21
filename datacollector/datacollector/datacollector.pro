#-------------------------------------------------
#
# Project created by QtCreator 2016-04-13T17:57:40
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = datacollector
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

INCLUDEPATH += ../../src
INCLUDEPATH += ../../ext/eigen


SOURCES += main.cpp \
    ../../src/midedge.cpp \
    ../../src/simulationmesh.cpp \
    ../../src/mesh.cpp

HEADERS += \
    ../../src/midedge.h \
    ../../src/simulationmesh.h \
    ../../src/mesh.h
