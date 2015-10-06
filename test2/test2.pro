#-------------------------------------------------
#
# Project created by QtCreator 2015-07-14T18:46:49
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = test2
TEMPLATE = app


SOURCES += main.cpp\
        widget.cpp \
    plotter.cpp \
    myemd.cpp \
    splineclass.cpp \
    boundaryclass.cpp \
    extrzeroclass.cpp \
    emdclass.cpp \
    fileclass.cpp \
    hilbertclass.cpp \
    fourierclass.cpp

HEADERS  += widget.h \
    plotter.h \
    splineclass.h \
    boundaryclass.h \
    extrzeroclass.h \
    emdclass.h \
    fileclass.h \
    hilbertclass.h \
    fourierclass.h

FORMS    += widget.ui
CONFIG += console
RESOURCES += plotter.qrc
