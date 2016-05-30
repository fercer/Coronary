#-------------------------------------------------
#
# Project created by QtCreator 2016-05-25T15:25:03
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = CoronaryGUI
TEMPLATE = app

INCLUDEPATH += /usr/local/include/vtk-6.3/

SOURCES += main.cpp\
        coronarygui.cpp

HEADERS  += coronarygui.h

FORMS    += coronarygui.ui
