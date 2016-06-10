#-------------------------------------------------
#
# Project created by QtCreator 2016-05-25T15:25:03
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = CoronaryGUI
TEMPLATE = app

INCLUDEPATH += /usr/local/include/vtk-6.3/ \
               /usr/include/gdcm-2.2/ \
               /usr/local/include/ \
               ../src/   

SOURCES += main.cpp\
           coronarygui.cpp\
           ../src/IMGVTK.cpp\
           ../src/reconstructor_3D.cpp\
           ../src/filtros.cpp

HEADERS  += coronarygui.h\
            ../src/IMGVTK.h\
            ../src/reconstructor_3D.h\
            ../src/filtros.h

FORMS    += coronarygui.ui
