/********************************************************************************
** Form generated from reading UI file 'coronarygui.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CORONARYGUI_H
#define UI_CORONARYGUI_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QProgressBar>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QTextEdit>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>
#include "QVTKWidget.h"

QT_BEGIN_NAMESPACE

class Ui_coronaryGUI
{
public:
    QAction *action_Open_file;
    QAction *action_Quit;
    QAction *actionGMF;
    QAction *actionSkeletonize;
    QAction *actionRidler_and_Calvard;
    QAction *actionOpen_file_as_ground_truth;
    QWidget *centralWidget;
    QGridLayout *gridLayout_2;
    QProgressBar *pbarLog;
    QGridLayout *gridLayout;
    QLabel *lblVP1;
    QLabel *lblVP2;
    QLabel *lblVP4;
    QTabWidget *tabWidget;
    QWidget *tab;
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout;
    QLabel *lblVP3;
    QWidget *tab_2;
    QWidget *gridLayoutWidget;
    QGridLayout *gridLayout_4;
    QVTKWidget *qvtkHist;
    QTextEdit *txtLog;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuDetection;
    QMenu *menuSegmentation;
    QMenu *menu3D;
    QMenu *menuStenosis_detection;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *coronaryGUI)
    {
        if (coronaryGUI->objectName().isEmpty())
            coronaryGUI->setObjectName(QString::fromUtf8("coronaryGUI"));
        coronaryGUI->resize(1185, 687);
        coronaryGUI->setContextMenuPolicy(Qt::CustomContextMenu);
        action_Open_file = new QAction(coronaryGUI);
        action_Open_file->setObjectName(QString::fromUtf8("action_Open_file"));
        action_Quit = new QAction(coronaryGUI);
        action_Quit->setObjectName(QString::fromUtf8("action_Quit"));
        actionGMF = new QAction(coronaryGUI);
        actionGMF->setObjectName(QString::fromUtf8("actionGMF"));
        actionSkeletonize = new QAction(coronaryGUI);
        actionSkeletonize->setObjectName(QString::fromUtf8("actionSkeletonize"));
        actionRidler_and_Calvard = new QAction(coronaryGUI);
        actionRidler_and_Calvard->setObjectName(QString::fromUtf8("actionRidler_and_Calvard"));
        actionOpen_file_as_ground_truth = new QAction(coronaryGUI);
        actionOpen_file_as_ground_truth->setObjectName(QString::fromUtf8("actionOpen_file_as_ground_truth"));
        centralWidget = new QWidget(coronaryGUI);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        gridLayout_2 = new QGridLayout(centralWidget);
        gridLayout_2->setSpacing(6);
        gridLayout_2->setContentsMargins(11, 11, 11, 11);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        pbarLog = new QProgressBar(centralWidget);
        pbarLog->setObjectName(QString::fromUtf8("pbarLog"));
        pbarLog->setMaximumSize(QSize(16777215, 10));
        QFont font;
        font.setPointSize(6);
        pbarLog->setFont(font);
        pbarLog->setValue(0);

        gridLayout_2->addWidget(pbarLog, 11, 0, 1, 1);

        gridLayout = new QGridLayout();
        gridLayout->setSpacing(6);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        lblVP1 = new QLabel(centralWidget);
        lblVP1->setObjectName(QString::fromUtf8("lblVP1"));
        lblVP1->setMinimumSize(QSize(256, 256));
        lblVP1->setMaximumSize(QSize(512, 512));

        gridLayout->addWidget(lblVP1, 0, 0, 1, 1);

        lblVP2 = new QLabel(centralWidget);
        lblVP2->setObjectName(QString::fromUtf8("lblVP2"));
        lblVP2->setMinimumSize(QSize(256, 256));

        gridLayout->addWidget(lblVP2, 0, 1, 1, 1);

        lblVP4 = new QLabel(centralWidget);
        lblVP4->setObjectName(QString::fromUtf8("lblVP4"));
        lblVP4->setMinimumSize(QSize(256, 256));
        lblVP4->setMaximumSize(QSize(512, 512));

        gridLayout->addWidget(lblVP4, 2, 1, 1, 1);

        tabWidget = new QTabWidget(centralWidget);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabWidget->setMaximumSize(QSize(572, 512));
        tabWidget->setTabPosition(QTabWidget::West);
        tabWidget->setTabShape(QTabWidget::Triangular);
        tab = new QWidget();
        tab->setObjectName(QString::fromUtf8("tab"));
        horizontalLayoutWidget = new QWidget(tab);
        horizontalLayoutWidget->setObjectName(QString::fromUtf8("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(10, 10, 531, 241));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        lblVP3 = new QLabel(horizontalLayoutWidget);
        lblVP3->setObjectName(QString::fromUtf8("lblVP3"));
        lblVP3->setMinimumSize(QSize(256, 256));
        lblVP3->setMaximumSize(QSize(512, 512));

        horizontalLayout->addWidget(lblVP3);

        tabWidget->addTab(tab, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QString::fromUtf8("tab_2"));
        gridLayoutWidget = new QWidget(tab_2);
        gridLayoutWidget->setObjectName(QString::fromUtf8("gridLayoutWidget"));
        gridLayoutWidget->setGeometry(QRect(-1, -1, 531, 461));
        gridLayout_4 = new QGridLayout(gridLayoutWidget);
        gridLayout_4->setSpacing(6);
        gridLayout_4->setContentsMargins(11, 11, 11, 11);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        gridLayout_4->setContentsMargins(0, 0, 0, 0);
        qvtkHist = new QVTKWidget(gridLayoutWidget);
        qvtkHist->setObjectName(QString::fromUtf8("qvtkHist"));
        qvtkHist->setMinimumSize(QSize(256, 256));

        gridLayout_4->addWidget(qvtkHist, 0, 0, 1, 1);

        tabWidget->addTab(tab_2, QString());

        gridLayout->addWidget(tabWidget, 2, 0, 1, 1);


        gridLayout_2->addLayout(gridLayout, 9, 0, 1, 1);

        txtLog = new QTextEdit(centralWidget);
        txtLog->setObjectName(QString::fromUtf8("txtLog"));
        txtLog->setMaximumSize(QSize(16777215, 100));
        txtLog->setContextMenuPolicy(Qt::CustomContextMenu);

        gridLayout_2->addWidget(txtLog, 10, 0, 1, 1);

        coronaryGUI->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(coronaryGUI);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1185, 20));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuDetection = new QMenu(menuBar);
        menuDetection->setObjectName(QString::fromUtf8("menuDetection"));
        menuSegmentation = new QMenu(menuBar);
        menuSegmentation->setObjectName(QString::fromUtf8("menuSegmentation"));
        menu3D = new QMenu(menuBar);
        menu3D->setObjectName(QString::fromUtf8("menu3D"));
        menuStenosis_detection = new QMenu(menuBar);
        menuStenosis_detection->setObjectName(QString::fromUtf8("menuStenosis_detection"));
        coronaryGUI->setMenuBar(menuBar);
        mainToolBar = new QToolBar(coronaryGUI);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        coronaryGUI->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(coronaryGUI);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        coronaryGUI->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuDetection->menuAction());
        menuBar->addAction(menuSegmentation->menuAction());
        menuBar->addAction(menuStenosis_detection->menuAction());
        menuBar->addAction(menu3D->menuAction());
        menuFile->addSeparator();
        menuFile->addAction(action_Open_file);
        menuFile->addAction(actionOpen_file_as_ground_truth);
        menuFile->addSeparator();
        menuFile->addAction(action_Quit);
        menuDetection->addSeparator();
        menuDetection->addAction(actionGMF);
        menuDetection->addSeparator();
        menuSegmentation->addAction(actionRidler_and_Calvard);
        menu3D->addAction(actionSkeletonize);

        retranslateUi(coronaryGUI);

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(coronaryGUI);
    } // setupUi

    void retranslateUi(QMainWindow *coronaryGUI)
    {
        coronaryGUI->setWindowTitle(QApplication::translate("coronaryGUI", "coronaryGUI", 0, QApplication::UnicodeUTF8));
        action_Open_file->setText(QApplication::translate("coronaryGUI", "&Open file as input", 0, QApplication::UnicodeUTF8));
        action_Open_file->setShortcut(QApplication::translate("coronaryGUI", "Ctrl+O", 0, QApplication::UnicodeUTF8));
        action_Quit->setText(QApplication::translate("coronaryGUI", "&Quit", 0, QApplication::UnicodeUTF8));
        action_Quit->setShortcut(QApplication::translate("coronaryGUI", "Ctrl+Q", 0, QApplication::UnicodeUTF8));
        actionGMF->setText(QApplication::translate("coronaryGUI", "GMF", 0, QApplication::UnicodeUTF8));
        actionSkeletonize->setText(QApplication::translate("coronaryGUI", "Skeletonize", 0, QApplication::UnicodeUTF8));
        actionRidler_and_Calvard->setText(QApplication::translate("coronaryGUI", "Ridler and Calvard", 0, QApplication::UnicodeUTF8));
        actionOpen_file_as_ground_truth->setText(QApplication::translate("coronaryGUI", "Open file as ground-truth", 0, QApplication::UnicodeUTF8));
        lblVP1->setText(QApplication::translate("coronaryGUI", "TextLabel", 0, QApplication::UnicodeUTF8));
        lblVP2->setText(QApplication::translate("coronaryGUI", "TextLabel", 0, QApplication::UnicodeUTF8));
        lblVP4->setText(QApplication::translate("coronaryGUI", "TextLabel", 0, QApplication::UnicodeUTF8));
        lblVP3->setText(QApplication::translate("coronaryGUI", "TextLabel", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("coronaryGUI", "Tab 1", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_2), QApplication::translate("coronaryGUI", "Tab 2", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("coronaryGUI", "File", 0, QApplication::UnicodeUTF8));
        menuDetection->setTitle(QApplication::translate("coronaryGUI", "Detection", 0, QApplication::UnicodeUTF8));
        menuSegmentation->setTitle(QApplication::translate("coronaryGUI", "Segmentation", 0, QApplication::UnicodeUTF8));
        menu3D->setTitle(QApplication::translate("coronaryGUI", "3D", 0, QApplication::UnicodeUTF8));
        menuStenosis_detection->setTitle(QApplication::translate("coronaryGUI", "Lesion detection", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class coronaryGUI: public Ui_coronaryGUI {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CORONARYGUI_H
