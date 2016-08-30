/********************************************************************************
** Form generated from reading UI file 'dialogbayes.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DIALOGBAYES_H
#define UI_DIALOGBAYES_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QHeaderView>

QT_BEGIN_NAMESPACE

class Ui_DialogBayes
{
public:

    void setupUi(QDialog *DialogBayes)
    {
        if (DialogBayes->objectName().isEmpty())
            DialogBayes->setObjectName(QString::fromUtf8("DialogBayes"));
        DialogBayes->resize(480, 640);

        retranslateUi(DialogBayes);

        QMetaObject::connectSlotsByName(DialogBayes);
    } // setupUi

    void retranslateUi(QDialog *DialogBayes)
    {
        DialogBayes->setWindowTitle(QApplication::translate("DialogBayes", "Dialog", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogBayes: public Ui_DialogBayes {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DIALOGBAYES_H
