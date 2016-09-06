/********************************************************************************
** Form generated from reading UI file 'filterpars.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_FILTERPARS_H
#define UI_FILTERPARS_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>

QT_BEGIN_NAMESPACE

class Ui_filterpars
{
public:
    QGridLayout *gridLayout;
    QPushButton *pushButton_3;
    QComboBox *cmbFilter;
    QFrame *line;
    QPushButton *pushButton;
    QGridLayout *gridLayout_5;
    QLabel *label_2;
    QCheckBox *chkLFix;
    QGridLayout *gridLayout_8;
    QDoubleSpinBox *dspnL;
    QGridLayout *gridLayout_6;
    QDoubleSpinBox *dspnLUpp;
    QLabel *lblLLow;
    QDoubleSpinBox *dspnLLow;
    QLabel *lblLUpp;
    QLabel *lblLDel;
    QDoubleSpinBox *dspnLDel;
    QLabel *lblLStat;
    QCheckBox *chkTune;
    QGridLayout *gridLayout_2;
    QLabel *lblTemplate;
    QLabel *label;
    QPushButton *pushButton_2;
    QFrame *line_2;
    QGridLayout *gridLayout_4;
    QPushButton *btnSetOMet;
    QSpacerItem *horizontalSpacer;
    QComboBox *cmbOMethod;
    QLabel *lblOMethod;
    QGridLayout *gridLayout_12;
    QLabel *label_4;
    QCheckBox *chkSFix;
    QGridLayout *gridLayout_13;
    QDoubleSpinBox *dspnS;
    QGridLayout *gridLayout_14;
    QDoubleSpinBox *dspnSUpp;
    QLabel *lblSLow;
    QDoubleSpinBox *dspnSLow;
    QLabel *lblSUpp;
    QLabel *lblSDel;
    QDoubleSpinBox *dspnSDel;
    QLabel *lblSStat;
    QGridLayout *gridLayout_9;
    QLabel *label_3;
    QCheckBox *chkTFix;
    QGridLayout *gridLayout_10;
    QDoubleSpinBox *dspnT;
    QGridLayout *gridLayout_11;
    QDoubleSpinBox *dspnTUpp;
    QLabel *lblTLow;
    QDoubleSpinBox *dspnTLow;
    QLabel *lblTUpp;
    QLabel *lblTDel;
    QDoubleSpinBox *dspnTDel;
    QLabel *lblTStat;
    QSpacerItem *verticalSpacer;

    void setupUi(QDialog *filterpars)
    {
        if (filterpars->objectName().isEmpty())
            filterpars->setObjectName(QString::fromUtf8("filterpars"));
        filterpars->resize(695, 546);
        gridLayout = new QGridLayout(filterpars);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        pushButton_3 = new QPushButton(filterpars);
        pushButton_3->setObjectName(QString::fromUtf8("pushButton_3"));

        gridLayout->addWidget(pushButton_3, 9, 1, 1, 1);

        cmbFilter = new QComboBox(filterpars);
        cmbFilter->setObjectName(QString::fromUtf8("cmbFilter"));

        gridLayout->addWidget(cmbFilter, 2, 0, 1, 2);

        line = new QFrame(filterpars);
        line->setObjectName(QString::fromUtf8("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        gridLayout->addWidget(line, 3, 0, 1, 3);

        pushButton = new QPushButton(filterpars);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));

        gridLayout->addWidget(pushButton, 9, 3, 1, 1);

        gridLayout_5 = new QGridLayout();
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        label_2 = new QLabel(filterpars);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        QFont font;
        font.setPointSize(12);
        label_2->setFont(font);

        gridLayout_5->addWidget(label_2, 2, 0, 1, 1);

        chkLFix = new QCheckBox(filterpars);
        chkLFix->setObjectName(QString::fromUtf8("chkLFix"));
        chkLFix->setChecked(true);

        gridLayout_5->addWidget(chkLFix, 3, 0, 1, 1);

        gridLayout_8 = new QGridLayout();
        gridLayout_8->setObjectName(QString::fromUtf8("gridLayout_8"));
        dspnL = new QDoubleSpinBox(filterpars);
        dspnL->setObjectName(QString::fromUtf8("dspnL"));

        gridLayout_8->addWidget(dspnL, 0, 0, 2, 1);


        gridLayout_5->addLayout(gridLayout_8, 4, 1, 2, 1);

        gridLayout_6 = new QGridLayout();
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        dspnLUpp = new QDoubleSpinBox(filterpars);
        dspnLUpp->setObjectName(QString::fromUtf8("dspnLUpp"));

        gridLayout_6->addWidget(dspnLUpp, 2, 1, 1, 1);

        lblLLow = new QLabel(filterpars);
        lblLLow->setObjectName(QString::fromUtf8("lblLLow"));
        lblLLow->setMaximumSize(QSize(40, 16777215));

        gridLayout_6->addWidget(lblLLow, 1, 0, 1, 1);

        dspnLLow = new QDoubleSpinBox(filterpars);
        dspnLLow->setObjectName(QString::fromUtf8("dspnLLow"));

        gridLayout_6->addWidget(dspnLLow, 1, 1, 1, 1);

        lblLUpp = new QLabel(filterpars);
        lblLUpp->setObjectName(QString::fromUtf8("lblLUpp"));

        gridLayout_6->addWidget(lblLUpp, 2, 0, 1, 1);

        lblLDel = new QLabel(filterpars);
        lblLDel->setObjectName(QString::fromUtf8("lblLDel"));
        lblLDel->setMaximumSize(QSize(60, 16777215));

        gridLayout_6->addWidget(lblLDel, 1, 2, 1, 1);

        dspnLDel = new QDoubleSpinBox(filterpars);
        dspnLDel->setObjectName(QString::fromUtf8("dspnLDel"));

        gridLayout_6->addWidget(dspnLDel, 2, 2, 1, 1);


        gridLayout_5->addLayout(gridLayout_6, 4, 2, 2, 1);

        lblLStat = new QLabel(filterpars);
        lblLStat->setObjectName(QString::fromUtf8("lblLStat"));

        gridLayout_5->addWidget(lblLStat, 4, 0, 2, 1);


        gridLayout->addLayout(gridLayout_5, 6, 0, 1, 3);

        chkTune = new QCheckBox(filterpars);
        chkTune->setObjectName(QString::fromUtf8("chkTune"));
        chkTune->setMaximumSize(QSize(80, 16777215));

        gridLayout->addWidget(chkTune, 2, 2, 1, 1);

        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        lblTemplate = new QLabel(filterpars);
        lblTemplate->setObjectName(QString::fromUtf8("lblTemplate"));
        lblTemplate->setFrameShape(QFrame::Box);
        lblTemplate->setFrameShadow(QFrame::Sunken);
        lblTemplate->setAlignment(Qt::AlignCenter);

        gridLayout_2->addWidget(lblTemplate, 0, 0, 1, 1);


        gridLayout->addLayout(gridLayout_2, 5, 3, 3, 1);

        label = new QLabel(filterpars);
        label->setObjectName(QString::fromUtf8("label"));
        QFont font1;
        font1.setPointSize(14);
        label->setFont(font1);

        gridLayout->addWidget(label, 1, 0, 1, 2);

        pushButton_2 = new QPushButton(filterpars);
        pushButton_2->setObjectName(QString::fromUtf8("pushButton_2"));

        gridLayout->addWidget(pushButton_2, 9, 0, 1, 1);

        line_2 = new QFrame(filterpars);
        line_2->setObjectName(QString::fromUtf8("line_2"));
        line_2->setFrameShape(QFrame::HLine);
        line_2->setFrameShadow(QFrame::Sunken);

        gridLayout->addWidget(line_2, 5, 0, 1, 3);

        gridLayout_4 = new QGridLayout();
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        btnSetOMet = new QPushButton(filterpars);
        btnSetOMet->setObjectName(QString::fromUtf8("btnSetOMet"));
        btnSetOMet->setMaximumSize(QSize(16777215, 20));

        gridLayout_4->addWidget(btnSetOMet, 4, 1, 1, 1);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer, 4, 0, 1, 1);

        cmbOMethod = new QComboBox(filterpars);
        cmbOMethod->setObjectName(QString::fromUtf8("cmbOMethod"));
        cmbOMethod->setMaximumSize(QSize(16777215, 20));

        gridLayout_4->addWidget(cmbOMethod, 2, 0, 1, 2);

        lblOMethod = new QLabel(filterpars);
        lblOMethod->setObjectName(QString::fromUtf8("lblOMethod"));
        lblOMethod->setMaximumSize(QSize(16777215, 20));
        lblOMethod->setFont(font);

        gridLayout_4->addWidget(lblOMethod, 1, 0, 1, 1);


        gridLayout->addLayout(gridLayout_4, 4, 0, 1, 3);

        gridLayout_12 = new QGridLayout();
        gridLayout_12->setObjectName(QString::fromUtf8("gridLayout_12"));
        label_4 = new QLabel(filterpars);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setFont(font);

        gridLayout_12->addWidget(label_4, 2, 0, 1, 1);

        chkSFix = new QCheckBox(filterpars);
        chkSFix->setObjectName(QString::fromUtf8("chkSFix"));
        chkSFix->setChecked(true);

        gridLayout_12->addWidget(chkSFix, 3, 0, 1, 1);

        gridLayout_13 = new QGridLayout();
        gridLayout_13->setObjectName(QString::fromUtf8("gridLayout_13"));
        dspnS = new QDoubleSpinBox(filterpars);
        dspnS->setObjectName(QString::fromUtf8("dspnS"));

        gridLayout_13->addWidget(dspnS, 0, 0, 2, 1);


        gridLayout_12->addLayout(gridLayout_13, 4, 1, 2, 1);

        gridLayout_14 = new QGridLayout();
        gridLayout_14->setObjectName(QString::fromUtf8("gridLayout_14"));
        dspnSUpp = new QDoubleSpinBox(filterpars);
        dspnSUpp->setObjectName(QString::fromUtf8("dspnSUpp"));

        gridLayout_14->addWidget(dspnSUpp, 2, 1, 1, 1);

        lblSLow = new QLabel(filterpars);
        lblSLow->setObjectName(QString::fromUtf8("lblSLow"));

        gridLayout_14->addWidget(lblSLow, 1, 0, 1, 1);

        dspnSLow = new QDoubleSpinBox(filterpars);
        dspnSLow->setObjectName(QString::fromUtf8("dspnSLow"));

        gridLayout_14->addWidget(dspnSLow, 1, 1, 1, 1);

        lblSUpp = new QLabel(filterpars);
        lblSUpp->setObjectName(QString::fromUtf8("lblSUpp"));

        gridLayout_14->addWidget(lblSUpp, 2, 0, 1, 1);

        lblSDel = new QLabel(filterpars);
        lblSDel->setObjectName(QString::fromUtf8("lblSDel"));

        gridLayout_14->addWidget(lblSDel, 1, 2, 1, 1);

        dspnSDel = new QDoubleSpinBox(filterpars);
        dspnSDel->setObjectName(QString::fromUtf8("dspnSDel"));

        gridLayout_14->addWidget(dspnSDel, 2, 2, 1, 1);


        gridLayout_12->addLayout(gridLayout_14, 4, 2, 2, 1);

        lblSStat = new QLabel(filterpars);
        lblSStat->setObjectName(QString::fromUtf8("lblSStat"));

        gridLayout_12->addWidget(lblSStat, 4, 0, 2, 1);


        gridLayout->addLayout(gridLayout_12, 8, 0, 1, 3);

        gridLayout_9 = new QGridLayout();
        gridLayout_9->setObjectName(QString::fromUtf8("gridLayout_9"));
        label_3 = new QLabel(filterpars);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setFont(font);

        gridLayout_9->addWidget(label_3, 2, 0, 1, 1);

        chkTFix = new QCheckBox(filterpars);
        chkTFix->setObjectName(QString::fromUtf8("chkTFix"));
        chkTFix->setChecked(true);

        gridLayout_9->addWidget(chkTFix, 3, 0, 1, 1);

        gridLayout_10 = new QGridLayout();
        gridLayout_10->setObjectName(QString::fromUtf8("gridLayout_10"));
        dspnT = new QDoubleSpinBox(filterpars);
        dspnT->setObjectName(QString::fromUtf8("dspnT"));

        gridLayout_10->addWidget(dspnT, 0, 0, 2, 1);


        gridLayout_9->addLayout(gridLayout_10, 4, 1, 2, 1);

        gridLayout_11 = new QGridLayout();
        gridLayout_11->setObjectName(QString::fromUtf8("gridLayout_11"));
        dspnTUpp = new QDoubleSpinBox(filterpars);
        dspnTUpp->setObjectName(QString::fromUtf8("dspnTUpp"));

        gridLayout_11->addWidget(dspnTUpp, 2, 1, 1, 1);

        lblTLow = new QLabel(filterpars);
        lblTLow->setObjectName(QString::fromUtf8("lblTLow"));

        gridLayout_11->addWidget(lblTLow, 1, 0, 1, 1);

        dspnTLow = new QDoubleSpinBox(filterpars);
        dspnTLow->setObjectName(QString::fromUtf8("dspnTLow"));

        gridLayout_11->addWidget(dspnTLow, 1, 1, 1, 1);

        lblTUpp = new QLabel(filterpars);
        lblTUpp->setObjectName(QString::fromUtf8("lblTUpp"));

        gridLayout_11->addWidget(lblTUpp, 2, 0, 1, 1);

        lblTDel = new QLabel(filterpars);
        lblTDel->setObjectName(QString::fromUtf8("lblTDel"));

        gridLayout_11->addWidget(lblTDel, 1, 2, 1, 1);

        dspnTDel = new QDoubleSpinBox(filterpars);
        dspnTDel->setObjectName(QString::fromUtf8("dspnTDel"));

        gridLayout_11->addWidget(dspnTDel, 2, 2, 1, 1);


        gridLayout_9->addLayout(gridLayout_11, 4, 2, 2, 1);

        lblTStat = new QLabel(filterpars);
        lblTStat->setObjectName(QString::fromUtf8("lblTStat"));

        gridLayout_9->addWidget(lblTStat, 4, 0, 2, 1);


        gridLayout->addLayout(gridLayout_9, 7, 0, 1, 3);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout->addItem(verticalSpacer, 0, 0, 1, 1);


        retranslateUi(filterpars);

        QMetaObject::connectSlotsByName(filterpars);
    } // setupUi

    void retranslateUi(QDialog *filterpars)
    {
        filterpars->setWindowTitle(QApplication::translate("filterpars", "Dialog", 0, QApplication::UnicodeUTF8));
        pushButton_3->setText(QApplication::translate("filterpars", "Export Configuration", 0, QApplication::UnicodeUTF8));
        cmbFilter->clear();
        cmbFilter->insertItems(0, QStringList()
         << QApplication::translate("filterpars", "Gaussian Matched Filter", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("filterpars", "Gabor Filter", 0, QApplication::UnicodeUTF8)
        );
        pushButton->setText(QApplication::translate("filterpars", "Set Configuration", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("filterpars", "Length (L)", 0, QApplication::UnicodeUTF8));
        chkLFix->setText(QApplication::translate("filterpars", "Fixed", 0, QApplication::UnicodeUTF8));
        lblLLow->setText(QApplication::translate("filterpars", "Lower", 0, QApplication::UnicodeUTF8));
        lblLUpp->setText(QApplication::translate("filterpars", "Upper", 0, QApplication::UnicodeUTF8));
        lblLDel->setText(QApplication::translate("filterpars", "Delta", 0, QApplication::UnicodeUTF8));
        lblLStat->setText(QApplication::translate("filterpars", "Fixed Value", 0, QApplication::UnicodeUTF8));
        chkTune->setText(QApplication::translate("filterpars", "Tune", 0, QApplication::UnicodeUTF8));
        lblTemplate->setText(QApplication::translate("filterpars", "Template", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("filterpars", "Filter", 0, QApplication::UnicodeUTF8));
        pushButton_2->setText(QApplication::translate("filterpars", "Import Configuration", 0, QApplication::UnicodeUTF8));
        btnSetOMet->setText(QApplication::translate("filterpars", "Set Parameters", 0, QApplication::UnicodeUTF8));
        cmbOMethod->clear();
        cmbOMethod->insertItems(0, QStringList()
         << QApplication::translate("filterpars", "Exhaustive search", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("filterpars", "Genetic Algorithm", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("filterpars", "Differential Evolution", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("filterpars", "Univariate Marginal Distribution Algorithm", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("filterpars", "Boltzmann Univariate Marginal Distribution Algorithm", 0, QApplication::UnicodeUTF8)
        );
        lblOMethod->setText(QApplication::translate("filterpars", "Optimization Method", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("filterpars", "Standard Deviation (\317\203)", 0, QApplication::UnicodeUTF8));
        chkSFix->setText(QApplication::translate("filterpars", "Fixed", 0, QApplication::UnicodeUTF8));
        lblSLow->setText(QApplication::translate("filterpars", "Lower", 0, QApplication::UnicodeUTF8));
        lblSUpp->setText(QApplication::translate("filterpars", "Upper", 0, QApplication::UnicodeUTF8));
        lblSDel->setText(QApplication::translate("filterpars", "Delta", 0, QApplication::UnicodeUTF8));
        lblSStat->setText(QApplication::translate("filterpars", "Fixed Value", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("filterpars", "Width (T)", 0, QApplication::UnicodeUTF8));
        chkTFix->setText(QApplication::translate("filterpars", "Fixed", 0, QApplication::UnicodeUTF8));
        lblTLow->setText(QApplication::translate("filterpars", "Lower", 0, QApplication::UnicodeUTF8));
        lblTUpp->setText(QApplication::translate("filterpars", "Upper", 0, QApplication::UnicodeUTF8));
        lblTDel->setText(QApplication::translate("filterpars", "Delta", 0, QApplication::UnicodeUTF8));
        lblTStat->setText(QApplication::translate("filterpars", "Fixed Value", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class filterpars: public Ui_filterpars {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_FILTERPARS_H
