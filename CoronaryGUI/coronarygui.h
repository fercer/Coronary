#ifndef CORONARYGUI_H
#define CORONARYGUI_H

#include <QMainWindow>
#include <QMenu>
#include <QAction>
#include <QPlainTextEdit>
#include <QFileDialog>
#include <QString>
#include <QImage>
#include <QPixmap>
#include <QVector>
#include <QRgb>
#include <QSlider>
#include <QCheckBox>
#include <QComboBox>
#include <QLabel>

#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../src/reconstructor_3D.h"
#include "filterpars.h"
#include "dialogbayes.h"

namespace Ui {
class coronaryGUI;
}

class coronaryGUI : public QMainWindow
{
    Q_OBJECT

public:
    explicit coronaryGUI(QWidget *parent = 0);
    ~coronaryGUI();

private slots:
    void defineColors();

    void showImage(IMGVTK::IMG_IDX, const int viewport_ID);

    void loadBase();

    void loadGroundtruth();

    void on_action_Open_file_triggered();

    void on_action_Quit_triggered();

    void on_actionSkeletonize_triggered();
    void on_actionOpen_file_as_ground_truth_triggered();

    void on_actionSet_filter_parameters_triggered();

private:
    QVector<QRgb> colors;
    Ui::coronaryGUI *ui;

    QImage *imgVP1, *imgVP2, *imgVP3, *imgVP4;

    RECONS3D mi_rec3D;

    filterpars mis_filterpars;

    FILE *fp_filters;
};

#endif // CORONARYGUI_H
