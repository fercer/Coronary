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

#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../src/reconstructor_3D.h"

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

    void showImage(IMGVTK::IMG_IDX, const int viewport_ID , const int img_i);

    void writeConfiguration(const char *path);
    void readConfiguration(const char *path);
    void loadBase();
    void loadGroundtruth();

    void minimizeVP4();
    void maximizeVP4();

    void on_action_Open_file_triggered();

    void on_action_Quit_triggered();

    void on_actionSkeletonize_triggered();

    void on_btnLoadGT_clicked();

    void on_btnMaximizeVP4_clicked();

    void on_hsldImages_valueChanged(int value);

    void on_actionSet_filter_parameters_triggered();

    void on_actionImport_parameters_triggered();

    void on_chkFixT_toggled(bool checked);

    void on_actionExport_parameters_triggered();

    void on_chkFixK_toggled(bool checked);

    void on_chkFixL_toggled(bool checked);

    void on_chkFixSigma_toggled(bool checked);

    void on_cmbAlgorithm_currentIndexChanged(int index);

    void on_cmbFilter_currentIndexChanged(int index);

    void on_btnRunConfiguration_clicked();

private:
    QVector<QRgb> colors;
    Ui::coronaryGUI *ui;

    QImage *imgBase;

    RECONS3D mi_rec3D;
    bool loaded, showing_gt;
    bool maximized;
};

#endif // CORONARYGUI_H
