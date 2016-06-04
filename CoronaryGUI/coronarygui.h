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

#include <vector>

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

    void showImage(IMGVTK::IMG_IDX, const int viewport_ID );

    void loadBase();
    void loadGroundtruth();

    void minimizeVP4();
    void maximizeVP4();

    void on_action_Open_file_triggered();

    void on_action_Quit_triggered();

    void on_actionSkeletonize_triggered();

    void on_btnLoadGT_clicked();

    void on_btnMaximizeVP4_clicked();

private:
    QVector<QRgb> colors;
    Ui::coronaryGUI *ui;

    QImage *imgBase;

    RECONS3D mi_rec3D;
    bool loaded, show_gt;
    bool maximized;
};

#endif // CORONARYGUI_H
