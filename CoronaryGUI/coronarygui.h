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

    void showImageBase(IMGVTK::IMG_IDX, const int viewport_ID );

    void loadBase();
    void loadGroundtruth();
    void showGroundtruth();
    void showBase();

    void minimizeVP4();
    void maximizeVP4();

    void on_action_Open_file_triggered();

    void on_action_Quit_triggered();

    void on_actionGMF_ROC_triggered();

    void on_actionChaudhuri_triggered();

    void on_actionSkeletonize_triggered();

    void on_qvtkVP4_customContextMenuRequested(const QPoint &pos);

    void on_lblVP1_customContextMenuRequested(const QPoint &pos);

    private:
    QVector<QRgb> colors;
    Ui::coronaryGUI *ui;
    QMenu *lblVP1_CM;
    QMenu *qvtkVP4_CM;
    QAction *lblVP1_CM_LoadGroundtruth, *lblVP1_CM_ShowGroundtruth, *lblVP1_CM_ShowBase;
    QAction *qvtkVP4_CM_Maximize, *qvtkVP4_CM_Minimize;

    QImage *imgBase;

    RECONS3D mi_rec3D;
    bool loaded, show_gt;
};

#endif // CORONARYGUI_H
