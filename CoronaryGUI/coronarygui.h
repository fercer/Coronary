#ifndef CORONARYGUI_H
#define CORONARYGUI_H

#include <QMainWindow>
#include <QMenu>
#include <QAction>
#include <QPlainTextEdit>
#include <QFileDialog>
#include <QString>

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
    void on_qvtkVP1_customContextMenuRequested(const QPoint &pos);
    void loadGroundtruth();
    void showGroundtruth();
    void showBase();

    void on_action_Open_file_triggered();

    void on_action_Quit_triggered();

    void on_actionGMF_ROC_triggered();

    void on_actionChaudhuri_triggered();

    void on_actionSkeletonize_triggered();

private:
    Ui::coronaryGUI *ui;
    QMenu* qvtkVP1_CM;
    QAction *qvtkVP1_CM_LoadGroundtruth, *qvtkVP1_CM_ShowGroundtruth, *qvtkVP1_CM_ShowBase;

    RECONS3D mi_rec3D;
    bool loaded, show_gt;
};

#endif // CORONARYGUI_H
