#ifndef CORONARYGUI_H
#define CORONARYGUI_H

#include <QMainWindow>

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
    void on_actionQuit_triggered();

    void on_actionBase_triggered();

    void on_actionGround_truth_triggered();

    void on_action_Skeletonize_triggered();

private:
    Ui::coronaryGUI *ui;
    RECONS3D mi_reconstructor;
};

#endif // CORONARYGUI_H
