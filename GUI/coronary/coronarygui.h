#ifndef CORONARYGUI_H
#define CORONARYGUI_H

#include <QMainWindow>

namespace Ui {
class coronaryGUI;
}

class coronaryGUI : public QMainWindow
{
    Q_OBJECT

public:
    explicit coronaryGUI(QWidget *parent = 0);
    ~coronaryGUI();

private:
    Ui::coronaryGUI *ui;
};

#endif // CORONARYGUI_H
