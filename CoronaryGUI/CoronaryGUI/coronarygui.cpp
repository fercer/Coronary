#include "coronarygui.h"
#include "ui_coronarygui.h"

coronaryGUI::coronaryGUI(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::coronaryGUI)
{
    ui->setupUi(this);
}

coronaryGUI::~coronaryGUI()
{
    delete ui;
}
