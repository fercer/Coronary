#include "filterpars.h"
#include "ui_filterpars.h"

filterpars::filterpars(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::filterpars)
{
    ui->setupUi(this);
}

filterpars::~filterpars()
{
    delete ui;
}
