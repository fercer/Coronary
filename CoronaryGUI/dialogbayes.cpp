#include "dialogbayes.h"
#include "ui_dialogbayes.h"

DialogBayes::DialogBayes(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogBayes)
{
    ui->setupUi(this);
}

DialogBayes::~DialogBayes()
{
    delete ui;
}
