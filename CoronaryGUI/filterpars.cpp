#include "filterpars.h"
#include "ui_filterpars.h"

filterpars::filterpars(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::filterpars)
{
    ui->setupUi(this);

    ui->lblOMethod->setMaximumHeight( 0 );
    ui->cmbOMethod->setMaximumHeight( 0 );
    ui->btnSetOMet->setMaximumHeight( 0 );

    ui->dspnLLow->setMaximumWidth( 0 );
    ui->dspnLUpp->setMaximumWidth( 0 );
    ui->dspnLDel->setMaximumWidth( 0 );
    ui->lblLLow->setMaximumWidth( 0 );
    ui->lblLUpp->setMaximumWidth( 0 );
    ui->lblLDel->setMaximumWidth( 0 );
}

//

filterpars::~filterpars()
{
    delete ui;
}

void filterpars::on_chkTune_toggled(bool checked)
{
    if( checked ){
        ui->lblOMethod->setMaximumHeight( 20 );
        ui->cmbOMethod->setMaximumHeight( 20 );
        ui->btnSetOMet->setMaximumHeight( 20 );
    }else{
        ui->lblOMethod->setMaximumHeight( 0 );
        ui->cmbOMethod->setMaximumHeight( 0 );
        ui->btnSetOMet->setMaximumHeight( 0 );
    }
}

void filterpars::on_chkLFix_toggled(bool checked)
{

    unsigned int width_fix, width_ran;

    if( checked ){
        width_fix = 16777215;
        width_ran = 0;

        ui->lblSStat->setText( "Fixed value" );

    }else{
        width_fix = 0;
        width_ran = 16777215;

        ui->lblSStat->setText( "Search space" );
    }

    ui->dspnLLow->setMaximumWidth( width_ran );
    ui->dspnLUpp->setMaximumWidth( width_ran );
    ui->dspnLDel->setMaximumWidth( width_ran );
    ui->lblLLow->setMaximumWidth( width_ran );
    ui->lblLUpp->setMaximumWidth( width_ran );
    ui->lblLDel->setMaximumWidth( width_ran );

    ui->dspnL->setMaximumWidth( width_fix );

}




void filterpars::on_cmbFilter_currentIndexChanged(int index)
{
    // 0: GMF
    // 1: Gabor Filter

    unsigned int heigth = 0;

    switch( index ){
    case 0:
        heigth = 0;
        break;
    case 1:
        height = 16777215;
        break;
    }

    ui->lblSBox->setMaximumWidth( height );
    ui->lblSDel->setMaximumWidth( height );
    ui->lblSLow->setMaximumWidth( height );
    ui->lblSUpp->setMaximumWidth( height );
    ui->lblLStat->setMaximumWidth( height );
    ui->chkSFix->setMaximumWidth( height );
    ui->dspnS->setMaximumWidth( height );
    ui->dspnSLow->setMaximumWidth( height );
    ui->dspnSUpp->setMaximumWidth( height );
    ui->dspnSDel->setMaximumWidth( height );
}
