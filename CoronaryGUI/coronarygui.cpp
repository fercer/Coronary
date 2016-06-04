#include "coronarygui.h"
#include "ui_coronarygui.h"


coronaryGUI::coronaryGUI(QWidget *parent) : QMainWindow(parent), ui(new Ui::coronaryGUI){
    ui->setupUi(this);


    //// Define characteristics for show images in QLabels
    defineColors();
    imgBase = NULL;

    //// Define status for load grond-truth:
    loaded = false;
    show_gt = true;
    ui->btnLoadGT->setEnabled(false);



    //// Connect Reconstructor to viewports and log:
    maximized = false;
    ui->btnMaximizeVP4->setEnabled(false);
    ui->qvtkVP4->GetRenderWindow()->AddRenderer( mi_rec3D.getRenderer() );
    mi_rec3D.setLog( ui->ptxtLog );
}

coronaryGUI::~coronaryGUI(){
    delete ui;
    if( imgBase ){
        delete imgBase;
    }
}




void coronaryGUI::loadBase(){
    QString filename = QFileDialog::getOpenFileName(
                this,
                tr("Open file as image base"),
                "~",
                "Png (*.png);;Bitmap (*.bmp);;Jpeg (*.jpg,*.jpeg);;DICOM (*.dcm,*.*)"
                );
    QByteArray filename_ba = filename.toLatin1();

    DEB_MSG( "Abriendo archivo: '" << filename_ba.data() << "'");

    if( strlen(filename_ba.data())  ){
        mi_rec3D.agregarInput( filename_ba.data() );
        ui->action_Open_file->setEnabled(false);
        ui->btnLoadGT->setEnabled(true);
        ui->btnMaximizeVP4->setEnabled(true);
        showImage( IMGVTK::BASE, 1);
    }
}



void coronaryGUI::defineColors(){

    for( double c = 0.0; c <= 255.0; c++){
        colors.append( qRgb( c , c , c ) );
    }
}



void coronaryGUI::showImage( IMGVTK::IMG_IDX img_idx, const int viewport_ID ){

    const int cols = mi_rec3D.getCols(0);
    const int rows = mi_rec3D.getRows(0);

    if( imgBase ){
        delete imgBase;
    }
    imgBase = new QImage( cols, rows, QImage::Format_Indexed8 );
    imgBase->setColorTable( colors );

    const double *img_ptr = mi_rec3D.get_pixelData( 0, img_idx );

    double min = INF;
    double max =-INF;

    for( int xy = 0; xy < rows*cols; xy++){
        if( min > *(img_ptr + xy) ){
            min = *(img_ptr + xy );
        }
        if( max < *(img_ptr + xy) ){
            max = *(img_ptr + xy );
        }
    }

    const double range = 1.0 / (max - min);

    for(int y = 0; y < cols; y++){
        for( int x = 0; x < rows; x++){
            const int idx = (int) (255.0 * (*(img_ptr + y*cols + x) - min) * range);
            imgBase->setPixel(x, y, idx);
        }
    }

    switch( viewport_ID ){
        case 1:
            ui->lblVP1->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
            ui->lblVP1->setAlignment(Qt::AlignCenter);
            ui->lblVP1->setMinimumSize(cols, rows);
            ui->lblVP1->setPixmap(QPixmap::fromImage(*imgBase));
            //ui->lblVP1->show();
            break;
    }

}


void coronaryGUI::loadGroundtruth(){
    QString filename = QFileDialog::getOpenFileName(
                this,
                tr("Open file as image base"),
                "~",
                "Png (*.png);;Bitmap (*.bmp);;Jpeg (*.jpg,*.jpeg);;DICOM (*.dcm,*.*)"
                );

    QByteArray filename_ba = filename.toLatin1();

    if( strlen(filename_ba.data())  ){
        mi_rec3D.agregarGroundtruth( filename_ba.data(), 0 );
        loaded = true;
    }
}


void coronaryGUI::minimizeVP4()
{
    ui->lblVP1->setVisible( true );
    ui->lblVP2->setVisible( true );
    ui->lblVP3->setVisible( true );

    maximized = false;
}

void coronaryGUI::maximizeVP4()
{
    ui->lblVP1->setVisible( false );
    ui->lblVP2->setVisible( false );
    ui->lblVP3->setVisible( false );

    maximized = true;
}


void coronaryGUI::on_action_Open_file_triggered()
{
    loadBase();
}

void coronaryGUI::on_action_Quit_triggered()
{
    close();
}

void coronaryGUI::on_actionGMF_ROC_triggered()
{

}

void coronaryGUI::on_actionChaudhuri_triggered()
{

}

void coronaryGUI::on_actionSkeletonize_triggered()
{
    mi_rec3D.skeletonize( 0 );
}


void coronaryGUI::on_btnLoadGT_clicked()
{
    if( !loaded ){
        loadGroundtruth();
        ui->btnLoadGT->setText( "Show Ground-truth" );
    }else{
        if( show_gt ){
            showImage( 0, IMGVTK::BASE );
            ui->btnLoadGT->setText( "Show Base" );
        }else{
            showImage( 0, IMGVTK::BASE );
            ui->btnLoadGT->setText( "Show Ground-truth" );
        }
    }
}

void coronaryGUI::on_btnMaximizeVP4_clicked()
{
    if( maximized ){
        minimizeVP4();
    }else{
        maximizeVP4();
    }
}
