#include "coronarygui.h"
#include "ui_coronarygui.h"


coronaryGUI::coronaryGUI(QWidget *parent) : QMainWindow(parent), ui(new Ui::coronaryGUI){
    ui->setupUi(this);
    defineColors();

    imgBase = NULL;

    //// Create context menu for vtk viewports
    lblVP1_CM = new QMenu(ui->lblVP1);
    qvtkVP4_CM = new QMenu(ui->qvtkVP4);


    //// Connect context menu actions
    connect(ui->lblVP1, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(on_lblVP1_customContextMenuRequested(QPoint)));
    connect(ui->qvtkVP4, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(on_qvtkVP4_customContextMenuRequested(QPoint)));

    qvtkVP4_CM_Maximize = qvtkVP4_CM->addAction("Maximize Viewport 3D");
    connect(qvtkVP4_CM_Maximize,SIGNAL(triggered()),this,SLOT(maximizeVP4()));

    //// Connect Reconstructor to viewports and log:
    ui->qvtkVP4->GetRenderWindow()->AddRenderer( mi_rec3D.getRenderer() );
    mi_rec3D.setLog( ui->ptxtLog );
}

coronaryGUI::~coronaryGUI(){
    delete ui;

    if( lblVP1_CM ){
        delete lblVP1_CM;
    }

    if( qvtkVP4_CM ){
        delete qvtkVP4_CM;
    }

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

        lblVP1_CM_LoadGroundtruth = lblVP1_CM->addAction("Load file as groundtruth");
        connect(lblVP1_CM_LoadGroundtruth,SIGNAL(triggered()),this,SLOT(loadGroundtruth()));
        ui->action_Open_file->setEnabled(false);

        showImageBase( IMGVTK::BASE, 1);
    }
}



void coronaryGUI::defineColors(){

    for( double c = 0.0; c <= 255.0; c++){
        colors.append( qRgb( c / 255.0 , c / 255.0, c / 255.0 ) );
    }
}



void coronaryGUI::showImageBase( IMGVTK::IMG_IDX img_idx, const int viewport_ID ){

    const int cols = mi_rec3D.getCols(0);
    const int rows = mi_rec3D.getRows(0);

    if( imgBase ){
        delete imgBase;
    }
    imgBase = new QImage( cols, rows, QImage::Format_Indexed8 );
    imgBase->setColorTable( colors );

    const double *img_ptr = mi_rec3D.get_pixelData( 0, IMGVTK::BASE );

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
        lblVP1_CM->removeAction( lblVP1_CM_LoadGroundtruth );
        lblVP1_CM_LoadGroundtruth = NULL;
        lblVP1_CM_ShowGroundtruth = lblVP1_CM->addAction("Show ground-truth");
        connect(lblVP1_CM_ShowGroundtruth,SIGNAL(triggered()),this,SLOT(showGroundtruth()));
    }
}

void coronaryGUI::showGroundtruth()
{
    lblVP1_CM->removeAction( lblVP1_CM_ShowGroundtruth );
    lblVP1_CM_ShowGroundtruth = NULL;
    lblVP1_CM_ShowBase = lblVP1_CM->addAction("Show base");
    connect(lblVP1_CM_ShowBase,SIGNAL(triggered()),this,SLOT(showBase()));

    mi_rec3D.mostrarGroundtruth( 0 );
}

void coronaryGUI::showBase()
{
    lblVP1_CM->removeAction( lblVP1_CM_ShowBase );
    lblVP1_CM_ShowBase = NULL;
    lblVP1_CM_ShowGroundtruth = lblVP1_CM->addAction("Show ground-truth");
    connect(lblVP1_CM_ShowGroundtruth,SIGNAL(triggered()),this,SLOT(showGroundtruth()));

    mi_rec3D.mostrarBase( 0 );
}

void coronaryGUI::minimizeVP4()
{
    ui->lblVP1->setVisible( true );
    ui->lblVP2->setVisible( true );
    ui->lblVP3->setVisible( true );

    qvtkVP4_CM->removeAction( qvtkVP4_CM_Minimize );
    qvtkVP4_CM_Minimize = NULL;
    qvtkVP4_CM_Maximize = qvtkVP4_CM->addAction("Maximize Viewport 3D");
    connect(qvtkVP4_CM_Maximize,SIGNAL(triggered()),this,SLOT(maximizeVP4()));
}

void coronaryGUI::maximizeVP4()
{
    ui->lblVP1->setVisible( false );
    ui->lblVP2->setVisible( false );
    ui->lblVP3->setVisible( false );

    qvtkVP4_CM->removeAction( qvtkVP4_CM_Maximize );
    qvtkVP4_CM_Maximize = NULL;
    qvtkVP4_CM_Minimize = qvtkVP4_CM->addAction("Minimize Viewport 3D");
    connect(qvtkVP4_CM_Minimize, SIGNAL(triggered()),this,SLOT(minimizeVP4()));
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

void coronaryGUI::on_qvtkVP4_customContextMenuRequested(const QPoint &pos)
{
    ui->ptxtLog->appendPlainText("Trying to pop-up the menu...");
    qvtkVP4_CM->popup(mapToGlobal(pos));
}

void coronaryGUI::on_lblVP1_customContextMenuRequested(const QPoint &pos)
{
    ui->ptxtLog->appendPlainText("Trying to pop-up form VP1");
    lblVP1_CM->popup(mapToGlobal(pos));
}
