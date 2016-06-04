#include "coronarygui.h"
#include "ui_coronarygui.h"

coronaryGUI::coronaryGUI(QWidget *parent) : QMainWindow(parent), ui(new Ui::coronaryGUI){
    ui->setupUi(this);


    //// Define characteristics for show images in QLabels
    defineColors();
    imgBase = NULL;

    //// Define status for load grond-truth:
    loaded = false;
    showing_gt = false;
    ui->btnLoadGT->setEnabled(false);



    //// Connect Reconstructor to viewports and log:
    maximized = false;
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
    QFileDialog dlgOpen(this);
    dlgOpen.setDirectory("~");
    dlgOpen.setFileMode(QFileDialog::ExistingFiles);
    dlgOpen.setNameFilter(trUtf8("Png (*.png);;Bitmap (*.bmp);;Jpeg (*.jpg,*.jpeg);;DICOM (*.dcm,*.*)"));
    dlgOpen.setWindowTitle("Open file as image base");
    if( dlgOpen.exec() ){
        QStringList filenames = dlgOpen.selectedFiles();
        const int n_paths = filenames.length();

        if( n_paths == 1 ){
            QByteArray filename_ba = filenames.at(0).toLatin1();

            if( strlen(filename_ba.data())  ){
                mi_rec3D.agregarInput( filename_ba.data() );
                ui->action_Open_file->setEnabled(false);
                ui->btnLoadGT->setEnabled(true);
                ui->btnMaximizeVP4->setEnabled(true);
                showImage( IMGVTK::BASE, 1, 0);
            }
        }else if( n_paths > 1){
            char **paths = new char* [ n_paths ];

            for( int i = 0; i < n_paths; i++){
                QByteArray filename_ba = filenames.at(i).toLatin1();
                const int len_path_i = strlen( filename_ba.data() );
                paths[i] = new char [len_path_i + 1];
                memcpy( paths[i], filename_ba.data(), len_path_i*sizeof(char));
                paths[i][len_path_i] = '\0';
            }

            mi_rec3D.agregarInput(paths, n_paths);

            ui->action_Open_file->setEnabled(false);
            ui->btnLoadGT->setEnabled(true);
            ui->btnMaximizeVP4->setEnabled(true);
            showImage( IMGVTK::BASE, 1, ui->hsldImages->value() );

            for( int i = 0; i < n_paths; i++ ){
                delete [] paths[i];
            }
            delete [] paths;

            ui->hsldImages->setMaximum( n_paths );
        }
    }
}



void coronaryGUI::defineColors(){

    for( double c = 0.0; c <= 255.0; c++){
        colors.append( qRgb( c , c , c ) );
    }
}



void coronaryGUI::showImage( IMGVTK::IMG_IDX img_idx, const int viewport_ID, const int img_i ){

    DEB_MSG(COLOR_BLUE "image id: " COLOR_WHITE << img_i << COLOR_NORMAL);

    const int cols_org = mi_rec3D.getCols(0);
    const int rows_org = mi_rec3D.getRows(0);
    int cols, rows;

    if( (cols_org % 512) == 0 ){
        cols = 512;
    }else if( (cols_org % 300) == 0 ){
        cols = 300;
    }

    if( (rows_org % 512) == 0 ){
        rows = 512;
    }else if( (rows_org % 300) == 0 ){
        rows = 300;
    }

    DEB_MSG( COLOR_BLUE "cols: " COLOR_GREEN << cols_org << COLOR_BLUE ", rows: " COLOR_GREEN << rows_org << COLOR_NORMAL);

    if( !imgBase ){
        imgBase = new QImage( cols, rows, QImage::Format_Indexed8 );
        imgBase->setColorTable( colors );
    }

    DEB_MSG(COLOR_BLUE "img cols: " COLOR_GREEN << imgBase->height() << COLOR_BLUE ", rows: " COLOR_GREEN << imgBase->height() << COLOR_NORMAL);

    const double *img_ptr = mi_rec3D.get_pixelData( 0, img_idx );

    double min = INF;
    double max =-INF;

    DEB_MSG("img_idx BASE?" << (img_idx == IMGVTK::BASE) << " img_ptr: " << img_ptr);

    for( int y = 0; y < rows; y++){
        for( int x = cols*img_i; x < cols*(img_i + 1); x++){

            if( min > *(img_ptr + x + y*cols_org) ){
                min = *(img_ptr + x + y*cols_org );
            }
            if( max < *(img_ptr + x + y*cols_org) ){
                max = *(img_ptr + x + y*cols_org );
            }
        }
    }

    const double range = 1.0 / (max - min);

    for(int y = 0; y < rows; y++){
        for( int x = 0; x < cols; x++){
            const int idx = (int) (255.0 * (*(img_ptr + y*cols_org + (x + cols*img_i)) - min) * range);
            imgBase->setPixel(x, y, idx);
        }
    }

    switch( viewport_ID ){
        case 1:
            ui->lblVP1->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
            ui->lblVP1->setAlignment(Qt::AlignCenter);
            ui->lblVP1->setMinimumSize(cols, rows);
            ui->lblVP1->setPixmap(QPixmap::fromImage(*imgBase));
            ui->lblVP1->show();
            break;
    }



}


void coronaryGUI::loadGroundtruth(){
    QFileDialog dlgOpen(this);
    dlgOpen.setDirectory("~");
    dlgOpen.setFileMode(QFileDialog::ExistingFiles);
    dlgOpen.setNameFilter(trUtf8("Png (*.png);;Bitmap (*.bmp);;Jpeg (*.jpg,*.jpeg)"));
    dlgOpen.setWindowTitle("Open file as ground-truth");

    if( dlgOpen.exec() ){
        QStringList filenames = dlgOpen.selectedFiles();
        const int n_paths = filenames.length();

        if( n_paths == 1 ){
            QByteArray filename_ba = filenames.at(0).toLatin1();

            if( strlen(filename_ba.data())  ){
                mi_rec3D.agregarGroundtruth( filename_ba.data(), 0 );
                loaded = true;
            }
        }else if( n_paths > 1 ){
            char **paths = new char* [ n_paths ];

            for( int i = 0; i < n_paths; i++){
                QByteArray filename_ba = filenames.at(i).toLatin1();
                const int len_path_i = strlen( filename_ba.data() );
                paths[i] = new char [len_path_i + 1];
                memcpy( paths[i], filename_ba.data(), len_path_i*sizeof(char));
                paths[i][len_path_i] = '\0';
            }

            mi_rec3D.agregarGroundtruth(paths, n_paths, 0);
            loaded = true;

            for( int i = 0; i < n_paths; i++ ){
                delete [] paths[i];
            }
            delete [] paths;
        }
    }
}


void coronaryGUI::minimizeVP4()
{
    ui->lblVP1->setVisible( true );
    ui->lblVP2->setVisible( true );
    ui->lblVP3->setVisible( true );
    ui->btnLoadGT->setVisible( true );
    ui->hsldImages->setVisible( true );

    maximized = false;
}

void coronaryGUI::maximizeVP4()
{
    ui->lblVP1->setVisible( false );
    ui->lblVP2->setVisible( false );
    ui->lblVP3->setVisible( false );
    ui->btnLoadGT->setVisible( false );
    ui->hsldImages->setVisible( false );

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
        if( !showing_gt ){
            showImage( IMGVTK::GROUNDTRUTH, 1, ui->hsldImages->value() );
            ui->btnLoadGT->setText( "Show Base" );
            showing_gt = true;
        }else{
            showImage( IMGVTK::BASE, 1, ui->hsldImages->value() );
            ui->btnLoadGT->setText( "Show Ground-truth" );
            showing_gt = false;
        }
    }
}

void coronaryGUI::on_btnMaximizeVP4_clicked()
{
    if( maximized ){
        minimizeVP4();
        ui->btnMaximizeVP4->setText( "Maximize 3D canvas" );
    }else{
        maximizeVP4();
        ui->btnMaximizeVP4->setText( "Minimize 3D canvas" );
    }
}

void coronaryGUI::on_hsldImages_valueChanged(int value)
{
    IMGVTK::IMG_IDX img_idx;
    if( showing_gt ){
        DEB_MSG("SLIDE: GROUNDTRUTH");
        img_idx = IMGVTK::GROUNDTRUTH;
    }else{
        DEB_MSG("SLIDE: BASE");
        img_idx = IMGVTK::BASE;
    }
    showImage( img_idx, 1, value-1 );
}
