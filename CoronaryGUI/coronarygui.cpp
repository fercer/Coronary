#include "coronarygui.h"
#include "ui_coronarygui.h"

coronaryGUI::coronaryGUI(QWidget *parent) : QMainWindow(parent), ui(new Ui::coronaryGUI){
    ui->setupUi(this);
    fp_filters = NULL;

    //// Define characteristics for show images in QLabels
    imgVP1 = NULL;
    imgVP2 = NULL;
    imgVP3 = NULL;
    imgVP4 = NULL;
    defineColors();


    //// Connect Reconstructor to viewports and log:
    mi_rec3D.setLog( ui->txtLog );
    mi_rec3D.setProgressBar( ui->pbarLog );
}

coronaryGUI::~coronaryGUI(){
    delete ui;
    if( imgVP1 ){
        delete imgVP1;
    }
    if( imgVP2 ){
        delete imgVP2;
    }
    if( imgVP3 ){
        delete imgVP3;
    }
    if( imgVP4 ){
        delete imgVP4;
    }
    if( fp_filters ){
        fclose( fp_filters );
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
                mi_rec3D.agregarInput( filename_ba.data(), true );
                ui->action_Open_file->setEnabled(false);
                showImage( IMGVTK::BASE, 1);
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

            mi_rec3D.agregarInput(paths, n_paths, true);

            ui->action_Open_file->setEnabled(false);

            for( int i = 0; i < n_paths; i++ ){
                delete [] paths[i];
            }
            delete [] paths;

        }
    }
}



void coronaryGUI::defineColors(){

    for( double c = 0.0; c <= 255.0; c++){
        colors.append( qRgb( c , c , c ) );
    }
}



void coronaryGUI::showImage( IMGVTK::IMG_IDX img_idx, const int viewport_ID ){
    const int cols_org = mi_rec3D.getCols(0);
    const int rows_org = mi_rec3D.getRows(0);

    QImage *imgBase = NULL;

    switch (viewport_ID){
    case 1:
        if( !imgVP1 ){
            imgVP1 = new QImage( cols_org, rows_org, QImage::Format_Indexed8 );
            imgVP1->setColorTable( colors );
        }
        imgBase = imgVP1;
        break;

    case 2:
        if( !imgVP2 ){
            imgVP2 = new QImage( cols_org, rows_org, QImage::Format_Indexed8 );
            imgVP2->setColorTable( colors );
        }
        imgBase = imgVP2;
        break;

    case 3:
        if( !imgVP3 ){
            imgVP3 = new QImage( cols_org, rows_org, QImage::Format_Indexed8 );
            imgVP3->setColorTable( colors );
        }
        imgBase = imgVP3;
        break;

    case 4:
        if( !imgVP4 ){
            imgVP4 = new QImage( cols_org, rows_org, QImage::Format_Indexed8 );
            imgVP4->setColorTable( colors );
        }
        imgBase = imgVP4;
        break;
    }


    const double *img_ptr = mi_rec3D.get_pixelData( 0, img_idx );

    int offset_x = 0, offset_y = 0;

    switch (img_idx){
    case IMGVTK::SKELETON:
        offset_x = 1;
        offset_y = 1;
        break;
    }

    double min = INF;
    double max =-INF;

    for( int y = 0; y < (2*offset_y + rows_org); y++){
        for( int x = 0; x < (2*offset_x + cols_org); x++){

            if( min > *(img_ptr + (x+offset_x) + (y+offset_y)*(2*offset_x + cols_org)) ){
                min = *(img_ptr + (x+offset_x) + (y+offset_y)*(2*offset_x + cols_org) );
            }
            if( max < *(img_ptr + (x+offset_x) + (y+offset_y)*(2*offset_x + cols_org)) ){
                max = *(img_ptr + (x+offset_x) + (y+offset_y)*(2*offset_x + cols_org) );
            }
        }
    }


    char mensaje[] = COLOR_BACK_BLACK COLOR_WHITE "min: X.XXX , max: X.XXX " COLOR_BACK_WHITE COLOR_BLACK;
    sprintf(mensaje, COLOR_BACK_BLACK COLOR_WHITE "min: %1.3f , max: %1.3f" COLOR_BACK_WHITE COLOR_BLACK, min, max);
    ui->txtLog->append( mensaje );

    const double range = 1.0 / (max - min);

    DEB_MSG("rows: " << rows_org << ", cols: " << cols_org);
    DEB_MSG("IMGBASE: " << imgBase);

    for(int y = 0; y < rows_org; y++){
        for( int x = 0; x < cols_org; x++){
            const int idx = (int) (255.0 * (*(img_ptr + (x+offset_x) + ( offset_y + rows_org-1 - y)*(2*offset_x + cols_org)) - min) * range);
            imgBase->setPixel(x, y, idx);
        }
    }

    QLabel *lbltemp;

    switch( viewport_ID ){
    case 1:
        lbltemp = ui->lblVP1;
        break;
    case 2:
        lbltemp = ui->lblVP2;
        break;
    case 3:
        lbltemp = ui->lblVP3;
        break;
    case 4:
        lbltemp = ui->lblVP4;
        break;
    }

    lbltemp->setAlignment(Qt::AlignCenter);
    lbltemp->setPixmap(QPixmap::fromImage(*imgBase));
    lbltemp->show();
}


void coronaryGUI::loadGroundtruth(){
    QFileDialog dlgOpen(this);
    dlgOpen.setDirectory("~/test_data/");
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
                showImage(IMGVTK::GROUNDTRUTH, 3);
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

            for( int i = 0; i < n_paths; i++ ){
                delete [] paths[i];
            }
            delete [] paths;
        }
    }
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

void coronaryGUI::on_actionOpen_file_as_ground_truth_triggered()
{
    loadGroundtruth();
}
