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


    /// Set visibility for fixed values of parameters:
    ui->chkFixT->toggle();
    ui->chkFixK->toggle();
    ui->chkFixL->toggle();
    ui->chkFixSigma->toggle();

    ui->cmbAlgorithm->setCurrentIndex( 0 );
    ui->lblPop->setEnabled( false );
    ui->lblGen->setEnabled( false );
    ui->lblFitness->setEnabled( false );

    ui->lblCR->setEnabled( false );
    ui->lblMR->setEnabled( false );

    ui->spbPop->setEnabled( false );
    ui->spbGen->setEnabled( false );
    ui->cmbFitness->setEnabled( false );

    ui->dspbCR->setEnabled( false );
    ui->dspbMR->setEnabled( false );

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



void coronaryGUI::writeConfiguration(const char *path)
{
    int len_path = strlen(path);
    char *my_path = NULL;
    if( strcmp(path + len_path - 4, ".dat") == 0 ){
        my_path = new char [len_path + 1];
        memcpy( my_path, path, len_path*sizeof(char) );
        *(my_path + len_path ) = '\0';
    }else{
        my_path = new char [len_path + 5];
        sprintf( my_path, "%s.dat", path );
    }


    FILE *fp_config = fopen(my_path, "w");

    if( !fp_config ){
        char message[513] = "\n<< Error: Cannot write file:\'%s\'\n\n";
        sprintf(message, "\n<< Error: Cannot write file:\'%s\'\n\n", my_path);
        ui->ptxtLog->appendPlainText( message );
        delete [] my_path;
        return;
    }
    DEB_MSG( COLOR_BACK_CYAN COLOR_BLACK "\nFile '" << my_path << "'[" << fp_config << "] opened for writing successfuly\n" COLOR_NORMAL);

    fprintf( fp_config, "#@ Filter\n");

    if( ui->cmbFilter->currentIndex() == 0 ){
        fprintf(fp_config, "GMF\t1\nSSG\t0\n");
    }else{
        fprintf(fp_config, "GMF\t0\nSSG\t1\n");
    }

    fprintf(fp_config, "#@ Parameters\n");
    if( ui->chkFixT->isChecked() ){
        fprintf(fp_config, "T\t%i\n", ui->spbFixedT->value());
    }else{
        fprintf(fp_config, "T\to\n");
    }

    if( ui->chkFixL->isChecked() ){
        fprintf(fp_config, "L\t%f\n", ui->dspbFixedL->value());
    }else{
        fprintf(fp_config, "L\to\n");
    }

    if( ui->chkFixK->isChecked() ){
        fprintf(fp_config, "K\t%i\n", ui->spbFixedK->value());
    }else{
        fprintf(fp_config, "K\to\n");
    }

    if( ui->chkFixSigma->isChecked() ){
        fprintf(fp_config, "sig\t%f\n", ui->dspbFixedSigma->value());
    }else{
        fprintf(fp_config, "sig\to\n");
    }


    fprintf(fp_config, "#@ Limits\n");
    if( ui->chkFixT->isChecked() ){
        fprintf(fp_config, "Tinf\t0\nTsup\t0\nTdel\t0.0\n");
    }else{
        fprintf(fp_config, "Tinf\t%i\nTsup\t%i\nTdel\t%f\n", ui->spbInfT->value(), ui->spbSupT->value(), ui->dspbDeltaT->value());
    }

    if( ui->chkFixL->isChecked() ){
        fprintf(fp_config, "Linf\t0.0\nLsup\t0.0\nLdel\t0.0\n");
    }else{
        fprintf(fp_config, "Linf\t%f\nLsup\t%f\nLdel\t%f\n", ui->dspbInfL->value(), ui->dspbSupL->value(), ui->dspbDeltaL->value());
    }

    if( ui->chkFixK->isChecked() ){
        fprintf(fp_config, "Kinf\t0\nKsup\t0\nKdel\t0.0\n");
    }else{
        fprintf(fp_config, "Kinf\t%i\nKsup\t%i\nKdel\t%f\n", ui->spbInfK->value(), ui->spbSupK->value(), ui->dspbDeltaK->value());
    }

    if( ui->chkFixSigma->isChecked() ){
        fprintf(fp_config, "Sinf\t0.0\nSsup\t0.0\nSdel\t0.0\n");
    }else{
        fprintf(fp_config, "Sinf\t%f\nSsup\t%f\nSdel\t%f\n", ui->dspbInfSigma->value(), ui->dspbSupSigma->value(), ui->dspbDeltaSigma->value());
    }


    fprintf(fp_config, "#@ Fitness\n");
    if( ui->cmbFitness->currentIndex() == 0 ){
        fprintf(fp_config, "ROC\t1\nCC\t0\n");
    }else{
        fprintf(fp_config, "ROC\t0\nCC\t1\n");
    }


    fprintf(fp_config, "#@ Algorithm\n");
    int algorithm[] = {0,0,0,0,0};
    algorithm[ ui->cmbAlgorithm->currentIndex() ] = 1;

    fprintf(fp_config, "Unset\t%i\nExh\t%i\nGA\t%i\nUMDA\t%i\nBUMDA\t%i\n", algorithm[0], algorithm[1], algorithm[2], algorithm[3], algorithm[4] );

    fprintf(fp_config, "#@ Algorithm_Parameters\n");
    if( algorithm[0] ){
        fprintf( fp_config, "pop_size\t0\nmax_gen\t0\nCR\t0.0\nMR\t0.0");
    }else{
        fprintf( fp_config, "pop_size\t%i\nmax_gen\t%i\nCR\t%f\nMR\t%f", ui->spbPop->value(), ui->spbGen->value(), ui->dspbCR->value(), ui->dspbMR->value());
    }

    DEB_MSG( COLOR_BACK_CYAN COLOR_BLACK "\nConfiguration exported successfuly" COLOR_NORMAL);

    fclose( fp_config );
    delete [] my_path;
}




void coronaryGUI::readConfiguration( const char *path )
{
    FILE *fp_config = fopen(path, "r");

    if( !fp_config ){
        char message[513] = "\n<< Error: Cannot open file:\'%s\'\n\n";
        sprintf(message, "\n<< Error: Cannot open file:\'%s\'\n\n", path);
        ui->ptxtLog->appendPlainText( message );
        return;
    }

    char tmp = fgetc( fp_config );
    char tmp_par[128] = "", tmp_val[128] = "";
    char *tmp_str_ptr = tmp_par;

    DEB_MSG( COLOR_BACK_CYAN COLOR_BLACK "\nFile '" << path << "'[" << fp_config << "] opened successfuly\n" COLOR_NORMAL);

    do{
        tmp = fgetc( fp_config );
        if( tmp == '\n' || tmp == EOF ){
            *(tmp_str_ptr) = '\0';
            if( strcmp( tmp_par, "GMF" ) == 0 ){
                if( atoi(tmp_val) ){
                    ui->cmbFilter->setCurrentIndex( 0 );
                }
            }else if( strcmp( tmp_par, "SSG") == 0 ){
                if( atoi(tmp_val) ){
                    ui->cmbFilter->setCurrentIndex( 1 );
                }
            }else if( strcmp( tmp_par, "T") == 0 ){
                if( tmp_val[0] != 'o' ){
                    const int par_t = atoi(tmp_val);
                    ui->spbFixedT->setValue( par_t );
                    ui->chkFixT->setChecked( true );
                }else{
                    ui->chkFixT->setChecked( false );
                }
            }else if( strcmp( tmp_par, "L") == 0 ){
                if( tmp_val[0] != 'o' ){
                    const double par_l = atof(tmp_val);
                    ui->dspbFixedL->setValue( par_l );
                    ui->chkFixL->setChecked( true );
                }else{
                    ui->chkFixL->setChecked( false );
                }
            }else if( strcmp( tmp_par, "K") == 0 ){
                if( tmp_val[0] != 'o' ){
                    const int par_k = atoi(tmp_val);
                    ui->spbFixedK->setValue( par_k );
                    ui->chkFixK->setChecked( true );
                }else{
                    ui->chkFixK->setChecked( false );
                }
            }else if( strcmp( tmp_par, "sig") == 0 ){
                if( tmp_val[0] != 'o' ){
                    const double par_s = atof(tmp_val);
                    ui->dspbFixedSigma->setValue( par_s );
                    ui->chkFixSigma->setChecked( true );
                }else{
                    ui->chkFixSigma->setChecked( false );
                }
            }else if( strcmp( tmp_par, "Tinf") == 0 ){
                const int inf_t = atoi(tmp_val);
                ui->spbInfT->setValue( inf_t );
            }else if( strcmp( tmp_par, "Tsup") == 0 ){
                const int sup_t = atoi(tmp_val);
                ui->spbSupT->setValue( sup_t );
            }else if( strcmp( tmp_par, "Tdel") == 0 ){
                const double del_t = atof(tmp_val);
                ui->dspbDeltaT->setValue( del_t );
            }else if( strcmp( tmp_par, "Linf") == 0 ){
                const double inf_l = atof(tmp_val);
                ui->dspbInfL->setValue( inf_l );
            }else if( strcmp( tmp_par, "Lsup") == 0 ){
                const double sup_l = atof(tmp_val);
                ui->dspbSupL->setValue( sup_l );
            }else if( strcmp( tmp_par, "Ldel") == 0 ){
                const double del_l = atof(tmp_val);
                ui->dspbDeltaL->setValue( del_l );
            }else if( strcmp( tmp_par, "Kinf") == 0 ){
                const int inf_k = atoi(tmp_val);
                ui->spbInfK->setValue( inf_k );
            }else if( strcmp( tmp_par, "Ksup") == 0 ){
                const int sup_k = atoi(tmp_val);
                ui->spbSupK->setValue( sup_k );
            }else if( strcmp( tmp_par, "Kdel") == 0 ){
                const double del_k = atof(tmp_val);
                ui->dspbDeltaK->setValue( del_k );
            }else if( strcmp( tmp_par, "Sinf") == 0 ){
                const double inf_s = atof(tmp_val);
                ui->dspbInfSigma->setValue(  inf_s );
            }else if( strcmp( tmp_par, "Ssup") == 0 ){
                const double sup_s = atof(tmp_val);
                ui->dspbSupSigma->setValue( sup_s );
            }else if( strcmp( tmp_par, "Sdel") == 0 ){
                const double del_s = atof(tmp_val);
                ui->dspbDeltaSigma->setValue( del_s );
            }else if( strcmp( tmp_par, "ROC") == 0 ){
                if( atoi(tmp_val) ){
                    ui->cmbFitness->setCurrentIndex( 0 );
                }
            }else if( strcmp( tmp_par, "CC") == 0 ){
                if( atoi(tmp_val) ){
                    ui->cmbFitness->setCurrentIndex( 1 );
                }
            }else if( strcmp( tmp_par, "Unset") == 0 ){
                if( atoi(tmp_val) ){
                    ui->cmbAlgorithm->setCurrentIndex( 0 );
                }
            }else if( strcmp( tmp_par, "Exh") == 0 ){
                if( atoi(tmp_val) ){
                    ui->cmbAlgorithm->setCurrentIndex( 1 );
                }
            }else if( strcmp( tmp_par, "GA") == 0 ){
                if( atoi(tmp_val) ){
                    ui->cmbAlgorithm->setCurrentIndex( 2 );
                }
            }else if( strcmp( tmp_par, "UMDA") == 0 ){
                if( atoi(tmp_val) ){
                    ui->cmbAlgorithm->setCurrentIndex( 3 );
                }
            }else if( strcmp( tmp_par, "BUMDA") == 0 ){
                if( atoi(tmp_val) ){
                    ui->cmbAlgorithm->setCurrentIndex( 4 );
                }
            }else if( strcmp( tmp_par, "pop_size") == 0 ){
                const int pop_size = atoi(tmp_val);
                ui->spbPop->setValue( pop_size );
            }else if( strcmp( tmp_par, "max_gen") == 0 ){
                const int max_gen = atoi(tmp_val);
                ui->spbGen->setValue( max_gen );
            }else if( strcmp( tmp_par, "CR") == 0 ){
                const double cr = atof(tmp_val);
                ui->dspbCR->setValue( cr );
            }else if( strcmp( tmp_par, "MR") == 0 ){
                const double mr = atof(tmp_val);
                ui->dspbMR->setValue( mr );
            }

            tmp_str_ptr = tmp_par;

        }else if( tmp == '\t'){
            *(tmp_str_ptr) = '\0';
            tmp_str_ptr = tmp_val;
        }else{
            *(tmp_str_ptr) = tmp;
            tmp_str_ptr++;
        }
    }while ( tmp != EOF);

    DEB_MSG( COLOR_BACK_CYAN COLOR_BLACK "\nConfiguration imported successfuly" COLOR_NORMAL);

    fclose( fp_config );
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




void coronaryGUI::on_actionSet_filter_parameters_triggered()
{
    ui->tabParameters->setCurrentIndex(1);
}

void coronaryGUI::on_actionImport_parameters_triggered()
{
    QFileDialog dlgImport(this);
    dlgImport.setDirectory("~");
    dlgImport.setFileMode(QFileDialog::ExistingFile);
    dlgImport.setNameFilter(trUtf8("Configuration File (*.dat)"));
    dlgImport.setWindowTitle("Import configuration for Filtering");

    if( dlgImport.exec() ){
        QStringList filename = dlgImport.selectedFiles();
        QByteArray filename_ba = filename.at(0).toLatin1();

        ui->ptxtLog->appendPlainText( filename_ba );

        if( strlen(filename_ba.data())  ){
            readConfiguration( filename_ba.data() );
        }
    }
}

void coronaryGUI::on_actionExport_parameters_triggered()
{
     QFileDialog dlgExport(this);
     dlgExport.setDirectory("~");

     dlgExport.setFileMode(QFileDialog::AnyFile);
     dlgExport.setAcceptMode(QFileDialog::AcceptSave);
     dlgExport.setNameFilter(trUtf8("Configuration File (*.dat)"));
     dlgExport.setWindowTitle("Export configuration for Filtering");

     if( dlgExport.exec() ){
         QStringList filename = dlgExport.selectedFiles();
         QByteArray filename_ba = filename.at(0).toLatin1();

         ui->ptxtLog->appendPlainText( filename_ba );
         if( strlen(filename_ba.data())  ){
             writeConfiguration( filename_ba.data() );
         }
     }
}


void coronaryGUI::on_chkFixT_toggled(bool checked)
{
    ui->spbFixedT->setVisible( checked );

    ui->lblLeftT->setVisible( !checked );
    ui->lblComaT->setVisible( !checked );
    ui->lblRightT->setVisible( !checked );
    ui->spbInfT->setVisible( !checked );
    ui->spbSupT->setVisible( !checked );


    ui->dspbDeltaT->setEnabled( !checked );
    ui->lblParT->setEnabled( !checked );
}



void coronaryGUI::on_chkFixK_toggled(bool checked)
{
    ui->spbFixedK->setVisible( checked );

    ui->lblLeftK->setVisible( !checked );
    ui->lblComaK->setVisible( !checked );
    ui->lblRightK->setVisible( !checked );
    ui->spbInfK->setVisible( !checked );
    ui->spbSupK->setVisible( !checked );

    ui->dspbDeltaK->setEnabled( !checked );
    ui->lblParK->setEnabled( !checked );
}

void coronaryGUI::on_chkFixL_toggled(bool checked)
{
    ui->dspbFixedL->setVisible( checked );

    ui->lblLeftL->setVisible( !checked );
    ui->lblComaL->setVisible( !checked );
    ui->lblRightL->setVisible( !checked );
    ui->dspbInfL->setVisible( !checked );
    ui->dspbSupL->setVisible( !checked );

    ui->dspbDeltaL->setEnabled( !checked );
    ui->lblParL->setEnabled( !checked );
}

void coronaryGUI::on_chkFixSigma_toggled(bool checked)
{
    ui->dspbFixedSigma->setVisible( checked );

    ui->lblLeftSigma->setVisible( !checked );
    ui->lblComaSigma->setVisible( !checked );
    ui->lblRightSigma->setVisible( !checked );
    ui->dspbInfSigma->setVisible( !checked );
    ui->dspbSupSigma->setVisible( !checked );

    ui->dspbDeltaSigma->setEnabled( !checked );
    ui->lblParSigma->setEnabled( !checked );
}

void coronaryGUI::on_cmbAlgorithm_currentIndexChanged(int index)
{
    ui->lblPop->setEnabled( (index > 1) );
    ui->lblGen->setEnabled( (index > 1) );
    ui->lblFitness->setEnabled( (index != 0) );

    ui->lblCR->setEnabled( (index == 2 || index == 3) );
    ui->lblMR->setEnabled( (index == 2) );

    ui->spbPop->setEnabled( (index > 1) );
    ui->spbGen->setEnabled( (index > 1) );
    ui->cmbFitness->setEnabled( (index != 0) );
    ui->dspbCR->setEnabled( (index == 2 || index == 3) );
    ui->dspbMR->setEnabled( (index == 2) );

    char label[] = "delta";
    if( index == 2 || index == 3 ){
        sprintf( label, "bits" );
    }

    ui->lblParT->setText( label );
    ui->lblParL->setText( label );
    ui->lblParK->setText( label );
    ui->lblParSigma->setText( label );

}

void coronaryGUI::on_cmbFilter_currentIndexChanged(int index)
{
    ui->chkFixSigma->setVisible( (index==0) );

    const bool fixed = ui->chkFixSigma->isChecked();

    ui->lblSigmaFixed->setVisible( (index == 0) );
    ui->lblLeftSigma->setVisible( (index==0) && !fixed );
    ui->lblRightSigma->setVisible( (index==0) && !fixed );
    ui->lblComaSigma->setVisible( (index==0) && !fixed );
    ui->lblParSigma->setVisible( (index==0) );
    ui->dspbDeltaSigma->setVisible( (index==0) );
    ui->dspbFixedSigma->setVisible( (index==0) && fixed  );
    ui->dspbInfSigma->setVisible( (index==0) && !fixed  );
    ui->dspbSupSigma->setVisible( (index==0) && !fixed  );
}

void coronaryGUI::on_btnRunConfiguration_clicked()
{
    switch( ui->cmbFilter->currentIndex() ){
    case 0:
        mi_rec3D.setFiltroMetodo( FILTROS::GMF );
        break;
    case 1:
        mi_rec3D.setFiltroMetodo( FILTROS::SS_GABOR );
        break;
    }

    switch( ui->cmbFitness->currentIndex() ){
    case 0:
        mi_rec3D.setFiltroEval( FILTROS::ROC );
        break;
    case 1:
        mi_rec3D.setFiltroEval( FILTROS::CORCON );
        break;
    }

    switch( ui->cmbAlgorithm->currentIndex() ) {
    case 0:
        mi_rec3D.setFiltroEntrenamiento( FILTROS::EVO_UNSET );
        break;
    case 1:
        mi_rec3D.setFiltroEntrenamiento( FILTROS::EXHAUSTIVA );
        break;
    case 2:
        mi_rec3D.setFiltroEntrenamiento( FILTROS::EA_GA );
        break;
    case 3:
        mi_rec3D.setFiltroEntrenamiento( FILTROS::EDA_UMDA );
        break;
    case 4:
        mi_rec3D.setFiltroEntrenamiento( FILTROS::EDA_BUMDA );
        break;
    }

    mi_rec3D.setFiltroEntrenamientoPars(FILTROS::POPSIZE, (double)ui->spbPop->value() );
    mi_rec3D.setFiltroEntrenamientoPars(FILTROS::MAXGEN, (double)ui->spbGen->value() );
    mi_rec3D.setFiltroEntrenamientoPars(FILTROS::CR, ui->dspbCR->value() );
    mi_rec3D.setFiltroEntrenamientoPars(FILTROS::MR, ui->dspbMR->value() );

    if( ui->chkFixT->isChecked() ){
        mi_rec3D.setFiltroParametros(FILTROS::PAR_T, (double)ui->spbFixedT->value() );
    }else{
        mi_rec3D.setFiltroParametros(FILTROS::PAR_T, FILTROS::INFERIOR, (double)ui->spbInfT->value() );
        mi_rec3D.setFiltroParametros(FILTROS::PAR_T, FILTROS::SUPERIOR, (double)ui->spbSupT->value() );
        mi_rec3D.setFiltroParametros(FILTROS::PAR_T, FILTROS::DELTA, ui->dspbDeltaT->value() );
    }

    if( ui->chkFixL->isChecked() ){
        mi_rec3D.setFiltroParametros(FILTROS::PAR_L, ui->dspbFixedL->value() );
    }else{
        mi_rec3D.setFiltroParametros(FILTROS::PAR_L, FILTROS::INFERIOR, ui->dspbInfL->value() );
        mi_rec3D.setFiltroParametros(FILTROS::PAR_L, FILTROS::SUPERIOR, ui->dspbSupL->value() );
        mi_rec3D.setFiltroParametros(FILTROS::PAR_L, FILTROS::DELTA, ui->dspbDeltaL->value() );
    }

    if( ui->chkFixT->isChecked() ){
        mi_rec3D.setFiltroParametros(FILTROS::PAR_K, (double)ui->spbFixedK->value() );
    }else{
        mi_rec3D.setFiltroParametros(FILTROS::PAR_K, FILTROS::INFERIOR, (double)ui->spbInfK->value() );
        mi_rec3D.setFiltroParametros(FILTROS::PAR_K, FILTROS::SUPERIOR, (double)ui->spbSupK->value() );
        mi_rec3D.setFiltroParametros(FILTROS::PAR_K, FILTROS::DELTA, ui->dspbDeltaK->value() );
    }

    if( ui->cmbFilter->currentIndex() == 0 ){
        if( ui->chkFixSigma->isChecked() ){
            mi_rec3D.setFiltroParametros(FILTROS::PAR_SIGMA, ui->dspbFixedSigma->value() );
        }else{
            mi_rec3D.setFiltroParametros(FILTROS::PAR_SIGMA, FILTROS::INFERIOR, ui->dspbInfSigma->value() );
            mi_rec3D.setFiltroParametros(FILTROS::PAR_SIGMA, FILTROS::SUPERIOR, ui->dspbSupSigma->value() );
            mi_rec3D.setFiltroParametros(FILTROS::PAR_SIGMA, FILTROS::DELTA, ui->dspbDeltaSigma->value() );
        }
    }

    mi_rec3D.segmentarImagenBase( 0 );
}
