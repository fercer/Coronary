#include "coronarygui.h"
#include "ui_coronarygui.h"

coronaryGUI::coronaryGUI(QWidget *parent) : QMainWindow(parent), ui(new Ui::coronaryGUI){
    ui->setupUi(this);

    //// Create context menu for vtk viewports
    qvtkVP1_CM = new QMenu(ui->qvtkVP1);

    //// Connect context menu actions
    connect(ui->qvtkVP1, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(on_qvtkVP1_customContextMenuRequested(QPoint)));


    //// Connect Reconstructor to viewports and log:
    ui->qvtkVP4->GetRenderWindow()->AddRenderer( mi_rec3D.getRenderer() );
    mi_rec3D.setLog( ui->ptxtLog );
}

coronaryGUI::~coronaryGUI(){
    delete ui;
}

void coronaryGUI::loadGroundtruth(){
    QString filename = QFileDialog::getOpenFileName(
                this,
                tr("Open file as image base"),
                "~",
                "Png (*.png);;Bitmap (*.bmp);;Jpeg (*.jpg,*.jpeg);;DICOM (*.dcm,*.*)"
                );

    QByteArray filename_ba = filename.toLatin1();
    mi_rec3D.agregarGroundtruth( filename_ba.data(), 0 );

    qvtkVP1_CM->removeAction( qvtkVP1_CM_LoadGroundtruth );
    qvtkVP1_CM_LoadGroundtruth = NULL;
    qvtkVP1_CM_ShowGroundtruth = qvtkVP1_CM->addAction("Show ground-truth");
    connect(qvtkVP1_CM_ShowGroundtruth,SIGNAL(triggered()),this,SLOT(showGroundtruth()));
}

void coronaryGUI::showGroundtruth()
{
    qvtkVP1_CM->removeAction( qvtkVP1_CM_ShowGroundtruth );
    qvtkVP1_CM_ShowGroundtruth = NULL;
    qvtkVP1_CM_ShowBase = qvtkVP1_CM->addAction("Show base");
    connect(qvtkVP1_CM_ShowBase,SIGNAL(triggered()),this,SLOT(showBase()));

    mi_rec3D.mostrarGroundtruth( 0 );
}

void coronaryGUI::showBase()
{
    qvtkVP1_CM->removeAction( qvtkVP1_CM_ShowBase );
    qvtkVP1_CM_ShowBase = NULL;
    qvtkVP1_CM_ShowGroundtruth = qvtkVP1_CM->addAction("Show ground-truth");
    connect(qvtkVP1_CM_ShowGroundtruth,SIGNAL(triggered()),this,SLOT(showGroundtruth()));

    mi_rec3D.mostrarBase( 0 );
}

void coronaryGUI::on_qvtkVP1_customContextMenuRequested(const QPoint &pos){
    qvtkVP1_CM->popup(mapToGlobal(pos));
}


void coronaryGUI::on_action_Open_file_triggered()
{
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
        vtkSmartPointer< vtkRenderer > nuevo_renderer = mi_rec3D.getRenderer( 0 );
        if( nuevo_renderer ){
            ui->qvtkVP1->GetRenderWindow()->AddRenderer( nuevo_renderer );
            qvtkVP1_CM_LoadGroundtruth = qvtkVP1_CM->addAction("Load file as groundtruth");
            connect(qvtkVP1_CM_LoadGroundtruth,SIGNAL(triggered()),this,SLOT(loadGroundtruth()));
            ui->action_Open_file->setEnabled(false);
            mi_rec3D.mostrarBase( 0 );
        }
    }
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
    mi_rec3D.skeletonize();
}
