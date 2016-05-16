#include "coronarygui.h"
#include "ui_coronarygui.h"
#include <QFileDialog>
#include <QMessageBox>

coronaryGUI::coronaryGUI(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::coronaryGUI)
{
    ui->setupUi(this);
    mi_reconstructor.setLog( ui->ptxtLog );
    ui->vtkVP4->GetRenderWindow()->AddRenderer( mi_reconstructor.getRenderer() );
}

coronaryGUI::~coronaryGUI()
{
    delete ui;

}

void coronaryGUI::on_actionQuit_triggered()
{
    close();
}

void coronaryGUI::on_actionBase_triggered()
{
    QString filename = QFileDialog::getOpenFileName(
                    this,
                    tr("Open File"),
                    "/home/",
                    "Portable Network Graphics (*.png);;Bitmap Image File (*.bmp);;Joint Photographic Experts Group (*.jpeg);;DICOM File (*.*)"
                );
    if( filename.length() > 0 ){
        QByteArray filename_c = filename.toLatin1();
        mi_reconstructor.agregarInput( filename_c.data() );

        vtkSmartPointer< vtkRenderer > nuevo_renderer = mi_reconstructor.getRenderer(0);
        if( nuevo_renderer ){
            ui->vtkVP1->GetRenderWindow()->AddRenderer( nuevo_renderer );
        }
    }
}

void coronaryGUI::on_actionGround_truth_triggered()
{
    QString filename = QFileDialog::getOpenFileName(
                    this,
                    tr("Open File as Ground-truth"),
                    "/home/",
                    "Portable Network Graphics (*.png);;Bitmap Image File (*.bmp);;Joint Photographic Experts Group (*.jpeg);;DICOM File (*.*)"
                );
    if( filename.length() > 0 ){
        QByteArray filename_c = filename.toLatin1();
        mi_reconstructor.agregarGroundtruth(filename_c.data(), 0);

        vtkSmartPointer< vtkRenderer > nuevo_renderer = mi_reconstructor.getRenderer(0);
        if( nuevo_renderer ){
            ui->vtkVP1->GetRenderWindow()->AddRenderer( nuevo_renderer );
        }
    }
}

void coronaryGUI::on_action_Skeletonize_triggered()
{
    mi_reconstructor.skeletonize();
}
