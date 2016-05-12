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
}

coronaryGUI::~coronaryGUI()
{
    delete ui;

}

void coronaryGUI::on_actionOpen_File_triggered()
{
    QString filename = QFileDialog::getOpenFileName(
                    this,
                    tr("Open File"),
                    "/home/",
                    "Portable Network Graphics (*.png);;Bitmap Image File (*.bmp);;Joint Photographic Experts Group (*.jpeg);;DICOM File (*.*)"
                );
    if( filename.length() > 0 ){
        QByteArray filename_c = filename.toLatin1();
        QString filename_message = "Opening input: ";
        filename_message.append( filename_c );
        QMessageBox::information( this, "Open input image", filename_message );
        mi_reconstructor.agregarInput( filename_c.data() );

        vtkSmartPointer< vtkRenderer > nuevo_renderer = mi_reconstructor.getRenderer(0);
        if( nuevo_renderer ){
            ui->vtkVP1->GetRenderWindow()->AddRenderer( nuevo_renderer );
        }
    }
}
