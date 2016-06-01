#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pars_fercer.c"
#include "IMGVTK.h"
#include "reconstructor_3D.h"

#include <iostream>

#include "gdcmReader.h"
#include "gdcmWriter.h"
#include "gdcmAttribute.h"
#include "gdcmImage.h"
#include "gdcmImageWriter.h"
#include "gdcmFileDerivation.h"
#include "gdcmUIDGenerator.h"

#include <iostream>

/*  Funcion: definirParametros
    Descripcion: Define los parametros de entrada del programa (Hard-coded).
*/
void definirParametros(PARS_ENTRADA *parametros){
    // Parametro input base:
    parametros[0].mi_tipo = CHAR;
    sprintf(parametros[0].short_tag, "-b1");
    sprintf(parametros[0].long_tag, "--base1");
    sprintf(parametros[0].mi_default.par_s, "ang_1.png");
    sprintf(parametros[0].pregunta, "Imagen angiografica BASE 1 de entrada (.PNG, .BMP, .JPEG/.JPG, o archivo DICOM)");
    parametros[0].opcional = 0;
    // Parametro nivel extraido del archivo DICOM:
    parametros[1].mi_tipo = INT;
    sprintf(parametros[1].short_tag, "-l1");
    sprintf(parametros[1].long_tag, "--desde1");
    parametros[1].mi_default.par_i = 0;
    sprintf(parametros[1].pregunta, "Nivel desde el que se extraen las imagenes del archivo DICOM como imagenes base 1");
    parametros[1].opcional = 1;
    // Parametro nivel extraido del archivo DICOM:
    parametros[2].mi_tipo = INT;
    sprintf(parametros[2].short_tag, "-u1");
    sprintf(parametros[2].long_tag, "--hasta1");
    parametros[2].mi_default.par_i = 0;
    sprintf(parametros[2].pregunta, "Nivel hasta el que se extraen las iamgenes del archivo DICOM como imagenes base 1");
    parametros[2].opcional = 1;


    // Parametro input base:
    parametros[3].mi_tipo = CHAR;
    sprintf(parametros[3].short_tag, "-b2");
    sprintf(parametros[3].long_tag, "--base2");
    sprintf(parametros[3].mi_default.par_s, "NULL");
    sprintf(parametros[3].pregunta, "Imagen angiografica BASE 2 de entrada (.PNG, .BMP, .JPEG/.JPG, o archivo DICOM)");
    parametros[3].opcional = 1;
    // Parametro nivel extraido del archivo DICOM:
    parametros[4].mi_tipo = INT;
    sprintf(parametros[4].short_tag, "-l2");
    sprintf(parametros[4].long_tag, "--desde2");
    parametros[4].mi_default.par_i = 0;
    sprintf(parametros[4].pregunta, "Nivel desde el que se extraen las imagenes del archivo DICOM como imagenes base 2");
    parametros[4].opcional = 1;
    // Parametro nivel extraido del archivo DICOM:
    parametros[5].mi_tipo = INT;
    sprintf(parametros[5].short_tag, "-u2");
    sprintf(parametros[5].long_tag, "--hasta2");
    parametros[5].mi_default.par_i = 0;
    sprintf(parametros[5].pregunta, "Nivel hasta el que se extraen las iamgenes del archivo DICOM como imagenes base 2");
    parametros[5].opcional = 1;


    // Parametro input base:
    parametros[6].mi_tipo = CHAR;
    sprintf(parametros[6].short_tag, "-b3");
    sprintf(parametros[6].long_tag, "--base3");
    sprintf(parametros[6].mi_default.par_s, "NULL");
    sprintf(parametros[6].pregunta, "Imagen angiografica BASE 3 de entrada (.PNG, .BMP, .JPEG/.JPG, o archivo DICOM)");
    parametros[6].opcional = 1;
    // Parametro nivel extraido del archivo DICOM:
    parametros[7].mi_tipo = INT;
    sprintf(parametros[7].short_tag, "-l3");
    sprintf(parametros[7].long_tag, "--desde3");
    parametros[7].mi_default.par_i = 0;
    sprintf(parametros[7].pregunta, "Nivel desde el que se extraen las imagenes del archivo DICOM como imagenes base 3");
    parametros[7].opcional = 1;
    // Parametro nivel extraido del archivo DICOM:
    parametros[8].mi_tipo = INT;
    sprintf(parametros[8].short_tag, "-u3");
    sprintf(parametros[8].long_tag, "--hasta3");
    parametros[8].mi_default.par_i = 0;
    sprintf(parametros[8].pregunta, "Nivel hasta el que se extraen las iamgenes del archivo DICOM como imagenes base 3 ");
    parametros[8].opcional = 1;




    // Parametro input ground:
    parametros[9].mi_tipo = CHAR;
    sprintf(parametros[9].short_tag, "-g1");
    sprintf(parametros[9].long_tag, "--ground1");
    sprintf(parametros[9].mi_default.par_s, "NULL");
    sprintf(parametros[9].pregunta, "Imagen usada como ground truth (.PNG, .BMP, .JPEG/.JPG) para la imagen base 1");
    parametros[9].opcional = 1;

    // Parametro input ground:
    parametros[10].mi_tipo = CHAR;
    sprintf(parametros[10].short_tag, "-g2");
    sprintf(parametros[10].long_tag, "--ground2");
    sprintf(parametros[10].mi_default.par_s, "NULL");
    sprintf(parametros[10].pregunta, "Imagen usada como ground truth (.PNG, .BMP, .JPEG/.JPG) para la imagen base 2");
    parametros[10].opcional = 1;

    // Parametro input ground:
    parametros[11].mi_tipo = CHAR;
    sprintf(parametros[11].short_tag, "-g3");
    sprintf(parametros[11].long_tag, "--ground3");
    sprintf(parametros[11].mi_default.par_s, "NULL");
    sprintf(parametros[11].pregunta, "Imagen usada como ground truth (.PNG, .BMP, .JPEG/.JPG) para la imagen base 3");
    parametros[11].opcional = 1;



    ////---------------------------------------------------------------------- PARA GENERAR IMAGENES PHANTOM EN ARCHIVO DICOM
    // Parametro output phantom
    parametros[12].mi_tipo = CHAR;
    sprintf(parametros[12].short_tag, "-ph");
    sprintf(parametros[12].long_tag, "--phantom");
    sprintf(parametros[12].mi_default.par_s, "NULL");
    sprintf(parametros[12].pregunta, "[GENERAR ARCHIVO DICOM PARA UN PHANTOM] Ruta del arhcivo de la imagen phantom");
    parametros[12].opcional = 1;

    // Parametro output phantom
    parametros[13].mi_tipo = CHAR;
    sprintf(parametros[13].short_tag, "-op");
    sprintf(parametros[13].long_tag, "--output");
    sprintf(parametros[13].mi_default.par_s, "NULL");
    sprintf(parametros[13].pregunta, "[GENERAR ARCHIVO DICOM PARA UN PHANTOM] Ruta de salida del archivo DICOM generado en base al phantom");
    parametros[13].opcional = 1;
}



/*  Funcion: genDICOM
    Descripcion: Genera un archivo DICOM en base a una imagen Phantom.
*/
void genDICOM(const char *ruta_origen, const char *ruta_img_phantom, const char *ruta_salida){
    // Abrir archivo DICOM base:
    gdcm::ImageReader DICOMreader;
    DICOMreader.SetFileName( ruta_origen );
    DICOMreader.Read();

    gdcm::File &file = DICOMreader.GetFile();
    gdcm::DataSet &ds = file.GetDataSet();


    // Abrir la imagen:
    IMGVTK img_phantom( ruta_img_phantom, false, 0);
    const int mis_cols = img_phantom.cols;
    const int mis_rens = img_phantom.rens;
    const int mis_rens_cols = img_phantom.rens_cols;

    // Alojar memoria para la imagen dentro del archivo DICOM:
    gdcm::SmartPointer<gdcm::Image> im = new gdcm::Image;

    char * buffer = new char[mis_rens_cols];

    im->SetNumberOfDimensions( 2 );
    im->SetDimension(0, mis_cols );
    im->SetDimension(1, mis_rens);

    // Definir el espaciado entre pixeles
    double pixX = 1.0, pixY = 1.0;

    const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1164) );
    const gdcm::ByteValue *bv = de.GetByteValue();
    if( bv ){
        std::string strm(bv->GetPointer(), bv->GetLength());
DEB_MSG("pixXY: " << strm);
        char *pixXYstr = new char [bv->GetLength()];
        memcpy(pixXYstr, strm.c_str(), bv->GetLength() * sizeof(char ));
        char *tmp = strchr(pixXYstr,'\\');
        pixXYstr[ tmp - pixXYstr ] = '\0';
        pixY = atof(pixXYstr);// / Magnification;
        pixX = atof(tmp+1);// / Magnification;
DEB_MSG("pixY: " << pixY << ", pixX: " << pixX);
        delete [] pixXYstr;
    }
    im->SetSpacing(0, pixY);
    im->SetSpacing(1, pixX);

    // Almacenar la imagen phantom al archivo DICOM:
    for( int y = 0; y < mis_rens; y++ ){
        for( int x = 0; x < mis_cols; x++ ){
            *(buffer + x + y*mis_cols) = (char) (255.0 * *(img_phantom.base_ptr + x + (mis_rens - y - 1)*mis_cols));
        }
    }

    im->GetPixelFormat().SetSamplesPerPixel(1);
    im->SetPhotometricInterpretation( gdcm::PhotometricInterpretation::MONOCHROME2 );

    unsigned long l = im->GetBufferLength();

    if( l != mis_rens_cols ){
        std::cout << "\33[44m" << "<<Error al generar la imagen para el archivo DICOM>>" << "\33[0m" << std::endl;
        delete[] buffer;
        return;
    }

    gdcm::DataElement pixeldata( gdcm::Tag(0x7fe0,0x0010) );
    pixeldata.SetByteValue( buffer, (uint32_t)l );
    im->SetDataElement( pixeldata );

    gdcm::Attribute<0x0020,0x4000> imagecomments;
    imagecomments.SetValue( "Phantom" );
    ds.Replace( imagecomments.GetAsDataElement() );

    // Incluir los datos al archivo DICOM
    gdcm::ImageWriter w;
    w.CheckFileMetaInformationOff();
    w.SetImage( *im );
    w.SetFile( file );

    // Set the filename:
    w.SetFileName( ruta_salida );
    if( !w.Write() ){
        std::cout << "\33[14m" << "<<Error al guardar el archivo DICOM>>" << "\33[0m" << std::endl;
        return;
    }

    return;
}


int main(int argc, char** argv ){
    // Definir los parametros de entrada:
    PARS_ENTRADA *parametros = new PARS_ENTRADA [14];
    definirParametros( parametros );

    if( argc < 2 ){
        mostrar_ayuda(parametros, 14, "Coronary");
        delete [] parametros;
        return EXIT_FAILURE;
    }
    // Revisar los parametros de entrada:
    revisar_pars(parametros, 14, &argc, argv);


    // Si se va a generar un archivo DICOM par aun phantom, no se genera el reconstructor 3D:
    if( strcmp( parametros[12].mi_valor.par_s , "NULL" ) && strcmp( parametros[13].mi_valor.par_s , "NULL" )){
        /// Generar archivo DICOM
        genDICOM( parametros[0].mi_valor.par_s, parametros[12].mi_valor.par_s, parametros[13].mi_valor.par_s);
    }else{
        /// Reconstruir arteria:
        RECONS3D reconstructor;

        DEB_MSG("Base: " << parametros[0].mi_valor.par_s <<  ", Ground truth: " << parametros[9].mi_valor.par_s);
        reconstructor.agregarInput(parametros[0].mi_valor.par_s, parametros[1].mi_valor.par_i, parametros[2].mi_valor.par_i, parametros[9].mi_valor.par_s);

        if( strcmp(parametros[3].mi_valor.par_s, "NULL") ){
            DEB_MSG("Base (" << strcmp(parametros[3].mi_valor.par_s, "NULL") << " " << (strcmp(parametros[3].mi_valor.par_s, "NULL") == 0) << "): " << parametros[3].mi_valor.par_s <<  ", Ground truth: " << parametros[10].mi_valor.par_s);
            reconstructor.agregarInput(parametros[3].mi_valor.par_s, parametros[4].mi_valor.par_i, parametros[5].mi_valor.par_i, parametros[10].mi_valor.par_s);
        }
        if( strcmp(parametros[6].mi_valor.par_s, "NULL")){
            DEB_MSG("Base: " << parametros[6].mi_valor.par_s <<  ", Ground truth: " << parametros[11].mi_valor.par_s);
            reconstructor.agregarInput(parametros[6].mi_valor.par_s, parametros[7].mi_valor.par_i, parametros[8].mi_valor.par_i, parametros[11].mi_valor.par_s);
        }

    ///*********************************************************************************************
    ///             BUSQUEDA EXHAUSTIVA DE LOS PARAMETROS:
    ///*********************************************************************************************

    //    char *rutas_base[20];
    //    char *rutas_grd[20];
    //    char ruta_dir[] = "/home/fercer/test_data/40Gabor/ang_";

    //    for(int i = 0; i < 20; i++){
    //        rutas_base[i] = new char[256];
    //        sprintf(rutas_base[i], "%s%i.png", ruta_dir, i+1);
    //        rutas_grd[i] = new char[256];
    //        sprintf(rutas_grd[i], "%s%i_gt.png", ruta_dir, i+1);;
    //    }

    //    RECONS3D reconstructor( rutas_base, rutas_grd, 20 );

    //    for(int i = 0; i < 20; i++){
    //        delete [] rutas_base[i];
    //        delete [] rutas_grd[i];
    //    }



        reconstructor.setFiltroEntrenamiento( FILTROS::EXHAUSTIVA, 0, 0);
        reconstructor.setFiltroEval( FILTROS::CORCON );
        reconstructor.setFiltroMetodo( FILTROS::SS_GABOR );
        reconstructor.setFiltroParametros( FILTROS::PAR_L, 2.9);//, 5.0, 0.1);
        reconstructor.setFiltroParametros( FILTROS::PAR_T, 11.0 );
        reconstructor.setFiltroParametros( FILTROS::PAR_K, 45.0);
        //reconstructor.setFiltroParametros( FILTROS::PAR_SIGMA, 1.5, 2.3, 0.001);

        /*
        reconstructor.setFiltroEntrenamiento( FILTROS::EDA_BUMDA, 2, 4);
        reconstructor.setFiltroEval( FILTROS::ROC );
        reconstructor.setFiltroMetodo( FILTROS::GMF );
        reconstructor.setFiltroLimites( FILTROS::PAR_L, 8.0, 15.0, 1e-4);
        reconstructor.setFiltroLimites( FILTROS::PAR_T, 8.0, 15.0, 1e-4);
        reconstructor.setFiltroLimites( FILTROS::PAR_SIGMA, 1.0, 5.0, 1e-4);
        reconstructor.setFiltroParametros( FILTROS::PAR_K, 12.0);
        */


        reconstructor.segmentarImagenBase();
        //reconstructor.skeletonize();
    }
    delete [] parametros;
    return EXIT_SUCCESS;
}
