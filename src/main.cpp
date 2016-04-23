
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
    sprintf(parametros[12].short_tag, "-op");
    sprintf(parametros[12].long_tag, "--output-phantom");
    sprintf(parametros[12].mi_default.par_s, "NULL");
    sprintf(parametros[12].pregunta, "[GENERAR ARHCIVO DICOM PARA UN PHANTOM] Ruta de salida del archivo DICOM generado en base al phantom");
    parametros[12].opcional = 1;

    // Parametro SID nuevo
    parametros[13].mi_tipo = DOUBLE;
    sprintf(parametros[13].short_tag, "-SI");
    sprintf(parametros[13].long_tag, "--SID");
    parametros[13].mi_default.par_d = 1000.0;
    sprintf(parametros[13].pregunta, "[GENERAR ARHCIVO DICOM PARA UN PHANTOM] Distancia de la fuente al detector para la imagen phantom generada");
    parametros[13].opcional = 1;

    // Parametro SOD nuevo
    parametros[14].mi_tipo = DOUBLE;
    sprintf(parametros[14].short_tag, "-SO");
    sprintf(parametros[14].long_tag, "--SOD");
    parametros[14].mi_default.par_d = 400.0;
    sprintf(parametros[14].pregunta, "[GENERAR ARHCIVO DICOM PARA UN PHANTOM] Distancia de la fuente al paciente para la imagen phantom generada");
    parametros[14].opcional = 1;

    // Parametro SISO nuevo
    parametros[15].mi_tipo = DOUBLE;
    sprintf(parametros[15].short_tag, "-ISO");
    sprintf(parametros[15].long_tag, "--S-ISO");
    parametros[15].mi_default.par_d = 400.0;
    sprintf(parametros[15].pregunta, "[GENERAR ARHCIVO DICOM PARA UN PHANTOM] Distancia de la fuente al iso-centro para la imagen phantom generada");
    parametros[15].opcional = 1;

    // Parametro LAO/RAO nuevo
    parametros[16].mi_tipo = DOUBLE;
    sprintf(parametros[16].short_tag, "-LR");
    sprintf(parametros[16].long_tag, "--LAO-RAO");
    parametros[16].mi_default.par_d = 0.0;
    sprintf(parametros[16].pregunta, "[GENERAR ARHCIVO DICOM PARA UN PHANTOM] Angulo LAO(-)/RAO(+)");
    parametros[16].opcional = 1;

    // Parametro CRA/CAU nuevo
    parametros[17].mi_tipo = DOUBLE;
    sprintf(parametros[17].short_tag, "-CC");
    sprintf(parametros[17].long_tag, "--CRA-CUA");
    parametros[17].mi_default.par_d = 0.0;
    sprintf(parametros[17].pregunta, "[GENERAR ARHCIVO DICOM PARA UN PHANTOM] Angulo CRA(-)/CAU(+)");
    parametros[17].opcional = 1;

    // Parametro pixX nuevo
    parametros[18].mi_tipo = DOUBLE;
    sprintf(parametros[18].short_tag, "-px");
    sprintf(parametros[18].long_tag, "--pix-X");
    parametros[18].mi_default.par_d = 0.2;
    sprintf(parametros[18].pregunta, "[GENERAR ARHCIVO DICOM PARA UN PHANTOM] Espaciado en mm entre el centro de un pixel a otro en el centro de la imagen");
    parametros[18].opcional = 1;

    // Parametro pixY nuevo
    parametros[19].mi_tipo = DOUBLE;
    sprintf(parametros[19].short_tag, "-py");
    sprintf(parametros[19].long_tag, "--pix-Y");
    parametros[19].mi_default.par_d = 0.2;
    sprintf(parametros[19].pregunta, "[GENERAR ARHCIVO DICOM PARA UN PHANTOM] Espaciado en mm entre el centro de un pixel a otro en el centro de la imagen");
    parametros[19].opcional = 1;
}



/*  Funcion: genDICOM
    Descripcion: Genera un archivo DICOM en base a una imagen Phantom.
*/
void genDICOM(const char *ruta_origen, const char *ruta_salida, const double SID, const double SOD, const double ISO, const double LAORAO, const double CRACAU, const double pixX, const double pixY ){

    // Abrir la imagen:
    IMGVTK img_phantom( ruta_origen, false, 0);
    const int mis_cols = img_phantom.cols;
    const int mis_rens = img_phantom.rens;
    const int mis_rens_cols = img_phantom.rens_cols;

    // Alojar memoria para la imagen dentro del archivo DICOM:
    gdcm::SmartPointer<gdcm::Image> im = new gdcm::Image;

    char * buffer = new char[mis_rens_cols];

    im->SetNumberOfDimensions( 2 );
    im->SetDimension(0, mis_cols );
    im->SetDimension(1, mis_rens);

    // Almacenar la imagen phantom al archivo DICOM:
    for( int y = 0; y < mis_rens; y++ ){
        for( int x = 0; x < mis_cols; x++ ){
            *(buffer + x + y*mis_cols) = (char) (255.0 * *(img_phantom.base_ptr + x + (mis_rens - y - 1)*mis_cols));
        }
    }

    im->GetPixelFormat().SetSamplesPerPixel(1);
    im->SetPhotometricInterpretation( gdcm::PhotometricInterpretation::MONOCHROME1 );

    unsigned long l = im->GetBufferLength();

    if( l != mis_rens_cols ){
        std::cout << "\33[44m" << "<<Error al generar la imagen para el archivo DICOM>>" << "\33[0m" << std::endl;
        delete[] buffer;
        return;
    }

    gdcm::DataElement pixeldata( gdcm::Tag(0x7fe0,0x0010) );
    pixeldata.SetByteValue( buffer, (uint32_t)l );
    im->SetDataElement( pixeldata );

    // Agregar la informacion:
    gdcm::UIDGenerator uid; // helper for uid generation

    gdcm::SmartPointer<gdcm::File> file = new gdcm::File; // empty file

    // Step 2: DERIVED object
    gdcm::FileDerivation fd;

    // For the pupose of this execise we will pretend that this image is referencing
    // two source image (we need to generate fake UID for that).
    const char ReferencedSOPClassUID[] = "1.2.840.10008.5.1.4.1.1.12.1"; // Secondary Capture
    fd.AddReference( ReferencedSOPClassUID, uid.Generate() );

    // Again for the purpose of the exercise we will pretend that the image is a
    // multiplanar reformat (MPR):
    // CID 7202 Source Image Purposes of Reference
    // {"DCM",121322,"Source image for image processing operation"},
    fd.SetPurposeOfReferenceCodeSequenceCodeValue( 121322 );

    // CID 7203 Image Derivation
    // { "DCM",113072,"Multiplanar reformatting" },
    fd.SetDerivationCodeSequenceCodeValue( 113072 );
    fd.SetFile( *file );



    // If all Code Value are ok the filter will execute properly
    if( !fd.Derive() ){
      std::cerr << "Sorry could not derive using input info" << std::endl;
      return;
    }
    gdcm::DataSet &ds = fd.GetFile().GetDataSet();

    // Decir que la imagen es PHANTOM:
    gdcm::Attribute<0x0020,0x4000> imagecomments;
    imagecomments.SetValue( "Phantom" );

    gdcm::Attribute<0x0018,0x1110> nuevo_SID;
    nuevo_SID.SetValue( SID );

    gdcm::Attribute<0x0018,0x1111> nuevo_SOD;
    nuevo_SOD.SetValue( SOD );

    gdcm::Attribute<0x0018,0x1111> nuevo_ISO;
    nuevo_ISO.SetValue( ISO );

    gdcm::Attribute<0x0018,0x1111> nuevo_LAORAO;
    nuevo_LAORAO.SetValue( LAORAO );

    gdcm::Attribute<0x0018,0x1111> nuevo_CRACAU;
    nuevo_CRACAU.SetValue( CRACAU );

    gdcm::Attribute<0x0018,0x1111> nuevo_pix;
    char pixXY[256];
    sprintf(pixXY, "%f\\%f", pixX, pixY);
    nuevo_pix.SetValue( pixXY );

    gdcm::Attribute<0x0018,0x1111> nuevo_winCenter;
    nuevo_winCenter.SetValue( 127.5 );

    gdcm::Attribute<0x0018,0x1111> nuevo_winWidth;
    nuevo_winWidth.SetValue( 255.0 );

    // Now replace the Image Comments from the dataset with our:
    ds.Replace( imagecomments.GetAsDataElement() );
    ds.Replace( nuevo_SID.GetAsDataElement() );
    ds.Replace( nuevo_SOD.GetAsDataElement() );
    ds.Replace( nuevo_ISO.GetAsDataElement() );
    ds.Replace( nuevo_LAORAO.GetAsDataElement() );
    ds.Replace( nuevo_CRACAU.GetAsDataElement() );
    ds.Replace( nuevo_pix.GetAsDataElement() );
    ds.Replace( nuevo_winCenter.GetAsDataElement() );
    ds.Replace( nuevo_winWidth.GetAsDataElement() );

    // Incluir los datos al archivo DICOM
    gdcm::ImageWriter w;
    w.CheckFileMetaInformationOff();
    w.SetImage( *im );
    w.SetFile( fd.GetFile() );

    // Set the filename:
    w.SetFileName( ruta_salida );
    if( !w.Write() ){
        std::cout << "\33[14m" << "<<Error al guardar el archivo DICOM>>" << "\33[0m" << std::endl;
        return;
    }

    return;

/*
         if( argc < 3 )
          {
          std::cerr << argv[0] << " input.dcm output.dcm" << std::endl;
          return 1;
          }
        const char *filename = argv[1];
        const char *outfilename = argv[2];

        // Instanciate the reader:
        gdcm::Reader reader;
        reader.SetFileName( filename );
        if( !reader.Read() )
          {
          std::cerr << "Could not read: " << filename << std::endl;
          return 1;
          }

        // If we reach here, we know for sure only 1 thing:
        // It is a valid DICOM file (potentially an old ACR-NEMA 1.0/2.0 file)
        // (Maybe, it's NOT a Dicom image -could be a DICOMDIR, a RTSTRUCT, etc-)

        // The output of gdcm::Reader is a gdcm::File
        gdcm::File &file = reader.GetFile();

        // the dataset is the the set of element we are interested in:
        gdcm::DataSet &ds = file.GetDataSet();

        // Contruct a static(*) type for Image Comments :
        gdcm::Attribute<0x0020,0x4000> imagecomments;
        imagecomments.SetValue( "Phantom" );

        gdcm::Attribute<0x0018,0x1110> nuevo_SID;
        nuevo_SID.SetValue( 1118.0 );

        gdcm::Attribute<0x0018,0x1110> nuevo_SID;
        nuevo_SID.SetValue( 1118.0 );

        gdcm::Attribute<0x0019,0x1000> nuevo_ReviewMode;
        nuevo_ReviewMode.SetValue( NULL );

        // Now replace the Image Comments from the dataset with our:
        ds.Replace( imagecomments.GetAsDataElement() );
        ds.Replace( nuevo_DSD.GetAsDataElement() );


        // Write the modified DataSet back to disk
        gdcm::Writer writer;
        writer.CheckFileMetaInformationOff(); // Do not attempt to reconstruct the file meta to preserve the file
                                              // as close to the original as possible.
        writer.SetFileName( outfilename );
        writer.SetFile( file );
        if( !writer.Write() )
          {
          std::cerr << "Could not write: " << outfilename << std::endl;
          return 1;
          }

        return 0;
*/
}


int main(int argc, char** argv ){
    // Definir los parametros de entrada:
    PARS_ENTRADA *parametros = new PARS_ENTRADA [20];
    definirParametros( parametros );

    if( argc < 2 ){
        mostrar_ayuda(parametros, 20, "Coronary");
        delete [] parametros;
        return EXIT_FAILURE;
    }
    // Revisar los parametros de entrada:
    revisar_pars(parametros, 20, &argc, argv);


    // Si se va a generar un archivo DICOM par aun phantom, no se genera el reconstructor 3D:
    if( strcmp( parametros[12].mi_valor.par_s , "NULL" ) ){ /// Generar archivo DICOM
        genDICOM( parametros[0].mi_valor.par_s, parametros[12].mi_valor.par_s, parametros[13].mi_valor.par_d, parametros[14].mi_valor.par_d, parametros[15].mi_valor.par_d, parametros[16].mi_valor.par_d, parametros[17].mi_valor.par_d, parametros[18].mi_valor.par_d, parametros[19].mi_valor.par_d);

    }else{ /// Reconstruir arteria:
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


    //    reconstructor.segmentarImagenBase();
    //    reconstructor.skeletonize();



    }
    delete [] parametros;
    return EXIT_SUCCESS;
}
