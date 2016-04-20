
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pars_fercer.c"
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
}





int main(int argc, char** argv ){
/*
    // Step 1: Fake Image
      gdcm::SmartPointer<gdcm::Image> im = new gdcm::Image;

      char * buffer = new char[ 256 * 256 * 3];
      char * p = buffer;
      int b = 128;
      //int ybr[3];
      int rgb[3];
      //int rgb[3];

      for(int y = 0; y < 256; ++y){
        for(int x = 0; x < 256; ++x){
          *p++ = (char)x;
          *p++ = (char)y;
          *p++ = (char)128;
        }
      }

      im->SetNumberOfDimensions( 2 );
      im->SetDimension(0, 256 );
      im->SetDimension(1, 256 );

      im->GetPixelFormat().SetSamplesPerPixel(3);
      im->SetPhotometricInterpretation( gdcm::PhotometricInterpretation::RGB );

      unsigned long l = im->GetBufferLength();
      if( l != 256 * 256 * 3 )
        {
        return 1;
        }

      gdcm::DataElement pixeldata( gdcm::Tag(0x7fe0,0x0010) );
      pixeldata.SetByteValue( buffer, (uint32_t)l );
      delete[] buffer;
      im->SetDataElement( pixeldata );

      gdcm::UIDGenerator uid; // helper for uid generation

      gdcm::SmartPointer<gdcm::File> file = new gdcm::File; // empty file

      // Step 2: DERIVED object
      gdcm::FileDerivation fd;
      // For the pupose of this execise we will pretend that this image is referencing
      // two source image (we need to generate fake UID for that).
      const char ReferencedSOPClassUID[] = "1.2.840.10008.5.1.4.1.1.7"; // Secondary Capture
      fd.AddReference( ReferencedSOPClassUID, uid.Generate() );
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
      if( !fd.Derive() )
        {
        std::cerr << "Sorry could not derive using input info" << std::endl;
        return 1;
        }

      // We pass both :
      // 1. the fake generated image
      // 2. the 'DERIVED' dataset object
      // to the writer.
      gdcm::ImageWriter w;
      w.SetImage( *im );
      w.SetFile( fd.GetFile() );

      // Set the filename:
      w.SetFileName( "rgb_new.dcm" );
      if( !w.Write() )
        {
        return 1;
        }

      return 0;
*/
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
    imagecomments.SetValue( "Hello, World !" );

    gdcm::Attribute<0x0018,0x1110> nuevo_DSD;
    nuevo_DSD.SetValue( 1118.0 );

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

    // Definir los parametros de entrada:
    PARS_ENTRADA *parametros = new PARS_ENTRADA [12];
    definirParametros( parametros );

    if( argc < 2 ){
        mostrar_ayuda(parametros, 12, "Coronary");
        delete [] parametros;
        return EXIT_FAILURE;
    }
    // Revisar los parametros de entrada:
    revisar_pars(parametros, 12, &argc, argv);

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
    reconstructor.skeletonize();

    delete [] parametros;

    return EXIT_SUCCESS;

}
