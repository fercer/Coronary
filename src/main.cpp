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
    sprintf(parametros[0].mi_default.par_s, "NULL");
    sprintf(parametros[0].pregunta, "Imagen angiografica BASE 1 de entrada (.PNG, .BMP, .JPEG/.JPG, o archivo DICOM)");
    parametros[0].opcional = 1;
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


    ////---------------------------------------------------------------------- PARA CARGAR VARIAS IMAGENES Y REALIZAR UN ENTRENAMIENTO DE LOS PARAMETROS
    // Parametro input data.set
    parametros[14].mi_tipo = CHAR;
    sprintf(parametros[14].short_tag, "-dsb");
    sprintf(parametros[14].long_tag, "--datasetbase");
    sprintf(parametros[14].mi_default.par_s, "NULL");
    sprintf(parametros[14].pregunta, "[ENTRENAR PARAMETROS] Ruta del archivo data.set con la lista de rutas de imagenes a ser usadas como base");
    parametros[14].opcional = 1;

    // Parametro input data.set
    parametros[15].mi_tipo = CHAR;
    sprintf(parametros[15].short_tag, "-dsg");
    sprintf(parametros[15].long_tag, "--datasetground");
    sprintf(parametros[15].mi_default.par_s, "NULL");
    sprintf(parametros[15].pregunta, "[ENTRENAR PARAMETROS] Ruta del archivo data.set con la lista de rutas de imagenes a ser usadas como ground-truth");
    parametros[15].opcional = 1;

    // Parametro input configuracion.dat
    parametros[16].mi_tipo = CHAR;
    sprintf(parametros[16].short_tag, "-c");
    sprintf(parametros[16].long_tag, "--condfig");
    sprintf(parametros[16].mi_default.par_s, "NULL");
    sprintf(parametros[16].pregunta, "[ENTRENAR PARAMETROS] Ruta del archivo configuracion.dat con los parametros del filtro utilizado en la deteccion");
    parametros[16].opcional = 1;

    // Parametro input filtrado.log
    parametros[17].mi_tipo = CHAR;
    sprintf(parametros[17].short_tag, "-lf");
    sprintf(parametros[17].long_tag, "--logfiltrado");
    sprintf(parametros[17].mi_default.par_s, "NULL");
    sprintf(parametros[17].pregunta, "[ENTRENAR PARAMETROS] Ruta donde se guarda el log del proceso de filtrado");
    parametros[17].opcional = 1;

    // Parametro output.pgm
    parametros[18].mi_tipo = CHAR;
    sprintf(parametros[18].short_tag, "-pgm");
    sprintf(parametros[18].long_tag, "--rutaPGM");
    sprintf(parametros[18].mi_default.par_s, "NULL");
    sprintf(parametros[18].pregunta, "[EXPORTAR COMO PGM] Ruta donde se guarda la imagen PGM para la ejecucion del No-GUI (si se considera como )");
    parametros[18].opcional = 1;

    // Parametro dataset as single image
    parametros[19].mi_tipo = CHAR;
    sprintf(parametros[19].short_tag, "-con");
    sprintf(parametros[19].long_tag, "--concatenar");
    sprintf(parametros[19].mi_default.par_s, "yes");
    sprintf(parametros[19].pregunta, "[ENTRENAR PARAMETROS] Concatenar el dataseten una sola imagen");
    parametros[19].opcional = 1;
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
    const int mis_rows = img_phantom.rows;
    const int mis_rows_cols = img_phantom.rows_cols;

    // Alojar memoria para la imagen dentro del archivo DICOM:
    gdcm::SmartPointer<gdcm::Image> im = new gdcm::Image;

    char * buffer = new char[mis_rows_cols];

    im->SetNumberOfDimensions( 2 );
    im->SetDimension(0, mis_cols );
    im->SetDimension(1, mis_rows);

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
    for( int y = 0; y < mis_rows; y++ ){
        for( int x = 0; x < mis_cols; x++ ){
            *(buffer + x + y*mis_cols) = (char) (255.0 * *(img_phantom.base_ptr + x + (mis_rows - y - 1)*mis_cols));
        }
    }

    im->GetPixelFormat().SetSamplesPerPixel(1);
    im->SetPhotometricInterpretation( gdcm::PhotometricInterpretation::MONOCHROME2 );

    unsigned long l = im->GetBufferLength();

    if( l != mis_rows_cols ){
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
    if( strcmp( parametros[12].mi_valor.par_s , "NULL" ) && strcmp( parametros[13].mi_valor.par_s , "NULL" )){
        /// Generar archivo DICOM
        genDICOM( parametros[0].mi_valor.par_s, parametros[12].mi_valor.par_s, parametros[13].mi_valor.par_s);
    }else if( strcmp(parametros[0].mi_valor.par_s, "NULL") && strcmp(parametros[16].mi_valor.par_s, "NULL") ){
        /// Reconstruir arteria:
        RECONS3D reconstructor;

        reconstructor.agregarInput(parametros[0].mi_valor.par_s, parametros[1].mi_valor.par_i, parametros[2].mi_valor.par_i, parametros[9].mi_valor.par_s, true);
        reconstructor.leerConfiguracion( parametros[16].mi_valor.par_s );

        if( strcmp( parametros[17].mi_valor.par_s, "NULL") ){
            reconstructor.setFiltroLog( parametros[17].mi_valor.par_s );
        }

        reconstructor.segmentarImagenBase( 0 );
        reconstructor.umbralizar(0, IMGVTK::OTSU, 0);
        reconstructor.lengthFilter(0, IMGVTK::THRESHOLD, 2000);
        reconstructor.Guardar("otsu.png", IMGVTK::THRESHOLD, IMGVTK::PNG, 0);
        reconstructor.skeletonize( 0, 1);

        reconstructor.Guardar("skel.png", IMGVTK::SKELETON, IMGVTK::PNG, 0);
        reconstructor.Guardar("boundaries.png", IMGVTK::BORDERS, IMGVTK::PNG, 0);

        reconstructor.skeletonize( 0, 4);
        reconstructor.skeletonize( 0, 20);
        reconstructor.skeletonize( 0, 100);
    }else if( strcmp(parametros[14].mi_valor.par_s, "NULL") && strcmp(parametros[15].mi_valor.par_s, "NULL") && strcmp(parametros[16].mi_valor.par_s, "NULL") ){

        FILE *fp_dataset = fopen(parametros[14].mi_valor.par_s, "r");

        char tmp;
        char tmp_str[512] = "";
        char *tmp_ptr = tmp_str;
        do{
            tmp = fgetc(fp_dataset);
            if( tmp == '\n' || tmp == EOF ){
                break;
            }
            *(tmp_ptr) = tmp;
            tmp_ptr++;
        }while( 1 );
        *(tmp_ptr) = '\0';
        tmp_ptr = tmp_str;
        const int n_imgs = atoi(tmp_str);

        char **rutas = new char* [ n_imgs ];
        for( int i = 0; i < n_imgs; i++){
            do{
                tmp = fgetc(fp_dataset);
                if(tmp == EOF || tmp == '\n'){
                    break;
                }
                *(tmp_ptr) = tmp;
                tmp_ptr++;
            }while( 1 );

            *(tmp_ptr) = '\0';
            tmp_ptr = tmp_str;

            rutas[i] = new char [(int)strlen(tmp_str)+1];
            sprintf(rutas[i], "%s", tmp_str);
            DEB_MSG("[" << i << "]: '" << rutas[i]);
        }

        fclose( fp_dataset );

        /// Leer el ground truth del dataset:
        fp_dataset = fopen( parametros[15].mi_valor.par_s, "r" );

        do{
            tmp = fgetc(fp_dataset);
            if( tmp == '\n' || tmp == EOF ){
                break;
            }
            *(tmp_ptr) = tmp;
            tmp_ptr++;
        }while( 1 );
        *(tmp_ptr) = '\0';
        tmp_ptr = tmp_str;

        char **rutas_gt = new char* [ n_imgs ];
        for( int i = 0; i < n_imgs; i++){
            do{
                tmp = fgetc(fp_dataset);
                if(tmp == EOF || tmp == '\n'){
                    break;
                }
                *(tmp_ptr) = tmp;
                tmp_ptr++;
            }while( 1 );

            *(tmp_ptr) = '\0';
            tmp_ptr = tmp_str;

            rutas_gt[i] = new char [(int)strlen(tmp_str)+1];
            sprintf(rutas_gt[i], "%s", tmp_str);
            DEB_MSG("[" << i << "]: '" << rutas_gt[i]);
        }

        fclose( fp_dataset );

        RECONS3D reconstructor;
        if( parametros[ 19 ].mi_valor.par_s[0] == 'y' ){
            if( strcmp( parametros[17].mi_valor.par_s, "NULL") ){
                reconstructor.setFiltroLog( parametros[17].mi_valor.par_s );
            }
            reconstructor.agregarInput( rutas, n_imgs, true );
            reconstructor.agregarGroundtruth( rutas_gt, n_imgs, 0);
            reconstructor.leerConfiguracion( parametros[16].mi_valor.par_s );
            reconstructor.segmentarImagenBase( 0 );
        }else{
            char ruta_log[512] = "";
            int tam_ruta_log = strlen(parametros[17].mi_valor.par_s);

            for(int i = 0; i < n_imgs; i++){
                if( strcmp( parametros[17].mi_valor.par_s, "NULL") ){
                    memcpy(ruta_log, parametros[17].mi_valor.par_s, (tam_ruta_log - 4)*sizeof(char));
                    sprintf(ruta_log, "%s_%i.log", ruta_log, i);
                    reconstructor.setFiltroLog( ruta_log );
                }

                reconstructor.agregarInput( rutas[i], true );
                reconstructor.agregarGroundtruth( rutas_gt[i], i);

                reconstructor.segmentarImagenBase( i );
            }
        }

        for( int i = 0; i < n_imgs; i++){
            delete [] rutas[i];
            delete [] rutas_gt[i];
        }
        delete [] rutas;
        delete [] rutas_gt;

    }else if( strcmp(parametros[14].mi_valor.par_s, "NULL") && strcmp(parametros[18].mi_valor.par_s, "NULL") ){

        FILE *fp_dataset = fopen(parametros[14].mi_valor.par_s, "r");

        int n_imgs;
        fscanf(fp_dataset, "%i" , &n_imgs);

        char tmp_str[512];
        char **rutas = new char* [ n_imgs ];
        for( int i = 0; i < n_imgs; i++){
            fscanf(fp_dataset, "%s", tmp_str);

            rutas[i] = new char [(int)strlen(tmp_str)+1];
            sprintf(rutas[i], "%s", tmp_str);
            DEB_MSG("[" << i << "]: '" << rutas[i]);
        }

        fclose( fp_dataset );

        RECONS3D reconstructor;

        if( parametros[19].mi_valor.par_s[0] == 'y' ){
            reconstructor.agregarInput( rutas, n_imgs, false );
            reconstructor.Guardar( parametros[18].mi_valor.par_s, IMGVTK::BASE, IMGVTK::PGM, 0);
        }else{
            fp_dataset = fopen(parametros[18].mi_valor.par_s, "r");
            fscanf(fp_dataset, "%i" , &n_imgs);
            char ruta_salida[512] = "";
            for( int i = 0; i < n_imgs; i++){
                fscanf(fp_dataset, "%s", ruta_salida);

                reconstructor.agregarInput( rutas[i], false );
                reconstructor.Guardar( ruta_salida, IMGVTK::BASE, IMGVTK::PGM, i);
            }
            fclose( fp_dataset );
        }

        for( int i = 0; i < n_imgs; i++){
            delete [] rutas[i];
        }
        delete [] rutas;
    }

    delete [] parametros;
    return EXIT_SUCCESS;
}
