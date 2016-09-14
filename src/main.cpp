#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pars_fercer.c"
#include "IMGVTK.h"
#include "reconstructor_3D.h"

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
    sprintf(parametros[16].long_tag, "--config");
    sprintf(parametros[16].mi_default.par_s, "NULL");
    sprintf(parametros[16].pregunta, "[ENTRENAR PARAMETROS] Ruta del archivo configuracion.dat con los parametros del filtro utilizado en la deteccion");
    parametros[16].opcional = 1;

    // Parametro input filtrado.log
    parametros[17].mi_tipo = CHAR;
    sprintf(parametros[17].short_tag, "-lf");
    sprintf(parametros[17].long_tag, "--log-filtrado");
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

    // Parametro log de reconstruccion
    parametros[20].mi_tipo = CHAR;
    sprintf(parametros[20].short_tag, "-lr");
    sprintf(parametros[20].long_tag, "--log-reconstruccion");
    sprintf(parametros[20].mi_default.par_s, "NULL");
    sprintf(parametros[20].pregunta, "[ENTRENAR PARAMETROS] Ruta donde se guarda el log del proceso de reconstruccion");
    parametros[20].opcional = 1;

    // Parametro umbralizar
    parametros[21].mi_tipo = CHAR;
    sprintf(parametros[21].short_tag, "-umb");
    sprintf(parametros[21].long_tag, "--umbralizar");
    sprintf(parametros[21].mi_default.par_s, "no");
    sprintf(parametros[21].pregunta, "[ENTRENAR PARAMETROS] Se umbraliza la respuesta");
    parametros[21].opcional = 1;


    // Parametro input dir_threshold
    parametros[22].mi_tipo = CHAR;
    sprintf(parametros[22].short_tag, "-dt");
    sprintf(parametros[22].long_tag, "--directoriothreshold");
    sprintf(parametros[22].mi_default.par_s, ".");
    sprintf(parametros[22].pregunta, "[ENTRENAR PARAMETROS] Ruta donde se guardan las imagenes resultantes del proceso de umbralizado");
    parametros[22].opcional = 1;
}



int main(int argc, char** argv ){
    // Definir los parametros de entrada:
    PARS_ENTRADA *parametros = new PARS_ENTRADA [23];
    definirParametros( parametros );

    if( argc < 2 ){
#ifdef BUILD_VTK_VERSION
        printf("\n" COLOR_INVERSE "Version VTK" COLOR_RESET "\n");
#else
        printf("\n" COLOR_INVERSE "Version No-VTK" COLOR_RESET "\n");
#endif
        mostrar_ayuda(parametros, 23, "Coronary");
        delete [] parametros;
        return EXIT_FAILURE;
    }
    // Revisar los parametros de entrada:
    revisar_pars(parametros, 23, &argc, argv);

    // Si se va a generar un archivo DICOM par aun phantom, no se genera el reconstructor 3D:
    if( strcmp( parametros[12].mi_valor.par_s , "NULL" ) && strcmp( parametros[13].mi_valor.par_s , "NULL" )){
        /// Generar archivo DICOM
        return 0;
    }else if( strcmp(parametros[0].mi_valor.par_s, "NULL") && strcmp(parametros[16].mi_valor.par_s, "NULL") ){
        /// Reconstruir arteria:
        RECONS3D reconstructor;
        reconstructor.leerConfiguracion( parametros[16].mi_valor.par_s );

        reconstructor.agregarInput( parametros[0].mi_valor.par_s, true );
        reconstructor.agregarGroundtruth( parametros[9].mi_valor.par_s, 0);

        if( strcmp( parametros[17].mi_valor.par_s, "NULL") ){
            reconstructor.setFiltroLog( parametros[17].mi_valor.par_s );
        }

        if( strcmp( parametros[20].mi_valor.par_s, "NULL") ){
            reconstructor.setLog( parametros[20].mi_valor.par_s );
        }

        reconstructor.segmentarImagenBase( 0 );
        reconstructor.Guardar("resp.pgm", IMGVTK::SEGMENT, IMGVTK::PGM, 0);

        reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
        reconstructor.Guardar("ryc.png", IMGVTK::THRESHOLD, IMGVTK::PNG, 0);

        reconstructor.lengthFilter(IMGVTK::THRESHOLD, 600, 0);
        reconstructor.Guardar("ryc_600.png", IMGVTK::THRESHOLD, IMGVTK::PNG, 0);

        reconstructor.skeletonize(0, 100);
        reconstructor.skeletonize(0, 50);
        reconstructor.skeletonize(0, 20);
        reconstructor.skeletonize(0, 5);

        if( parametros[21].mi_valor.par_s[0] == 'y' ){

            double acc;

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.Guardar( "otsu.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );
            acc = reconstructor.medirExactitud( 0 );

            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 100, 0);
            acc = reconstructor.medirExactitud( 0 );
            reconstructor.Guardar( "otsu_0100.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 200, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_0200.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 300, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_0300.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 400, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_0400.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 500, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_0500.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 600, 0);
            acc = reconstructor.medirExactitud( 0 );
            //reconstructor.Guardar( "otsu_0600.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 700, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_0700.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 800, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_0800.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 900, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_0900.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1000, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_1000.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1100, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_1100.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1200, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_1200.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1300, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_1300.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1400, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_1400.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1500, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_1500.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1600, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_1600.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1700, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_1700.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1800, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_1800.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1900, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_1900.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 2000, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "otsu_2000.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
//            reconstructor.Guardar( "ryc.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );
            acc = reconstructor.medirExactitud( 0 );

            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 100, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_0100.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 200, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_0200.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 300, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_0300.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 400, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_0400.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 500, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_0500.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 600, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_0600.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 700, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_0700.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 800, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_0800.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 900, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_0900.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1000, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_1000.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1100, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_1100.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1200, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_1200.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1300, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_1300.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1400, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_1400.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1500, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_1500.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1600, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_1600.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1700, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_1700.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1800, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_1800.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1900, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_1900.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0, 0);
            reconstructor.lengthFilter(IMGVTK::THRESHOLD, 2000, 0);
            acc = reconstructor.medirExactitud( 0 );
//            reconstructor.Guardar( "ryc_2000.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );
        }

        //reconstructor.skeletonize( 0 );

    }else if( strcmp(parametros[14].mi_valor.par_s, "NULL") && strcmp(parametros[15].mi_valor.par_s, "NULL") && strcmp(parametros[16].mi_valor.par_s, "NULL") ){

        FILE *fp_dataset = fopen(parametros[14].mi_valor.par_s, "r");
        char tmp_str[512] = "";

        int n_imgs;
        fscanf(fp_dataset, "%i", &n_imgs);

        char **rutas = new char* [ n_imgs ];
        for( int i = 0; i < n_imgs; i++){

            fscanf(fp_dataset, "%s", tmp_str);

            rutas[i] = new char [(int)strlen(tmp_str)+1];
            sprintf(rutas[i], "%s", tmp_str);
            DEB_MSG("[" << i << "]: '" << rutas[i]);
        }

        fclose( fp_dataset );

        /// Leer el ground truth del dataset:
        fp_dataset = fopen( parametros[15].mi_valor.par_s, "r" );

        char **rutas_gt = new char* [ n_imgs ];

        fscanf( fp_dataset, "%i", &n_imgs);

        for( int i = 0; i < n_imgs; i++){

            fscanf(fp_dataset, "%s", tmp_str);

            rutas_gt[i] = new char [(int)strlen(tmp_str)+1];
            sprintf(rutas_gt[i], "%s", tmp_str);
            DEB_MSG("[" << i << "]: '" << rutas_gt[i]);
        }

        fclose( fp_dataset );

        RECONS3D reconstructor;
        reconstructor.leerConfiguracion( parametros[16].mi_valor.par_s );
        if( parametros[ 19 ].mi_valor.par_s[0] == 'y' ){
            DEB_MSG("Concatenando");
            if( strcmp( parametros[17].mi_valor.par_s, "NULL") ){
                reconstructor.setFiltroLog( parametros[17].mi_valor.par_s );
            }
            reconstructor.agregarInput( rutas, n_imgs, true );
            reconstructor.agregarGroundtruth( rutas_gt, n_imgs, 0);
            reconstructor.segmentarImagenBase( 0 );
            reconstructor.Guardar("filter_response.pgm", IMGVTK::SEGMENT, IMGVTK::PGM, 0 );
        }else{
            DEB_MSG("No Concatenado");
            char ruta_log[512] = "";
            char ruta_logrec[512] = "";
            int tam_ruta_log = strlen(parametros[17].mi_valor.par_s);
            int tam_ruta_logrec = strlen(parametros[20].mi_valor.par_s);

            for(int i = 0; i < n_imgs; i++){
                if( strcmp(parametros[17].mi_valor.par_s, "NULL") ){
                    memset(ruta_log, 0, 512*sizeof(char));
                    memcpy(ruta_log, parametros[17].mi_valor.par_s, (tam_ruta_log - 4)*sizeof(char));
                    sprintf(ruta_log, "%s_%i.log", ruta_log, i);
                    reconstructor.setFiltroLog( ruta_log );
                }

                if( strcmp(parametros[20].mi_valor.par_s, "NULL") ){
                    memset(ruta_logrec, 0, 512*sizeof(char));
                    memcpy(ruta_logrec, parametros[20].mi_valor.par_s, (tam_ruta_logrec - 4)*sizeof(char));
                    sprintf(ruta_logrec, "%s_%i.log", ruta_logrec, i);
                    reconstructor.setLog( ruta_logrec );
                }

                reconstructor.agregarInput( rutas[i], true );
                reconstructor.agregarGroundtruth( rutas_gt[i], i);
                reconstructor.segmentarImagenBase( i );



                if( parametros[21].mi_valor.par_s[0] != 'n' ){
                    char ruta_thresholds[512];
                    double accuracy;

                    sprintf(ruta_thresholds, "%s/tout_otsu_100_%i.pgm", parametros[22].mi_valor.par_s, i );
                    reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0.0, i);
                    reconstructor.lengthFilter(IMGVTK::THRESHOLD, 100, i);
                    reconstructor.Guardar(ruta_thresholds, IMGVTK::THRESHOLD, IMGVTK::PGM, i);
                    accuracy = reconstructor.medirExactitud(i);

                    sprintf(ruta_thresholds, "%s/tout_otsu_500_%i.pgm", parametros[22].mi_valor.par_s, i );
                    reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0.0, i);
                    reconstructor.lengthFilter(IMGVTK::THRESHOLD, 500, i);
                    reconstructor.Guardar(ruta_thresholds, IMGVTK::THRESHOLD, IMGVTK::PGM, i);
                    accuracy = reconstructor.medirExactitud(i);

                    sprintf(ruta_thresholds, "%s/tout_otsu_1000_%i.pgm", parametros[22].mi_valor.par_s, i );
                    reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0.0, i);
                    reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1000, i);
                    reconstructor.Guardar(ruta_thresholds, IMGVTK::THRESHOLD, IMGVTK::PGM, i);
                    accuracy = reconstructor.medirExactitud(i);

                    sprintf(ruta_thresholds, "%s/tout_otsu_2000_%i.pgm", parametros[22].mi_valor.par_s, i );
                    reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0.0, i);
                    reconstructor.lengthFilter(IMGVTK::THRESHOLD, 2000, i);
                    reconstructor.Guardar(ruta_thresholds, IMGVTK::THRESHOLD, IMGVTK::PGM, i);
                    accuracy = reconstructor.medirExactitud(i);

                    sprintf(ruta_thresholds, "%s/tout_rc_100_%i.pgm", parametros[22].mi_valor.par_s, i );
                    reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0.0, i);
                    reconstructor.lengthFilter(IMGVTK::THRESHOLD, 100, i);
                    reconstructor.Guardar(ruta_thresholds, IMGVTK::THRESHOLD, IMGVTK::PGM, i);
                    accuracy = reconstructor.medirExactitud(i);

                    sprintf(ruta_thresholds, "%s/tout_rc_500_%i.pgm", parametros[22].mi_valor.par_s, i );
                    reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0.0, i);
                    reconstructor.lengthFilter(IMGVTK::THRESHOLD, 500, i);
                    reconstructor.Guardar(ruta_thresholds, IMGVTK::THRESHOLD, IMGVTK::PGM, i);
                    accuracy = reconstructor.medirExactitud(i);

                    sprintf(ruta_thresholds, "%s/tout_rc_1000_%i.pgm", parametros[22].mi_valor.par_s, i );
                    reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0.0, i);
                    reconstructor.lengthFilter(IMGVTK::THRESHOLD, 1000, i);
                    reconstructor.Guardar(ruta_thresholds, IMGVTK::THRESHOLD, IMGVTK::PGM, i);
                    accuracy = reconstructor.medirExactitud(i);

                    sprintf(ruta_thresholds, "%s/tout_rc_2000_%i.pgm", parametros[22].mi_valor.par_s, i );
                    reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::RIDLER_CALVARD, 0.0, i);
                    reconstructor.lengthFilter(IMGVTK::THRESHOLD, 2000, i);
                    reconstructor.Guardar(ruta_thresholds, IMGVTK::THRESHOLD, IMGVTK::PGM, i);
                    accuracy = reconstructor.medirExactitud(i);

                }
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
    return 0;
}
