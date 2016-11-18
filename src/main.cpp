#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "args_fercer.c"
#include "IMGVTK.h"
#include "reconstructor_3D.h"

#include <iostream>

/*  Funcion: definirParametros
    Descripcion: Define los parametros de entrada del programa (Hard-coded).
*/
INPUT_ARGS *definirParametros(){
    // Parametro input base:
	int n_pars = 0;
    (parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-b1");
    sprintf((parametros + n_pars)->long_tag, "--base1");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "Imagen angiografica BASE 1 de entrada (.PNG, .BMP, .JPEG/.JPG, o archivo DICOM)");
	(parametros + n_pars)->opcional = 1;

    // Parametro nivel extraido del archivo DICOM:
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_INT;
    sprintf((parametros + n_pars)->short_tag, "-l1");
    sprintf((parametros + n_pars)->long_tag, "--desde1");
	(parametros + n_pars)->mi_default.par_i = 0;
    sprintf((parametros + n_pars)->pregunta, "Nivel desde el que se extraen las imagenes del archivo DICOM como imagenes base 1");
	(parametros + n_pars)->opcional = 1;

    // Parametro nivel extraido del archivo DICOM:
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_INT;
    sprintf((parametros + n_pars)->short_tag, "-u1");
    sprintf((parametros + n_pars)->long_tag, "--hasta1");
	(parametros + n_pars)->mi_default.par_i = 0;
    sprintf((parametros + n_pars)->pregunta, "Nivel hasta el que se extraen las iamgenes del archivo DICOM como imagenes base 1");
	(parametros + n_pars)->opcional = 1;


    // Parametro input base:
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-b2");
    sprintf((parametros + n_pars)->long_tag, "--base2");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "Imagen angiografica BASE 2 de entrada (.PNG, .BMP, .JPEG/.JPG, o archivo DICOM)");
	(parametros + n_pars)->opcional = 1;

    // Parametro nivel extraido del archivo DICOM:
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_INT;
    sprintf((parametros + n_pars)->short_tag, "-l2");
    sprintf((parametros + n_pars)->long_tag, "--desde2");
	(parametros + n_pars)->mi_default.par_i = 0;
    sprintf((parametros + n_pars)->pregunta, "Nivel desde el que se extraen las imagenes del archivo DICOM como imagenes base 2");
	(parametros + n_pars)->opcional = 1;

    // Parametro nivel extraido del archivo DICOM:
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_INT;
    sprintf((parametros + n_pars)->short_tag, "-u2");
    sprintf((parametros + n_pars)->long_tag, "--hasta2");
	(parametros + n_pars)->mi_default.par_i = 0;
    sprintf((parametros + n_pars)->pregunta, "Nivel hasta el que se extraen las iamgenes del archivo DICOM como imagenes base 2");
	(parametros + n_pars)->opcional = 1;


    // Parametro input base:
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-b3");
    sprintf((parametros + n_pars)->long_tag, "--base3");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "Imagen angiografica BASE 3 de entrada (.PNG, .BMP, .JPEG/.JPG, o archivo DICOM)");
	(parametros + n_pars)->opcional = 1;

    // Parametro nivel extraido del archivo DICOM:
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_INT;
    sprintf((parametros + n_pars)->short_tag, "-l3");
    sprintf((parametros + n_pars)->long_tag, "--desde3");
	(parametros + n_pars)->mi_default.par_i = 0;
    sprintf((parametros + n_pars)->pregunta, "Nivel desde el que se extraen las imagenes del archivo DICOM como imagenes base 3");
	(parametros + n_pars)->opcional = 1;

    // Parametro nivel extraido del archivo DICOM:
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_INT;
    sprintf((parametros + n_pars)->short_tag, "-u3");
    sprintf((parametros + n_pars)->long_tag, "--hasta3");
	(parametros + n_pars)->mi_default.par_i = 0;
    sprintf((parametros + n_pars)->pregunta, "Nivel hasta el que se extraen las iamgenes del archivo DICOM como imagenes base 3 ");
	(parametros + n_pars)->opcional = 1;




    // Parametro input ground:
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-g1");
    sprintf((parametros + n_pars)->long_tag, "--ground1");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "Imagen usada como ground truth (.PNG, .BMP, .JPEG/.JPG) para la imagen base 1");
	(parametros + n_pars)->opcional = 1;

    // Parametro input ground:
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-g2");
    sprintf((parametros + n_pars)->long_tag, "--ground2");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "Imagen usada como ground truth (.PNG, .BMP, .JPEG/.JPG) para la imagen base 2");
	(parametros + n_pars)->opcional = 1;

    // Parametro input ground:
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-g3");
    sprintf((parametros + n_pars)->long_tag, "--ground3");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "Imagen usada como ground truth (.PNG, .BMP, .JPEG/.JPG) para la imagen base 3");
	(parametros + n_pars)->opcional = 1;



    ////---------------------------------------------------------------------- PARA GENERAR IMAGENES PHANTOM EN ARCHIVO DICOM
    // Parametro output phantom
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-ph");
    sprintf((parametros + n_pars)->long_tag, "--phantom");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "[GENERAR ARCHIVO DICOM PARA UN PHANTOM] Ruta del arhcivo de la imagen phantom");
	(parametros + n_pars)->opcional = 1;

    // Parametro output phantom
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-op");
    sprintf((parametros + n_pars)->long_tag, "--output");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "[GENERAR ARCHIVO DICOM PARA UN PHANTOM] Ruta de salida del archivo DICOM generado en base al phantom");
	(parametros + n_pars)->opcional = 1;


    ////---------------------------------------------------------------------- PARA CARGAR VARIAS IMAGENES Y REALIZAR UN ENTRENAMIENTO DE LOS PARAMETROS
    // Parametro input data.set
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-dsb");
    sprintf((parametros + n_pars)->long_tag, "--datasetbase");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "[ENTRENAR PARAMETROS] Ruta del archivo data.set con la lista de rutas de imagenes a ser usadas como base");
	(parametros + n_pars)->opcional = 1;

    // Parametro input data.set
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-dsg");
    sprintf((parametros + n_pars)->long_tag, "--datasetground");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "[ENTRENAR PARAMETROS] Ruta del archivo data.set con la lista de rutas de imagenes a ser usadas como ground-truth");
	(parametros + n_pars)->opcional = 1;

    // Parametro input configuracion.dat
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-c");
    sprintf((parametros + n_pars)->long_tag, "--config");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "[ENTRENAR PARAMETROS] Ruta del archivo configuracion.dat con los parametros del filtro utilizado en la deteccion");
	(parametros + n_pars)->opcional = 1;

    // Parametro input filtrado.log
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-lf");
    sprintf((parametros + n_pars)->long_tag, "--log-filtrado");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "[ENTRENAR PARAMETROS] Ruta donde se guarda el log del proceso de filtrado");
	(parametros + n_pars)->opcional = 1;

    // Parametro output.pgm
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-pgm");
    sprintf((parametros + n_pars)->long_tag, "--rutaPGM");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "[EXPORTAR COMO PGM] Ruta donde se guarda la imagen PGM para la ejecucion del No-GUI (si se considera como )");
	(parametros + n_pars)->opcional = 1;
	
    // Parametro log de reconstruccion
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-lr");
    sprintf((parametros + n_pars)->long_tag, "--log-reconstruccion");
    sprintf((parametros + n_pars)->mi_default.par_s, "NULL");
    sprintf((parametros + n_pars)->pregunta, "[ENTRENAR PARAMETROS] Ruta donde se guarda el log del proceso de reconstruccion");
	(parametros + n_pars)->opcional = 1;

    // Parametro umbralizar
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-umb");
    sprintf((parametros + n_pars)->long_tag, "--umbralizar");
    sprintf((parametros + n_pars)->mi_default.par_s, "no");
    sprintf((parametros + n_pars)->pregunta, "[ENTRENAR PARAMETROS] Se umbraliza la respuesta");
	(parametros + n_pars)->opcional = 1;


    // Parametro input dir_threshold
	n_pars++;
	(parametros + n_pars)->mi_tipo = MI_CHAR;
    sprintf((parametros + n_pars)->short_tag, "-dt");
    sprintf((parametros + n_pars)->long_tag, "--directoriothreshold");
    sprintf((parametros + n_pars)->mi_default.par_s, ".");
    sprintf((parametros + n_pars)->pregunta, "[ENTRENAR PARAMETROS] Ruta donde se guardan las imagenes resultantes del proceso de umbralizado");
	(parametros + n_pars)->opcional = 1;
}



int main(int argc, char** argv ){
    // Definir los parametros de entrada:
    INPUT_ARGS *parametros = definirParametros();

    if( argc < 2 ){
#ifdef BUILD_VTK_VERSION
        printf("\n" COLOR_INVERSE "Version VTK" COLOR_RESET "\n");
#else
        printf("\n" COLOR_INVERSE "Version No-VTK" COLOR_RESET "\n");
#endif
        showHelp(parametros, "Coronary");
        delete [] parametros;
        return EXIT_FAILURE;
    }
    // Revisar los parametros de entrada:
    checkArgs(parametros, &argc, argv);

    // Si se va a generar un archivo DICOM para un phantom, no se genera el reconstructor 3D:
    if( strcmp( parametros[12].mi_valor.par_s , "NULL" ) && strcmp( parametros[13].mi_valor.par_s , "NULL" )){
        /// Generar archivo DICOM
        return 0;
    }/* Segmentar una unica imagen: */
	else if( strcmp(parametros[0].mi_valor.par_s, "NULL") && strcmp(parametros[16].mi_valor.par_s, "NULL") ){
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
        reconstructor.Guardar("ryc.png", IMGVTK::THRESHOLD, IMGVTK::PGM, 0);

        reconstructor.lengthFilter(IMGVTK::THRESHOLD, 600, 0);
        reconstructor.Guardar("ryc_600.png", IMGVTK::THRESHOLD, IMGVTK::PGM, 0);

        reconstructor.skeletonize(0, 100);
        reconstructor.skeletonize(0, 50);
        reconstructor.skeletonize(0, 20);
        reconstructor.skeletonize(0, 5);

        if( parametros[21].mi_valor.par_s[0] == 'y' ){

            double acc;

            reconstructor.umbralizar(IMGVTK::SEGMENT, IMGVTK::OTSU, 0, 0);
            reconstructor.Guardar( "otsu.pgm", IMGVTK::THRESHOLD, IMGVTK::PGM, 0 );
            acc = reconstructor.medirExactitud( 0 );
        }

    }/* Segmentar varias imagenes: */
	else if( strcmp(parametros[14].mi_valor.par_s, "NULL") &&
		strcmp(parametros[15].mi_valor.par_s, "NULL") &&
		strcmp(parametros[16].mi_valor.par_s, "NULL") )
	{
        FILE *fp_dataset = fopen(parametros[14].mi_valor.par_s, "r");
        char tmp_str[512] = "";

        int n_imgs;
        fscanf(fp_dataset, "%i", &n_imgs);

        char **rutas = new char* [ n_imgs ];
        for( int i = 0; i < n_imgs; i++){

            fscanf(fp_dataset, "%s", tmp_str);

            rutas[i] = new char [(int)strlen(tmp_str)+1];
            sprintf(rutas[i], "%s", tmp_str);
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
            if( strcmp( parametros[17].mi_valor.par_s, "NULL") ){
                reconstructor.setFiltroLog( parametros[17].mi_valor.par_s );
            }
            reconstructor.agregarInput( rutas, n_imgs, true );
            reconstructor.agregarGroundtruth( rutas_gt, n_imgs, 0);
            reconstructor.segmentarImagenBase( 0 );
            reconstructor.Guardar("filter_response.pgm", IMGVTK::SEGMENT, IMGVTK::PGM, 0 );
        }else{
            char ruta_log[512] = "";
            char ruta_logrec[512] = "";
            int tam_ruta_log = (int)strlen(parametros[17].mi_valor.par_s);
            int tam_ruta_logrec = (int)strlen(parametros[20].mi_valor.par_s);

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
