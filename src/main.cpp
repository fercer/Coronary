
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pars_fercer.c"
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
    /*
    *********************************************************************************************
                BUSQUEDA EXHAUSTIVA DE LOS PARAMETROS:
    *********************************************************************************************
    */
/*
    char *rutas_base[20];
    char *rutas_grd[20];
    char ruta_dir[] = "/home/fercer/test_data/40Gabor/ang_";

    for(int i = 0; i < 20; i++){
        rutas_base[i] = new char[256];
        sprintf(rutas_base[i], "%s%i.png", ruta_dir, i+1);
        rutas_grd[i] = new char[256];
        sprintf(rutas_grd[i], "%s%i_gt.png", ruta_dir, i+1);;
    }

    RECONS3D reconstructor( rutas_base, rutas_grd, 20 );

    for(int i = 0; i < 20; i++){
        delete [] rutas_base[i];
        delete [] rutas_grd[i];
    }
*/

    reconstructor.segmentarImagenBase();
//    reconstructor.skeletonize( 0 );

    delete [] parametros;

    return EXIT_SUCCESS;
}
