
#include <stdio.h>
#include <stdlib.h>

#include "pars_fercer.c"
#include "reconstructor_3D.h"

#include <iostream>




/*  Funcion: definirParametros
    Descripcion: Define los parametros de entrada del programa (Hard-coded).
*/
void definirParametros(PARS_ENTRADA *parametros){

    // Parametro input base:
    parametros[0].mi_tipo = CHAR;
    sprintf(parametros[0].short_tag, "-b");
    sprintf(parametros[0].long_tag, "--base");
    sprintf(parametros[0].mi_default.par_s, "ang_1.png");
    sprintf(parametros[0].pregunta, "Imagen angiografica BASE de entrada (.PNG, .BMP, .JPEG/.JPG, o archivo DICOM)");
    parametros[0].opcional = 0;

    // Parametro input delin:
    parametros[1].mi_tipo = CHAR;
    sprintf(parametros[1].short_tag, "-d");
    sprintf(parametros[1].long_tag, "--delineada");
    sprintf(parametros[1].mi_default.par_s, "ang_1_gt.png");
    sprintf(parametros[1].pregunta, "Imagen angiografica DELINEADA de entrada (.PNG, .BMP, .JPEG/.JPG)");
    parametros[1].opcional = 0;

    // Parametro nivel extraido del archivo DICOM:
    parametros[2].mi_tipo = INT;
    sprintf(parametros[2].short_tag, "-n");
    sprintf(parametros[2].long_tag, "--nivel");
    parametros[2].mi_default.par_i = 0;
    sprintf(parametros[2].pregunta, "Nivel que se extrae del archivo DICOM como imagen base");
    parametros[2].opcional = 1;
}


int main(int argc, char** argv ){
    // Definir los parametros de entrada:
    PARS_ENTRADA *parametros = new PARS_ENTRADA [3];
    definirParametros( parametros );

    if( argc < 2 ){
        mostrar_ayuda(parametros, 3, "Coronary");
        delete [] parametros;
        return EXIT_FAILURE;
    }
    // Revisar los parametros de entrada:
    revisar_pars(parametros, 3, &argc, argv);

DEB_MSG("Parametros revisados ...");
    RECONS3D reconstructor;
    reconstructor.agregarInput(parametros[0].mi_valor.par_s, parametros[1].mi_valor.par_s, parametros[2].mi_valor.par_i);
    reconstructor.agregarInput(parametros[0].mi_valor.par_s, parametros[1].mi_valor.par_s, parametros[2].mi_valor.par_i );

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

    //reconstructor.skeletonize();
    reconstructor.segmentarImagenBase( 0 );

    delete [] parametros;

    return EXIT_SUCCESS;
}
