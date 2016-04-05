
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


    // Parametro Distancia desde la fuente al detector (SID)
    parametros[3].mi_tipo = DOUBLE;
    sprintf(parametros[3].short_tag, "-i");
    sprintf(parametros[3].long_tag, "--sid");
    parametros[3].mi_default.par_d = 1000.0;
    sprintf(parametros[3].pregunta, "Distancia desde la fuente de Rayos-X al detector (en mm)");
    parametros[3].opcional = 0;


    // Parametro Distancia desde la fuente al paciente (SOD)
    parametros[4].mi_tipo = DOUBLE;
    sprintf(parametros[4].short_tag, "-o");
    sprintf(parametros[4].long_tag, "--sod");
    parametros[4].mi_default.par_d = 500.0;
    sprintf(parametros[4].pregunta, "Distancia desde la fuente de Rayos-X al paciente (en mm)");
    parametros[4].opcional = 0;


    // Parametro Angulo de deteccion CAU/CRA
    parametros[5].mi_tipo = DOUBLE;
    sprintf(parametros[5].short_tag, "-c");
    sprintf(parametros[5].long_tag, "--cra");
    parametros[5].mi_default.par_d = 0.0;
    sprintf(parametros[5].pregunta, "Angulo de deteccion CAU/CRA en direccion CRA (en grados)");
    parametros[5].opcional = 0;


    // Parametro Angulo de deteccion CAU/CRA
    parametros[6].mi_tipo = DOUBLE;
    sprintf(parametros[6].short_tag, "-l");
    sprintf(parametros[6].long_tag, "--lao");
    parametros[6].mi_default.par_d = 0.0;
    sprintf(parametros[6].pregunta, "Angulo de deteccion RAO/LAO en direccion LAO (en grados)");
    parametros[6].opcional = 0;
}


int main(int argc, char** argv ){
    // Definir los parametros de entrada:
    PARS_ENTRADA *parametros = new PARS_ENTRADA [7];
    definirParametros( parametros );
    if( argc < 2 ){
        mostrar_ayuda(parametros, 7, "Coronary");
        return EXIT_FAILURE;
    }
    // Revisar los parametros de entrada:
    revisar_pars(parametros, 7, &argc, argv);


    RECONS3D reconstructor( parametros[0].mi_valor.par_s, parametros[1].mi_valor.par_s, parametros[2].mi_valor.par_i );
    reconstructor.moverPosicion( 0, parametros[6].mi_valor.par_d, parametros[5].mi_valor.par_d, parametros[4].mi_valor.par_d, parametros[3].mi_valor.par_d);

    reconstructor.agregarInput(parametros[0].mi_valor.par_s, parametros[1].mi_valor.par_s, parametros[2].mi_valor.par_i );
    reconstructor.moverPosicion( 1, -parametros[6].mi_valor.par_d, -parametros[5].mi_valor.par_d, parametros[4].mi_valor.par_d, parametros[3].mi_valor.par_d);

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
    //reconstructor.segmentarImagenBase();

    return EXIT_SUCCESS;
}
