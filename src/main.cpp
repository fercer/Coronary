
#include <stdio.h>
#include <stdlib.h>

#include "reconstructor_3D.h"


#include <iostream>

int main(int argc, char** argv ){

    if( argc < 4 ){
        printf("usage: Segmentador <ruta_imagen_base> <ruta_ground_truth> <DICOM iamge level>\n");
        return EXIT_FAILURE;
    }
    RECONS3D reconstructor( argv[1] , argv[2], atoi(argv[3]) );

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
    reconstructor.segmentarImagenBase();

    return EXIT_SUCCESS;
}
