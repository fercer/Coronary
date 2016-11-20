#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "args_fercer.h"
#include "reconstructor_3D.h"

#include <iostream>


int main(int argc, char** argv ){

    // Definir los parametros de entrada:
    ARGUMENTS parametros( argc, argv);
	RECONS3D mi_reconstructor(&parametros);

    if( argc < 2 || parametros.getArgumentINT("-h") ){
		printf("\n" COLOR_INVERSE);
#ifndef BUILD_VTK_VERSION
		printf("No-");
#endif
		printf("VTK/");
#ifndef BUILD_GUI_VERSION
		printf("No-");
#endif
		printf("GUI/");
#ifndef BUILD_GCM_VERSION
		printf("No-");
#endif
		printf("GDCM-");
		printf("Version" COLOR_RESET "\n");

		parametros.setArgumentINT("-h", 1);
		parametros.showHelp();
        return EXIT_FAILURE;
    }

	/*
    // Si se va a generar un archivo DICOM para un phantom, no se genera el reconstructor 3D:
    if( strcmp( parametros.getArgumentCHAR("-ph") , "NULL" ) &&
		strcmp( parametros.getArgumentCHAR("-oph") , "NULL" )){
        /// Generar archivo DICOM
        return EXIT_SUCCESS;
    }/* Segmentar una unica imagen: *
	else if (strcmp(parametros.getArgumentCHAR("-b"), "NULL") &&
		strcmp(parametros.getArgumentCHAR("g"), "NULL")) {
		/// Reconstruir arteria:
	}
	*/

    return EXIT_SUCCESS;
}
