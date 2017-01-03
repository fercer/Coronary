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
		printf("GDCM-Version" COLOR_RESET "\n");

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
    }/* Segmentar una unica imagen: */

	DEB_MSG("Input: " << parametros.getArgumentCHAR("-b"));
	DEB_MSG("Groundtruth: " << parametros.getArgumentCHAR("-g"));
	DEB_MSG("Data base dir: " << parametros.getArgumentCHAR("-dsb"));
	DEB_MSG("Data Groundtruth dir: " << parametros.getArgumentCHAR("-dsg"));

	if ((strcmp(parametros.getArgumentCHAR("-b"), "NULL") != 0) &&
		(strcmp(parametros.getArgumentCHAR("-g"), "NULL") != 0)) {
		/// Reconstruir arteria:
		mi_reconstructor.agregarInput(parametros.getArgumentCHAR("-b"),parametros.getArgumentINT("-bl"), parametros.getArgumentINT("-bu"), parametros.getArgumentCHAR("-g"), false);

		DEB_MSG("Saving image...");
		mi_reconstructor.Guardar("img.pgm", RECONS3D::IMG_BASE, IMGCONT::IMGPGM, 0);
	}

	if ((strcmp(parametros.getArgumentCHAR("-dsb"), "NULL") != 0) &&
		(strcmp(parametros.getArgumentCHAR("-dsg"), "NULL") != 0)) {

		FILE *fp_base = fopen(parametros.getArgumentCHAR("-dsb"), "r");
		FILE *fp_gt   = fopen(parametros.getArgumentCHAR("-dsg"), "r");

		int n_imgs;

		fscanf(fp_base, "%i\n", &n_imgs);
		fscanf(fp_gt, "%i\n", &n_imgs);

		DEB_MSG("n imgs in directory: " << n_imgs);

		char *ruta_base = new char[512];
		char *ruta_gt   = new char[512];
		int ruta_len;
		char *resp;
		for (int i = 0; i < n_imgs; i++) {
			resp = fgets(ruta_base, 512, fp_base);
			resp = fgets(ruta_gt, 512, fp_gt);

			ruta_len = (int)strlen(ruta_base);
			if ((int)*(ruta_base + ruta_len - 1) == 10) {
				*(ruta_base + ruta_len - 1) = '\0';
			}

			ruta_len = (int)strlen(ruta_gt);
			if ((int)*(ruta_gt + ruta_len - 1) == 10) {
				*(ruta_gt + ruta_len - 1) = '\0';
			}

			DEB_MSG("[" << i << "] base: " << ruta_base <<", gt: " << ruta_gt);
			mi_reconstructor.agregarInput(ruta_base, 0, 0, ruta_gt, true);

		}

		fclose(fp_base);
		fclose(fp_gt);
	}

	getchar();

    return EXIT_SUCCESS;
}
