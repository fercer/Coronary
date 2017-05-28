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

	
	if ((strcmp(parametros.getArgumentCHAR("-b"), "NULL") != 0) &&
		(strcmp(parametros.getArgumentCHAR("-g"), "NULL") != 0)) {

		DEB_MSG("Input: " << parametros.getArgumentCHAR("-b"));
		DEB_MSG("Groundtruth: " << parametros.getArgumentCHAR("-g"));

		mi_reconstructor.agregarInput(parametros.getArgumentCHAR("-b"), parametros.getArgumentINT("-bl"), parametros.getArgumentINT("-bu"), parametros.getArgumentCHAR("-g"), true);

		mi_reconstructor.setFiltroMetodo(FILTROS::SS_GABOR);

		mi_reconstructor.setFiltroParametros(OPTI_PARS::PAR_K, 180.0);
		mi_reconstructor.setFiltroParametros(OPTI_PARS::PAR_T, 12.0);
		mi_reconstructor.setFiltroParametros(OPTI_PARS::PAR_L, 2.5);

		mi_reconstructor.segmentar();

		mi_reconstructor.Guardar("img_resp.png", RECONS3D::IMG_RESPONSE, IMGCONT::IMGPNG, 0);
	}


	if ((strcmp(parametros.getArgumentCHAR("-b"), "NULL") != 0) &&
		(strcmp(parametros.getArgumentCHAR("-g"), "NULL") != 0) &&
		(strcmp(parametros.getArgumentCHAR("-c"), "NULL") != 0)) {
		/// Reconstruir arteria:
		mi_reconstructor.agregarInput(parametros.getArgumentCHAR("-b"),parametros.getArgumentINT("-bl"), parametros.getArgumentINT("-bu"), parametros.getArgumentCHAR("-g"), true);
		
#ifndef NDEBUG
		mi_reconstructor.Guardar("img_base.pgm", RECONS3D::IMG_BASE, IMGCONT::IMGPGM, 0);
		mi_reconstructor.Guardar("img_mask.pgm", RECONS3D::IMG_MASK, IMGCONT::IMGPGM, 0);
#endif // !NDEBUG

		if (strcmp(parametros.getArgumentCHAR("-lf"), "NULL") != 0) {
			mi_reconstructor.setFiltroLog(parametros.getArgumentCHAR("-lf"));
		}
		mi_reconstructor.leerConfiguracion(parametros.getArgumentCHAR("-c"));
		mi_reconstructor.segmentar();
		mi_reconstructor.Guardar("img_resp.pgm", RECONS3D::IMG_RESPONSE, IMGCONT::IMGPGM, 0);
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

		for (int i = 0; i < n_imgs; i++) {
			fscanf(fp_base, "%s", ruta_base);
			fscanf(fp_gt, "%s", ruta_gt);
			
			DEB_MSG("[" << i << "] base: " << ruta_base <<", gt: " << ruta_gt);
			mi_reconstructor.agregarInput(ruta_base, 0, 0, ruta_gt, true);
		}

		if (strcmp(parametros.getArgumentCHAR("-lf"), "NULL") != 0) {
			mi_reconstructor.setFiltroLog(parametros.getArgumentCHAR("-lf"));
		}

		if (strcmp(parametros.getArgumentCHAR("-lr"), "NULL") != 0) {
			mi_reconstructor.setLog(parametros.getArgumentCHAR("-lr"));
		}

		mi_reconstructor.leerConfiguracion(parametros.getArgumentCHAR("-c"));
		mi_reconstructor.segmentar();

		char out_path[512] = "000_res.pgm";
		for (int i = 0; i < n_imgs; i++) {
#if defined(_WIN32) || defined(_WIN64)
			sprintf_s(out_path, 512, "%s/%i_res.pgm", parametros.getArgumentCHAR("-odir"), i);
#else
	sprintf(out_path, "%s/%i_res.pgm", parametros.getArgumentCHAR("-odir"), i);
#endif
			printf("\n<<%s>>\n", out_path);	
			mi_reconstructor.Guardar(out_path, RECONS3D::IMG_RESPONSE, IMGCONT::IMGPGM, i);
		}

		fclose(fp_base);
		fclose(fp_gt);

		delete[] ruta_base;
		delete[] ruta_gt;
	}
    return EXIT_SUCCESS;
}
