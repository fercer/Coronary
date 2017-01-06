/******
    CENTRO DE INVESTIGACION EN MATEMATICAS
    MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

    FERNANDO CERVANTES SANCHEZ
    NOV - 2015
*****/

#ifndef FILTROS_H_INCLUDED
#define FILTROS_H_INCLUDED



#include <assert.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#if defined(_WIN32)
    #include MY_FFTW_PATH
#else
    #include <fftw3.h>
#endif

#include <math.h>
#include <omp.h>

#include <iostream>
#include <vector>

#ifdef BUILD_GUI_VERSION
    #include <QPlainTextEdit>
    #include <QProgressBar>
#endif

//#include <vtkSmartPointer.h>
//#include <vtkImageData.h>

#include "args_fercer.h"
#include "reconstructor_3D.h"


#ifndef MY_PI
#define Mi_PI 3.1415926535897932384626433832795
#endif



// C L A S E: FILTRO  ---------------------------------------------------------------------------------------------- v
class FILTROS {

public: //----------------------------------------------------------------------------- PUBLIC ------- v
	typedef enum SEG_FILTER { SEG_UNSET, GMF, SS_GABOR } SEG_FILTER;

	FILTROS();
	~FILTROS();

	void setFilterMethod(const SEG_FILTER new_filtering_method);

	void filtrar();

	void setLog(FILE *fplog);
	void setLog(const char *ruta_log);

#ifdef BUILD_GUI_VERSION
	void setLog(QTextEdit *txtLog);
	void setProgressBar(QProgressBar *pBar);
#endif

	void setInputBase(std::vector<IMGCONT>* new_img_base);
	void setInputBaseMask(std::vector<IMGCONT>* new_img_base_mask);
	void setInputGroundtruth(std::vector<IMGCONT>* new_img_groundtruth);
	void setInputResponse(std::vector<IMGCONT>* new_img_response);

	void setParL(const double new_L_value);
	void setParT(const unsigned int new_T_value);
	void setParK(const unsigned int new_K_value);
	void setParSigma(const double new_sigma_value);

private:
	bool input_already_set;
	bool already_transformed;

	fftw_complex **Img_fft;
	fftw_complex **Img_fft_HPF;

	void fftImgOrigen();

	ARGUMENTS *my_args;
	SEG_FILTER chosen_filter;

	double my_L;
	unsigned int my_T;
	unsigned int my_K;
	double my_sigma;

	std::vector< IMGCONT > * my_img_response;
	std::vector< IMGCONT > * my_img_base;
	std::vector< IMGCONT > * my_img_groundtruth;
	std::vector< IMGCONT > * my_img_base_mask;

	unsigned int my_filters_imgs_count;

	inline double interpolacion(const double *pix, const int j, const int i, const double x, const double y, const int mis_rens, const int mis_cols);

	void rotateImg(const double *org, double *rot, const double ctheta, const double stheta, const int mis_rens, const int mis_cols, const int org_rens, const int org_cols);

	void respGMF();
	void respGabor();

	void escribirLog(const char *mensaje);
	void barraProgreso(const int avance, const int max_progress);

	char *mi_ruta_log;
	FILE *mi_fplog;

#ifdef BUILD_GUI_VERSION
	QTextEdit *mi_txtLog;
	QProgressBar *mi_pBar;
#endif

};
// C L A S E: FILTROS  ----------------------------------------------------------------------------------------- ^

#endif //FILTROS_H_INCLUDED
