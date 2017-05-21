/******
    CENTRO DE INVESTIGACION EN MATEMATICAS
    MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

    FERNANDO CERVANTES SANCHEZ
   NOV - 2015
*****/


#include "filtros.h"

/*  Metodo: escribirLog

    Funcion: Escribe un mensaje en el log.
*/
void FILTROS::escribirLog( const char *mensaje )
{
    if(mi_fplog){
        fprintf(mi_fplog, "%s", mensaje);
#ifdef BUILD_GUI_VERSION
    }else if(mi_txtLog){
        mi_txtLog->append( mensaje );
#endif
    }else{
        std::cout << mensaje << std::endl;
        fflush(stdout);
    }
}




/*  Metodo: barraProgreso

    Funcion: Actualiza la barra de progreso.
*/
void FILTROS::barraProgreso( const int avance, const int max_progress ){

#ifdef BUILD_GUI_VERSION
    if( mi_pBar ){
        mi_pBar->setMaximum( max_progress );
        mi_pBar->setValue( avance );
        return;
    }
#endif

    /// Limpiar el resto de la linea:
    int max_ancho = 100;
    for( int i = 0; i <= max_ancho; i++){
        printf("\r");
    }
    printf(COLOR_BACK_RED "[");
    int avance_progreso = (int)((double)max_ancho * (double)avance / (double)max_progress);
    for( int i = 0; i < avance_progreso; i++){
#if defined(_WIN32) || defined(_WIN64)
		printf("=");
#else
        printf(COLOR_BACK_GREEN " ");
#endif
    }
    for( int i = avance_progreso; i < max_ancho; i++){
        printf(COLOR_BACK_CYAN " ");
    }
    printf(COLOR_BACK_RED "]" COLOR_NORMAL);
    fflush(stdout);
}







/************************************************************************************************************
* FILTROS::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: setFilterMethod                                                                            *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_filtering_method      const SEG_FILTER           I   The filtering method defined.                    *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void FILTROS::setFilterMethod( const SEG_FILTER new_filtering_method)
{
    chosen_filter = new_filtering_method;
}




/*  Metodo: FILTROS

    Funcion: Inicializa los parametros en valores por defecto.
*/
FILTROS::FILTROS(){
#ifdef BUILD_GUI_VERSION
    mi_pBar = NULL;
    mi_txtLog = NULL;
#endif

    mi_ruta_log = NULL;
    mi_fplog = NULL;

	input_already_set = false;
	already_transformed = false;

	Img_fft = NULL;
	Img_fft_HPF = NULL;
	
    chosen_filter = SEG_UNSET;

	my_filters_imgs_count = 0;
}




/*  Metodo: ~FILTROS

    Funcion: Libera la memoria requerida por el objeto.
*/
FILTROS::~FILTROS() {
	if (mi_fplog) {
		fclose(mi_fplog);
		mi_fplog = NULL;
	}

	if (mi_ruta_log) {
		delete[] mi_ruta_log;
		mi_ruta_log = NULL;
	}

	if (already_transformed) {

		for (unsigned int i = 0; i < my_filters_imgs_count; i++) {
			fftwl_free(Img_fft[i]);
			fftwl_free(Img_fft_HPF[i]);
		}

		delete[] Img_fft;
		delete[] Img_fft_HPF;
	}
}




/*  Metodo: filtrar
    Funcion: Aplica el filtro sobre la imagen de origen y la almacena sobre la imagen destino.
*/
void FILTROS::filter()
{
	switch (chosen_filter) {
	case GMF:
		respGMF();
		break;
	case SS_GABOR:
		respGabor();
		break;
	}
}




/*  Metodo: setLog
    Funcion: Define el objeto tipo QPlainTextEdit donde se imprimen los logs del sistema.
*/
#ifdef BUILD_GUI_VERSION
void FILTROS::setLog(QTextEdit *txtLog)
{
    mi_txtLog = txtLog;
}
#endif


/*  Metodo: setLog
    Funcion: Define el stream donde se imprimen los logs del sistema.
*/
void FILTROS::setLog(FILE *fplog){
    mi_fplog = fplog;
}



/*  Metodo: setLog
    Funcion: Define el stream donde se imprimen los logs del sistema.
*/
void FILTROS::setLog( const char *ruta_log){
    mi_ruta_log = new char [512];
#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(mi_ruta_log, 512 * sizeof(char), "%s", ruta_log);
#else
    sprintf(mi_ruta_log, "%s", ruta_log);
#endif
    if(mi_fplog){
        fclose( mi_fplog );
    }
#if defined(_WIN32) || defined(_WIN64)
	fopen_s(&mi_fplog, mi_ruta_log, "a");
#else
    mi_fplog = fopen( mi_ruta_log, "a" );
#endif
}





/*  Metodo: setProgressBar
    Funcion: Define el objeto para mostrar el avance del proceso.
*/
#ifdef BUILD_GUI_VERSION
void FILTROS::setProgressBar( QProgressBar *pBar ){
    mi_pBar = pBar;
}
#endif





/*  respGMF

    Function: Applies the Gaussian Matched Filter to the input image using the parameters defined in the
    INDIV structure 'test', and saves the output in a pointer 'resp'.

    ******** Code used for the computational experiments presented in the chapter:

    "Automatic Detection of Coronary Artery Stenosis using Bayesian Classification and
            Gaussian Filters based on Differential Evolution"

    for the book: "Hybrid Intelligence for Image Analysis and Understanding" published by John Wiley, UK.

    Authors: Ivan Cruz-Aceves, Fernando Cervantes-Sanchez and Arturo Hernandez-Aguirre

    *************************************
*/
void FILTROS::respGMF()
{
	const unsigned int temp_dims = (int)(1.5*((my_T > (unsigned int)my_L) ? my_T : (unsigned int)my_L));

	IMGCONT*GMF_templates = new IMGCONT[my_K];

	/* Calculate a single Gaussian curve along the template width, at 0 degrees: */
	double *gauss_0 = (double*)malloc(my_T * sizeof(double));
	double *gauss_ptr = gauss_0;
	const double sig_2 = 2.0 * my_sigma * my_sigma;
	double sum = 0.0;
	for (double x = -floor((double)my_T / 2.0) + (double)(1 - my_T % 2) / 2.0;
		x <= floor((double)my_T / 2.0) - (double)(1 - my_T % 2) / 2.0;
		x += 1.0)
	{
		*(gauss_ptr) = 1.0 - exp(-(x*x / sig_2));
		sum += *(gauss_ptr);
		gauss_ptr++;
	}

	sum *= (unsigned int)my_L;
	gauss_ptr = gauss_0;
	const double mean = sum / ((double)my_T * (double)(unsigned int)my_L);

	for (double x = -floor((double)my_T / 2.0) + (double)(1 - my_T % 2) / 2.0;
		x <= floor((double)my_T / 2.0) - (double)(1 - my_T % 2) / 2.0;
		x += 1.0)
	{
		*(gauss_ptr) = (*(gauss_ptr)-mean) / sum;
		gauss_ptr++;
	}

	/* Copy the Gaussian curve along the template height: */
	(GMF_templates)->setDimensions(temp_dims, temp_dims);

	for (unsigned int y = 0; y < (unsigned int)my_L; y++)
	{
		for (unsigned int x = 0; x < my_T; x++) {
			(GMF_templates)->setPix( (temp_dims / 2 - (unsigned int)my_L / 2 + y), temp_dims / 2 - my_T / 2 + x, *(gauss_0 + x));
		}
	}

	/* Rotate template at each 'my_K' orientation: */
	const double theta_inc = 180.0 / (double)my_K;
	double theta = 0.0;

	for (unsigned int k = 1; k < my_K; k++)
	{
		theta += theta_inc;
		const double ctheta = cos(-theta * MY_PI / 180.0);
		const double stheta = sin(-theta * MY_PI / 180.0);

		*(GMF_templates + k) = *(GMF_templates);

		(GMF_templates + k)->Rotate(theta);
	}

	free(gauss_0);
	
	/* Initialize the filter response: */
	IMGCONT*resp_tmp = new IMGCONT[my_filters_imgs_count];

	const unsigned int offset = (unsigned int)(temp_dims / 2);
	for (unsigned int i = 0; i < my_filters_imgs_count; i++)
	{
		(resp_tmp + i)->setDimensions((my_img_base->at(i)).getHeight() + temp_dims - 1,
			(my_img_base->at(i)).getWidth() + temp_dims - 1, -MY_INF);
	}


	/* Convolve the image with the templates at each 'my_K' orientation: */
	for (unsigned int k = 0; k < my_K; k++)
	{

		DEB_MSG("rotation: " << k);

		for (unsigned int i = 0; i < my_filters_imgs_count; i++)
		{
			for (unsigned int yR = 0; yR < ((my_img_base->at(i)).getHeight() + temp_dims - 1); yR++)
			{
				/* Determine the limits in the 'y' axis to keep the convolution inside the image rows: */
				const unsigned int min_y = (yR > (temp_dims - 1)) ?
					(yR - temp_dims + 1) : 0;

				const unsigned int max_y = (yR < (my_img_base->at(i)).getHeight()) ?
					(yR + 1) : (my_img_base->at(i)).getHeight();


				for (unsigned int xR = 0; xR < ((my_img_base->at(i)).getWidth() + temp_dims - 1); xR++)
				{
					/* Determine the limits in the 'x' axis to keep the convolution inside the image columns:*/
					const unsigned int min_x = (xR > (temp_dims - 1)) ? (xR - temp_dims + 1) : 0;
					const unsigned int max_x = (xR < (my_img_base->at(i)).getWidth()) ? (xR + 1) : (my_img_base->at(i)).getWidth();

					double resp_k = 0.0;
					/* Convolve the template with the input image: */
					for (unsigned int y = min_y; y < max_y; y++)
					{
						for (unsigned int x = min_x; x < max_x; x++)
						{
							resp_k += (GMF_templates + k)->getPix(max_y - y - 1, max_x - x - 1) *
								(my_img_base->at(i)).getPix(y, x);
						}
					}

					/* Save the response for each pixel among the 'my_K' orientations: */
					if (resp_k > (resp_tmp + i)->getPix(yR, xR))
					{
						(resp_tmp + i)->setPix(yR, xR, resp_k);
					}
				}
			}
		}
	}

	/* Keep the highest response for each pixel among the whole orientations as the final filter response:*/
	for (unsigned int i = 0; i < my_filters_imgs_count; i++) {
		(my_img_response->at(i)).setDimensions((my_img_base->at(i)).getHeight(),
			(my_img_base->at(i)).getWidth());

		for (unsigned int y = 0; y < (my_img_base->at(i)).getHeight(); y++)
		{
			for (unsigned int x = 0; x < (my_img_base->at(i)).getWidth(); x++)
			{
				(my_img_response->at(i)).setPix(y, x, (resp_tmp + i)->getPix(y + offset, x + offset));

			}
		}
	}

	/* Free memory */
	delete[] GMF_templates;
	delete[] resp_tmp;
}







/************************************************************************************************************
* FILTROS::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: fftImgOrigen                                                                               *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* -------                   --------                   -   -------------------------------------------      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Computes the Fourier transformation of the base images.                                                   *
*                                                                                                           *
************************************************************************************************************/
void FILTROS::fftImgOrigen()
{
	DEB_MSG("Transform image to the Fourier space");

	if (!already_transformed) {
		if (input_already_set) {
			if (!Img_fft) {
				Img_fft = new fftwl_complex*[my_filters_imgs_count];
			}
			if (!Img_fft_HPF) {
				Img_fft_HPF = new fftwl_complex*[my_filters_imgs_count];
			}

			for (unsigned int i = 0; i < my_filters_imgs_count; i++) {
				*(Img_fft + i) = (fftwl_complex*)fftwl_malloc((my_img_base->at(i)).getHeight() *
					(my_img_base->at(i)).getWidth() * sizeof(fftwl_complex));

				*(Img_fft_HPF + i) = (fftwl_complex*)fftwl_malloc((my_img_base->at(i)).getHeight() *
					(my_img_base->at(i)).getWidth() * sizeof(fftwl_complex));
			}
			already_transformed = true;

			fftwl_complex *Img_org = (fftwl_complex*) fftwl_malloc((my_img_base->at(0)).getHeight() * (my_img_base->at(0)).getWidth() * sizeof(fftwl_complex));
			
			DEB_MSG("Height: " << (my_img_base->at(0)).getHeight());
			DEB_MSG("Width: " << (my_img_base->at(0)).getWidth());

			for (unsigned int i = 0; i < my_filters_imgs_count; i++) {
				const unsigned int my_img_width = (my_img_base->at(i)).getWidth();
				const unsigned int my_img_height = (my_img_base->at(i)).getHeight();

				/* Copy data from image 'i' to temporal Img_org array */
				for (unsigned int xy = 0; xy < (my_img_height*my_img_width); xy++) {
					*(*(Img_org + xy)) = 1.0L - (long double)(my_img_base->at(i)).getPix(xy);
					*(*(Img_org + xy)+1) = 0.0L;
				}

				DEB_MSG("img: " << *(*(Img_org)));

				fftwl_plan p_r2c = fftwl_plan_dft_2d(my_img_width, my_img_height, Img_org, *(Img_fft + i), FFTW_FORWARD, FFTW_ESTIMATE);
				fftwl_execute(p_r2c);

				fftwl_destroy_plan(p_r2c);
				
				DEB_MSG("fft 0: " << *(*(*(Img_fft + i) + 0)) << " + " << *(*(*(Img_fft + i) + 0) + 1) << "i");
				DEB_MSG("fft 1: " << *(*(*(Img_fft + i) + 1)) << " + " << *(*(*(Img_fft + i) + 1) + 1) << "i");
				DEB_MSG("fft 2: " << *(*(*(Img_fft + i) + 2)) << " + " << *(*(*(Img_fft + i) + 2) + 1) << "i");
				DEB_MSG("fft 3: " << *(*(*(Img_fft + i) + 3)) << " + " << *(*(*(Img_fft + i) + 3) + 1) << "i");
				DEB_MSG("fft 4: " << *(*(*(Img_fft + i) + 4)) << " + " << *(*(*(Img_fft + i) + 4) + 1) << "i");
			}
			fftwl_free(Img_org);
		}

		fftwl_cleanup();

	}
}



/*  Metodo: respGabor

    Funcion: Obtiene la respuesta del filtro de escala simple de Gabor, largo del template: my_L, ancho del template: my_T, y numero de rotaciones que se hacen al filtro entre 0 y 180°: my_K.
*/
void FILTROS::respGabor() {

	fftImgOrigen();

	// Calculate sx y sy:
	long double sx2 = (long double)my_T / (2.0*sqrt(2.0 * log(2.0)));
	sx2 *= sx2;
	const long double sy2 = my_L*my_L * sx2;

	const long double fx = 1.0 / (long double)my_T;
	const long double fy = 0.0;

	/* Generate the shifted high-pass template HPF */
	IMGCONT HPF((my_img_base->at(0)).getHeight(), (my_img_base->at(0)).getWidth()/2 + 1);

	long double *u = (long double*)malloc((my_img_base->at(0)).getWidth() * sizeof(long double));
	long double *v = (long double*)malloc((my_img_base->at(0)).getHeight() * sizeof(long double));

	const long double my_half_kernel_width = (long double)(my_img_base->at(0)).getWidth() / 2.0;
	const long double my_kernel_width_denominator = 2.0 * MY_PI / (long double)(my_img_base->at(0)).getWidth();
	const long double my_half_kernel_height = (long double)(my_img_base->at(0)).getHeight() / 2.0;
	const long double my_kernel_height_denominator = 2.0 * MY_PI / (long double)(my_img_base->at(0)).getHeight();

	for (unsigned int x = 0; x < (my_img_base->at(0)).getWidth(); x++) {
		*(u + x) = (x - my_half_kernel_width) / my_kernel_width_denominator;
	}
	for (unsigned int y = 0; y < (my_img_base->at(0)).getHeight(); y++) {
		*(v + y) = (y - my_half_kernel_height) / my_kernel_height_denominator;
	}

	for (unsigned int y = 0; y < (my_img_base->at(0)).getHeight() / 2; y++) {
		const long double v_y2 = *(v + y) * *(v + y);
		{
			const long double u_x2 = *(u) * *(u);
			HPF.setPix(y + (my_img_base->at(0)).getHeight() / 2, (my_img_base->at(0)).getWidth() / 2, 1.0 - exp(-(sy2*(u_x2 + v_y2)) / 2.0));
		}
		for (unsigned int x = (my_img_base->at(0)).getWidth() / 2; x < (my_img_base->at(0)).getWidth(); x++) {
			const long double u_x2 = *(u + x) * *(u + x);
			HPF.setPix(y + (my_img_base->at(0)).getHeight() / 2, x - (my_img_base->at(0)).getWidth() / 2, 1.0 - exp(-(sy2*(u_x2 + v_y2)) / 2.0));
		}
	}

	for (unsigned int y = (my_img_base->at(0)).getHeight() / 2; y < (my_img_base->at(0)).getHeight(); y++) {
		const long double v_y2 = *(v + y) * *(v + y);
		{
			const long double u_x2 = *(u) * *(u);
			HPF.setPix(y - (my_img_base->at(0)).getHeight() / 2, (my_img_base->at(0)).getWidth() / 2, 1.0 - exp(-(sy2*(u_x2 + v_y2)) / 2.0));
		}
		for (unsigned int x = (my_img_base->at(0)).getWidth() / 2; x < (my_img_base->at(0)).getWidth(); x++) {
			const long double u_x2 = *(u + x) * *(u + x);
			HPF.setPix(y - (my_img_base->at(0)).getHeight() / 2,  x - (my_img_base->at(0)).getWidth() / 2, 1.0 - exp(-(sy2*(u_x2 + v_y2)) / 2.0));
		}
	}

	/* 'Img_filter' is the temporal filtered image in the frequencies domain: */
	fftwl_complex *Img_filter = (fftwl_complex*)fftwl_malloc((my_img_base->at(0)).getHeight()*((my_img_base->at(0)).getWidth() / 2 + 1) * sizeof(fftwl_complex));

	/* 'Img_resp' contains the response of the gabor filter at certain orientation 'theta': */
	fftwl_complex *Img_resp = (fftwl_complex*)fftwl_malloc((my_img_base->at(0)).getHeight() * (my_img_base->at(0)).getWidth() * sizeof(fftwl_complex));

	/* 'max_resp' saves the highest response for every pixel of the input: */
	long double **max_resp = (long double**)malloc(my_filters_imgs_count * sizeof(long double*));

	fftwl_plan p_c2r;
	/* Apply the high-pass filter to the image in the frequencies domain: */
	for (unsigned int i = 0; i < my_filters_imgs_count; i++) {
		for (unsigned int y = 0; y < (my_img_base->at(i)).getHeight(); y++) {
			for (unsigned int x = 0; x <= (my_img_base->at(i)).getWidth(); x++) {
				*(*(*(Img_fft_HPF + i) + x + y*((my_img_base->at(i)).getWidth()))) =
					*(*(*(Img_fft + i) + x + y*((my_img_base->at(i)).getWidth()))) * HPF.getPix(y, x);
			}
		}

		*(max_resp + i) = (long double*)malloc((my_img_base->at(i)).getHeight() * (my_img_base->at(i)).getWidth() * sizeof(long double));
		for (unsigned int xy = 0; xy < (my_img_base->at(i)).getHeight() * (my_img_base->at(i)).getWidth(); xy++) {
			*(*(max_resp + i) + xy) = -MY_INF;
		}
	}

	const long double theta_increment = 180.0 / (long double)my_K;

	long double Vr, Ur;

	long double * Gabor_filter = (long double *)malloc((my_img_base->at(0)).getHeight() * (my_img_base->at(0)).getWidth() * sizeof(long double));
	const unsigned int gabor_kernel_width = (my_img_base->at(0)).getWidth();

	unsigned int k = 0;
	for (long double theta = 0.0; theta < 180.0; theta += theta_increment, k++) {
		const long double stheta = sin(theta*MY_PI / 180.0);
		const long double ctheta = cos(theta*MY_PI / 180.0);

		/* The Gabor filter is calculated for the rotated base at 'theta' degrees: */
		for (unsigned int y = 0; y < (my_img_base->at(0)).getHeight() / 2; y++) {
			const double v_y = *(v + y + (my_img_base->at(0)).getHeight() / 2);

			unsigned int x;
			for (x = 0; x < (my_img_base->at(0)).getWidth() / 2; x++) {
				Ur = *(u + x + (my_img_base->at(0)).getWidth() / 2)*ctheta + v_y*stheta;
				Vr = -*(u + x + (my_img_base->at(0)).getWidth() / 2)*stheta + v_y*ctheta;

				*(Gabor_filter + y*gabor_kernel_width + x) = exp(-(0.5)*(sx2*(Ur*Ur + 4.0*MY_PI*fx*MY_PI*fx) + sy2*(Vr*Vr + 4.0*MY_PI*fy*MY_PI*fy))) * cosh(2.0*MY_PI*(sx2*fx*Ur + sy2*fy*Vr));
			}

			Ur = *(u)*ctheta + v_y*stheta;
			Vr = -*(u)*stheta + v_y*ctheta;

			*(Gabor_filter + y*gabor_kernel_width + x) = exp(-(0.5)*(sx2*(Ur*Ur + 4.0*MY_PI*fx*MY_PI*fx) + sy2*(Vr*Vr + 4.0*MY_PI*fy*MY_PI*fy))) * cosh(2.0*MY_PI*(sx2*fx*Ur + sy2*fy*Vr));
		}

		for (unsigned int y = (my_img_base->at(0)).getHeight() / 2; y < (my_img_base->at(0)).getHeight(); y++) {
			const long double v_y = *(v + y - (my_img_base->at(0)).getHeight() / 2);

			unsigned int x;
			for (x = 0; x < (my_img_base->at(0)).getWidth() / 2; x++) {
				Ur = *(u + x + (my_img_base->at(0)).getWidth() / 2)*ctheta + v_y*stheta;
				Vr = -*(u + x + (my_img_base->at(0)).getWidth() / 2)*stheta + v_y*ctheta;

				*(Gabor_filter + y*gabor_kernel_width + x) = exp(-(0.5)*(sx2*(Ur*Ur + 4.0*MY_PI*fx*MY_PI*fx) + sy2*(Vr*Vr + 4.0*MY_PI*fy*MY_PI*fy)))	* cosh(2.0*MY_PI*(sx2*fx*Ur + sy2*fy*Vr));
			}

			Ur = *(u)*ctheta + v_y*stheta;
			Vr = -*(u)*stheta + v_y*ctheta;

			*(Gabor_filter + y*gabor_kernel_width + x) = exp(-(0.5)*(sx2*(Ur*Ur + 4.0*MY_PI*fx*MY_PI*fx) + sy2*(Vr*Vr + 4.0*MY_PI*fy*MY_PI*fy)))
				* cosh(2.0*MY_PI*(sx2*fx*Ur + sy2*fy*Vr));
		}

		/* Apply the filter to the n images: */
		for (unsigned int i = 0; i < my_filters_imgs_count; i++) {
			for (unsigned int y = 0; y < (my_img_base->at(i)).getHeight(); y++) {
				for (unsigned int x = 0; x <= (my_img_base->at(i)).getWidth() / 2; x++) {
					*(*(Img_filter + y*((my_img_base->at(i)).getWidth() / 2 + 1) + x)) =
						*(*(*(Img_fft_HPF + i) + y*((my_img_base->at(i)).getWidth() / 2 + 1) + x)) *
						*(Gabor_filter + y*gabor_kernel_width + x);

					*(*(Img_filter + y*((my_img_base->at(i)).getWidth() / 2 + 1) + x) + 1) =
						*(*(*(Img_fft_HPF + i) + y*((my_img_base->at(i)).getWidth() / 2 + 1) + x) + 1) *
						*(Gabor_filter + y*gabor_kernel_width + x);
				}
			}

			/* Transform the response to the original domain: */
			p_c2r = fftwl_plan_dft_2d((my_img_base->at(i)).getWidth(), (my_img_base->at(i)).getHeight(), Img_filter, Img_resp, FFTW_FORWARD, FFTW_ESTIMATE);
			fftwl_execute(p_c2r);

			/* Update the highest response for every pixel: */
			for (unsigned int xy = 0; xy < (my_img_base->at(i)).getHeight() * (my_img_base->at(i)).getWidth(); xy++) {
				if (*(*(max_resp + i) + xy) < *(*(Img_resp + xy))) {
					*(*(max_resp + i) + xy) = *(*(Img_resp + xy));
				} 
			}

			fftwl_destroy_plan(p_c2r);
		}
	}


	for (unsigned int i = 0; i < my_filters_imgs_count; i++) {
		(my_img_response->at(i)).setDimensions((my_img_base->at(i)).getHeight(), (my_img_base->at(i)).getWidth());

		for (unsigned int y = 0; y < (my_img_base->at(i)).getHeight(); y++) {
			for (unsigned int x = 0; x < (my_img_base->at(i)).getWidth(); x++) {
				if ((my_img_base_mask->at(i)).getPix(y, x) > 0.5) {
					(my_img_response->at(i)).setPix(y, x, (*(*(max_resp + i) + x + y*(my_img_base->at(i)).getWidth()) / (my_img_base->at(i)).getHeight() * (my_img_base->at(i)).getWidth()));
				}
			}
		}

		free(*(max_resp + i));
	}

	fftwl_free(Img_filter);
	fftwl_free(Img_resp);
	free(u);
	free(v);
	free(Gabor_filter);
	free(max_resp);
}













/************************************************************************************************************
* FILTROS::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: setInputFilterBase                                                                         *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_base              std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void FILTROS::setInputFilterBase(std::vector<IMGCONT>* new_img_base)
{
	my_img_base = new_img_base;
	input_already_set = true;
	my_filters_imgs_count = (unsigned int)new_img_base->size();
}













/************************************************************************************************************
* FILTROS::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: setInputFilterMask                                                                         *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_mask              std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void FILTROS::setInputFilterMask(std::vector<IMGCONT>* new_img_mask)
{
	my_img_base_mask = new_img_mask;

	my_filters_imgs_count = (unsigned int)new_img_mask->size();
}













/************************************************************************************************************
* FILTROS::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: setInputFilterGroundtruth                                                                  *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_groundtruth       std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void FILTROS::setInputFilterGroundtruth(std::vector<IMGCONT>* new_img_groundtruth)
{
	my_img_groundtruth = new_img_groundtruth;

	my_filters_imgs_count = (unsigned int)new_img_groundtruth->size();
}













/************************************************************************************************************
* FILTROS::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: setInputFilterResponse                                                                     *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_response          std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void FILTROS::setInputFilterResponse(std::vector<IMGCONT>* new_img_response)
{
	my_img_response = new_img_response;

	my_filters_imgs_count = (unsigned int)new_img_response->size();
}













/************************************************************************************************************
* FILTROS::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: setParL                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_L_value               const double               I   Fixed L value for GMF or SSG filter.             *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void FILTROS::setParL(const double new_L_value)
{
	my_L = new_L_value;
}













/************************************************************************************************************
* FILTROS::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: setParT                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_T_value               const unsigned int         I   Fixed T value for GMF or SSG filter.             *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void FILTROS::setParT(const unsigned int new_T_value)
{
	my_T = new_T_value;
}













/************************************************************************************************************
* FILTROS::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: setParK                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_K_value               const usnigned int         I   Fixed K value for GMF or SSG filter.             *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void FILTROS::setParK(const unsigned int new_K_value)
{
	my_K = new_K_value;
}













/************************************************************************************************************
* FILTROS::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: setParSigma                                                                                *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_sigma_value           const double               I   Fixed sigma value for GMF filter.                *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void FILTROS::setParSigma(const double new_sigma_value)
{
	my_sigma = new_sigma_value;
}