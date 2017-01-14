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


    //DEB_MSG("No se definio el pBar");

    /// Limpiar el resto de la linea:
    int max_ancho = 100;
    for( int i = 0; i <= max_ancho; i++){
        printf("\r");
    }
    printf(COLOR_BACK_RED "[");
    int avance_progreso = (int)((double)max_ancho * (double)avance / (double)max_progress);
    for( int i = 0; i < avance_progreso; i++){
        printf(COLOR_BACK_GREEN " ");
    }
    for( int i = avance_progreso; i < max_ancho; i++){
        printf(COLOR_BACK_CYAN " ");
    }
    printf(COLOR_BACK_RED "]" COLOR_NORMAL);
    fflush(stdout);
}




/*  Metodo: interpolacion (Lineal)
    Funcion: Interpola el valor del pixel destino en base a los 4 pixeles origen.
*/
inline double FILTROS::interpolacion(const double *pix, const int j, const int i, const double x, const double y, const int mis_rens, const int mis_cols){
    double intensidad = 0.0;

    if( j >= 0 && j < mis_cols ){
        if( i >= 0 && i < mis_rens ){
            intensidad += *(pix +   i  *mis_cols +   j  ) * (1.0 - x) * (1.0 - y);
        }
        if( i >= -1  && i < (mis_rens-1) ){
            intensidad += *(pix + (i+1)*mis_cols +   j  ) * (1.0 - x) * (   y   );
        }
    }
    if( j >= -1 && j < (mis_cols-1)){
        if( i >=  0 && i < mis_rens ){
            intensidad += *(pix +   i  *mis_cols + (j+1)) * (   x   ) * (1.0 - y);
        }
        if( i >= -1 && i < (mis_rens-1) ){
            intensidad += *(pix + (i+1)*mis_cols + (j+1)) * (   x   ) * (   y   );
        }
    }

    return intensidad;
}




/*  Metodo: rotarImg
    Funcion: Rota una imagen almacenada como intensidad de 0 a 1 en un angulo theta.
*/
void FILTROS::rotateImg( const double *org, double *rot, const double ctheta, const double stheta, const int mis_rens, const int mis_cols, const int org_rens, const int org_cols){

    const double mitad_x = (double)(mis_cols - 1) / 2.0;
    const double mitad_y = (double)(mis_rens - 1) / 2.0;
    const double mitad_orgx = (double)(org_cols - 1) / 2.0;
    const double mitad_orgy = (double)(org_rens - 1) / 2.0;


    double x, y;
    for(int i = 0; i < (mis_rens-1); i++){
        for(int j = 0; j < (mis_cols-1); j++){
            x = ((double)j-mitad_x)*ctheta + ((double)i-mitad_y)*stheta + mitad_orgx;
            y =-((double)j-mitad_x)*stheta + ((double)i-mitad_y)*ctheta + mitad_orgy;
            const int flx = (int)x;
            const int fly = (int)y;
            const double delta_x = x - (double)flx;
            const double delta_y = y - (double)fly;

            *(rot + i*mis_cols + j) = interpolacion(org, flx, fly, delta_x, delta_y, org_rens, org_cols);
        }
    }
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
			fftw_free(Img_fft[i]);
			fftw_free(Img_fft_HPF[i]);
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
void FILTROS::respGMF() {

	const unsigned int temp_dims = (int)(1.5*((my_T > (unsigned int)my_L) ? my_T : (unsigned int)my_L));

	double **templates = (double**)malloc(my_K * sizeof(double*));

	////// Calculate a single Gaussian curve along the template width, at 0 degrees:
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
	*(templates) = (double*)malloc(temp_dims * temp_dims * sizeof(double));

	//// Copy the Gaussian curve along the template height:
	memset(*(templates), 0, temp_dims*temp_dims * sizeof(double));
	for (unsigned int y = 0; y < (unsigned int)my_L; y++)
	{
		memcpy(*(templates) + (temp_dims / 2 - (unsigned int)my_L / 2 + y)*temp_dims + (temp_dims / 2 - my_T / 2), gauss_0, my_T * sizeof(double));
	}


	// Rotate template at each 'my_K' orientation:
	const double theta_inc = 180.0 / (double)my_K;
	double theta = 0.0;

	for (unsigned int k = 1; k < my_K; k++)
	{
		theta += theta_inc;
		const double ctheta = cos(-theta * MY_PI / 180.0);
		const double stheta = sin(-theta * MY_PI / 180.0);

		*(templates + k) = (double*)malloc(temp_dims * temp_dims * sizeof(double));
		memset(*(templates + k), 0, temp_dims*temp_dims * sizeof(double));
		rotateImg(templates[0], templates[k], ctheta, stheta, temp_dims, temp_dims, temp_dims, temp_dims);
	}

	free(gauss_0);


	//// Initialize the filter response:
	double *resp_tmp = (double*)malloc(((my_img_base->at(0)).getWidth() + temp_dims - 1) *
		((my_img_base->at(0)).getHeight() + temp_dims - 1) * sizeof(double));

	const unsigned int offset = (unsigned int)(temp_dims / 2);

	//// Convolve the image with the templates at each 'my_K' orientation:
	for (unsigned int i = 0; i < my_filters_imgs_count; i++)
	{
		for (unsigned int y = 0; y < ((my_img_base->at(i)).getHeight() + temp_dims - 1); y++) {
			for (unsigned int x = 0; x < ((my_img_base->at(i)).getWidth() + temp_dims - 1); x++) {
				*(resp_tmp + x + y * ((my_img_base->at(i)).getWidth() + temp_dims - 1)) = -MY_INF;
			}
		}

		for (unsigned int k = 0; k < my_K; k++)
		{
			for (unsigned int yR = 0; yR < ((my_img_base->at(i)).getHeight() + temp_dims - 1); yR++)
			{
				//// Determine the limits in the 'y' axis to keep the convolution inside the image rows:
				const unsigned int min_y = (yR > (temp_dims - 1)) ? (yR - temp_dims + 1) : 0;
				const unsigned int max_y = (yR < (my_img_base->at(i)).getHeight()) ? (yR + 1) : (my_img_base->at(i)).getHeight();
				for (unsigned int xR = 0; xR < ((my_img_base->at(i)).getWidth() + temp_dims - 1); xR++)
				{
					//// Determine the limits in the 'x' axis to keep the convolution inside the image columns:
					const unsigned int min_x = (xR > (temp_dims - 1)) ? (xR - temp_dims + 1) : 0;
					const unsigned int max_x = (xR < (my_img_base->at(i)).getWidth()) ? (xR + 1) : (my_img_base->at(i)).getWidth();

						double resp_k = 0.0;
						//// Convolve the template with the input image:
						for (unsigned int y = min_y; y < max_y; y++)
						{
							for (unsigned int x = min_x; x < max_x; x++)
							{
								resp_k += *(*(templates + k) + (max_y - y - 1)*temp_dims + (max_x - x - 1)) *
									(my_img_base->at(i)).getPix(y, x);
							}
						}

					//// Save the response for each pixel among the 'my_K' orientations:
					if (resp_k > *(resp_tmp + yR*((my_img_base->at(i)).getWidth() + temp_dims - 1) + xR))
					{
						*(resp_tmp +
							yR*((my_img_base->at(i)).getWidth() + temp_dims - 1) + xR +
							((my_img_base->at(i)).getWidth() + temp_dims - 1) *
								((my_img_base->at(i)).getHeight() + temp_dims - 1)) = resp_k;
					}
				}
			}
		}

		//// Keep the highest response for each pixel among the whole orientations as the final filter response:
		for (unsigned int y = 0; y < (my_img_base->at(i)).getHeight(); y++)
		{
			for (unsigned int x = 0; x < (my_img_base->at(i)).getWidth(); x++)
			{
				(my_img_response->at(i)).setPix(y, x,
					*(resp_tmp + (x + offset) + (y + offset)*((my_img_base->at(i)).getWidth() + temp_dims - 1)));

			}
		}
	}

	// Free memory
	for (unsigned int k = 0; k < my_K; k++)
	{
		free(*(templates + k));
	}
	free(templates);
	free(resp_tmp);
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

	DEB_MSG("my img base ptr: " << my_img_base);

    if(input_already_set){
		if (!Img_fft) {
			Img_fft = new fftw_complex*[my_filters_imgs_count];
		}
		if (!Img_fft_HPF) {
			Img_fft_HPF = new fftw_complex*[my_filters_imgs_count];
		}

		if (!already_transformed) {
			for (unsigned int i = 0; i < my_filters_imgs_count; i++) {
				*(Img_fft + i) = (fftw_complex*)fftw_malloc((my_img_base->at(i)).getHeight() *
					((my_img_base->at(i)).getWidth() / 2 + 1) * sizeof(fftw_complex));
				*(Img_fft_HPF + i) = (fftw_complex*)fftw_malloc((my_img_base->at(i)).getHeight() *
					((my_img_base->at(i)).getWidth() / 2 + 1) * sizeof(fftw_complex));
			}
			already_transformed = true;

			double *Img_org = new double[(my_img_base->at(0)).getHeight() * (my_img_base->at(0)).getWidth()];

			for (unsigned int i = 0; i < my_filters_imgs_count; i++) {
				for (unsigned int y = 0; y < (my_img_base->at(i)).getHeight(); y++) {
					for (unsigned int x = 0; x < (my_img_base->at(i)).getWidth(); x++) {
						*(Img_org + x + y * (my_img_base->at(i)).getWidth()) = 1.0 - (my_img_base->at(i)).getPix(y, x);
					}
				}

				fftw_plan p_r2c = fftw_plan_dft_r2c_2d((my_img_base->at(i)).getHeight(), (my_img_base->at(i)).getWidth(),
					Img_org, *(Img_fft + i), FFTW_ESTIMATE);
				fftw_execute(p_r2c);
				fftw_destroy_plan(p_r2c);
			}

			delete[] Img_org;
		}
    }
}



/*  Metodo: respGabor

    Funcion: Obtiene la respuesta del filtro de escala simple de Gabor, largo del template: my_L, ancho del template: my_T, y numero de rotaciones que se hacen al filtro entre 0 y 180°: my_K.
*/
void FILTROS::respGabor() {

	fftImgOrigen();

	DEB_MSG("L: " << my_L << ", T: " << my_T << ", K: " << my_K);
	DEB_MSG("height: " << (my_img_base->at(0)).getHeight() << ", width: " << (my_img_base->at(0)).getWidth());

	// Calculate sx y sy:
	double sx2 = (double)my_T / (2.0*sqrt(2.0 * log(2.0)));
	sx2 *= sx2;
	const double sy2 = my_L*my_L * sx2;

	const double fx = 1.0 / (double)my_T;
	const double fy = 0.0;

	double *u = (double*)malloc((my_img_base->at(0)).getWidth() * sizeof(double));
	double *v = (double*)malloc((my_img_base->at(0)).getHeight() * sizeof(double));

	for (unsigned int x = 0; x < (my_img_base->at(0)).getWidth(); x++) {
		*(u + x) = (x - (double)(my_img_base->at(0)).getWidth() / 2.0) * (2.0 * MY_PI / (double)(my_img_base->at(0)).getWidth());
	}
	for (unsigned int y = 0; y < (my_img_base->at(0)).getHeight(); y++) {
		*(v + y) = (y - (double)(my_img_base->at(0)).getHeight() / 2.0) * (2.0 * MY_PI / (double)(my_img_base->at(0)).getHeight());
	}

	/* Generate the high-pass template HPF */
	IMGCONT HPF((my_img_base->at(0)).getHeight(), (my_img_base->at(0)).getWidth()/2 + 1);

	for (unsigned int y = 0; y < (my_img_base->at(0)).getHeight() / 2; y++) {
		const double v_y2 = *(v + y) * *(v + y);
		{
			const double u_x2 = *(u) * *(u);
			HPF.setPix(y + (my_img_base->at(0)).getHeight() / 2, (my_img_base->at(0)).getWidth() / 2, 1.0 - exp(-(sy2*(u_x2 + v_y2)) / 2.0));
		}
		for (unsigned int x = (my_img_base->at(0)).getWidth() / 2; x < (my_img_base->at(0)).getWidth(); x++) {
			const double u_x2 = *(u + x) * *(u + x);
			HPF.setPix(y + (my_img_base->at(0)).getHeight() / 2, x - (my_img_base->at(0)).getWidth() / 2, 1.0 - exp(-(sy2*(u_x2 + v_y2)) / 2.0));
		}
	}

	for (unsigned int y = (my_img_base->at(0)).getHeight() / 2; y < (my_img_base->at(0)).getHeight(); y++) {
		const double v_y2 = *(v + y) * *(v + y);
		{
			const double u_x2 = *(u) * *(u);
			HPF.setPix(y - (my_img_base->at(0)).getHeight() / 2, (my_img_base->at(0)).getWidth() / 2, 1.0 - exp(-(sy2*(u_x2 + v_y2)) / 2.0));
		}
		for (unsigned int x = (my_img_base->at(0)).getWidth() / 2; x < (my_img_base->at(0)).getWidth(); x++) {
			const double u_x2 = *(u + x) * *(u + x);
			HPF.setPix(y - (my_img_base->at(0)).getHeight() / 2,  x - (my_img_base->at(0)).getWidth() / 2, 1.0 - exp(-(sy2*(u_x2 + v_y2)) / 2.0));
		}
	}

	/* 'Img_filter' is the temporal filtered image in the frequencies domain: */
	fftw_complex *Img_filter = (fftw_complex*)fftw_malloc((my_img_base->at(0)).getHeight()*((my_img_base->at(0)).getWidth() / 2 + 1) * sizeof(fftw_complex));

	/* 'Img_resp' contains the response of the gabor filter at certain orientation 'theta': */
	double *Img_resp = (double*)calloc((my_img_base->at(0)).getHeight() * (my_img_base->at(0)).getWidth(), sizeof(double));

	/* 'max_resp' saves the highest response for every pixel of the input: */
	double **max_resp = (double**)malloc(my_filters_imgs_count * sizeof(double*));

	fftw_plan p_c2r;
	/* Apply the high-pass filter to the image in the frequencies domain: */
	for (unsigned int i = 0; i < my_filters_imgs_count; i++) {
		for (unsigned int y = 0; y < (my_img_base->at(i)).getHeight(); y++) {
			for (unsigned int x = 0; x <= (my_img_base->at(i)).getWidth() / 2; x++) {
				*(*(*(Img_fft_HPF + i) + x + y*((my_img_base->at(i)).getWidth() / 2 + 1))) =
					*(*(*(Img_fft + i) + x + y*((my_img_base->at(i)).getWidth() / 2 + 1))) * HPF.getPix(y, x);
				*(*(*(Img_fft_HPF + i) + x + y*((my_img_base->at(i)).getWidth() / 2 + 1)) + 1) =
					*(*(*(Img_fft + i) + x + y*((my_img_base->at(i)).getWidth() / 2 + 1)) + 1) * HPF.getPix(y, x);
			}
		}

		*(max_resp + i) = (double*)malloc((my_img_base->at(i)).getHeight() * (my_img_base->at(i)).getWidth() * sizeof(double));
		for (unsigned int xy = 0; xy < (my_img_base->at(i)).getHeight() * (my_img_base->at(i)).getWidth(); xy++) {
			*(*(max_resp + i) + xy) = -MY_INF;
		}
	}

	const double theta_increment = 180.0 / (double)my_K;

	double Vr, Ur;

	unsigned int k = 0;
	for (double theta = 0.0; theta < 180.0; theta += theta_increment, k++) {

		DEB_MSG("theta: " << theta);

		const double stheta = sin(theta*MY_PI / 180.0);
		const double ctheta = cos(theta*MY_PI / 180.0);

		IMGCONT Gabor_filter((my_img_base->at(0)).getHeight(), ((my_img_base->at(0)).getWidth() / 2 + 1));

		/* The Gabor filter is calculated for the rotated base at 'theta' degrees: */
		for (unsigned int y = 0; y < (my_img_base->at(0)).getHeight() / 2; y++) {
			const double v_y = *(v + y + (my_img_base->at(0)).getHeight() / 2);

			unsigned int x;
			for (x = 0; x < (my_img_base->at(0)).getWidth() / 2; x++) {
				Ur = *(u + x + (my_img_base->at(0)).getWidth() / 2)*ctheta + v_y*stheta;
				Vr = -*(u + x + (my_img_base->at(0)).getWidth() / 2)*stheta + v_y*ctheta;

				Gabor_filter.setPix(y, x, exp(-(0.5)*(sx2*(Ur*Ur + 4.0*MY_PI*fx*MY_PI*fx) + sy2*(Vr*Vr + 4.0*MY_PI*fy*MY_PI*fy))) * cosh(2.0*MY_PI*(sx2*fx*Ur + sy2*fy*Vr)));
			}

			Ur = *(u)*ctheta + v_y*stheta;
			Vr = -*(u)*stheta + v_y*ctheta;

			Gabor_filter.setPix(y, x, exp(-(0.5)*(sx2*(Ur*Ur + 4.0*MY_PI*fx*MY_PI*fx) + sy2*(Vr*Vr + 4.0*MY_PI*fy*MY_PI*fy))) * cosh(2.0*MY_PI*(sx2*fx*Ur + sy2*fy*Vr)));
		}

		for (unsigned int y = (my_img_base->at(0)).getHeight() / 2; y < (my_img_base->at(0)).getHeight(); y++) {
			const double v_y = *(v + y - (my_img_base->at(0)).getHeight() / 2);

			unsigned int x;
			for (x = 0; x < (my_img_base->at(0)).getWidth() / 2; x++) {
				Ur = *(u + x + (my_img_base->at(0)).getWidth() / 2)*ctheta + v_y*stheta;
				Vr = -*(u + x + (my_img_base->at(0)).getWidth() / 2)*stheta + v_y*ctheta;

				Gabor_filter.setPix(y, x, exp(-(0.5)*(sx2*(Ur*Ur + 4.0*MY_PI*fx*MY_PI*fx) + sy2*(Vr*Vr + 4.0*MY_PI*fy*MY_PI*fy)))	* cosh(2.0*MY_PI*(sx2*fx*Ur + sy2*fy*Vr)));

			}

			Ur = *(u)*ctheta + v_y*stheta;
			Vr = -*(u)*stheta + v_y*ctheta;

			Gabor_filter.setPix(y, x, exp(-(0.5)*(sx2*(Ur*Ur + 4.0*MY_PI*fx*MY_PI*fx) + sy2*(Vr*Vr + 4.0*MY_PI*fy*MY_PI*fy)))
				* cosh(2.0*MY_PI*(sx2*fx*Ur + sy2*fy*Vr)));
		}

		/* Apply the filter to the n images: */
		for (unsigned int i = 0; i < my_filters_imgs_count; i++) {
			for (unsigned int y = 0; y < (my_img_base->at(i)).getHeight(); y++) {
				for (unsigned int x = 0; x <= (my_img_base->at(i)).getWidth() / 2; x++) {
					*(*(Img_filter + y*((my_img_base->at(i)).getWidth() / 2 + 1) + x)) =
						*(*(*(Img_fft_HPF + i) + y*((my_img_base->at(i)).getWidth() / 2 + 1) + x)) *
						Gabor_filter.getPix(y, x);

					*(*(Img_filter + y*((my_img_base->at(i)).getWidth() / 2 + 1) + x) + 1) =
						*(*(*(Img_fft_HPF + i) + y*((my_img_base->at(i)).getWidth() / 2 + 1) + x) + 1) *
						Gabor_filter.getPix(y, x);
				}
			}

			/* Transform the response to the original domain: */
			p_c2r = fftw_plan_dft_c2r_2d((my_img_base->at(i)).getHeight(), (my_img_base->at(i)).getWidth(), Img_filter, Img_resp, FFTW_ESTIMATE);
			fftw_execute(p_c2r);

			/* Update the highest response for every pixel: */
			for (unsigned int xy = 0; xy < (my_img_base->at(i)).getHeight() * (my_img_base->at(i)).getWidth(); xy++) {
				if (*(*(max_resp + i) + xy) < *(Img_resp + xy)) {
					*(*(max_resp + i) + xy) = *(Img_resp + xy);
				} 
			}

			fftw_destroy_plan(p_c2r);
		}
	}


	for (unsigned int i = 0; i < my_filters_imgs_count; i++) {

		DEB_MSG("Image " << i << " filtered ...");

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

	DEB_MSG("All images filtered ...................................");

	fftw_free(Img_filter);
	free(Img_resp);
	free(u);
	free(v);
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