/******
    CENTRO DE INVESTIGACION EN MATEMATICAS
    MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

    FERNANDO CERVANTES SANCHEZ
   NOV - 2015
*****/


#include "filtros.h"

#include <assert.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <fftw3.h>

#include <math.h>
#include <omp.h>

#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>


#ifdef _OPENMP
    #define TIMERS double t_ini, t_fin
    #define GETTIME_INI t_ini = omp_get_wtime()
    #define GETTIME_FIN t_fin = omp_get_wtime()
    #define DIFTIME (t_fin - t_ini)
    #define OMP_ENABLED true
#else
    #include <sys/time.h>
    #define TIMERS struct timeval t_ini, t_fin
    #define GETTIME_INI gettimeofday( &t_ini, NULL)
    #define GETTIME_FIN gettimeofday( &t_fin, NULL)
    #define DIFTIME ((t_fin.tv_sec*1e6 + t_fin.tv_usec) - (t_ini.tv_sec*1e6 + t_ini.tv_usec) )/ 1e6
    #define omp_get_num_threads() 1
    #define omp_set_num_threads(cores)
    #define omp_get_thread_num() 0
    #define OMP_ENABLED false
#endif


#define PI 3.14159265


#ifndef NDEBUG
    #define DEB_MSG(MENSAJE) using namespace std;\
                             cout << MENSAJE << endl;
#else
    #define DEB_MSG(MENSAJE)
#endif

#define PI 3.14159265




/*  Metodo: interpolacion(const double *pix, const int x, const int y, const double delta_x, const double delta_y)
    Funcion: Interpola el valor del pixel destino en base a los 4 pixeles origen.
*/
inline double interpolacion(const double *pix, const int x, const int y, const double delta_x, const double delta_y, const int cols){
    return ( (1.0-delta_x)*(1.0-delta_y) * *(pix + y*cols + x) +
             (  delta_x  )*(1.0-delta_y) * *(pix + y*cols + x+1) +
             (1.0-delta_x)*(  delta_y  ) * *(pix + (y+1)*cols + x) +
             (  delta_x  )*(  delta_y  ) * *(pix + (y+1)*cols + x+1)
           );
}




/*  Metodo: rotarImg( const double *org, double *rot, const double theta, const int rens, const int cols)
    Funcion: Rota una imagen almacenada como intensidad de 0 a 1 en un angulo theta.
*/
void rotarImg( const double *org, double *rot, const double theta, const int rens, const int cols){
    int x_i, y_i;
    double stheta = sin( theta * PI/180.0 );
    double ctheta = cos( theta * PI/180.0 );

    double delta_x, delta_y, x, y;

    for(register int i = 0; i < (rens-1); i++){
        for(register int j = 0; j < (cols-1); j++){
            x = (double)(j - cols/2)*ctheta - (double)(i - rens/2)*stheta + cols/2;
            y = (double)(j - cols/2)*stheta + (double)(i - rens/2)*ctheta + rens/2;

            x_i = (int)x;
            y_i = (int)y;

            // Si el pixel rotado sale de los limites del template, se descarta:
            if( x_i >= cols || x_i <= 0 || y_i >= rens || y_i <= 0){
                continue;
            }

            delta_x = x - (double)x_i;
            delta_y = y - (double)y_i;

            *(rot + y_i*cols + x_i) = interpolacion(org, j, i, delta_x, delta_y, cols);
        }
    }
}




/*//        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
 *
 *              M E T O D O S       P U B L I C O S     C O N S T R U C T O R E S / D E S T R U C T O R E S / S E T / G E T
 *
 * //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
*/
/*  Metodo: setEvoMet
    Funcion: Define el metodo evolutivo que se usara para establecer los mejores aprametros del filtro.
*/
void FILTROS::setEvoMet( const EVO_MET evo_met, const int m_iters, const int pob){
    metodo_elegido = evo_met;
    max_iters = m_iters;
    n_pob = pob;
}



/*  Metodo: setFiltro
    Funcion: Define el metodo de fitrado que se usara sobre la imagen de origen.
*/
void FILTROS::setFiltro( const SEG_FILTRO seg_fil){
    filtro_elegido = seg_fil;
    switch( filtro_elegido ){
        case SS_GABOR:
            pars_optim[3] = false;
            mi_elite->vars[3] = 0.0;
            max_bits[3] = 0;
            transformada = false;
            break;
    }
}


/*  Metodo: setFitness
    Funcion: Define que funcion de evaluacion se usara en el metodo de optimizacion como funcion objetivo.
*/
void FILTROS::setFitness( const FITNESS fit_fun){
    fitness_elegido = fit_fun;
}



/*  Metodo: FILTROS

    Funcion: Inicializa los parametros en valores por defecto.
*/
FILTROS::FILTROS(){
    n_pars = 5;
    for( int p = 0; p < 5; p++){
        // Por defecto, se optimizan todos los parametros.
        pars_optim[p] = true;
        max_bits[p] = 0;
    }

    mi_elite = new INDIV;
    semilla = NULL;

    setFiltro(GMF);
    setEvoMet(EDA_BUMDA, 15, 20);
    setFitness(ROC);

    transformada = false;
}


/*  Metodo: ~FILTROS

    Funcion: Libera la memoria requerida por el objeto.
*/
FILTROS::~FILTROS(){
    delete mi_elite;
    if( semilla ){
        delete semilla;
    }
    if(transformada){
        delete [] Img_org;
        fftw_free(Img_fft);
    }
}



/*  Metodo: setInput
    Funcion: Establece las imagenes origen yu ground truth.
*/
void FILTROS::setInput( IMGVTK &img_org, IMGVTK &img_ground){
    org = img_org.base_ptr;
    ground_truth = img_ground.base_ptr;
    mask = img_org.mask_ptr;

    // Obtener las dimensiones de la imagen:
    cols = img_org.cols;
    rens = img_org.rens;
    rens_cols = img_org.rens_cols;

    switch( filtro_elegido ){
        case SS_GABOR:
            fftImgOrigen();
            break;
    }
}


/*  Metodo: setOutput
    Funcion: Establece sobre que objeto se escribira la imagen resultante del filtro.
*/
void FILTROS::setOutput( IMGVTK &img_dest){
    dest = img_dest.base_ptr;
}


/*  Metodo: setPar
    Funcion: Ejecuta el metodo de optimizacion seleccionado para encontrar los mejores parametros para el filtro sobre la imagen de origen.
*/
void FILTROS::setPar(){
    switch( metodo_elegido ){
        case EA_GA:
            GA();
            break;
        case EDA_BUMDA:
            BUMDA();
            break;
        case EDA_UMDA:
            UMDA();
            break;
        case EXHAUSTIVA:
            TIMERS;
            for( double L = 1.3; L <= 18.1; L+= 0.4){
                mi_elite->vars[0] = L;
                for( double T = 1.0; T <= 16.00; T+=1.0 ){
                    mi_elite->vars[1] = T;
                    GETTIME_INI;
                    mi_elite->eval = fitnessROC( mi_elite );
                    GETTIME_FIN;
                    cout << "Evaluado en " << DIFTIME << " segundos. eval: " << mi_elite->eval << " [L: " << L << ", T: " << T << "]." << endl;
                }
            }
            break;
    }
}


/*  Metodo: setPar
    Funcion: Establece los parametros del filtro.
*/
void FILTROS::setPar( const PARAMETRO par, const double val ){
    mi_elite->vars[ par ] = val;
    pars_optim[ par ] = false;
    max_bits[ par ] = 0;
}


/*  Metodo: setLim
    Funcion: Establece los limites de busqueda de los parametros, asi como la varianza minima que se espera obtener (para BUMDA)
*/
void FILTROS::setLim( const PARAMETRO par, const double inf, const double sup, const unsigned char bits){
    lim_inf[ par ] = inf;
    lim_sup[ par ] = sup;
    pars_optim[ par ] = true; // Si se dan los limites de busqueda, se reactiva el parametro.
    max_bits[ par ] = bits;
}


/*  Metodo: setLim
    Funcion: Establece los limites de busqueda de los parametros, asi como la varianza minima que se espera obtener (para BUMDA)
*/
void FILTROS::setLim( const PARAMETRO par, const double inf, double sup, const double min_var){
    lim_inf[ par ] = inf;
    lim_sup[ par ] = sup;
    min_vars[ par ] = min_var;
    pars_optim[ par ] = true; // Si se dan los limites de busqueda, se reactiva el parametro.
}




/*  Metodo: filtrar
    Funcion: Aplica el filtro sobre la imagen de origen y la almacena sobre la imagen destino.
*/
void FILTROS::filtrar(){
    double *resp = new double [rens_cols];

    // Ejecutar el filtro con los parametros optimos:
    switch( filtro_elegido ){
        case GMF:
            respGMF(mi_elite, resp);
            break;
        case SS_GABOR:
            respGabor(mi_elite, resp);
            break;
    }

    double min_resp =  1e100;
    double max_resp = -1e100;

    for(int xy = 0; xy < rens_cols; xy++){
        if( resp[xy] < min_resp ){
            min_resp = resp[xy];
        }
        if( resp[xy] > max_resp ){
            max_resp = resp[xy];
        }
    }
    const double rango_resp = max_resp - min_resp;
    for( int xy = 0; xy < rens_cols; xy++){
        dest[xy] = (unsigned char)(255.0 * (resp[xy] - min_resp) / rango_resp);
    }

    DEB_MSG("delta: " << mi_elite->vars[4])

    calcROC(mi_elite, resp);

    delete [] resp;
}




/*//        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
 *
 *              M E T O D O S       D E     F I L T R A D O     D E     I M A G E N E S
 *
 * //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
*/


/*  Metodo: respGMF

    Funcion: Obtiene la respuesta del filtro Gaussiano Multiescala, con desviacion estandar: sigma, largo del template: L, ancho del template: T, y numero de rotaciones que se hacen al filtro entre 0 y 180°: K.
*/
void FILTROS::respGMF(INDIV *test, double *resp){

    TIMERS;

    GETTIME_INI;
    const int L = (int)test->vars[0];
    const int T = (int)test->vars[1];
    const int K = (int)test->vars[2];
    const double sigma = test->vars[3];

    DEB_MSG("L: " << L << ", T: " << T << ", K: " << K << ", Sigma: " << sigma);

    // Se calculan los templates para las rotaciones:
    double **templates = new double*[K];

    //// Se calcula el template Gaussiano en rotacion 0°.
    const int ancho_tmp = T + 3; // El ancho del template es el ancho de la gaussiana + 2 pixeles a cada lado para dejar espacio a la rotacion.
    const int alto_tmp = L + 6; // El alto del template es el largo de la gaussiana dado por L + 3 pixeles abajo y otros 3 arriba para dar espacio a la rotacion.

    const int ex_alto = alto_tmp%2;
    const int ex_ancho = ancho_tmp%2;

    templates[0] = new double[alto_tmp*ancho_tmp];
    memset( templates[0], 0, alto_tmp*ancho_tmp*sizeof(double));

    const double sig_2 = 2.0 * sigma * sigma;
    double *gaussiana_org = templates[0] + ancho_tmp*3 + ancho_tmp/2;

    ////// Se calcula una linea de la gaussiana para el template, luego se copia hacia abajo:
    double sum = 0.0;
    for( int x = -(int)(T/2); x <= (int)(T/2); x++){
        *(gaussiana_org + x) = 1.0 - exp( -((double)(x*x) / sig_2) );
        sum += *(gaussiana_org + x);
    }

    const double media = sum / (double)T;
    //// Se resta la media a todo el template, y se divide entre la suma:
    for( int x = -(int)(T/2); x <= (int)(T/2); x++){
        *(gaussiana_org + x) = (*(gaussiana_org + x) - media) / (sum*L);
    }

    //// Se termina de construir el template a 0°:
    gaussiana_org = templates[0] + 3*ancho_tmp;
    double *gaussiana_cpy = templates[0] + 3*ancho_tmp;
    for( int y = 1; y < L; y++ ){
        memcpy( gaussiana_cpy + y*ancho_tmp, gaussiana_org, ancho_tmp*sizeof(double) );
    }

    // Rotar el template segun el numero de rotaciones 'K':
    const double theta_inc = 180.0 / (double)K;
    double theta = 180.0;
    for( int k = 1; k < K; k++){
        theta -= theta_inc;
        templates[k] = new double[alto_tmp*ancho_tmp];
        memset(templates[k], 0, alto_tmp*ancho_tmp*sizeof(double));
        rotarImg( templates[0], templates[k], theta, alto_tmp, ancho_tmp);
    }
    DEB_MSG("Templates generados: " << K);

    ////--------------------------------------------------------- Aplicacion del filtro:
    // Imagen auxiliar:
    unsigned char **Img_aux = new unsigned char* [rens + 2*alto_tmp];
    for( int y = 0; y < alto_tmp; y++){
        Img_aux[y] = new unsigned char[cols + 2*ancho_tmp];
        Img_aux[y+alto_tmp+rens] = new unsigned char[cols + 2*ancho_tmp];
        memset( Img_aux[y], 0, (cols + 2*ancho_tmp)*sizeof(unsigned char));
        memset( Img_aux[y+alto_tmp+rens], 0, (cols + 2*ancho_tmp)*sizeof(unsigned char));
    }
    for(int y = 0; y < rens; y++){
        Img_aux[y+alto_tmp] = new unsigned char[cols + 2*ancho_tmp];
        memset( Img_aux[y+alto_tmp], 0, (cols + 2*ancho_tmp)*sizeof(unsigned char));
        memcpy( Img_aux[y+alto_tmp] + ancho_tmp, org + y*cols, cols*sizeof(unsigned char));
    }

    DEB_MSG("Imagen auxiliar generada");
    //    #pragma omp parallel for shared(resp, Img_amp, templates) firstprivate(rens, cols, ancho_tmp, alto_tmp, K) reduction(min: min_img) reduction(max: max_img)
    for( int yI = 0; yI < rens; yI++){
        for( int xI = 0; xI < cols; xI++){
            // Para quedarme con la mayor respuesta de entre los filtros:
            double max_resp = -1e10;

            // Aplicar los templates a la imagen:
            for( int k = 0; k < K; k++ ){
                double resp_k = 0.0;

                // Recorrer todo el template:
                for( int y = -alto_tmp/2; y < alto_tmp/2 - (1 - ex_alto); y++){
                    for( int x = -ancho_tmp/2; x < ancho_tmp/2 - (1 - ex_ancho); x++){
                        resp_k += *(templates[k] + (y+alto_tmp/2)*ancho_tmp + (x+ancho_tmp/2)) * (double)Img_aux[yI + y + alto_tmp][xI + x + ancho_tmp] / 255.0;
                        if( resp_k > 1000000  || resp_k < -1000000){
                            DEB_MSG("Woa woa algo anda mal aqui [" << x << "," <<  y << "] {" << xI << ", " << yI << "}: " << resp_k);
                            DEB_MSG("con el template[" << k << "]: " << (*(templates[k] + (y+alto_tmp/2)*ancho_tmp + (x+ancho_tmp/2))) << " en: " << (y-alto_tmp/2)*ancho_tmp + (x+ancho_tmp/2) << " de: " << alto_tmp*ancho_tmp);
                            DEB_MSG("O con la imagen " << ((double)Img_aux[yI + y + alto_tmp][xI + x + ancho_tmp] / 255.0));
                            exit(0);
                        }
                    }
                }
                if( resp_k > max_resp ){
                    max_resp = resp_k;
                }
            }

            *(resp + yI*cols + xI) = max_resp;
        }
    }

    DEB_MSG("filtrado");

    // Liberar la matriz de templates:
    for( int k = 0; k < K; k++){
        delete [] templates[k];
    }
    delete [] templates;

    //Liberar la imagen auxiliar:
    for(int y = 0; y < (rens+alto_tmp); y++){
        delete [] Img_aux[y];
    }
    delete [] Img_aux;
    GETTIME_FIN;
}



/*  Metodo: fftImgOrigen

    Funcion: Obtiene la transformada de Fourier de la imagen original.
*/
void FILTROS::fftImgOrigen(){
    if(!transformada){
        transformada = true;
        Img_org = new double[rens_cols];
        Img_fft = (fftw_complex*) fftw_malloc(rens*(cols/2+1)*sizeof(fftw_complex));
        for(int xy = 0; xy < rens_cols; xy++){
            Img_org[xy] = 1.0 - (double)org[xy] / 255.0;
        }
        fftw_plan p_r2c = fftw_plan_dft_r2c_2d( rens, cols, Img_org, Img_fft, FFTW_ESTIMATE);
        fftw_execute(p_r2c);
        fftw_destroy_plan(p_r2c);
    }
}



/*  Metodo: respGabor

    Funcion: Obtiene la respuesta del filtro de escala simple de Gabor, largo del template: L, ancho del template: T, y numero de rotaciones que se hacen al filtro entre 0 y 180°: K.
*/
void FILTROS::respGabor(INDIV *test, double *resp){

    TIMERS;

    GETTIME_INI;

    const double L = test->vars[0];
    const int T = (int)test->vars[1];
    const int K = (int)test->vars[2];

    // Calcular sx y sy:
    const double sx = (double)T / (2.0*sqrt(2.0 * log(2.0)));
    const double sy = L * sx;

    const double fx = 1.0 / (double)T;
    const double fy = 0.0;

    // Para generar los 'u' y 'v':
    double *u = new double[cols*sizeof(double)];
    double *v = new double[rens*sizeof(double)];

    for( int x = 0; x < cols; x++){
        u[x] = (x - (double)cols/2.0) * (2.0 * PI / (double)cols);
    }
    for( int y = 0; y < rens; y++){
        v[y] = (y - (double)rens/2.0) * (2.0 * PI / (double)rens);
    }

    //// Generar HPF como la gaussiana, se hace por partes para cambiar los cuadrantes (como fftshift en Matlab)
    //// Convertir la imagen del rango de 0 a 255 a uno continuo de 0 a 1:
    double *HPF = new double[rens_cols];

    for( int y = 0; y < rens/2; y++ ){
        const double v_y = v[y];
        for( int x = 0; x < cols/2; x++){
            HPF[(y+rens/2)*cols + x+cols/2] = 1.0 - exp( -( sy*sy*(u[x]*u[x] + v_y*v_y) )/2.0 );
        }
        for( int x = cols/2; x < cols; x++){
            HPF[(y+rens/2)*cols + x-cols/2] = 1.0 - exp( -( sy*sy*(u[x]*u[x] + v_y*v_y) )/2.0 );
        }
    }
    for( int y = rens/2; y < rens; y++ ){
        const double v_y = v[y];
        for( int x = 0; x < cols/2; x++){
            HPF[(y-rens/2)*cols + x+cols/2] = 1.0 - exp( -( sy*sy*(u[x]*u[x] + v_y*v_y) )/2.0 );
        }
        for( int x = cols/2; x < cols; x++){
            HPF[(y-rens/2)*cols + x-cols/2] = 1.0 - exp( -( sy*sy*(u[x]*u[x] + v_y*v_y) )/2.0 );
        }
    }

    // Multiplicar la imagen transformada por HPF:
    for(int y = 0; y < rens; y++){
        for( int x = 0; x <= cols/2; x++){
            Img_fft[x+y*(cols/2+1)][0] *= HPF[x+y*cols];
            Img_fft[x+y*(cols/2+1)][1] *= HPF[x+y*cols];
        }
    }
    delete [] HPF;

    // Aplicar los filtros de Gabor a diferentes angulos:
    fftw_complex *Img_filtro = (fftw_complex*) fftw_malloc(rens*(cols/2+1)*sizeof(fftw_complex));
    double *Img_resp = new double[rens_cols];
    double *max_resp = new double[rens_cols];
    memset(Img_resp, 0, rens_cols*sizeof(double));
    memset(max_resp, 0, rens_cols*sizeof(double));


    fftw_plan p_c2r;
    const double theta_inc = 180.0 / (double)K;

    double Gabor, Vr, Ur;

    for( double theta = 0.0; theta < 180.0; theta+=theta_inc){
        //// Se rota la mascara del filtro en theta grados:
        const double stheta = sin(theta*PI/180.0);
        const double ctheta = cos(theta*PI/180.0);

        // Se calcula el filtro de Gabor y se hace el shift similar a fftshift en Matlab:
        for( int y = 0; y < rens/2; y++){
            const double v_y = v[y+rens/2];
            for( int x = 0; x < cols/2; x++){
                Ur = u[x+cols/2]*ctheta + v_y*stheta;
                Vr =-u[x+cols/2]*stheta + v_y*ctheta;
                Gabor = exp(-(0.5)*(sx*sx*(Ur*Ur + 4.0*PI*fx*PI*fx) + sy*sy*(Vr*Vr + 4.0*PI*fy*PI*fy))) * cosh(2.0*PI*(sx*sx*fx*Ur + sy*sy*fy*Vr));
                Img_filtro[y*(cols/2+1) + x][0] = Img_fft[y*(cols/2+1) + x][0] * Gabor;
                Img_filtro[y*(cols/2+1) + x][1] = Img_fft[y*(cols/2+1) + x][1] * Gabor;
            }
            Ur = u[0]*ctheta + v_y*stheta;
            Vr =-u[0]*stheta + v_y*ctheta;
            Gabor = exp(-(0.5)*(sx*sx*(Ur*Ur + 4.0*PI*fx*PI*fx) + sy*sy*(Vr*Vr + 4.0*PI*fy*PI*fy))) * cosh(2.0*PI*(sx*sx*fx*Ur + sy*sy*fy*Vr));
            Img_filtro[y*(cols/2+1) + cols/2][0] = Img_fft[y*(cols/2+1) + cols/2][0] * Gabor;
            Img_filtro[y*(cols/2+1) + cols/2][1] = Img_fft[y*(cols/2+1) + cols/2][1] * Gabor;
        }
        for( int y = rens/2; y < rens; y++){
            const double v_y = v[y-rens/2];
            for( int x = 0; x < cols/2; x++){
                Ur = u[x+cols/2]*ctheta + v_y*stheta;
                Vr =-u[x+cols/2]*stheta + v_y*ctheta;
                Gabor = exp(-(0.5)*(sx*sx*(Ur*Ur + 4.0*PI*fx*PI*fx) + sy*sy*(Vr*Vr + 4.0*PI*fy*PI*fy))) * cosh(2.0*PI*(sx*sx*fx*Ur + sy*sy*fy*Vr));
                Img_filtro[y*(cols/2+1) + x][0] = Img_fft[y*(cols/2+1) + x][0] * Gabor;
                Img_filtro[y*(cols/2+1) + x][1] = Img_fft[y*(cols/2+1) + x][1] * Gabor;

            }
            Ur = u[0]*ctheta + v_y*stheta;
            Vr =-u[0]*stheta + v_y*ctheta;
            Gabor = exp(-(0.5)*(sx*sx*(Ur*Ur + 4.0*PI*fx*PI*fx) + sy*sy*(Vr*Vr + 4.0*PI*fy*PI*fy))) * cosh(2.0*PI*(sx*sx*fx*Ur + sy*sy*fy*Vr));
            Img_filtro[y*(cols/2+1) + cols/2][0] = Img_fft[y*(cols/2+1) + cols/2][0] * Gabor;
            Img_filtro[y*(cols/2+1) + cols/2][1] = Img_fft[y*(cols/2+1) + cols/2][1] * Gabor;
        }

        //// Ap[licar la transformada inversa de fourier:
        p_c2r = fftw_plan_dft_c2r_2d(rens, cols, Img_filtro, Img_resp, FFTW_ESTIMATE);
        fftw_execute(p_c2r);
        fftw_destroy_plan(p_c2r);

        //// Almacenar la respuesta maxima encontrada:
        for( int xy = 0; xy < rens_cols; xy++){
            if( max_resp[xy] < Img_resp[xy]){
                max_resp[xy] = Img_resp[xy];
            }
        }
    }


    fftw_free(Img_filtro);
    delete [] Img_resp;
    delete [] u;
    delete [] v;

    for( int xy = 0; xy < rens_cols; xy++){
        // La inversa de la FFT retorna los valores multiplicados por el tamaño del arreglo que se mete como argumento, por eso se divide entre el tamaño del arreglo.
        resp[xy] = (mask[xy] > 125) ? (max_resp[xy] / rens_cols) : 0.0;
        //resp[xy] = (max_resp[xy] / rens_cols);
    }

    delete [] max_resp;
    GETTIME_FIN;

    DEB_MSG("Tiempo de filtrado: " << DIFTIME << "s." );
}





/*  Metodo: comp_resp

    Funcion: Compara la respuesta de dos pixeles.
*/
int comp_resp( const void *a, const void *b){
    double resp_a = *(double*)a;
    double resp_b = *(double*)b;
    if( resp_a > resp_b ){
        return 1;
    }else if(resp_a < resp_b){
        return -1;
    }else{
        return 0;
    }
}





/*  Metodo: calcROC

    Funcion: Calcula la 'receiver operating characteristic curve' ROC curve y el area bajo la curva, para prueba en espacio continuo.
*/
double FILTROS::calcROC(INDIV *test, double *resp){

    TIMERS;

    GETTIME_INI;

    const double delta = test->vars[4];
    DEB_MSG("Calculando ROC");
    //// Contar la cantidad de pixeles blancos (unos) y negros (ceros)
    int n_unos = 0, n_ceros = 0;
    for( int xy = 0; xy < rens_cols; xy++){
        if( ground_truth[xy] > 125 ){
            n_unos++;
        }else{
            n_ceros++;
        }

    }


    //// Pasar las respuestas a los arreglos correspondientes a la respuesta en el Ground truth
    double *arr_unos = new double[n_unos];
    double *arr_ceros = new double[n_ceros];
    n_unos = 0;
    n_ceros = 0;
    double max_resp = 0.0;
    double min_resp = 1e100;

    for( int xy = 0; xy < rens_cols; xy++){
        if( ground_truth[xy] > 125 ){
            arr_unos[n_unos] = resp[xy];
            n_unos++;
        }else{
            arr_ceros[n_ceros] = resp[xy];
            n_ceros++;
        }

        if(max_resp < resp[xy]){
            max_resp = resp[xy];
        }
        if(min_resp > resp[xy]){
            min_resp = resp[xy];
        }
    }


    //// Ordenar ambos arreglos de acuerdo a la respuesta (de la mas alta a la mas baja):
    qsort(arr_unos, n_unos, sizeof(double), comp_resp);
    qsort(arr_ceros, n_ceros, sizeof(double), comp_resp);

    //// Comienzar con el threshold mas discriminativo:
    double threshold = max_resp;
    int i_unos = (n_unos-1), i_ceros = (n_ceros-1), i_threshold = 0;
    double TP = 0.0, TN = (double)n_ceros;

    double curva_TPR_old = 0.0;
    double curva_FPR_old = 0.0;

    double curva_TPR_new, curva_FPR_new;
    double Az = 0.0;

    // Continuar con el thresholding acumulando el area bajo cada punto en Az:
    min_resp -= delta;
    while( threshold > min_resp ){
        //// Contar los verdaderos positivos:
        while( i_unos > 0 ){
            if(arr_unos[i_unos] >= threshold){
                TP+=1.0;
                i_unos--;
            }else{
                break;
            }
        }

        while( i_ceros > 0){
            if( arr_ceros[i_ceros] > threshold){
                TN-=1.0;
                i_ceros--;
            }else{
                break;
            }
        }

        //// Almacenar la curva en este arreglo pasado como argumento:
        curva_TPR_new = TP/(double)n_unos;
        curva_FPR_new = 1.0 - TN/(double)n_ceros;

        // Calcular el area con la Regla del Trapecio:
        Az += (curva_FPR_new - curva_FPR_old)*(curva_TPR_new + curva_TPR_old)/2.0;

        threshold-=delta;
        i_threshold++;

        // Mover los valores de la curva:
        curva_TPR_old = curva_TPR_new;
        curva_FPR_old = curva_FPR_new;

        //DEB_MSG("[" << i_threshold << "]: " << threshold << " => " << curva_TPR_old << ", " << curva_FPR_old)

    }

    delete [] arr_unos;
    delete [] arr_ceros;

    GETTIME_FIN;

    DEB_MSG("Tiempo para obtener la curva de ROC: " << DIFTIME << " s.  " << Az)

    return Az;
}



/*  Metodo: fitnessROC
    Funcion: Evalua los parametros 'L', 'T' y //'K'// para el filtro de Gabor y la curva de ROC:
*/
double FILTROS::fitnessROC( INDIV *test ){
    //// Generar la respuesta del filtro de establecido para los parametros dados:
    double *resp = new double [rens_cols];

    switch(filtro_elegido){
        case SS_GABOR:
            respGabor(test, resp);
            break;
        case GMF:
            respGMF(test, resp);
            break;
    }
    double Az = calcROC(test, resp);
    delete [] resp;
    return Az;
}






/*//        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
 *
 *              G E N E R A D O R E S       D E     N U M E R O S       A L E A T O R I O S
 *
 * //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
*/


/*	Metodo:         lcg_s
    Funcion: Genera un numero pseudo-aleatorio por medio del metodo Linear Congruential Generator. El generador utiliza una semilla static y los parametros que utiliza son los mismos que utiliza gcc.
*/
unsigned int FILTROS::lcg_s(){
    static unsigned int st_semilla = 0;
    if(!st_semilla){
        st_semilla = rand();
    }
    const unsigned int m = pow(2,31);
    return (1103515245u*st_semilla + 12345) % m;
}




/* Metodo:        lcg_r
    Funcion: Genera un numero pseudo-aleatorio por medio del metodo Linear Congruential Generator. El generador utiliza una semilla enviada por referencia, la cual actualiza y los parametros que utiliza son los mismos que utiliza gcc.
*/
unsigned int FILTROS::lcg_r(unsigned int *mi_semilla){
    const unsigned int m = pow(2,31);
    return *mi_semilla= (1103515245u* *mi_semilla + 12345) % m;
}




/*	Metodo:        tausStep
    Funcion: Genera un numero pseudo-aleatorio por medio del metodo Taus Step.
*/
unsigned int FILTROS::tausStep(unsigned int *z, const int S1, const int S2, const int S3, const unsigned int M){
    const unsigned int b = (((*z << S1) ^ *z) >> S2);
    return *z = (((*z & M) << S3) ^ b);
}




/*	Metodo         anorm_est
    Funcion: Genera un número aleatorio de la distribución normal estándar, por el método de transformación Box-Muller.
*/
double FILTROS::anorm_est(){
    double x1, x2, w, y1;//, y2;
    do{
        x1 = HybTaus(-1.0, 1.0);
        x2 = HybTaus(-1.0, 1.0);
        w = x1*x1 + x2*x2;
    }while( w >= 1.0);
    w = sqrt( (-2.0 * log(w)) / w);
    y1 = x1 * w;
//	y2 = x2 * w;
    return y1;
}




/*//        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
 *
 *              A L G O R I T M O S     E V O L U T I V O S
 *
 * //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
*/







/*
    Metodo:        compIndiv
    Funcion:    Compara el indidividuo A contra el B.
        Si la evaluacion del individuo A es menor a la del B, se retorna +,
                             si es igual, se retorna 0,
                             si es mayor, se retorna -.
*/
int FILTROS::compIndiv(const void* A, const void* B){
    INDIV *indA = (INDIV*)A;
    INDIV *indB = (INDIV*)B;
    if(indA->eval < indB->eval) return  1;
    if(indA->eval == indB->eval) return  0;
    if(indA->eval > indB->eval) return -1;
}




/*
    Metodo:        generarPobInicial
    Funcion:    Genera una poblacion inicial dentro de los limites establecidos.

*/
void FILTROS::generarPobInicial(INDIV *poblacion){
    TIMERS;
    //Se generan n_pob - 1 individuos, y permanece el elite que se tiene.
    // VERIFICAR CUALES PARAMETROS SE BUSCAN:
    using namespace std;
    for(int i = 0; i < n_pob; i++){
        GETTIME_INI;
        memcpy( poblacion[i].vars, mi_elite->vars, 5*sizeof(double) );

        for( int j = 0; j < n_pars; j++){
            const int k = idx_pars[j];
            poblacion[i].vars[k] = HybTaus(lim_inf[k], lim_sup[k]);
        }

        //// Evaluar los parametros en el filtro de Gabor con el area bajo de la curva de ROC:
        switch( fitness_elegido ){
            case ROC:
                poblacion[i].eval = fitnessROC( &poblacion[i] );
                break;
        }
        DEB_MSG("[" << omp_get_thread_num() << "/" << omp_get_num_threads() << ":" << i << "] " << poblacion[i].eval << ": " << poblacion[i].vars[0] << ", " << poblacion[i].vars[1] << ", " << poblacion[i].vars[2] << ", " << poblacion[i].vars[3] << ", " << poblacion[i].vars[4])
        GETTIME_FIN;
        cout << "Individuo " << i << " generado y evaluado en " << DIFTIME << " segundos." << endl;
    }

    //Se reordena toda la poblacion, no solo la parte generada.
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv);
}





// BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA
/*
    Metodo:        generarPob
    Funcion:    Genera una poblacion de tamaño ('n_pob' - 'ini') individuos con una distribucion indicada con 'fun_aln'. Evalua cada individuo conforme es evaluado en la funcion objetivo. Al final se ordenen de menor a mayor (porque se busca un minimo).

*/
void FILTROS::generarPob(INDIV *poblacion, const int n_gen, double medias[5], double varianzas[5]){
    TIMERS;
    double val_gen;

    //Se generan n_pob - 1 individuos, y permanece el elite que se tiene.
    // VERIFICAR CUALES PARAMETROS SE BUSCAN:
    using namespace std;
    for(register int i = n_gen; i < n_pob; i++){
        GETTIME_INI;
        for( int j = 0; j < n_pars; j++){
            const int k = idx_pars[j];
            val_gen =  anorm( medias[k], varianzas[k]);
            val_gen = (val_gen <= lim_sup[k]) ? ((val_gen >= lim_inf[k]) ? val_gen : lim_inf[k]) : lim_sup[k];
            poblacion[i].vars[k] = val_gen;
        }
        //// Evaluar los parametros en el filtro de Gabor con el area bajo de la curva de ROC:
        switch( fitness_elegido ){
            case ROC:
                poblacion[i].eval = fitnessROC( &poblacion[i] );
                break;
        }

        DEB_MSG("[" << i << "] " << poblacion[i].eval << ": " << poblacion[i].vars[0] << ", " << poblacion[i].vars[1] << ", " << poblacion[i].vars[2] << ", " << poblacion[i].vars[3] << ", " << poblacion[i].vars[4])
        GETTIME_FIN;
        cout << "Individuo " << i << " generado y evaluado en " << DIFTIME << " segundos. f(x) = " << poblacion[i].eval << "[L = " << poblacion[i].vars[0] << ", T = " << poblacion[i].vars[1] << "]" << endl;
    }

    //Se reordena toda la poblacion, no solo la parte generada.
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv);
}



/*
    Metodo:        calcularPars
    Funcion:    Calula la media y la varianza para cada variable de la poblacion, utilizando solo a los individuos de la pobllacion.
*/
void FILTROS::calcularPars(INDIV *poblacion, const int truncamiento, double *medias, double *varianzas){

    double sum_evals = 0.0;
    double gx_sel = poblacion[truncamiento].eval - 1e-5;
        //Se calculan las medias primero:
        memset(medias, 0, 5*sizeof(double));

        for( int i = 0; i < truncamiento; i++){//Este es el ciclo para toda la seleccion.
            const double g_bar = poblacion[i].eval - gx_sel;
            sum_evals += g_bar;
            for( int j = 0; j < n_pars; j++){
                const int k = idx_pars[j];
                medias[k] += poblacion[i].vars[k] * g_bar;
            }
        }

        for( int j = 0; j < n_pars; j++){
            medias[ idx_pars[j] ] /= sum_evals;
        }

        //Se calculan las varianzas:
        double tmp;
        memset(varianzas, 0, 5*sizeof(double));
        for( int i = 0; i < truncamiento; i++){
            const double g_bar = poblacion[i].eval - gx_sel;
            for( int j = 0; j < n_pars; j++){
                const int k = idx_pars[j];
                tmp = (poblacion[i].vars[k] - medias[k]);
                varianzas[k] += tmp*tmp * g_bar;
            }
        }

        for( int j = 0; j < n_pars; j++){
            varianzas[ idx_pars[j] ] /= 1.0 + sum_evals;
        }
}



/*
    Metodo:        seleccionarPob
    Funcion:    Realiza la seleccion por truncamiento, redefiniendo el tetha para esta iteracion dentro del mismo apuntador de tetha.
*/
int FILTROS::seleccionarPob(double *tetha_t, INDIV *poblacion){
    int truncamiento = 0;

    for(register int i = n_pob/2; i >= 0; i--){
        // Se trunca hasta el individuo con evaluacion en la funcion objetivo encima del tetha anterior.
        if(poblacion[i].eval > (*tetha_t)){
            truncamiento = i;
            break;
        }
    }

//    //Si todos los individuos son superiores al theta anterior
//    if(truncamiento == (n_pob-1)){
//        truncamiento = n_pob / 2;
//    }

    //Se asigna el nuevo tetha como la evaluacion del individuo en la posicion 'truncamiento'
    *tetha_t = poblacion[truncamiento].eval;

    return truncamiento;
}



/*
    Metodo:        BUMDA
    Funcion:    Utiliza el algoritmo Boltzmann Univariate Marginal Distribution para encotnrar los parametros automaticamente.
*/
void FILTROS::BUMDA(){

    TIMERS;
    GETTIME_INI;

    INDIV *poblacion = new INDIV [n_pob];

    double *medias = new double [5];
    double *varianzas = new double [5];

    // Se cuantan los parametros activos:
    n_pars = 0;
    for( int p = 0; p < 5; p++){
        if( pars_optim[p] ){
            idx_pars[n_pars] = p;
            n_pars ++;
        }
    }

    if(semilla){
        delete semilla;
    }
    semilla = ini_semilla(0);

    // Inicia el algoritmo BUMDA
    int k = 0, truncamiento = n_pob-1, procesar = 1;

    //// Generar la primer poblacion:
    DEB_MSG("Generando " << n_pob << " individuos como poblacion inicial ...")
    generarPobInicial(poblacion);
    DEB_MSG("Poblacion inicial generada, comenzando BUMDA ...");

    //Se define el theta en el tiempo 0 como el minimo de la poblacion inicial.
    double tetha_t = poblacion[n_pob-1].eval;
    using namespace std;
    do{
        calcularPars(poblacion, truncamiento, medias, varianzas);
        generarPob(poblacion, 1, medias, varianzas);

        //En truncamiento se guarda el indice del individuo a partir del cual se trunca la poblacion.
        truncamiento = seleccionarPob(&tetha_t, poblacion);
        k++;

        //// Verificar condiciones de paro:
        int condicion = 0;
        for(int j = 0; j < n_pars; j++){
            const int k = idx_pars[j];
            condicion += (varianzas[ k ] <= min_vars[ k ]);
        }
        procesar = (k < max_iters) && (condicion < n_pars);

        cout << "[" << k << "/" << max_iters << "] Best fit: " << poblacion[0].eval << ", L: " << poblacion[0].vars[0] << ", T: " << poblacion[0].vars[1] << endl;

    }while(procesar);

    delete [] medias;
    delete [] varianzas;

    memcpy(mi_elite->vars, poblacion[0].vars, 5*sizeof(double));
    mi_elite->eval = poblacion[0].eval;
    delete [] poblacion;
    GETTIME_FIN;
    cout << " Mejor candidato encontrada con " << k << " iteraciones y " << DIFTIME << "segundos.\n Con evaluacion de la funcion objetivo: " << mi_elite->eval << endl;
}




// UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA

/*
    Metodo:     generar_pob
    Funcion:   Genera una nueva poblacion en bsae al muestreo de cada atributo.
*/
void FILTROS::generarPob(INDIV *poblacion, const int n_gen, double *probs, double **tabla){
    TIMERS;
    using namespace std;

    //Se generan n_pob - 1 individuos, y permanece el elite que se tiene.
    // VERIFICAR CUALES PARAMETROS SE BUSCAN:
    for(register int i = n_gen; i < n_pob; i++){
        GETTIME_INI;
        // Se muestrean los genes:
        for( int b = 0; b < n_bits; b++){
            poblacion[i].cadena[b] = (HybTaus(0.0, 1.0) <= probs[b]) ? true : false;
        }

        // Decodificar la cadena en los valores de los parametros:
        int idx_bit = 0;
        for( int j = 0; j < n_pars; j++){
            const int k = idx_pars[j];
            unsigned int bit = 1, cadena_tmp = 0;
            for( int b = 0; b < bits_var[k]; b++){
                cadena_tmp |= poblacion[i].cadena[idx_bit] ? bit : 0;
                bit = bit << 1;
                idx_bit++;
            }

            // Asignar el valor segun la cadena formada para el atributo 'k':
            poblacion[i].vars[k] = tabla[k][ (int)cadena_tmp ];
        }

        //// Evaluar los parametros en el filtro de Gabor con el area bajo de la curva de ROC:
        switch( fitness_elegido ){
            case ROC:
                poblacion[i].eval = fitnessROC( &poblacion[i] );
                break;
        }

        DEB_MSG("[" << i << "] " << poblacion[i].eval << ": " << poblacion[i].vars[0] << ", " << poblacion[i].vars[1] << ", " << poblacion[i].vars[2] << ", " << poblacion[i].vars[3] << ", " << poblacion[i].vars[4])
        GETTIME_FIN;
        cout << "Individuo " << i << " generado y evaluado en " << DIFTIME << " segundos." << endl;
    }

    DEB_MSG("Reordenando poblacion ...");
    //Se reordena toda la poblacion, no solo la parte generada.
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv);

    DEB_MSG("Poblacion rordenada");
}




/*
    Metodo:     calcular_pars
    Funcion:    Calcula las probabilidades de cada distribucion marginal.
*/
void FILTROS::calcularPars(INDIV *poblacion, const int truncamiento, double *probs){
    // Se estiman las probabilidades:
    memset(probs, 0, n_bits*sizeof(double));
    for(int i = 0; i < truncamiento; i++){//Este es el ciclo para toda la seleccion.
        for( int b = 0; b < n_bits; b++){
            probs[b] += poblacion[i].cadena[b] ? 1.0 : 0.0;
        }
    }

    for( int b = 0; b < n_bits; b++){
        probs[b] /= (double)truncamiento;
    }
}





/*
    Metodo:     UMDA
    Funcion:    Utiliza el algoritmo Univariate Marginal Distribution para encotnrar los parametros automaticamente.
*/
void FILTROS::UMDA(){

    TIMERS;
    GETTIME_INI;

    INDIV *poblacion = new INDIV [n_pob];
    // Poner los valores del elite por defecto:
    for( int i = 0; i < n_pob; i++){
        memcpy(poblacion[i].vars, mi_elite->vars, 5*sizeof(double));
    }

    // Se cuentan los parametros activos:
    n_pars = 0;
    for( int p = 0; p < 5; p++){
        if( pars_optim[p] ){
            idx_pars[n_pars] = p;
            n_pars ++;
        }
    }

    // Se cuentan el numero total de bits que se van a usar:
    n_bits = 0;
    for( int j = 0; j < n_pars; j++){
        const int k = idx_pars[j];
        bits_var[k] = (unsigned char)log2(max_bits[k]);
        n_bits += bits_var[k];
    }

    double *probs = new double [n_bits];
    // Inicializar las probabilidades en 0.5:
    for( int b = 0; b < n_bits; b++){
        probs[b] = 0.5;
    }

    // Generar la tabla de valores:
    double *tabla[5];
    for( int j = 0; j < n_pars; j++){
        const int k = idx_pars[j];
        tabla[k] = new double [(int)max_bits[k]];
        const double incremento = (lim_sup[k] - lim_inf[k]) / ((int)max_bits[k] - 1);
        for( int b = 0; b < (int)max_bits[k]; b++){
            tabla[k][b] = lim_inf[k] + b*incremento;
        }
    }

    if(semilla){
        delete semilla;
    }
    semilla = ini_semilla(0);

    // Inicia el algoritmo UMDA
    int k = 0, truncamiento = (int)( (double)n_pob * 0.6), procesar = 1;

    //// Generar la primer poblacion:
    DEB_MSG("Generando " << n_pob << " individuos como poblacion inicial ...")
    generarPob( poblacion, 0, probs, tabla);
    DEB_MSG("Poblacion inicial generada, comenzando UMDA ...");

    using namespace std;
    do{
        calcularPars(poblacion, truncamiento, probs);
        generarPob(poblacion, 1, probs, tabla);
        k++;

        //// Verificar condiciones de paro:
        procesar = (k < max_iters);

        cout << "[" << k << "/" << max_iters << "] Best fit: " << poblacion[0].eval << ", L: " << poblacion[0].vars[0] << ", T: " << poblacion[0].vars[1] << endl;

    }while(procesar);

    delete [] probs;

    for(int j = 0; j < n_pars; j++){
        delete [] tabla[ idx_pars[j] ];
    }

    memcpy(mi_elite->vars, poblacion[0].vars, 5*sizeof(double));
    mi_elite->eval = poblacion[0].eval;
    delete [] poblacion;
    GETTIME_FIN;
    cout << " Mejor candidato encontrado con " << k << " iteraciones. y " << DIFTIME << "segundos. \n Con evaluacion de la funcion objetivo: " << mi_elite->eval << endl;
}




// GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA
/*
    Metodo:     generar_pob
    Funcion:   Genera una nueva poblacion en bsae al muestreo de cada atributo.
*/
void FILTROS::generarPobInicial(INDIV *poblacion, double **tabla){
    using namespace std;

    //Se generan todos los individuos de la poblacion con genes elegidos aleatoriamente.
    // VERIFICAR CUALES PARAMETROS SE BUSCAN:
    for(int i = 0; i < n_pob; i++){
        // Se muestrean los genes aleatoriamente:
        for( int b = 0; b < n_bits; b++){
            poblacion[i].cadena[b] = (HybTaus(0.0, 1.0) <= 0.5) ? true : false;
            DEB_MSG("cadena [" << b << "]: " << poblacion[i].cadena[b]);

        }

        // Decodificar la cadena en los valores de los parametros:
        int idx_bit = 0;
        for( int j = 0; j < n_pars; j++){
            const int k = idx_pars[j];
            unsigned int bit = 1, cadena_tmp = 0;
            for( int b = 0; b < bits_var[k]; b++){
                cadena_tmp |= poblacion[i].cadena[idx_bit] ? bit : 0;
                bit = bit << 1;
                idx_bit++;
            }
            // Asignar el valor segun la cadena formada para el atributo 'k':
            poblacion[i].vars[k] = tabla[k][ (int)cadena_tmp ];
        }
    }

    // Evaluacion de los individuos:
    switch( fitness_elegido ){
        case ROC:
            //#pragma omp parallel for shared(poblacion)
            for(int i = 0; i < n_pob; i++){
                TIMERS;
                GETTIME_INI;
                //// Evaluar los parametros en el filtro de Gabor con el area bajo de la curva de ROC:
                poblacion[i].eval = fitnessROC( &poblacion[i] );

                DEB_MSG("[PROCESS_" << omp_get_thread_num() << "][" << i << "] " << poblacion[i].eval << ": " << poblacion[i].vars[0] << ", " << poblacion[i].vars[1] << ", " << poblacion[i].vars[2] << ", " << poblacion[i].vars[3] << ", " << poblacion[i].vars[4])
                GETTIME_FIN;
                cout << "Individuo " << i << " generado y evaluado en " << DIFTIME << " segundos." << endl;
            }
            break;
    }
    //Se reordena toda la poblacion, no solo la parte generada.
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv);

}




/*
    Metodo:    seleccionPob
    Funcion:   Selecciona individuos que pasan a ser los padres dela siguiente generacion mediante el metodo de la ruleta.
*/
void FILTROS::seleccionarPob(INDIV* poblacion, INDIV* probs, INDIV* pob_tmp, int *seleccion){
    // Sumar todos los fitnes, a la par, generar las probabilidades de la ruleta para la seleccion de individuos:
    double sum_fitness = 0.0;
    for(int i = 0; i < n_pob; i++){
        sum_fitness += poblacion[i].eval;
        probs[i].eval = HybTaus(0.0, 1.0); // Se genera la probabilidad de la ruleta uniformemente.
        probs[i].vars[0] = (double)i; // Se guarda el indice donde se genera la probabilidad
        memcpy(pob_tmp[i].vars, poblacion[i].vars, 5*sizeof(double));
        memcpy(pob_tmp[i].cadena, poblacion[i].cadena, n_bits*sizeof(bool));
        pob_tmp[i].eval = poblacion[i].eval;
    }

    DEB_MSG("Probabilidades de la ruleta calculadas, el fitness global es: " << sum_fitness);

    // Se requiere reordenar las probabilidades para elegir los padres a partir de las probabilidades acumuladas:
    qsort((void*)probs, n_pob, sizeof(INDIV), compIndiv);
    DEB_MSG("Probabilidades ordenadas: " << probs[0].eval << " === " << probs[n_pob-1].eval);

    int i_selec = n_pob-1;

    // Calcular la probabilidad relativa y acumulada de cada individuo, y determinar los individuos que pasan al grupo seleccionado.
    double acum_fitness = 0.0;
    double nuevo_acum;
    for(int i = n_pob-1; i >= 0; i--){
        nuevo_acum = acum_fitness + poblacion[i].eval / sum_fitness;
        DEB_MSG("Fitness acumulado: " << nuevo_acum);
        while( (i_selec >= 0) && (nuevo_acum > probs[i_selec].eval) && (probs[i_selec].eval > acum_fitness) ){
            seleccion[i_selec] = (int)probs[i_selec].vars[0];
            DEB_MSG("       Probabilidad: " << probs[i_selec].eval);
            i_selec--;
        }
        acum_fitness = nuevo_acum;
    }
    DEB_MSG("faltan por seleccionar: " << i_selec+1);

    // Pasar la seleccion a la poblacion:
    for(int i = 0; i < n_pob; i++){
        memcpy(poblacion[i].vars, pob_tmp[seleccion[i]].vars, 5*sizeof(double));
        memcpy(poblacion[i].cadena, pob_tmp[seleccion[i]].cadena, n_bits*sizeof(bool));
        poblacion[i].eval = pob_tmp[seleccion[i]].eval;
    }
}




/*
    Metodo:     generarPob
    Funcion:    Calcula la siguiente generacion a partir de la cruza multipunto de los padres seleccionados, ademas se realiza una mutacion sobre un porcentage de los individuos.
*/
void FILTROS::generarPob(INDIV* poblacion, const double prob_mutacion, double **tabla){

    TIMERS;

    // Generar dos nuevos individuos a partir de dos padres:
    for(int i = 0; i < n_pob; i+=2){
        INDIV hijo_1, hijo_2;
        int inicio = 0;
        int padre_1 = 0;
        // Realizar cortes hasta terminar con la secuencia de los individuos:

        while(inicio < n_bits){
            int corte = (int)HybTaus((double)inicio, (double)n_bits - 1);
            DEB_MSG("inicio: " << inicio << ", corte: " << corte);
            for(int b = inicio; b <= corte; b++){
                hijo_1.cadena[b] = poblacion[i + padre_1].cadena[b];
                hijo_2.cadena[b] = poblacion[i + 1 - padre_1].cadena[b];
                padre_1 = 1 - padre_1;
                DEB_MSG("hijo 1 [" << b << "]: " << hijo_1.cadena[b]);
                DEB_MSG("hijo 2 [" << b << "]: " << hijo_2.cadena[b]);
            }
            inicio = corte + 1;
        }

        GETTIME_INI;
        // Realizar mutacion sobre los hijos generados con cierta probabilidad:
        // Mutar el hijo 1:
        if( HybTaus(0.0, 1.0) <= prob_mutacion ){
            // Seleccionar un gen a mutar:
            const int gen = (int)(HybTaus(-1.0, (double)n_bits-1.0))+1;
            hijo_1.cadena[gen] != hijo_1.cadena[gen];
        }

        // Evaluar las cadenas de los nuevos individuos:
        //----------------------------------------------------------------- Hijo 1
        // Decodificar la cadena en los valores de los parametros:
        int idx_bit = 0;
        for( int j = 0; j < n_pars; j++){
            const int k = idx_pars[j];
            unsigned int bit = 1, cadena_tmp = 0;
            for( int b = 0; b < bits_var[k]; b++){
                cadena_tmp |= hijo_1.cadena[idx_bit] ? bit : 0;
                poblacion[i].cadena[idx_bit] = hijo_1.cadena[idx_bit];
                bit = bit << 1;
                idx_bit++;
            }

            // Asignar el valor segun la cadena formada para el atributo 'k':
            poblacion[i].vars[k] = tabla[k][ (int)cadena_tmp ];
        }

        //// Evaluar los parametros en el filtro de Gabor con el area bajo de la curva de ROC:
        switch( fitness_elegido ){
            case ROC:
                poblacion[i].eval = fitnessROC( &poblacion[i] );
                break;
        }

        DEB_MSG("[" << i << "] " << poblacion[i].eval << ": " << poblacion[i].vars[0] << ", " << poblacion[i].vars[1] << ", " << poblacion[i].vars[2] << ", " << poblacion[i].vars[3] << ", " << poblacion[i].vars[4])
        GETTIME_FIN;
        cout << "Individuo " << i << " generado y evaluado en " << DIFTIME << " segundos." << endl;


        GETTIME_INI;
        // Mutar el hijo 2:
        if( HybTaus(0.0, 1.0) <= prob_mutacion ){
            // Seleccionar un gen a mutar:
            const int gen = (int)(HybTaus(-1.0, (double)n_bits-1.0))+1;
            hijo_2.cadena[gen] != hijo_2.cadena[gen];
        }

        //----------------------------------------------------------------- Hijo 2
        // Decodificar la cadena en los valores de los parametros:
        idx_bit = 0;
        for( int j = 0; j < n_pars; j++){
            const int k = idx_pars[j];
            unsigned int bit = 1, cadena_tmp = 0;
            for( int b = 0; b < bits_var[k]; b++){
                cadena_tmp |= hijo_2.cadena[idx_bit] ? bit : 0;
                poblacion[i+1].cadena[idx_bit] = hijo_2.cadena[idx_bit];
                bit = bit << 1;
                idx_bit++;
            }

            // Asignar el valor segun la cadena formada para el atributo 'k':
            poblacion[i+1].vars[k] = tabla[k][ (int)cadena_tmp ];
        }

        //// Evaluar los parametros en el filtro de Gabor con el area bajo de la curva de ROC:
        switch( fitness_elegido ){
            case ROC:
                poblacion[i+1].eval = fitnessROC( &poblacion[i+1] );
                break;
        }

        DEB_MSG("[" << i+1 << "] " << poblacion[i+1].eval << ": " << poblacion[i+1].vars[0] << ", " << poblacion[i+1].vars[1] << ", " << poblacion[i+1].vars[2] << ", " << poblacion[i+1].vars[3] << ", " << poblacion[i+1].vars[4])
        GETTIME_FIN;
        cout << "Individuo " << i+1 << " generado y evaluado en " << DIFTIME << " segundos." << endl;
    }


    //Se reordena toda la poblacion, no solo la parte generada.
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv);
}




/*
    Metodo:     GA
    Funcion:    Utiliza el algoritmo genetico para encotnrar los parametros automaticamente.
*/
void FILTROS::GA(){

    TIMERS;
    GETTIME_INI;

    INDIV *poblacion = new INDIV [n_pob];
    INDIV *probs = new INDIV [n_pob];
    INDIV *pob_tmp = new INDIV [n_pob];

    // Poner los valores del elite por defecto:
    for( int i = 0; i < n_pob; i++){
        memcpy(poblacion[i].vars, mi_elite->vars, 5*sizeof(double));
    }

    // Se cuentan los parametros activos:
    n_pars = 0;
    for( int p = 0; p < 5; p++){
        if( pars_optim[p] ){
            idx_pars[n_pars] = p;
            n_pars ++;
        }
    }

    // Se cuentan el numero total de bits que se van a usar:
    n_bits = 0;
    for( int j = 0; j < n_pars; j++){
        const int k = idx_pars[j];
        bits_var[k] = (unsigned char)log2(max_bits[k]);
        n_bits += bits_var[k];
    }


    // Generar la tabla de valores:
    double *tabla[5];
    for( int j = 0; j < n_pars; j++){
        const int k = idx_pars[j];
        tabla[k] = new double [(int)max_bits[k]];
        const double incremento = (lim_sup[k] - lim_inf[k]) / ((int)max_bits[k] - 1);
        for( int b = 0; b < (int)max_bits[k]; b++){
            tabla[k][b] = lim_inf[k] + b*incremento;
        }
    }

    if(semilla){
        delete semilla;
    }
    semilla = ini_semilla(0);

    int *seleccion = new int [n_pob];

    // Inicia el algoritmo GA
    int k = 0, procesar = 1;

    //// Generar la primer poblacion:
    DEB_MSG("Generando " << n_pob << " individuos como poblacion inicial ...")
    generarPobInicial(poblacion, tabla);
    DEB_MSG("Poblacion inicial generada, comenzando GA ...");

    mi_elite->eval = 0.0;

    using namespace std;
    do{
        seleccionarPob(poblacion, probs, pob_tmp, seleccion);

        generarPob(poblacion, 0.05, tabla);

        k++;

        //// Verificar condiciones de paro:
        procesar = (k < max_iters);

        cout << "[" << k << "/" << max_iters << "] Best fit: " << poblacion[0].eval << ", L: " << poblacion[0].vars[0] << ", T: " << poblacion[0].vars[1] << endl;

        // Extraer el elite de la poblacion:
        if(mi_elite->eval < poblacion[0].eval){
            memcpy(mi_elite->vars, poblacion[0].vars, 5*sizeof(double));
            mi_elite->eval = poblacion[0].eval;
        }
    }while(procesar);

    for(int j = 0; j < n_pars; j++){
        delete [] tabla[ idx_pars[j] ];
    }

//    memcpy(mi_elite->vars, poblacion[0].vars, 5*sizeof(double));
//    mi_elite->eval = poblacion[0].eval;
    delete [] poblacion;
    delete [] pob_tmp;
    delete [] probs;
    delete [] seleccion;

    GETTIME_FIN;
    cout << " Mejor candidato encontrado con " << k << " iteraciones. y " << DIFTIME << "segundos. \n Con evaluacion de la funcion objetivo: " << mi_elite->eval << endl;
}
