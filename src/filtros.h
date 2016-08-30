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

#if _WIN32
    #include "D:\Apps\FFTW\include\fftw3.h"
#else
    #include <fftw3.h>
#endif
#include <math.h>
#include <omp.h>

#include <iostream>

#include <QPlainTextEdit>
#include <QProgressBar>

//#include <vtkSmartPointer.h>
//#include <vtkImageData.h>

#include "IMGVTK.h"


// C L A S E: FILTRO  ---------------------------------------------------------------------------------------------- v
class FILTROS{
    public: //----------------------------------------------------------------------------- PUBLIC ------- v
        // T I P O S        D E     D A T O S      P U B L I C O S
        typedef enum LIMITES { INFERIOR, SUPERIOR, DELTA } LIMITES;
        typedef enum SEG_FILTRO { SEG_UNSET, GMF, SS_GABOR } SEG_FILTRO;
        typedef enum EVO_MET { EVO_UNSET, EXHAUSTIVA, EDA_BUMDA, EDA_UMDA, EA_GA } EVO_MET;
        typedef enum EVO_MET_PAR { POPSIZE, MAXGEN, CR, MR} EVO_MET_PAR;
        typedef enum FITNESS { FIT_UNSET, ROC, CORCON } FITNESS;
        typedef enum PARAMETRO { PAR_L, PAR_T, PAR_K, PAR_SIGMA } PARAMETRO;

        /** INDIV:	Define la estructura que contiene los atributos del individuo, y el valor de la funcion para este.  **/
        typedef struct INDIV {
            double eval;
            double vars[4]; // 1: L, 2: T, 3: K, 4: sigma
            unsigned char cadena[128];
        } INDIV;

        typedef double (*FITNESS_PTR)( FILTROS::INDIV *test );

        // M E T O D O S      P U B L I C O S
        FILTROS();
        ~FILTROS();

        void setEvoMet( const EVO_MET evo_met);
        void setEvoMetPar( const EVO_MET_PAR evo_par, const double val);

        void setFiltro( const SEG_FILTRO seg_fil);
        void setFitness( const FITNESS fit_fun);

        void setInput(IMGVTK &img_org);

        void setPar();
        void setPar( const PARAMETRO par, const double val);
        INDIV getPars();
        int getParametrosOptimizar();

        void setLim(const PARAMETRO par, const LIMITES lim, const double val);

        void filtrar();

        void setLog( QTextEdit *txtLog );
        void setLog( FILE *fplog );
        void setLog(const char *ruta_log);

        void setProgressBar( QProgressBar *pBar );

    private: //----------------------------------------------------------------------------- PRIVATE ----- v
        // T I P O S        D E     D A T O S      P R I V A D O S
        /** STAUS:	Define una estructura que almacena las semillas requeridas por el generador de numeros pseudo-aleatorios Hybrid Taus    **/
        typedef struct STAUS{
            unsigned int z1, z2, z3;
        } STAUS;

        /** GEN_PNT:	Define el apuntador a funcion generadora de numeros aleatorios. **/
        typedef double (*GEN_PNT) (const double par1, const double par2);

        // M E T O D O S      P R I V A D O S

        void escribirLog( const char *mensaje );
        void barraProgreso(const int avance, const int max_progress );

        inline double interpolacion(const double *pix, const int j, const int i, const double x, const double y, const int mis_rens, const int mis_cols);
        void rotateImg(const double *org, double *rot, const double ctheta, const double stheta, const int mis_rens, const int mis_cols, const int org_rens, const int org_cols);

        static int compIndiv(const void* A, const void* B);

        //================================================================================== GENERADORES DE NUMEROS ALEATORIOS:

        STAUS* ini_semilla(unsigned int semilla_i);
        unsigned int lcg_s();
        unsigned int lcg_r(unsigned int *mi_semilla);
        unsigned int tausStep(unsigned int *z, const int S1, const int S2, const int S3, const unsigned int M);
        double HybTaus(const double par1, const double par2);

        double anorm_est();

        /*	Metodo:        anorm
            Funcion:: Genera un número aleatorio de la distribución normal con parámetros mu y sigma cuadrada, a partir de la transformación de una variable normal estandar.
        */
        double anorm(const double par1, const double par2){
            return sqrt(par2)*anorm_est() + par1;
        }

        //================================================================================== ALGORITMOS EVOLUTIVOS:
        //// FUNCIONES DE FITNESS:
        double fitnessROC(INDIV *test , double *mi_resp);
        double fitnessCorCon(INDIV *test , double *resp);

        void generarPobInicial(INDIV *poblacion);
        double generarPobInicial(INDIV *poblacion, const double *deltas_var);

        //// ALGORITMOS:
        //--------------------------------------------------------------------------------------------------------------------------------- BUMDA:
        void generarPob(double medias[], double varianzas[], INDIV *poblacion);
        void calcularPars(const INDIV *poblacion, const int truncamiento, double *medias, double *varianzas);
        int seleccionarPob(double *theta_t, const INDIV *poblacion);
        void BUMDA();


        //--------------------------------------------------------------------------------------------------------------------------------- UMDA:
        void generarPob(INDIV *poblacion, const double *probs, const double *deltas_var);
        void calcularPars(const INDIV *poblacion, const int n_bits, const int truncamiento, double *probs);
        void UMDA();


        //--------------------------------------------------------------------------------------------------------------------------------- GA:
        void selecPob(INDIV *sel_grp, const INDIV* poblacion, double *fitness_acum, const double suma_fitness);
        void cruzaPob(INDIV *cruza, const INDIV *sel_grp, const unsigned int n_bits);
        double generarPob(INDIV *poblacion, const INDIV *cruza, const INDIV *sel_grp, const double *deltas_var);
        void GA();

        //--------------------------------------------------------------------------------------------------------------------------------- GA:
        void busquedaExhaustiva();

        // M I E M B R O S      P R I V A D O S
        STAUS *semilla;
        unsigned int semilla_g;

        bool pars_optim[4]; // Indica cuales parametros se van a optimizar: 1: L, 2: T, 3: K, 4: sigma(GMF), 5: delta del umbralizado(Gabor).
        unsigned int idx_pars[4], n_pars;
        SEG_FILTRO filtro_elegido;
        EVO_MET metodo_elegido;
        FITNESS fitness_elegido;

        // Parametros para los filtros:
        INDIV *mi_elite;

        // Entradas comunes:
        double *resp;
        double *org, *dest;
        double *ground_truth, *mask;
        int rows, cols, rows_cols;
        int n_pob, max_iters, seleccion;
        double prob_mutacion;

        double min_vars[4], lim_inf[4], lim_sup[4];

        char *mi_ruta_log;
        FILE *mi_fplog;
        QTextEdit *mi_txtLog;
        QProgressBar *mi_pBar;

        //================================================================================== FILTROS:
        void respGMF(INDIV *test, double *resp);
        bool transformada;
        fftw_complex *Img_fft;
        fftw_complex *Img_fft_HPF;
        void fftImgOrigen();
        void respGabor(INDIV *test, double *resp);


        double calcROC( double *resp );
        double calcCorCon( double *resp );
};
// C L A S E: FILTROS  ----------------------------------------------------------------------------------------- ^

#endif //FILTROS_H_INCLUDED
