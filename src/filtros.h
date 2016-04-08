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

#include <fftw3.h>

#include <math.h>
#include <omp.h>

#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>

#include "IMGVTK.h"

// C L A S E: FILTRO  ---------------------------------------------------------------------------------------------- v
class FILTROS{
    public: //----------------------------------------------------------------------------- PUBLIC ------- v
        // T I P O S        D E     D A T O S      P U B L I C O S
        typedef enum{ GMF, SS_GABOR } SEG_FILTRO;
        typedef enum{ EXHAUSTIVA, EDA_BUMDA, EDA_UMDA, EA_GA } EVO_MET;
        typedef enum{ ROC, ENTROPIA } FITNESS;
        typedef enum{ PAR_L, PAR_T, PAR_K, PAR_SIGMA, PAR_DELTA} PARAMETRO;

        /** INDIV:	Define la estructura que contiene los atributos del individuo, y el valor de la funcion para este.  **/
        typedef struct{
            double eval;
            double vars[5]; // 1: L, 2: T, 3: K, 4: sigma, 5: delta.
            bool cadena[64];
        } INDIV;

        // M E T O D O S      P U B L I C O S
        FILTROS();
        ~FILTROS();

        void setEvoMet( const EVO_MET evo_met, const int m_iters, const int pob);
        void setFiltro( const SEG_FILTRO seg_fil);
        void setFitness( const FITNESS fit_fun);

        void setInputOriginal(IMGVTK &img_org);
        void setInputGround(IMGVTK &img_ground);

        void setOutput( IMGVTK &img_dest);

        void setPar();
        void setPar( const PARAMETRO par, const double val);
        void setLim( const PARAMETRO par, const double inf, const double sup, const double min_var);
        void setLim( const PARAMETRO par, const double inf, const double sup, const unsigned char bits);


        /*  Metodo: getPars
            Funcion: Retorna los parametros utilzados para el filtro.
        */
        INDIV getPars(){
            INDIV out_pars;
            memcpy(out_pars.vars, mi_elite->vars, 5*sizeof(double));
            out_pars.eval = mi_elite->eval;
            return out_pars;
        }

        void filtrar();

    private: //----------------------------------------------------------------------------- PRIVATE ----- v
        // T I P O S        D E     D A T O S      P R I V A D O S
        /** STAUS:	Define una estructura que almacena las semillas requeridas por el generador de numeros pseudo-aleatorios Hybrid Taus    **/
        typedef struct{
            unsigned int z1, z2, z3;
        } STAUS;

        /** GEN_PNT:	Define el apuntador a funcion generadora de numeros aleatorios. **/
        typedef double (*GEN_PNT) (const double par1, const double par2);

        // M E T O D O S      P R I V A D O S
        void deriv( const double *org, double *dX, double *dY, double *dXY, const int mis_rens, const int mis_cols );
        double **coeficientesBiCubico( const double *org, const double *dX, const double *dY, const double *dXY, const int mis_rens, const int mis_cols);
        inline double interpolacion(const double *pix, const int x, const int y, const double delta_x, const double delta_y, const int mis_cols, double **coeficientes);
        inline double interpolacion(const double *pix, const int x, const int y, const double delta_x, const double delta_y, const int mis_cols);
        void rotarImg( const double *org, double *rot, const double theta, const int mis_rens, const int mis_cols, double **coeficientes );

        static int compIndiv(const void* A, const void* B);

        //================================================================================== GENERADORES DE NUMEROS ALEATORIOS:

        /*	Metodo:        ini_semilla
            Funcion:: Inicializa una semilla para el generador HybTaus. Si la semilla dada por el usuario es 0, se genera una semilla aleatoriamente y se utliza, de otro modo se usa la semilla dada por el usuario.
        */
        STAUS* ini_semilla(unsigned int semilla_i){
            if(!semilla_i){
                srand(clock());
                semilla_i = rand();
            }

            STAUS *mi_semilla = new STAUS;
            mi_semilla->z1 = lcg_r(&semilla_i);
            mi_semilla->z2 = lcg_r(&semilla_i);
            mi_semilla->z3 = lcg_r(&semilla_i);

            return mi_semilla;
        }

        unsigned int lcg_s();
        unsigned int lcg_r(unsigned int *mi_semilla);
        unsigned int tausStep(unsigned int *z, const int S1, const int S2, const int S3, const unsigned int M);

        /*	Metodo:        HybTaus
            Funcion:: Genera un numero pseudo-aleatorio por medio del metodo Hybrid Taus Step entre par1 y par2.
        */
        double HybTaus(const double par1, const double par2){
        // Combined period is lcm(p1,p2,p3,p4)~ 2^121
            double num = 2.3283064365387e-10 * (
                // Periods
                tausStep(&semilla->z1, 13, 19, 12, 4294967294UL) ^ // p1=2^31-1
                tausStep(&semilla->z2, 2, 25, 4, 4294967288UL) ^ // p2=2^30-1
                tausStep(&semilla->z3, 3, 11, 17, 4294967280UL) ^ // p3=2^28-1
                lcg_s()	// p4=2^32
            );

            return (par2 - par1) * num + par1;
        }

        double anorm_est();

        /*	Metodo:        anorm
            Funcion:: Genera un número aleatorio de la distribución normal con parámetros mu y sigma cuadrada, a partir de la transformación de una variable normal estandar.
        */
        double anorm(const double par1, const double par2){
            return sqrt(par2)*anorm_est() + par1;
        }

        //================================================================================== ALGORITMOS EVOLUTIVOS:
        //// FUNCIONES DE FITNESS:
        double fitnessROC(INDIV *test);
        void generarPobInicial(INDIV *poblacion);
        void generarPobInicial(INDIV *poblacion, double **tabla);

        //// ALGORITMOS:
        //--------------------------------------------------------------------------------------------------------------------------------- BUMDA:
        void generarPob(INDIV *poblacion, const int n_gen, double medias[5], double varianzas[5]);
        void calcularPars(INDIV *poblacion, const int truncamiento, double *medias, double *varianzas);
        int seleccionarPob(double *tetha_t, INDIV *poblacion);
        void BUMDA();


        //--------------------------------------------------------------------------------------------------------------------------------- UMDA:
        void generarPob(INDIV *poblacion, const int n_gen, double *probs, double **tabla);
        void calcularPars(INDIV *poblacion, const int truncamiento, double *probs);
        void UMDA();


        //--------------------------------------------------------------------------------------------------------------------------------- GA:
        void seleccionarPob(INDIV *poblacion, INDIV *probs, INDIV *pob_tmp, int *seleccion);
        void generarPob(INDIV *poblacion, const double prob_mutacion, double **tabla);
        void GA();



        // M I E M B R O S      P R I V A D O S
        STAUS *semilla;
        unsigned int semilla_g;

        bool pars_optim[5]; // Indica cuales parametros se van a optimizar: 1: L, 2: T, 3: K, 4: sigma(GMF), 5: delta del umbralizado(Gabor).
        int idx_pars[5], n_pars;
        SEG_FILTRO filtro_elegido;
        EVO_MET metodo_elegido;
        FITNESS fitness_elegido;

        // Parametros para los filtros:
        INDIV *mi_elite;

        // Entradas comunes:
        unsigned char *org, *dest;
        unsigned char *ground_truth, *mask;
        int rens, cols, rens_cols, n_pob, max_iters;

        double min_vars[5], lim_inf[5], lim_sup[5];
        unsigned char bits_var[5], max_bits[5];
        int n_bits;

        //================================================================================== FILTROS:
        void respGMF(INDIV *test, double *resp);
        bool transformada;
        double *Img_org;
        fftw_complex *Img_fft;
        void fftImgOrigen();
        void respGabor(INDIV *test, double *resp);
        double calcROC(INDIV *test, double *resp);
};
// C L A S E: FILTROS  ----------------------------------------------------------------------------------------- ^

#endif //FILTROS_H_INCLUDED
