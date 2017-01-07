/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* NAME:    opti_pars.h                                                                                      *
*                                                                                                           *
* PURPOSE: Definition of the OPTI_PARS class.                                                               *
*                                                                                                           *
* GLOBAL VARIABLES:                                                                                         *
* Variable    Type    Description                                                                           *
* none        ---     ----------                                                                            *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date           Author        Change Id    Release        Description Of Change                            *
* 06/Jan/2017    Fernando C.   0            1.0            Creation                                         *
*                                                                                                           *
************************************************************************************************************/



#ifndef OPTI_PARS_H_INCLUDED
#define OPTI_PARS_H_INCLUDED



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

#ifdef BUILD_GUI_VERSION
    #include <QPlainTextEdit>
    #include <QProgressBar>
#endif

//#include <vtkSmartPointer.h>
//#include <vtkImageData.h>

#include "rand_fercer.h"
#include "filtros.h"
#include "performance_functions.h"


class OPTI_PARS : public FILTROS, public PERFORMANCE_FUNCTIONS {
    public: //----------------------------------------------------------------------------- PUBLIC ------- v
        // T I P O S        D E     D A T O S      P U B L I C O S
        typedef enum LIMITES { INFERIOR, SUPERIOR, DELTA } LIMITES;
        typedef enum EVO_MET { EVO_UNSET, EXHAUSTIVA, EDA_BUMDA, EDA_UMDA, EA_DE, EA_GA } EVO_MET;
        typedef enum EVO_MET_PAR { POPSIZE, MAXGEN, CR, MR} EVO_MET_PAR;
        typedef enum FITNESS { FIT_UNSET, ROC, CORCON } FITNESS;
        typedef enum PARAMETRO { PAR_L, PAR_T, PAR_K, PAR_SIGMA } PARAMETRO;

        /** INDIV:	Define la estructura que contiene los atributos del individuo, y el valor de la funcion para este.  **/
        typedef struct INDIV {
            double eval;
            double vars[4]; // 1: L, 2: T, 3: K, 4: sigma
            unsigned char cadena[128];
        } INDIV;

        typedef double (*FITNESS_PTR)(INDIV *test);

        // M E T O D O S      P U B L I C O S
		OPTI_PARS();
        ~OPTI_PARS();

        void setEvoMet( const EVO_MET evo_met);
        void setEvoMetPar( const EVO_MET_PAR evo_par, const double val);
        void setFitness( const FITNESS fit_fun);
		void setPar();
        void setPar( const PARAMETRO par, const double val);
        INDIV getPars();
        int getParametrosOptimizar();
        void setLim(const PARAMETRO par, const LIMITES lim, const double val);

		void setInputBase(std::vector<IMGCONT>* new_img_base);
		void setInputMask(std::vector<IMGCONT>* new_img_mask);
		void setInputGroundtruth(std::vector<IMGCONT>* new_img_groundtruth);
		void setInputResponse(std::vector<IMGCONT>* new_img_response);
		void setInputThreshold(std::vector<IMGCONT> * new_img_response_threshold);

    private: //----------------------------------------------------------------------------- PRIVATE ----- v
        static int compIndiv(const void* A, const void* B);

        void generarPobInicial(INDIV *poblacion);
        double generarPobInicial(INDIV *poblacion, const double *deltas_var);

        //// ALGORITMOS:
        //--------------------------------------------------------------------------------------------------------------------------------- BUMDA:
        void generarPob(double * medias, double * varianzas, INDIV *poblacion);
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


        //--------------------------------------------------------------------------------------------------------------------------------- DE:
        void diferenciarPoblacion( INDIV* poblacion, INDIV* pob_base );
        void DE();

        //--------------------------------------------------------------------------------------------------------------------------------- Exhaustive Search:
        void busquedaExhaustiva();

        // M I E M B R O S      P R I V A D O S
		ARGUMENTS *my_args;

        STAUS *semilla;
        unsigned int semilla_g;

        bool pars_optim[4]; /* Indica cuales parametros se van a optimizar: 1: L, 2: T, 3: K, 4: sigma(GMF) */
        unsigned int idx_pars[4], n_pars;
        EVO_MET metodo_elegido;
        FITNESS fitness_elegido;

        // Parametros para los filtros:
        INDIV *mi_elite;

		int n_imgs;
        int rows, my_width, rows_cols;
        int n_pob, max_iters, seleccion;
        double prob_mutacion, prob_cruza;

        double min_vars[4], lim_inf[4], lim_sup[4];

		double fitnessROC(INDIV *test);
		double fitnessCorCon(INDIV *test);
};

#endif //OPTI_PARS_H_INCLUDED
