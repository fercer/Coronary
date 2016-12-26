/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* NAME:    IMGPROC.h                                                                                        *
*                                                                                                           *
* PURPOSE: Definition of the IMGPROC Class inherited from IMGCONT.                                          *
*                                                                                                           *
* GLOBAL VARIABLES:                                                                                         *
* Variable    Type    Description                                                                           *
* none        ---     ----------                                                                            *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date           Author        Change Id    Release        Description Of Change                            *
* 25/Dic/2016    Fernando C.   0            1.0            Creation                                         *
*                                                                                                           *
************************************************************************************************************/


#ifndef IMGPROC_H_INCLUDED
#define IMGPROC_H_INCLUDED

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <iostream>

#include <omp.h>

#include "IMGCONT.h"

class IMGPROC : public IMGCONT{
    public:
        /** PIX_TYPE: **/
        typedef enum PIX_TYPE { PIX_SKL, PIX_END, PIX_BRANCH, PIX_CROSS } PIX_TYPE;

		/** THRESHOLD_TYPE: **/
		typedef enum THRESHOLD_TYPE { TRESH_LEVEL, TRESH_OTSU, TRESH_RIDLER_CALVARD } THRESHOLD_TYPE;

		/** CONNECTED_ALG: **/
		typedef enum CONNECTED_ALG { CONN_DYN, CONN_ITER } CONNECTED_ALG;

        /** PIX_PAR:   **/
        typedef struct PIX_PAR {
            double x, y, x_r, y_r;
            double radio, alpha;
            int nivel;
            int n_hijos;
            PIX_PAR *hijos[3];
			PIX_TYPE pix_tipo;
        } PIX_PAR;

		IMGCONT * definirMask(IMGCONT *img_src);
        void skeletonization(IMG_IDX img_idx);
        void umbralizar(IMG_IDX img_idx, const TIPO_UMBRAL tipo_umb, const double nivel);

        void lengthFilter(IMG_IDX img_idx, const int min_length, CONNECTED_ALG mi_alg);
        void regionFill(IMG_IDX img_idx);

        IMGCONT & getDistancesMap();

        void detectarBorde(IMG_IDX img_idx);

    private: //----------------------------------------------------------------------------- PRIVATE ----- v
		void calcDistancesMap();

		IMGCONT* maskFOV(IMGCONT * img_src);
		void fillMask(IMGCONT *img_src, IMGCONT *mask_src);

        PIX_PAR *grafoSkeleton(double *skl_tmp, const int x, const int y, int *nivel, const unsigned char *lutabla, bool *visitados);
        void extraerCaract(IMG_IDX img_idx);
        void borrarSkeleton( PIX_PAR *raiz );

        // ---------------- Mascaras para region filling
		bool regionFilling9(IMGCONT *img_src, const int x, const int y);
        bool regionFilling7(IMGCONT *img_src, const int x, const int y);
        bool regionFilling5(IMGCONT *img_src, const int x, const int y);
        bool regionFilling3(IMGCONT *img_src, const int x, const int y);
		void regionFill(IMGCONT *img_src);

        double umbralizarOTSU(const double *img_ptr, const double min, const double max);
        double umbralizarRIDCAL( const double *img_ptr, const double min, const double max);

        inline unsigned char sklMask(IMGCONT *skl_ptr, const int x, const int y);

		inline unsigned char dilMask(IMGCONT *mask_dil, const int x, const int y);

        inline unsigned char erosionMask(IMGCONT *ptr_tmp, const int x, const int y);
        void erosionar(IMGCONT *img_src);

		void conexo(IMGCONT *img_src, const int x, const int y, int *conjuntos, unsigned int* n_etiquetados, bool* visitados, const int num_etiquetas);

		unsigned int *conjuntosConexosDinamico(IMGCONT *img_src, int *conjuntos);

        inline void ampliarConjunto(int *etiquetas, const int equiv_A, const int equiv_B, const int max_etiquetas);

        unsigned int *conjuntosConexos(IMGCONT *img_src, int *conjuntos);

		void lengthFilter(IMGCONT *img_src, const unsigned int min_length, ALG_CONJUNTOS mi_alg);

        int max_dist;
		PIX_PAR *pix_caract;

		double *my_dist_map;
};

#endif // IMGPROC_H_INCLUDED
