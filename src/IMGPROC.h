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


        /** PIX_PAIR:   **/
        typedef struct PIX_PAIR {
            double x, y, x_r, y_r;
            double radious, alpha;
            int deep_level;
            int n_childs;
            PIX_PAIR *childs[3];
			PIX_TYPE pix_type;
        } PIX_PAIR;


        void skeletonization(IMG_IDX img_idx);
        void umbralizar(IMG_IDX img_idx, const TIPO_UMBRAL tipo_umb, const double nivel);

        void lengthFilter(IMG_IDX img_idx, const int min_length, CONNECTED_ALG mi_alg);
        void regionFill(IMG_IDX img_idx);

        IMGCONT & getDistancesMap();

        void detectarBorde(IMG_IDX img_idx);

    private: //----------------------------------------------------------------------------- PRIVATE ----- v
		void calcDistancesMap();

		void maskFOV();
		void fillMask();

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

		IMGCONT auxiliary_img;
		IMGCONT my_FOV_mask;
};

#endif // IMGPROC_H_INCLUDED
