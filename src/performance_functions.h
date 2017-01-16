/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* FILE NAME: performance_functions.h                                                                        *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date        Author        Change Id    Release    Description Of Change                                   *
* Nov/2016    Fernando C.   0            1.0        Creation                                                *
*                                                                                                           *
************************************************************************************************************/

#ifndef PERFORMANCE_FUNCTIONS_H_INCLUDED
#define PERFORMANCE_FUNCTIONS_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#include "IMGCONT.h"

class PERFORMANCE_FUNCTIONS {
private:
	typedef struct PERFORMANCE_PAIR {
		unsigned int my_x;
		unsigned int my_y;
		unsigned int my_idx;
		double my_intensity;
	} PERFORMANCE_PAIR;

	static int comp_resp(const void *pix_A, const void * pix_B);
	void computeActiveGroundTruth();
	bool active_gt_already_computed;

	unsigned int n_true_positive;
	unsigned int n_true_negative;
	unsigned int total_active_response;

	PERFORMANCE_PAIR *all_groundtruth;
	PERFORMANCE_PAIR *all_active_response;

	double *true_positive_fraction_array;
	double *false_positive_fraction_array;


protected:
	std::vector< IMGCONT > * my_img_response;
	std::vector< IMGCONT > * my_img_response_mask;
	std::vector< IMGCONT > * my_img_thresholded_response;
	std::vector< IMGCONT > * my_img_groundtruth;

	unsigned int my_performance_imgs_count;

	double my_ROC_value;
	double my_ACC_value;

public:
	double calcROC();
	double calcCorCon();
	double calcAccuracy();

	void setInputPerformanceResponse(std::vector<IMGCONT> * new_img_response );
	void setInputPerformanceThreshold(std::vector<IMGCONT> * new_img_response_threshold);
	void setInputPerformanceGroundtruth(std::vector<IMGCONT> * new_img_groundtruth);
	void setInputPerformanceMask(std::vector<IMGCONT>* new_img_mask);

	PERFORMANCE_FUNCTIONS();
	~PERFORMANCE_FUNCTIONS();
};
#endif //PERFORMANCE_FUNCTIONS_H_INCLUDED
