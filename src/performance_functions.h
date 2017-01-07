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

	void setInputResponse(std::vector<IMGCONT> * new_img_response );
	void setInputThreshold(std::vector<IMGCONT> * new_img_response_threshold);
	void setInputGroundtruth(std::vector<IMGCONT> * new_img_groundtruth);
	void setInputResponseMask(std::vector<IMGCONT>* new_img_response_mask);
};
#endif //PERFORMANCE_FUNCTIONS_H_INCLUDED
