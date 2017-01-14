/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* FILE NAME: performance_functions.cpp                                                                      *
*                                                                                                           *
* PURPOSE: Implementation of performance functions for filters based methods.                               *
*                                                                                                           *
* FILE REFERENCES:                                                                                          *
* Name        I/O        Description                                                                        *
* None        ----       ----------                                                                         *
*                                                                                                           *
* ABNORMAL TERMINATION CONDITIONS, ERROR AND WARNING MESSAGES:                                              *
* None                                                                                                      *
*                                                                                                           *
* NOTES:                                                                                                    *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date        Author        Change Id    Release    Description Of Change                                   *
* 05/Ene/17   Fernando C.   0            1.0        Creation                                                *
*                                                                                                           *
************************************************************************************************************/

#include "performance_functions.h"



int PERFORMANCE_FUNCTIONS::comp_resp( const void *pix_A, const void * pix_B ) {
	if ( ((PERFORMANCE_PAIR*)pix_A)->my_intensity < ((PERFORMANCE_PAIR*)pix_B)->my_intensity) {
		return -1;
	}
	else if (((PERFORMANCE_PAIR*)pix_A)->my_intensity > ((PERFORMANCE_PAIR*)pix_B)->my_intensity) {
		return 1;
	}
	else {
		return 0;
	}
}














/************************************************************************************************************
* PERFORMANCE_FUNCTIONS::PUBLIC                                                                             *
*                                                                                                           *
* FUNCTION NAME: calcAccuracy                                                                               *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Calculates the area under the ROC curve constructed form the response and its ground-truth.               *
*                                                                                                           *
************************************************************************************************************/
double PERFORMANCE_FUNCTIONS::calcROC() {
	
	DEB_MSG("Calculating Az ...");

	unsigned int total_active_response = 0;
	PERFORMANCE_PAIR *all_active_response = (PERFORMANCE_PAIR*)malloc(my_performance_imgs_count *
		(my_img_groundtruth->at(0)).getHeight() * (my_img_groundtruth->at(0)).getWidth() * sizeof(PERFORMANCE_PAIR));

	PERFORMANCE_PAIR *all_groundtruth = (PERFORMANCE_PAIR*)malloc(my_performance_imgs_count *
		(my_img_groundtruth->at(0)).getHeight() * (my_img_groundtruth->at(0)).getWidth() * sizeof(PERFORMANCE_PAIR));

	/* Count the number of positive and negative individuals on the ground-truth that belong to the masked area: */
	unsigned int n_true_positive = 0;
	for (unsigned int i = 0; i < my_performance_imgs_count; i++) {
		for (unsigned int y = 0; y < (my_img_groundtruth->at(i)).getHeight(); y++) {
			for (unsigned int x = 0; x < (my_img_groundtruth->at(i)).getWidth(); x++) {
				if ((my_img_response_mask->at(i)).getPix(y, x) > 0.5) {

					if ((my_img_groundtruth->at(i)).getPix(y, x) > 0.5) {
						n_true_positive++;
					}

					(all_active_response + total_active_response)->my_idx = total_active_response;
					(all_active_response + total_active_response)->my_intensity = (my_img_response->at(i)).getPix(y, x);

					(all_groundtruth + total_active_response)->my_idx = total_active_response;
					(all_groundtruth + total_active_response)->my_intensity = (my_img_groundtruth->at(i)).getPix(y, x);

					total_active_response++;
				}
			}
		}
	}

	unsigned n_true_negative = total_active_response - n_true_positive;

	qsort(all_active_response, total_active_response, sizeof(PERFORMANCE_PAIR), comp_resp);

	DEB_MSG("n negative: " << n_true_negative << ", n positive: " << n_true_positive);
	double *true_positive_fraction_array = (double*)malloc(total_active_response * sizeof(double));
	double *false_positive_fraction_array = (double*)malloc(total_active_response * sizeof(double));

	unsigned int cumulative_false_negative = 0;
	unsigned int cumulative_true_negative = 0;


	for (unsigned int i = 0; i < total_active_response; i++) {
		const double curr_groundtruth = (all_groundtruth + (all_active_response + i)->my_idx)->my_intensity;

		if (curr_groundtruth > 0.5) {
			cumulative_false_negative++;
		}
		else {
			cumulative_true_negative++;
		}

		*(true_positive_fraction_array + i) = (double)(n_true_positive - cumulative_false_negative) /
			(double)n_true_positive;

		*(false_positive_fraction_array + i) = (double)(n_true_negative - cumulative_true_negative) /
			(double)n_true_negative;
	}
	free(all_active_response);
	free(all_groundtruth);

	/* Calculate the Area underd the ROC curve: */
	double Az = 0.0;

#if !defined(NDEBUG)
#if defined(_WIN32) || defined(_WIN64)
	FILE *fp_ROC;
	fopen_s(&fp_ROC, "ROC_curve.fcs", "w");
#else
	FILE *fp_ROC = fopen("ROC_curve.fcs", "w");
#endif
#endif

	for(unsigned int i = 0; i < (total_active_response  -  1); i++) {
		/* Approximate the area under the ROC curve with the Trapezoid Rule: */
		Az += 0.5 * (*(true_positive_fraction_array + i) + *(true_positive_fraction_array + i + 1)) *
			(*(false_positive_fraction_array + i) - *(false_positive_fraction_array + i + 1));

#ifndef NDEBUG
#if defined(_WIN32) || defined(_WIN64)
		fprintf_s(fp_ROC, "%1.12f %1.12f\n", *(true_positive_fraction_array + i), *(false_positive_fraction_array + i));
#else
		fprintf(fp_ROC, "%1.12f %1.12f\n", *(true_positive_fraction_array + i), *(false_positive_fraction_array + i));
#endif
#endif
	}

#ifndef NDEBUG
	fclose(fp_ROC);
#endif

	free(true_positive_fraction_array);
	free(false_positive_fraction_array);

	DEB_MSG("Overall area under the ROC curve: " << Az);

	return Az;
}









/************************************************************************************************************
* PERFORMANCE_FUNCTIONS::PUBLIC                                                                             *
*                                                                                                           *
* FUNCTION NAME: calcCorCon                                                                                 *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Calculates the Correlation and Contrast of the response image.                                            *
*                                                                                                           *
************************************************************************************************************/
double PERFORMANCE_FUNCTIONS::calcCorCon()
{
	/* Scale the response to 8 levels: */
	const int levels = 8;

	double max_resp = -MY_INF;
	double min_resp = MY_INF;

	for (unsigned int i = 0; i < my_performance_imgs_count; i++) {
		for (unsigned int y = 0; y < (my_img_response->at(i)).getHeight(); y++) {
			for (unsigned int x = 0; y < (my_img_response->at(i)).getWidth(); y++) {
				if (max_resp < (my_img_response->at(i)).getPix(y, x)) {
					max_resp = (my_img_response->at(i)).getPix(y, x);
				}
				if (min_resp > (my_img_response->at(i)).getPix(y, x)) {
					min_resp = (my_img_response->at(i)).getPix(y, x);
				}
			}
		}
	}

	const double range = max_resp - min_resp;
	const double scalar = (double)(levels - 1) / range;

	int *scaled_resp = (int*)malloc(my_performance_imgs_count * (my_img_response->at(0)).getHeight() * (my_img_response->at(0)).getWidth() * sizeof(int));

	for (unsigned int i = 0; i < my_performance_imgs_count; i++) {
		for (unsigned int y = 0; y < (my_img_response->at(i)).getHeight(); y++) {
			for (unsigned int x = 0; y < (my_img_response->at(i)).getWidth(); y++) {
				*(scaled_resp + i * (my_img_response->at(i)).getHeight() * (my_img_response->at(i)).getWidth() +
					y * (my_img_response->at(i)).getWidth() * my_performance_imgs_count + x) =
					(int)round(scalar * ((my_img_response->at(i)).getPix(y, x) - min_resp));
			}
		}
	}

	//// Compute the GLCM matrix using an intensity at 8 levels, for 4 orientations:
	double *GLCM_mat = (double*)calloc((levels * levels * 4), sizeof(double));

	const double fraction_1 = 1.0 / (double)((my_img_response->at(0)).getWidth() - 1) *
		((my_img_response->at(0)).getHeight());

	const double fraction_2 = 1.0 / (double)(((my_img_response->at(0)).getWidth() - 1) *
		((my_img_response->at(0)).getHeight() - 1));

	const double fraction_3 = 1.0 / (double)((my_img_response->at(0)).getWidth()   *
		((my_img_response->at(0)).getHeight() - 1));

	const double fraction_4 = 1.0 / (double)(((my_img_response->at(0)).getWidth() - 1) *
		((my_img_response->at(0)).getHeight() - 1));

	double contrast_1 = 0.0, contrast_2 = 0.0, contrast_3 = 0.0, contrast_4 = 0.0;

	double m_i[4] = { 0.0, 0.0, 0.0, 0.0 };
	double m_j[4] = { 0.0, 0.0, 0.0, 0.0 };

	for (unsigned int x = 0; x < my_performance_imgs_count*(my_img_response->at(0)).getWidth() - 1; x++) {
		// SI(x + 1, y    )
		const int i = *(scaled_resp + x);
		const int j = *(scaled_resp + (x + 1));

		*(GLCM_mat + (levels * 4)*i + 4 * j) = *(GLCM_mat + (levels * 4)*i + 4 * j) + fraction_1;
		contrast_1 += (double)((i - j)*(i - j)) * fraction_1;

		*(m_i) += (double)(i + 1) * fraction_1;
		*(m_j) += (double)(j + 1) * fraction_1;
	}

	for (unsigned int y = 1; y < (my_img_response->at(0)).getHeight(); y++) {
		// SI(x   , y - 1)
		const int i = *(scaled_resp + (y + 1)*(my_performance_imgs_count*(my_img_response->at(0)).getWidth()) - 1);
		const int j = *(scaled_resp + y  *(my_performance_imgs_count*(my_img_response->at(0)).getWidth()) - 1);
		*(GLCM_mat + (levels * 4)*i + 4 * j + 2) = *(GLCM_mat + (levels * 4)*i + 4 * j + 2) + fraction_3;
		contrast_3 += (double)((i - j)*(i - j)) * fraction_3;

		*(m_i + 2) += (double)(i + 1) * fraction_3;
		*(m_j + 2) += (double)(j + 1) * fraction_3;
	}

	for (unsigned int y = 1; y < (my_img_response->at(0)).getHeight(); y++) {
		for (unsigned int x = 0; x < my_performance_imgs_count*(my_img_response->at(0)).getWidth() - 1; x++) {
			{// SI(x + 1, y    )
				const int i = *(scaled_resp + (x)+(y)*(my_performance_imgs_count*(my_img_response->at(0)).getWidth()));
				const int j = *(scaled_resp + (x + 1) + (y)*(my_performance_imgs_count*(my_img_response->at(0)).getWidth()));
				*(GLCM_mat + (levels * 4)*i + 4 * j) = *(GLCM_mat + (levels * 4)*i + 4 * j) + fraction_1;
				contrast_1 += (double)((i - j)*(i - j)) * fraction_1;

				*(m_i) += (double)(i + 1) * fraction_1;
				*(m_j) += (double)(j + 1) * fraction_1;

			}
			{// SI(x + 1, y - 1)
				const int i = *(scaled_resp + (x)+(y)*(my_performance_imgs_count*(my_img_response->at(0)).getWidth()));
				const int j = *(scaled_resp + (x + 1) + (y - 1)*(my_performance_imgs_count*(my_img_response->at(0)).getWidth()));
				*(GLCM_mat + (levels * 4)*i + 4 * j + 1) = *(GLCM_mat + (levels * 4)*i + 4 * j + 1) + fraction_2;
				contrast_2 += (double)((i - j)*(i - j)) * fraction_2;

				*(m_i + 1) += (double)(i + 1) * fraction_2;
				*(m_j + 1) += (double)(j + 1) * fraction_2;
			}
			{// SI(x    , y - 1)
				const int i = *(scaled_resp + (x)+(y)*(my_performance_imgs_count*(my_img_response->at(0)).getWidth()));
				const int j = *(scaled_resp + (x)+(y - 1)*(my_performance_imgs_count*(my_img_response->at(0)).getWidth()));
				*(GLCM_mat + (levels * 4)*i + 4 * j + 2) = *(GLCM_mat + (levels * 4)*i + 4 * j + 2) + fraction_3;
				contrast_3 += (double)((i - j)*(i - j)) * fraction_3;

				*(m_i + 2) += (double)(i + 1) * fraction_3;
				*(m_j + 2) += (double)(j + 1) * fraction_3;
			}
			{// SI(x - 1, y - 1)
				const int i = *(scaled_resp + (x + 1) + (y)*(my_performance_imgs_count*(my_img_response->at(0)).getWidth()));
				const int j = *(scaled_resp + (x)+(y - 1)*(my_performance_imgs_count*(my_img_response->at(0)).getWidth()));
				*(GLCM_mat + (levels * 4)*i + 4 * j + 3) = *(GLCM_mat + (levels * 4)*i + 4 * j + 3) + fraction_4;
				contrast_4 += (double)((i - j)*(i - j)) * fraction_4;

				*(m_i + 3) += (double)(i + 1) * fraction_4;
				*(m_j + 3) += (double)(j + 1) * fraction_4;
			}
		}
	}
	free(scaled_resp);

	double s_i[4] = { 0.0, 0.0, 0.0, 0.0 };
	double s_j[4] = { 0.0, 0.0, 0.0, 0.0 };

	//// Calculate variance for i and j form GLCM:
	for (int i = 0; i < levels; i++) {
		const double i_m_1 = (double)(i + 1) - *(m_i);
		const double i_m_2 = (double)(i + 1) - *(m_i + 1);
		const double i_m_3 = (double)(i + 1) - *(m_i + 2);
		const double i_m_4 = (double)(i + 1) - *(m_i + 3);

		for (int j = 0; j < levels; j++) {
			const double j_m_1 = (double)(j + 1) - *(m_j);
			const double j_m_2 = (double)(j + 1) - *(m_j + 1);
			const double j_m_3 = (double)(j + 1) - *(m_j + 2);
			const double j_m_4 = (double)(j + 1) - *(m_j + 3);

			*(s_i) += i_m_1 * i_m_1 * *(GLCM_mat + (levels * 4)*i + 4 * j);
			*(s_j) += j_m_1 * j_m_1 * *(GLCM_mat + (levels * 4)*i + 4 * j);

			*(s_i + 1) += i_m_2 * i_m_2 * *(GLCM_mat + (levels * 4)*i + 4 * j + 1);
			*(s_j + 1) += j_m_2 * j_m_2 * *(GLCM_mat + (levels * 4)*i + 4 * j + 1);

			*(s_i + 2) += i_m_3 * i_m_3 * *(GLCM_mat + (levels * 4)*i + 4 * j + 2);
			*(s_j + 2) += j_m_3 * j_m_3 * *(GLCM_mat + (levels * 4)*i + 4 * j + 2);

			*(s_i + 3) += i_m_4 * i_m_4 * *(GLCM_mat + (levels * 4)*i + 4 * j + 3);
			*(s_j + 3) += j_m_4 * j_m_4 * *(GLCM_mat + (levels * 4)*i + 4 * j + 3);
		}
	}

	//// Calculate standard deviation for i and j form GLCM:
	*(s_i) = sqrt(*(s_i));
	*(s_j) = sqrt(*(s_j));

	*(s_i + 1) = sqrt(*(s_i + 1));
	*(s_j + 1) = sqrt(*(s_j + 1));

	*(s_i + 2) = sqrt(*(s_i + 2));
	*(s_j + 2) = sqrt(*(s_j + 2));

	*(s_i + 3) = sqrt(*(s_i + 3));
	*(s_j + 3) = sqrt(*(s_j + 3));

	//// Calculate correlation form GLCM:
	double correlation_1 = 0.0, correlation_2 = 0.0, correlation_3 = 0.0, correlation_4 = 0.0;
	for (int i = 0; i < levels; i++) {
		const double i_m_1 = (double)(i + 1) - *(m_i);
		const double i_m_2 = (double)(i + 1) - *(m_i + 1);
		const double i_m_3 = (double)(i + 1) - *(m_i + 2);
		const double i_m_4 = (double)(i + 1) - *(m_i + 3);

		const double i_s_1 = *(s_i);
		const double i_s_2 = *(s_i + 1);
		const double i_s_3 = *(s_i + 2);
		const double i_s_4 = *(s_i + 3);

		for (int j = 0; j < levels; j++) {
			const double j_m_1 = (double)(j + 1) - *(m_j);
			const double j_m_2 = (double)(j + 1) - *(m_j + 1);
			const double j_m_3 = (double)(j + 1) - *(m_j + 2);
			const double j_m_4 = (double)(j + 1) - *(m_j + 3);

			correlation_1 += *(GLCM_mat + (levels * 4)*i + 4 * j) * i_m_1*j_m_1 / (i_s_1 * *(s_j));
			correlation_2 += *(GLCM_mat + (levels * 4)*i + 4 * j + 1) * i_m_2*j_m_2 / (i_s_2 * *(s_j + 1));
			correlation_3 += *(GLCM_mat + (levels * 4)*i + 4 * j + 2) * i_m_3*j_m_3 / (i_s_3 * *(s_j + 2));
			correlation_4 += *(GLCM_mat + (levels * 4)*i + 4 * j + 3) * i_m_4*j_m_4 / (i_s_4 * *(s_j + 3));
		}
	}
	free(GLCM_mat);

	double corcon = (correlation_1 + correlation_2 + correlation_3 + correlation_4) / 4.0 + (double)((levels - 1)*(levels - 1)) - (contrast_1 + contrast_2 + contrast_3 + contrast_4) / 4.0;

	return corcon;
}








/************************************************************************************************************
* PERFORMANCE_FUNCTIONS::PUBLIC                                                                             *
*                                                                                                           *
* FUNCTION NAME: calcAccuracy                                                                               *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Calculates the accuracy between the trhesholded image and its respective ground-truth.                    *
*                                                                                                           *
************************************************************************************************************/
double PERFORMANCE_FUNCTIONS::calcAccuracy() 
{
	unsigned int FP = 0, FN = 0, TP = 0, TN = 0;

	for (unsigned int i = 0; i < my_performance_imgs_count; i++) {
		for (unsigned int y = 0; y < (my_img_thresholded_response->at(i)).getHeight(); y++) {
			for (unsigned int x = 0; x < (my_img_thresholded_response->at(i)).getWidth(); x++) {
				if ( (my_img_response_mask->at(i)).getPix(y, x) > 0.5) {
					if ( (my_img_groundtruth->at(i)).getPix(y, x) > 0.5) {
						if ( (my_img_thresholded_response->at(i)).getPix(y, x) > 0.5) {
							TP++;
						}
						else {
							FN++;
						}
					}
					else {
						if ((my_img_thresholded_response->at(i)).getPix(y, x) < 0.5) {
							TN++;
						}
						else {
							FP++;
						}
					}
				}
			}
		}
	}

	return (double)(TP + TN) / (double)(TP + TN + FP + FN);
}












/************************************************************************************************************
* PERFORMANCE_FUNCTIONS::PUBLIC                                                                             *
*                                                                                                           *
* FUNCTION NAME: setInputPerformanceResponse                                                                *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_response          std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void PERFORMANCE_FUNCTIONS::setInputPerformanceResponse(std::vector<IMGCONT>* new_img_response)
{
	my_img_response = new_img_response;

	my_performance_imgs_count = (unsigned int)new_img_response->size();
}












/************************************************************************************************************
* PERFORMANCE_FUNCTIONS::PUBLIC                                                                             *
*                                                                                                           *
* FUNCTION NAME: setInputPerformanceThreshold                                                               *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                   TYPE                      I/O  DESCRIPTION                                     *
* new_img_response_threshold std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.     *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void PERFORMANCE_FUNCTIONS::setInputPerformanceThreshold(std::vector<IMGCONT>* new_img_response_threshold)
{
	my_img_thresholded_response = new_img_response_threshold;

	my_performance_imgs_count = (unsigned int)new_img_response_threshold->size();
}












/************************************************************************************************************
* PERFORMANCE_FUNCTIONS::PUBLIC                                                                             *
*                                                                                                           *
* FUNCTION NAME: setInputPerformanceGroundtruth                                                             *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_groundtruth       std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void PERFORMANCE_FUNCTIONS::setInputPerformanceGroundtruth(std::vector<IMGCONT>* new_img_groundtruth)
{
	my_img_groundtruth = new_img_groundtruth;

	my_performance_imgs_count = (unsigned int)new_img_groundtruth->size();
}












/************************************************************************************************************
* PERFORMANCE_FUNCTIONS::PUBLIC                                                                             *
*                                                                                                           *
* FUNCTION NAME: setInputPerformanceMask                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_mask              std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void PERFORMANCE_FUNCTIONS::setInputPerformanceMask(std::vector<IMGCONT>* new_img_mask)
{
	my_img_response_mask = new_img_mask;

	my_performance_imgs_count = (unsigned int)new_img_mask->size();
}
