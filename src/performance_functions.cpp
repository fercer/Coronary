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



int comp_resp( const void *pix_A, const void * pix_B ) {
	if (*(double*)pix_A > *(double*)pix_B) {
		return -1;
	}
	else if (*(double*)pix_A < *(double*)pix_B) {
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
	
	const double default_delta = 1e-4;

	/* Count the number of positive and negative individuals on the ground-truth that belong to the masked area: */
	unsigned int n_positive = 0, n_negative = 0;
	for (unsigned int i = 0; i < my_performance_imgs_count; i++) {
		for (unsigned int y = 0; y < (my_img_groundtruth->at(i)).getHeight(); y++) {
			for (unsigned int x = 0; x < (my_img_groundtruth->at(i)).getWidth(); x++) {
				if (((my_img_groundtruth->at(i)).getPix(y, x) > 0.5) && ((my_img_response_mask->at(i)).getPix(y, x) > 0.5)) {
					n_positive++;
				}
				else {
					n_negative++;
				}
			}
		}
	}

	double *positive_array = (double*)malloc(n_positive * sizeof(double));
	double *negative_array = (double*)malloc(n_negative * sizeof(double));
	n_positive = 0;
	n_negative = 0;

	/* Store the values of the response in the corresponding array according to the ground-truth: */
	double max_resp = -MY_INF;
	double min_resp = MY_INF;

	for (unsigned int i = 0; i < my_performance_imgs_count; i++) {
		for (unsigned int y = 0; y < (my_img_response->at(i)).getHeight(); y++) {
			for (unsigned int x = 0; x < (my_img_response->at(i)).getWidth(); x++) {
				/* Store the value only if belongs to the masked area: */
				if ((my_img_response_mask->at(i)).getPix(y, x) > 0.5) {

					const double curr_pixel_intensity = (my_img_response->at(i)).getPix(y, x);

					if ((my_img_groundtruth->at(i)).getPix(y, x) >= 0.5) {
						*(positive_array + n_positive) = curr_pixel_intensity;
						n_positive++;
					}
					else {
						*(negative_array + n_negative) = curr_pixel_intensity;
						n_negative++;
					}

					/* Determine minima and maxima in the response: */
					if (max_resp < curr_pixel_intensity) {
						max_resp = curr_pixel_intensity;
					}
					if (min_resp > curr_pixel_intensity) {
						min_resp = curr_pixel_intensity;
					}
				}
			}
		}
	}

	/* Sort both arrays: */
	qsort(positive_array, n_positive, sizeof(double), comp_resp);
	qsort(negative_array, n_negative, sizeof(double), comp_resp);

	/* Start with the most discriminant threshold */
	double threshold = max_resp;
	int i_positive = (n_positive - 1), i_negative = (n_negative - 1), i_threshold = 0;
	double TP = 0.0, TN = (double)n_negative;

	double TPF_old = 0.0;
	double FPF_old = 0.0;

	double TPF_new, FPF_new;
	double Az = 0.0;

#if !defined(NDEBUG)
#if defined(_WIN32) || defined(_WIN64)
	FILE *fp_ROC;
	fopen_s(&fp_ROC, "ROC_curve.fcs", "w");
#else
	FILE *fp_ROC = fopen("ROC_curve.fcs", "w");
#endif
#endif

	min_resp -= default_delta;
	while (threshold > min_resp) {
		/* Count the positive and negative individuals at this threshold: */
		if (i_negative == 0 && i_positive == 0) {
			break;
		}

		while (i_positive > 0) {
			if (*(positive_array + i_positive) >= threshold) {
				TP += 1.0;
				i_positive--;
			}
			else {
				break;
			}
		}

		while (i_negative > 0) {
			if (*(negative_array + i_negative) >= threshold) {
				TN -= 1.0;
				i_negative--;
			}
			else {
				break;
			}
		}

		/* Compute the True Positive Fraction and the False Positive Fraction: */
		TPF_new = TP / (double)n_positive;
		FPF_new = 1.0 - TN / (double)n_negative;

		/* Approximate the area under the ROC curve with the Trapezoid Rule: */
		Az += (FPF_new - FPF_old)*(TPF_new + TPF_old) / 2.0;

		/* Update the threshold: */
		threshold -= default_delta;
		i_threshold++;

		TPF_old = TPF_new;
		FPF_old = FPF_new;
#ifndef NDEBUG
#if defined(_WIN32) || defined(_WIN64)
		fprintf_s(fp_ROC, "%1.12f %1.12f\n", TPF_old, FPF_old);
#else
		fprintf(fp_ROC, "%1.12f %1.12f\n", TPF_old, FPF_old);
#endif
#endif
	}

#ifndef NDEBUG
	fclose(fp_ROC);
#endif

	free(positive_array);
	free(negative_array);

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

	int *scaled_resp = (int*)malloc(my_performance_imgs_count * (my_img_response->at(i)).getHeight() * (my_img_response->at(i)).getWidth() * sizeof(int));

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

	const double fraction_1 = 1.0 / (double)((my_img_response->at(0)).getWidth() - 1) *   (my_img_response->at(i)).getHeight());
	const double fraction_2 = 1.0 / (double)(((my_img_response->at(0)).getWidth() - 1) * ((my_img_response->at(i)).getHeight() - 1));
	const double fraction_3 = 1.0 / (double)((my_img_response->at(0)).getWidth()   * ((my_img_response->at(i)).getHeight() - 1));
	const double fraction_4 = 1.0 / (double)(((my_img_response->at(0)).getWidth() - 1) * ((my_img_response->at(i)).getHeight() - 1));

	double contrast_1 = 0.0, contrast_2 = 0.0, contrast_3 = 0.0, contrast_4 = 0.0;

	double m_i[4] = { 0.0, 0.0, 0.0, 0.0 };
	double m_j[4] = { 0.0, 0.0, 0.0, 0.0 };

	for (int x = 0; x < my_performance_imgs_count*(my_img_response->at(0)).getWidth() - 1; x++) {
		// SI(x + 1, y    )
		const int i = *(scaled_resp + x);
		const int j = *(scaled_resp + (x + 1));

		*(GLCM_mat + (levels * 4)*i + 4 * j) = *(GLCM_mat + (levels * 4)*i + 4 * j) + fraction_1;
		contrast_1 += (double)((i - j)*(i - j)) * fraction_1;

		*(m_i) += (double)(i + 1) * fraction_1;
		*(m_j) += (double)(j + 1) * fraction_1;
	}

	for (int y = 1; y < (my_img_response->at(0)).getHeight(); y++) {
		// SI(x   , y - 1)
		const int i = *(scaled_resp + (y + 1)*(my_performance_imgs_count*(my_img_response->at(0)).getWidth()) - 1);
		const int j = *(scaled_resp + y  *(my_performance_imgs_count*(my_img_response->at(0)).getWidth()) - 1);
		*(GLCM_mat + (levels * 4)*i + 4 * j + 2) = *(GLCM_mat + (levels * 4)*i + 4 * j + 2) + fraction_3;
		contrast_3 += (double)((i - j)*(i - j)) * fraction_3;

		*(m_i + 2) += (double)(i + 1) * fraction_3;
		*(m_j + 2) += (double)(j + 1) * fraction_3;
	}

	for (int y = 1; y < (my_img_response->at(0)).getHeight(); y++) {
		for (int x = 0; x < my_performance_imgs_count*(my_img_response->at(0)).getWidth() - 1; x++) {
			{// SI(x + 1, y    )
				const int i = *(scaled_resp + (x)+(y)*(my_performance_imgs_count*(my_img_response->at(0)).getWidth()));
				const int j = *(scaled_resp + (x + 1) + (y)*(my_performance_imgs_count*(my_img_response->at(0)).getWidth()));
				*(GLCM_mat + (levels * 4)*i + 4 * j) = *(GLCM_mat + (levels * 4)*i + 4 * j) + fraction_1;
				contrast_1 += (double)((i - j)*(i - j)) * fraction_1;

				*(m_i) += (double)(i + 1) * fraction_1;
				*(m_j) += (double)(j + 1) * fraction_1;

			}
			{// SI(x + 1, y - 1)
				const int i = *(scaled_resp + (x)+(y)*(n_imgs*(my_img_response->at(0)).getWidth()));
				const int j = *(scaled_resp + (x + 1) + (y - 1)*(n_imgs*(my_img_response->at(0)).getWidth()));
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
* FUNCTION NAME: setInputResponse                                                                           *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_response          std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void PERFORMANCE_FUNCTIONS::setInputResponse(std::vector<IMGCONT>* new_img_response)
{
	my_img_response = new_img_response;

	my_performance_imgs_count = new_img_response->size();
}












/************************************************************************************************************
* PERFORMANCE_FUNCTIONS::PUBLIC                                                                             *
*                                                                                                           *
* FUNCTION NAME: setInputThreshold                                                                          *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                   TYPE                      I/O  DESCRIPTION                                     *
* new_img_response_threshold std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.     *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void PERFORMANCE_FUNCTIONS::setInputThreshold(std::vector<IMGCONT>* new_img_response_threshold)
{
	my_img_thresholded_response = new_img_response_threshold;

	my_performance_imgs_count = new_img_response_threshold->size();
}












/************************************************************************************************************
* PERFORMANCE_FUNCTIONS::PUBLIC                                                                             *
*                                                                                                           *
* FUNCTION NAME: setInputGroundtruth                                                                        *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_groundtruth       std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void PERFORMANCE_FUNCTIONS::setInputGroundtruth(std::vector<IMGCONT>* new_img_groundtruth)
{
	my_img_groundtruth = new_img_groundtruth;

	my_performance_imgs_count = new_img_groundtruth->size();
}












/************************************************************************************************************
* PERFORMANCE_FUNCTIONS::PUBLIC                                                                             *
*                                                                                                           *
* FUNCTION NAME: setInputResponseMask                                                                       *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_response_mask     std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void PERFORMANCE_FUNCTIONS::setInputResponseMask(std::vector<IMGCONT>* new_img_response_mask)
{
	my_img_response_mask = new_img_response_mask;

	my_performance_imgs_count = new_img_response_mask->size();
}
