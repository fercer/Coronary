/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* FILE NAME: IMGPROC.cpp                                                                                    *
*                                                                                                           *
* PURPOSE: Implementation of the image processing routines of the IMGCONT class.                            *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date           Author        Change Id    Release    Description Of Change                                *
* 25/Dic/2016    Fernando C.   0            1.0        Creation                                             *
*                                                                                                           *
************************************************************************************************************/


#include "IMGCONT.h"

/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: computeDistancesMap                                                                        *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* -------                   --------                   -   -------------------------------------------      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Calculates the distances map of the image data (a black and white image is required).                     *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::computeDistancesMap()
{
	my_dist_map = new double[my_height * my_width];

	for (int xy = 0; xy < (my_height * my_width); xy++) {
		*(my_dist_map + xy) = (*(my_img_data + xy) < 1.0) ? 0.0 : MY_INF;
	}

	double *f = new double[my_height > my_width ? my_height : my_width];
	double *dh = new double[my_height];
	int *vh = new int[my_height];
	double *zh = new double[my_height + 1];

	// transform along columns
	for (int x = 0; x < my_width; x++) {
		for (int y = 0; y < my_height; y++) {
			f[y] = *(my_dist_map + y*my_width + x);
		}

		int k = 0;
		vh[0] = 0;
		zh[0] = -MY_INF;
		zh[1] = MY_INF;

		for (int q = 1; q < my_height; q++) {
			double s = ((f[q] + (q*q)) - (f[vh[k]] + (vh[k] * vh[k]))) / (2 * q - 2 * vh[k]);

			while (s <= zh[k]) {
				k--;
				s = ((f[q] + (q*q)) - (f[vh[k]] + (vh[k] * vh[k]))) / (2 * q - 2 * vh[k]);
			}

			k++;
			vh[k] = q;
			zh[k] = s;
			zh[k + 1] = MY_INF;
		}

		k = 0;
		for (int y = 0; y < my_height; y++) {
			while (zh[k + 1] < y) {
				k++;
			}
			*(my_dist_map + y*my_width + x) = ((y - vh[k])*(y - vh[k])) + f[vh[k]];
		}
	}
	delete[] dh;
	delete[] vh;
	delete[] zh;

	double *dw = new double[my_width];
	int *vw = new int[my_width];
	double *zw = new double[my_width + 1];

	// transform along my_height
	for (int y = 0; y < my_height; y++) {
		for (int x = 0; x < my_width; x++) {
			f[x] = *(my_dist_map + y*my_width + x);
		}
		int k = 0;
		vw[0] = 0;
		zw[0] = -MY_INF;
		zw[1] = +MY_INF;

		for (int q = 1; q < my_width; q++) {
			double s = ((f[q] + (q*q)) - (f[vw[k]] + (vw[k] * vw[k]))) / (2 * q - 2 * vw[k]);
			while (s <= zw[k]) {
				k--;
				s = ((f[q] + (q*q)) - (f[vw[k]] + (vw[k] * vw[k]))) / (2 * q - 2 * vw[k]);
			}
			k++;
			vw[k] = q;
			zw[k] = s;
			zw[k + 1] = +MY_INF;
		}

		k = 0;
		for (int x = 0; x < my_width; x++) {
			while (zw[k + 1] < x) {
				k++;
			}
			*(my_dist_map + y*my_width + x) = sqrt(((x - vw[k])*(x - vw[k])) + f[vw[k]]);

			if (((int)*(my_dist_map + y*my_width + x) + 1) > my_max_distance) {
				my_max_distance = (int)*(my_dist_map + y*my_width + x) + 1;
			}
		}
	}

	delete[] dw;
	delete[] vw;
	delete[] zw;
	delete[] f;
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: getDistancesMap                                                                            *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* -------                   --------                   -   -------------------------------------------      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The pointer to the distances map array.                                                                   *
*                                                                                                           *
************************************************************************************************************/
double * IMGCONT::getDistancesMap()
{

	/* If the distances map has not been computed, it is calculated */
	if (!my_dist_map) {
		computeDistancesMap();
	}

	return my_dist_map;
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: computeBoundaries                                                                          *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* -------                   --------                   -   -------------------------------------------      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Computes the boudaries of the image data (a black and white image is required).                           *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::computeBoundaries()
{
	/* If the distances map has not been computed, it is calculated */
	if (!my_dist_map) {
		computeDistancesMap();
	}

	for (int xy = 0; xy < (my_height * my_width); xy++) {
		if ((*(my_img_data + xy) > 0.5) && (*(my_dist_map + xy) < 2.0)) {
			*(my_boundaries + xy) = 1.0;
		}
	}
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: getBoundaries                                                                              *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* -------                   --------                   -   -------------------------------------------      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The pointer to the boundaries array.                                                                      *
*                                                                                                           *
************************************************************************************************************/
double * IMGCONT::getBoundaries()
{
	/* If the boundaries have not been computed, their are calculated */
	if (!my_boundaries) {
		my_boundaries = new double[my_height * my_width];
		memset(my_boundaries, 0, my_height * my_width * sizeof(double));
	}

	return my_boundaries;
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: regionFilling9                                                                             *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* pos_x                     const unsigned int         I   Postion in the X-axis in the image.              *
* pos_y                     const unsigned int         I   Maximum intensity in the image.                  *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The fit, hit or miss response of a square shaped mask of 9x9 pixels.                                      *
*                                                                                                           *
************************************************************************************************************/
bool IMGCONT::regionFilling9(const unsigned int pos_x, const unsigned int pos_y)
{
	int n_hits = 0;

	for (int m = 0; m < 9; m++) {
		// Arriba:
		if ((pos_y - 4) > 0) {
			if (((pos_x - 4 + m) > 0) && ((pos_x - 4 + m) < my_width)) {
				n_hits += (*(my_img_data + (pos_x - 4 + m) + (pos_y - 4)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Abajo:
		if ((pos_y + 4) < my_height) {
			if (((pos_x - 4 + m) > 0) && ((pos_x - 4 + m) < my_width)) {
				n_hits += (*(my_img_data + (pos_x - 4 + m) + (pos_y + 4)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Izquierda:
		if ((pos_x - 4) > 0) {
			if (((pos_y - 4 + m) > 0) && ((pos_y - 4 + m) < my_height)) {
				n_hits += (*(my_img_data + (pos_x - 4) + (pos_y - 4 + m)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Derecha:
		if ((pos_x + 4) < my_width) {
			if (((pos_y - 4 + m) > 0) && ((pos_y - 4 + m) < my_height)) {
				n_hits += (*(my_img_data + (pos_x + 4) + (pos_y - 4 + m)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}
	}
	return (n_hits == 36);
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: regionFilling7                                                                             *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* pos_x                     const unsigned int         I   Postion in the X-axis in the image.              *
* pos_y                     const unsigned int         I   Maximum intensity in the image.                  *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The fit, hit or miss response of a square shaped mask of 7x7 pixels.                                      *
*                                                                                                           *
************************************************************************************************************/
bool IMGCONT::regionFilling7(const unsigned int pos_x, const unsigned int pos_y)
{

	int n_hits = 0;

	for (int m = 0; m < 7; m++) {
		// Arriba:
		if ((pos_y - 3) > 0) {
			if (((pos_x - 3 + m) > 0) && ((pos_x - 3 + m) < my_width)) {
				n_hits += (*(my_img_data + (pos_x - 3 + m) + (pos_y - 3)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Abajo:
		if ((pos_y + 3) < my_height) {
			if (((pos_x - 3 + m) > 0) && ((pos_x - 3 + m) < my_width)) {
				n_hits += (*(my_img_data + (pos_x - 3 + m) + (pos_y + 3)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Izquierda:
		if ((pos_x - 3) > 0) {
			if (((pos_y - 3 + m) > 0) && ((pos_y - 3 + m) < my_height)) {
				n_hits += (*(my_img_data + (pos_x - 3) + (pos_y - 3 + m)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Derecha:
		if ((pos_x + 3) < my_width) {
			if (((pos_y - 3 + m) > 0) && ((pos_y - 3 + m) < my_height)) {
				n_hits += (*(my_img_data + (pos_x + 3) + (pos_y - 3 + m)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}
	}

	return (n_hits == 28);
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: regionFilling5                                                                             *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* pos_x                     const unsigned int         I   Postion in the X-axis in the image.              *
* pos_y                     const unsigned int         I   Maximum intensity in the image.                  *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The fit, hit or miss response of a square shaped mask of 5x5 pixels.                                      *
*                                                                                                           *
************************************************************************************************************/
bool IMGCONT::regionFilling5(const unsigned int pos_x, const unsigned int pos_y)
{
	int n_hits = 0;

	for (int m = 0; m < 5; m++) {
		// Arriba:
		if ((pos_y - 2) > 0) {
			if (((pos_x - 2 + m) > 0) && ((pos_x - 2 + m) < my_width)) {
				n_hits += (*(my_img_data + (pos_x - 2 + m) + (pos_y - 2)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Abajo:
		if ((pos_y + 2) < my_height) {
			if (((pos_x - 2 + m) > 0) && ((pos_x - 2 + m) < my_width)) {
				n_hits += (*(my_img_data + (pos_x - 2 + m) + (pos_y + 2)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Izquierda:
		if ((pos_x - 2) > 0) {
			if (((pos_y - 2 + m) > 0) && ((pos_y - 2 + m) < my_height)) {
				n_hits += (*(my_img_data + (pos_x - 2) + (pos_y - 2 + m)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Derecha:
		if ((pos_x + 2) < my_width) {
			if (((pos_y - 2 + m) > 0) && ((pos_y - 2 + m) < my_height)) {
				n_hits += (*(my_img_data + (pos_x + 2) + (pos_y - 2 + m)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}
	}

	return (n_hits == 20);
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: regionFilling3                                                                             *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* pos_x                     const unsigned int         I   Postion in the X-axis in the image.              *
* pos_y                     const unsigned int         I   Maximum intensity in the image.                  *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The fit, hit or miss response of a square shaped mask of 3x3 pixels.                                      *
*                                                                                                           *
************************************************************************************************************/
bool IMGCONT::regionFilling3(const unsigned int pos_x, const unsigned int pos_y)
{
	int n_hits = 0;

	for (int m = 0; m < 3; m++) {
		// Arriba:
		if ((pos_y - 1) > 0) {
			if (((pos_x - 1 + m) > 0) && ((pos_x - 1 + m) < my_width)) {
				n_hits += (*(my_img_data + (pos_x - 1 + m) + (pos_y - 1)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Abajo:
		if ((pos_y + 1) < my_height) {
			if (((pos_x - 1 + m) > 0) && ((pos_x - 1 + m) < my_width)) {
				n_hits += (*(my_img_data + (pos_x - 1 + m) + (pos_y + 1)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Izquierda:
		if ((pos_x - 1) > 0) {
			if (((pos_y - 1 + m) > 0) && ((pos_y - 1 + m) < my_height)) {
				n_hits += (*(my_img_data + (pos_x - 1) + (pos_y - 1 + m)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Derecha:
		if ((pos_x + 1) < my_width) {
			if (((pos_y - 1 + m) > 0) && ((pos_y - 1 + m) < my_height)) {
				n_hits += (*(my_img_data + (pos_x + 1) + (pos_y - 1 + m)*my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}
	}

	return (n_hits == 12);
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: regionFill                                                                                 *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Fills the blank space inside the image data (a black and white image is required).                        *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::regionFill()
{
	// -------- Mascara 9x9
	for (unsigned int y = 0; y < my_height; y++) {
		for (unsigned int x = 0; x < my_width; x++) {
			if (*(my_img_data + x + y*my_width) < 1.0) {
				*(my_img_data + x + y*my_width) = (regionFilling9(x, y) ? 1.0 : 0.0);
			}
		}
	}
	// -------- Mascara 7x7
	for (unsigned int y = 0; y < my_height; y++) {
		for (unsigned int x = 0; x < my_width; x++) {
			if (*(my_img_data + x + y*my_width) < 1.0) {
				*(my_img_data + x + y*my_width) = (regionFilling7(x, y) ? 1.0 : 0.0);
			}
		}
	}
	// -------- Mascara 5x5
	for (unsigned int y = 0; y < my_height; y++) {
		for (unsigned int x = 0; x < my_width; x++) {
			if (*(my_img_data + x + y*my_width) < 1.0) {
				*(my_img_data + x + y*my_width) = (regionFilling5(x, y) ? 1.0 : 0.0);
			}
		}
	}
	// -------- Mascara 3x3
	for (unsigned int y = 0; y < my_height; y++) {
		for (unsigned int x = 0; x < my_width; x++) {
			if (*(my_img_data + x + y*my_width) < 1.0) {
				*(my_img_data + x + y*my_width) = (regionFilling3(x, y) ? 1.0 : 0.0);
			}
		}
	}
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: computeConnected                                                                           *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* img_ptr                   double *                   I   An image data array.                             *
* pos_x                     const unsigned int         I   The position in the X-axis                       *
* pos_y                     const unsigned int         I   Yhe position in the Y-axis                       *
* my_sets                   int *                      O   The connected sets in the image data             *
* number_of_labeled         unsigned int *             O   An array to store the number of pixels labeled   *
*                                                          in each set.                                     *
* was_visited               bool *                     O   Indicates if a pixel was visited already         *
* number_of_labels          conts int                  I   The current label identifier                     *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The label of the current pixel's set                                                                      *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::computeConnected(double * img_ptr, const int x, const int y, int *my_sets, unsigned int* number_of_labeled, bool* was_visited, const int number_of_labels)
{
	*(was_visited + x + y*my_width) = true;

	/* If the pixel is part of the background, it is ignored */
	if (*(img_ptr + x + y*my_width) < 1.0) {
		return;
	}
	else {
		*(my_sets + x + y*my_width) = number_of_labels;
		*(number_of_labeled + number_of_labels) = *(number_of_labeled + number_of_labels) + 1;
	}

	/* Otherwise, the neigthbors are checked */
	if (y > 0) {
		if (x > 0) {
			if (!*(was_visited + x - 1 + (y - 1)*my_width)) {
				computeConnected(img_ptr, x - 1, y - 1, my_sets, number_of_labeled, was_visited, number_of_labels);
			}
		}

		if (!*(was_visited + x + (y - 1)*my_width)) {
			computeConnected(img_ptr, x, y - 1, my_sets, number_of_labeled, was_visited, number_of_labels);
		}
		if (x < (my_width - 1)) {
			if (!*(was_visited + x + 1 + (y - 1)*my_width)) {
				computeConnected(img_ptr, x + 1, y - 1, my_sets, number_of_labeled, was_visited, number_of_labels);
			}
		}
	}

	if (x > 0) {
		if (!*(was_visited + x - 1 + y*my_width)) {
			computeConnected(img_ptr, x - 1, y, my_sets, number_of_labeled, was_visited, number_of_labels);
		}
	}

	if (x < (my_width - 1)) {
		if (!*(was_visited + x + 1 + y*my_width)) {
			computeConnected(img_ptr, x + 1, y, my_sets, number_of_labeled, was_visited, number_of_labels);
		}
	}

	if (y < (my_height - 1)) {
		if (x > 0) {
			if (!*(was_visited + x - 1 + (y + 1)*my_width)) {
				computeConnected(img_ptr, x - 1, y + 1, my_sets, number_of_labeled, was_visited, number_of_labels);
			}
		}
		
		if (!*(was_visited + x + (y + 1)*my_width)) {
			computeConnected(img_ptr, x, y + 1, my_sets, number_of_labeled, was_visited, number_of_labels);
		}

		if ((x < (my_width - 1))) {
			if (!*(was_visited + x + 1 + (y + 1)*my_width)) {
				computeConnected(img_ptr, x + 1, y + 1, my_sets, number_of_labeled, was_visited, number_of_labels);
			}
		}
	}
}








/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: connectedSets_Dynamic                                                                      *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* img_ptr                   double *                   I   An array to be filtered by its set's length      *
* my_sets                   int *                      O   An array with the different sets found           *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The number of pixels that belong to each set                                                              *
*                                                                                                           *
************************************************************************************************************/
unsigned int* IMGCONT::connectedSets_Dynamic(double * img_ptr, int * my_sets) {
	int number_of_labels = 0;

	bool *was_visited = new bool[my_height * my_width];
	memset(was_visited, 0, my_height * my_width * sizeof(bool));

	unsigned int *temp_labels = new unsigned int[my_height * my_width];
	memset(temp_labels, 0, sizeof(unsigned int)*my_height * my_width);

	memset(my_sets, -1, sizeof(int) * my_height * my_width);

	for (int y = 0; y < my_height; y++) {
		for (int x = 0; x < my_width; x++) {
			/// Se mueve la etiqueta si el pixel es de un nuevo conjunto:
			if (!was_visited[x + y*my_width]) {
				computeConnected(img_ptr, x, y, my_sets, temp_labels, was_visited, number_of_labels);
				if (*(img_ptr + x + y*my_width) > 0) {
					number_of_labels++;
				}
			}
		}
	}

	unsigned int *sets_length = new unsigned int[number_of_labels + 1];
	memcpy(sets_length + 1, temp_labels, number_of_labels * sizeof(unsigned int));
	sets_length[0] = number_of_labels;

	delete[] temp_labels;
	delete[] was_visited;

	return sets_length;
}





/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: increaseSetSize                                                                            *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* my_labels                 int *                      O   An array with the different labels defined       *
* equiv_A                   const int                  I   A label equivalent to the label B                *
* equiv_B                   const int                  I   A label equivalent to the label A                *
* max_number_of_labels      const int                  I   The maximum number of labels to be set as equiv. *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Creates a new equivalency between the label A and B.                                                      *
*                                                                                                           *
************************************************************************************************************/
inline void IMGCONT::increaseSetSize(int * my_labels, const int equiv_A, const int equiv_B, const int max_number_of_labels) {
	int default_label, previous_label;

	if (equiv_A > equiv_B) {
		default_label = equiv_A;
		previous_label = equiv_B;
	}
	else {
		default_label = equiv_B;
		previous_label = equiv_A;
	}

	for (int i = 0; i < max_number_of_labels; i++) {
		if (*(my_labels + i) == -1) {
			break;
		}
		else if (*(my_labels + i) == previous_label) {
			*(my_labels + i) = default_label;
		}
	}
}






/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: connectedSets_Iterative                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* img_ptr                   double *                   I   An array to be filtered by its set's length      *
* my_sets                   int *                      O   An array with the different sets found           *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The number of pixels that belong to each set                                                              *
*                                                                                                           *
************************************************************************************************************/
unsigned int* IMGCONT::connectedSets_Iterative(double * img_ptr, int * my_sets) {

	int max_number_of_labels = my_height * my_width;
	int * labels_values = new int[max_number_of_labels];
	memset(labels_values, -1, max_number_of_labels * sizeof(int));

	int *pix_labels = new int[(my_height + 2)*(my_width + 2)];
	memset(pix_labels, -1, (my_height + 2)*(my_width + 2) * sizeof(int));

	int number_of_sets = -1;

	for (int y = 1; y <= my_height; y++) {
		for (int x = 1; x <= my_width; x++) {
			if (*(img_ptr + (y - 1)*my_width + (x - 1)) > 0.0) {

				const int NorthWest = *(pix_labels + (y - 1)*(my_width + 2) + (x - 1));
				const int North = *(pix_labels + (y - 1)*(my_width + 2) + x);
				const int NorthEast = *(pix_labels + (y - 1)*(my_width + 2) + (x + 1));
				const int East = *(pix_labels + y  *(my_width + 2) + (x + 1));
				const int SouthEast = *(pix_labels + (y + 1)*(my_width + 2) + (x + 1));
				const int South = *(pix_labels + (y + 1)*(my_width + 2) + x);
				const int SouthWest = *(pix_labels + (y + 1)*(my_width + 2) + (x - 1));
				const int West = *(pix_labels + y  *(my_width + 2) + (x - 1));

				int default_label = -1;

				/* North West */
				if (NorthWest >= 0) {
					default_label = NorthWest;
				}

				/// North
				if (North >= 0) {
					if ((default_label >= 0) && (default_label != North)) {
						increaseSetSize(labels_values, default_label, North, max_number_of_labels);
					}
					else {
						default_label = North;
					}
				}

				/// NorthEast
				if (NorthEast >= 0) {
					if ((default_label >= 0) && (default_label != NorthEast)) {
						increaseSetSize(labels_values, default_label, NorthEast, max_number_of_labels);
					}
					else {
						default_label = NorthEast;
					}
				}

				/// East
				if (East >= 0) {
					if ((default_label >= 0) && (default_label != East)) {
						increaseSetSize(labels_values, default_label, East, max_number_of_labels);
					}
					else {
						default_label = East;
					}
				}

				/// SouthEast
				if (SouthEast >= 0) {
					if ((default_label >= 0) && (default_label != SouthEast)) {
						increaseSetSize(labels_values, default_label, SouthEast, max_number_of_labels);
					}
					else {
						default_label = SouthEast;
					}
				}

				/// South
				if (South >= 0) {
					if ((default_label >= 0) && (default_label != South)) {
						increaseSetSize(labels_values, default_label, South, max_number_of_labels);
					}
					else {
						default_label = South;
					}
				}

				/// SouthWest
				if (SouthWest >= 0) {
					if ((default_label >= 0) && (default_label != SouthWest)) {
						increaseSetSize(labels_values, default_label, SouthWest, max_number_of_labels);
					}
					else {
						default_label = SouthWest;
					}
				}

				/// West
				if (West >= 0) {
					if ((default_label >= 0) && (default_label != West)) {
						increaseSetSize(labels_values, default_label, West, max_number_of_labels);
					}
					else {
						default_label = West;
					}
				}

				if (default_label >= 0) {
					*(pix_labels + y*(my_width + 2) + x) = default_label;
				}
				else {
					number_of_sets++;
					*(labels_values + number_of_sets) = number_of_sets;
					*(pix_labels + y*(my_width + 2) + x) = number_of_sets;
				}
			}
		}
	}

	number_of_sets++;

	// Indicar a que conjunto pertenece cada pixel:
	memset(my_sets, -1, my_height * my_width * sizeof(int));
	for (int y = 0; y < my_height; y++) {
		for (int x = 0; x < my_width; x++) {
			const int pixel_label = *(pix_labels + (y + 1)*(my_width + 2) + (x + 1));
			if (pixel_label >= 0) {
				*(my_sets + x + y*my_width) = *(labels_values + pixel_label);
			}
		}
	}

	delete[] labels_values;
	delete[] pix_labels;


	/// Contar cuantos elementos hay en cada grupo:
	unsigned int *sets_length = new unsigned int[number_of_sets + 1];
	memset(sets_length + 1, 0, number_of_sets * sizeof(unsigned int));
	sets_length[0] = number_of_sets;

	for (int xy = 0; xy < my_height * my_width; xy++) {
		if (*(my_sets + xy) >= 0) {
			*(sets_length + 1 + *(my_sets + xy)) = *(sets_length + 1 + *(my_sets + xy)) + 1;
		}
	}

	return sets_length;
}





/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: lengthFilter                                                                               *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* img_ptr                   double *                   I   An array to be filtered by its set's length      *
* threshold_length          const unsigned int         I   The minimum length of a set must have to remain  *
* my_connected_algorithm    CONNECTED_ALG              I   The algorithm to be used in the length filtering *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The array in 'img_ptr' filtered by the length of its sets.                                                *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::lengthFilter(double *img_ptr, const unsigned int threshold_length, CONNECTED_ALG  my_connected_algorithm)
{
    int *my_sets = new int [my_height * my_width];

	unsigned int *my_sets_lengths = NULL;

    switch (my_connected_algorithm) {
    case CONN_DYN:
		my_sets_lengths = connectedSets_Dynamic(img_ptr, my_sets);
        break;
    case CONN_ITER:
        my_sets_lengths = connectedSets_Iterative(img_ptr, my_sets);
        break;
    }

    for( int xy = 0; xy < (my_height * my_width); xy++){
        if( ( *(my_sets + xy) >= 0) && ( *(my_sets_lengths + *(my_sets + xy ) + 1 ) < threshold_length)){
            *(img_ptr + xy) = 0.0;
        }
    }

    delete [] my_sets;
    delete [] my_sets_lengths;
}





/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: lengthFilter                                                                               *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* threshold_length          const unsigned int         I   The minimum length of a set must have to remain  *
* my_connected_algorithm    CONNECTED_ALG              I   The algorithm to be used in the length filtering *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The image contained in the 'my_img_data' array filtered by length.                                        *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::lengthFilter(const unsigned int threshold_length, CONNECTED_ALG  my_connected_algorithm)
{
    lengthFilter(my_img_data, threshold_length, my_connected_algorithm);
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: erosionMask                                                                                *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* erode_ptr                 double *                   I   An auxiliary array for the mask dilatation       *
* pos_x                     const unsigned int         I   Postion in the X-axis in the image.              *
* pos_y                     const unsigned int         I   Maximum intensity in the image.                  *
*                                                                                                           *
* RETURNS:                                                                                                  *
* A fit hit or miss response for the erosion of a 8x8 disk shaped mask.                                     *
*                                                                                                           *
************************************************************************************************************/
inline unsigned char IMGCONT::erosionMask(double * erode_ptr, const int pos_x, const int pos_y)
{
	return (*(erode_ptr + (pos_x - 2) + (pos_y - 4)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 1) + (pos_y - 4)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x)+(pos_y - 4)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 1) + (pos_y - 4)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 2) + (pos_y - 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 3) + (pos_y - 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 2) + (pos_y - 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 1) + (pos_y - 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x)+(pos_y - 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 1) + (pos_y - 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 2) + (pos_y - 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 3) + (pos_y - 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 4) + (pos_y - 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 3) + (pos_y - 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 2) + (pos_y - 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 1) + (pos_y - 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x)+(pos_y - 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 1) + (pos_y - 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 2) + (pos_y - 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 3) + (pos_y - 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 4) + (pos_y - 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 4) + (pos_y - 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 3) + (pos_y - 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 2) + (pos_y - 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 1) + (pos_y - 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x)+(pos_y - 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 1) + (pos_y - 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 2) + (pos_y - 4)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 3) + (pos_y - 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 4) + (pos_y - 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 4) + (pos_y)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 3) + (pos_y)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 2) + (pos_y)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 1) + (pos_y)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 1) + (pos_y)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 2) + (pos_y)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 3) + (pos_y)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 4) + (pos_y)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 4) + (pos_y + 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 3) + (pos_y + 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 2) + (pos_y + 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 1) + (pos_y + 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x)+(pos_y + 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 1) + (pos_y + 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 2) + (pos_y + 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 3) + (pos_y + 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 4) + (pos_y + 1)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 4) + (pos_y + 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 3) + (pos_y + 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 2) + (pos_y + 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 1) + (pos_y + 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x)+(pos_y + 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 1) + (pos_y + 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 2) + (pos_y + 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 3) + (pos_y + 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 4) + (pos_y + 2)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 3) + (pos_y + 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 2) + (pos_y + 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 1) + (pos_y + 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x)+(pos_y + 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 1) + (pos_y + 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 2) + (pos_y + 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 3) + (pos_y + 3)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 2) + (pos_y + 4)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x - 1) + (pos_y + 4)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x)+(pos_y + 4)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 1) + (pos_y + 4)*(my_width + 8)) > 0.0) +
		(*(erode_ptr + (pos_x + 2) + (pos_y + 4)*(my_width + 8)) > 0.0);
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: computeMaskFOV                                                                             *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Estimates the field of view of the image and stores it in the 'img_FOV_mask' array.                       *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::erode(double * img_ptr)
{
	double *erosion_temp = new double[(my_height + 8) *  (my_width + 8)];

	/* Copy the image in 'skl_temp' to 'ptr_tmp' */
	for (unsigned int x = 0; x < (my_width + 8); x++) {
		*(erosion_temp + x) = 1.0;
		*(erosion_temp + (my_width + 8) + x) = 1.0;
		*(erosion_temp + 2 * (my_width + 8) + x) = 1.0;
		*(erosion_temp + 3 * (my_width + 8) + x) = 1.0;
		*(erosion_temp + (my_height + 4)*(my_width + 8) + x) = 1.0;
		*(erosion_temp + (my_height + 5)*(my_width + 8) + x) = 1.0;
		*(erosion_temp + (my_height + 6)*(my_width + 8) + x) = 1.0;
		*(erosion_temp + (my_height + 7)*(my_width + 8) + x) = 1.0;
	}

	for (unsigned int y = 0; y < my_height; y++) {
		*(erosion_temp + (y+4) * (my_width + 8)) = 1.0;
		*(erosion_temp + (y+4) * (my_width + 8) + 1) = 1.0;
		*(erosion_temp + (y+4) * (my_width + 8) + 2) = 1.0;
		*(erosion_temp + (y+4) * (my_width + 8) + 3) = 1.0;

		memcpy(erosion_temp + (y+4) * (my_width + 8) + 4, img_ptr + y * my_width, my_width * sizeof(double));

		*(erosion_temp + (y+4) * (my_width + 8) + my_width + 4) = 1.0;
		*(erosion_temp + (y+4) * (my_width + 8) + my_width + 5) = 1.0;
		*(erosion_temp + (y+4) * (my_width + 8) + my_width + 6) = 1.0;
		*(erosion_temp + (y+4) * (my_width + 8) + my_width + 7) = 1.0;
	}

	for (int y = 0; y < my_height; y++) {
		for (int x = 0; x < my_width; x++) {
			if (*(erosion_temp + (x + 4) + (y + 4)*(my_width + 8)) > 0) {
				const unsigned char resp = erosionMask(erosion_temp, x + 4, y + 4);
				if ((0 < resp) && (resp < 68)) {
					*(img_ptr + x + y*my_width) = 0.0;
				}
			}
		}
	}

	delete [] erosion_temp;
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: computeMaskFOV                                                                             *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Estimates the field of view of the image and stores it in the 'img_FOV_mask' array.                       *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::computeMaskFOV()
{
    /* Threshold the image */
    for( unsigned int y = 0; y < my_height; y++){
		for (unsigned int x = 0; x < my_width; x++) {
			*(my_FOV_mask + x + y * my_width) = (this->getPix(y, x) < 0.1) ? 0.0 : 1.0;
		}
    }

    /* Remove the small connected sets */
    lengthFilter(my_FOV_mask, 1000);

    /*  Erode the FOV mask*/
    erode(my_FOV_mask);
	
	/* Identify the connected sets of the corners of the image */
    for( int xy = 0; xy < my_width * my_height; xy++){
        *(my_FOV_mask + xy) = 1.0 - *(my_FOV_mask + xy);
    }

    int *my_sets = new int [my_width * my_height];
	unsigned int *my_sets_lengths = connectedSets_Iterative(my_FOV_mask, my_sets);
    delete [] my_sets_lengths;


    //// etiqueta de los conjuntos donde existe una esquina
	const int NorthWest = *(my_sets);
	const int NorthEast = *(my_sets + my_width - 1);
	const int SouthWest = *(my_sets + (my_height - 1)*my_width);
	const int SouthEast = *(my_sets + my_height * my_width - 1);

    for( int xy = 0; xy < (my_height * my_width); xy++){
        if( (*(my_sets + xy) >= 0) && ((*(my_sets + xy) == NorthWest) || (*(my_sets + xy) == NorthEast) ||
			(*(my_sets + xy) == SouthWest) || (*(my_sets + xy) == SouthEast)) ){
            *(my_FOV_mask + xy) = 0.0;
        }else{
			*(my_FOV_mask + xy) = 1.0;
        }
    }

    delete [] my_sets;
}






/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: dilMask                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* mask_dil_ptr              double *                   I   An auxiliary array for the mask dilatation       *
* pos_x                     const unsigned int         I   Postion in the X-axis in the image.              *
* pos_y                     const unsigned int         I   Maximum intensity in the image.                  *
*                                                                                                           *
* RETURNS:                                                                                                  *
* A hit, fit or miss response for an 1x1 dilatation mask.                                                   *
*                                                                                                           *
************************************************************************************************************/
inline unsigned char IMGCONT::dilMask(double * mask_dil_ptr, const unsigned int pos_x, const unsigned int pos_y)
{
	return (*(mask_dil_ptr + (pos_x)+(pos_y - 1)*(my_width + 2)) > 0.0) +
		(*(mask_dil_ptr + (pos_x + 1) + (pos_y)*(my_width + 2)) > 0.0) +
		(*(mask_dil_ptr + (pos_x)+(pos_y + 1)*(my_width + 2)) > 0.0) +
		(*(mask_dil_ptr + (pos_x - 1) + (pos_y)*(my_width + 2)) > 0.0);
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: fillMask                                                                                   *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The FOVin the image data with the area out of the FOV filled with the average intensity of the neighbor-  *
* hood values in a mask of 21x21 pixels.                                                                    *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::fillMask()
{
	PIX_PAIR par_tmp;
	par_tmp.my_pix_type = PIX_CROSS;

	double * fill_temp = new double[(my_height + 2) * (my_width + 2)];

	int iter = 0;
	while (iter < my_width) {
		std::vector< PIX_PAIR > boundaries;

		/* Find the boundaries of the mask at the current iteration */
		for (unsigned int y = 0; y < my_height; y++) {
			for (unsigned int x = 0; x < my_width; x++) {
				if (*(fill_temp + (x + 1) + (y + 1)*(my_width + 2)) < 1.0) {
					if (dilMask(fill_temp, x + 1, y + 1) > 0) {
						par_tmp.my_pos_x = x;
						par_tmp.my_pos_y = y;
						boundaries.push_back(par_tmp);
					}
				}
			}
		}

		const int n_pixels_in_boundary = (int)boundaries.size();

		/* If there are no more pixels in the boundaries: Exit */
		if (n_pixels_in_boundary == 0) {
			break;
		}

		/* For each pixel in the boundary, its value in the FOV mask is defined as the average intensity of an 21x21 mask */
		for (int b = 0; b < n_pixels_in_boundary; b++) {
			const int curr_x = (int)boundaries[b].my_pos_x;
			const int curr_y = (int)boundaries[b].my_pos_y;

			const int offset_x_left = (curr_x < 10) ?
				0 : (curr_x - 10);

			const int offset_x_right = (curr_x >= (my_width - 10)) ?
				(my_width - 1) : (curr_x + 10);

			const int offset_y_upper = (curr_y < 10) ?
				0 : (curr_y - 10);

			const int offset_y_lower = (curr_y >= (my_height - 10)) ?
				(my_height - 1) : (curr_y + 10);

			double intensities_sum = 0.0;
			int n_in_neighborhood = 0;

			for (int y = offset_y_upper; y <= offset_y_lower; y++) {
				for (int x = offset_x_left; x <= offset_x_right; x++) {
					if (*(fill_temp + (x + 1) + (y + 1)*(my_width + 2)) > 0.0) {
						intensities_sum += *(my_img_data + x + y*my_width);
						n_in_neighborhood++;
					}
				}
			}

			*(my_img_data + curr_x + curr_y*my_width) = intensities_sum / n_in_neighborhood;
			*(fill_temp + (curr_x + 1) + (curr_y + 1)*(my_width + 2)) = 1.0;
		}
		iter++;
	}

	delete[] fill_temp;
}





/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: extractSkeletonFeatures                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* skl_temp                  double *                   I   An array whith the temporary skeleton            *
* pos_x                     const unsigned int         I   The position in the X-axis.                      *
* pos_y                     const unsigned int         I   The position in the Y-axis.                      *
* deep_level                int *                      O   The current deep level of the graph.             *
* lutable                   const unsigned char *      I   A look up table with the response codes          *
* was_visited               bool *                     I   An array of identifiers                          *
*                                                                                                           *
* RETURNS:                                                                                                  *
* A graph with the pixles of the skeleton characterized by its pixel type.                                  *
*                                                                                                           *
************************************************************************************************************/
IMGCONT::PIX_PAIR * IMGCONT::computeSkeletonGraph(double * skl_temp, const unsigned int pos_x, const unsigned int pos_y, int *deep_level, const unsigned char *lutable, bool *was_visited)
{
	if( *(was_visited + pos_x + pos_y*my_width) ){
		return NULL;
	}

	PIX_PAIR *pix_feratures_temp = new PIX_PAIR;

	const unsigned char resp = sklMask(skl_temp, pos_x, pos_y);

	pix_feratures_temp->my_pos_x = (pos_x-1 - (double)my_width/2)*pixX;
	pix_feratures_temp->my_pos_y = (pos_y-1 - (double)my_height/2)*pixY;

	pix_feratures_temp->my_n_children = 0;

	pix_feratures_temp->my_children[0] = NULL;
	pix_feratures_temp->my_children[1] = NULL;
	pix_feratures_temp->my_children[2] = NULL;

	pix_feratures_temp->my_deep_level = *deep_level;

	/// Calcular el curr_radious de la arteria en el pixel actual:
	const int min_x = ((pos_x - my_max_distance - 1) < 0) ? 0 : (pos_x - my_max_distance - 1);
	const int max_x = ((pos_x + my_max_distance - 1) > my_width) ? my_width : (pos_x + my_max_distance - 1);
	const int min_y = ((pos_y - my_max_distance - 1) < 0) ? 0 : (pos_y - my_max_distance - 1);
	const int max_y = ((pos_y + my_max_distance - 1) > my_height) ? my_height : (pos_y + my_max_distance - 1);

	double curr_distance;
	double curr_radious = MY_INF;
	bool already_visited = false;
	double curr_x_r;
	double curr_y_r;
	for( int yy = min_y; yy < max_y; yy++){
		for( int xx = min_x; xx < max_x; xx++){
			curr_distance = (double)(yy - pos_y + 0.5)*(double)(yy - pos_y + 0.5) + (double)(xx - pos_x + 0.5)*(double)(xx - pos_x + 0.5);
			if( (*(my_boundaries + xx + yy*my_width) > 0.0) && (curr_distance < curr_radious) ){
				curr_radious = curr_distance;
				curr_x_r = (double)xx;
				curr_y_r = (double)yy;
				already_visited = true;
			}
		}
	}
	
	pix_feratures_temp->my_radious = sqrt(curr_radious) * pixX;
	pix_feratures_temp->my_y_r = (curr_y_r - (double)my_height/2)*pixY;
	pix_feratures_temp->my_x_r = (curr_x_r - (double)my_width/2)*pixX;
	pix_feratures_temp->my_angle_alpha = atan2(pix_feratures_temp->my_y_r - pix_feratures_temp->my_pos_y, pix_feratures_temp->my_x_r - pix_feratures_temp->my_pos_x);// + MY_PI / 2.0;

	switch( lutable[ resp ] ){
		case (unsigned char)1:{ /* END point*/
			pix_feratures_temp->my_pix_type = PIX_END;
			break;
		}
		case (unsigned char)2:{ /* BRANCH point */
			pix_feratures_temp->my_pix_type = PIX_BRANCH;
			break;
		}
		case (unsigned char)3:{ /* CROSS point */
			pix_feratures_temp->my_pix_type = PIX_CROSS;
			break;
		}
		default:{
			pix_feratures_temp->my_pix_type = PIX_SKL;
			break;
		}
	}

	*(was_visited + pos_x + pos_y*my_width) = true;

	/* NorthWest */
	if( (resp & (unsigned char)1) && (*(skl_temp + (pos_x-1) + (pos_y-1)*(my_width+2)) < 2.0)){
		if( pix_feratures_temp->my_pix_type == PIX_CROSS || pix_feratures_temp->my_pix_type == PIX_BRANCH ){
			*deep_level = *deep_level + 1;
		}
		pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] = computeSkeletonGraph(skl_temp, pos_x-1, pos_y-1, deep_level, lutable, was_visited);

		if( pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] ){
			pix_feratures_temp->my_n_children++;
		}
	}

	/* North */
	if( (resp & (unsigned char)2) && (*(skl_temp + pos_x + (pos_y-1)*(my_width+2)) < 2.0) ){
		if( pix_feratures_temp->my_pix_type == PIX_CROSS || pix_feratures_temp->my_pix_type == PIX_BRANCH ){
			*deep_level = *deep_level + 1;
		}
		pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] = computeSkeletonGraph(skl_temp, pos_x, pos_y-1, deep_level, lutable, was_visited);
		if( pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] ){
			pix_feratures_temp->my_n_children++;
		}
	}

	/* NorthEast */
	if( (resp & (unsigned char)4) && (*(skl_temp + (pos_x+1) + (pos_y-1)*(my_width+2)) < 2.0) ){
		if( pix_feratures_temp->my_pix_type == PIX_CROSS || pix_feratures_temp->my_pix_type == PIX_BRANCH ){
			*deep_level = *deep_level + 1;
		}
		pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] = computeSkeletonGraph(skl_temp, pos_x+1, pos_y-1, deep_level, lutable, was_visited);

		if( pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] ){
			pix_feratures_temp->my_n_children++;
		}
	}

	/* East */
	if( (resp & (unsigned char)8) && (*(skl_temp + (pos_x+1) + pos_y*(my_width+2)) < 2.0) ){
		if( pix_feratures_temp->my_pix_type == PIX_CROSS || pix_feratures_temp->my_pix_type == PIX_BRANCH ){
			*deep_level = *deep_level + 1;
		}
		pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] = computeSkeletonGraph(skl_temp, pos_x+1, pos_y, deep_level, lutable, was_visited);

		if( pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] ){
			pix_feratures_temp->my_n_children++;
		}
	}

	/* SouthEast */
	if( (resp & (unsigned char)16) && (*(skl_temp + (pos_x+1) + (pos_y+1)*(my_width+2)) < 2.0) ){
		if( pix_feratures_temp->my_pix_type == PIX_CROSS || pix_feratures_temp->my_pix_type == PIX_BRANCH ){
			*deep_level = *deep_level + 1;
		}
		pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] = computeSkeletonGraph(skl_temp, pos_x+1, pos_y+1, deep_level, lutable, was_visited);

		if( pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] ){
			pix_feratures_temp->my_n_children++;
		}
	}

	/* South */
	if( (resp & (unsigned char)32) && (*(skl_temp + pos_x + (pos_y+1)*(my_width+2)) < 2.0) ){
		if( pix_feratures_temp->my_pix_type == PIX_CROSS || pix_feratures_temp->my_pix_type == PIX_BRANCH ){
			*deep_level = *deep_level + 1;
		}
		pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] = computeSkeletonGraph(skl_temp, pos_x, pos_y+1, deep_level, lutable, was_visited);

		if( pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] ){
			pix_feratures_temp->my_n_children++;
		}
	}

	/* SouthWest */
	if( (resp & (unsigned char)64) && (*(skl_temp + (pos_x-1) + (pos_y+1)*(my_width+2)) < 2.0) ){
		if( pix_feratures_temp->my_pix_type == PIX_CROSS || pix_feratures_temp->my_pix_type == PIX_BRANCH ){
			*deep_level = *deep_level + 1;
		}
		pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] = computeSkeletonGraph(skl_temp, pos_x-1, pos_y+1, deep_level, lutable, was_visited);

		if( pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] ){
			pix_feratures_temp->my_n_children++;
		}
	}

	/* West */
	if( (resp & (unsigned char)128) && (*(skl_temp + (pos_x-1) + pos_y*(my_width+2)) < 2.0)){
		if( pix_feratures_temp->my_pix_type == PIX_CROSS || pix_feratures_temp->my_pix_type == PIX_BRANCH ){
			*deep_level = *deep_level + 1;
		}
		pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] = computeSkeletonGraph(skl_temp, pos_x-1, pos_y, deep_level, lutable, was_visited);

		if( pix_feratures_temp->my_children[pix_feratures_temp->my_n_children] ){
			pix_feratures_temp->my_n_children++;
		}
	}

	/* Check the final status of the current pixel: */
	switch( pix_feratures_temp->my_n_children ){
		case 0: pix_feratures_temp->my_pix_type = PIX_END;
				break;
		case 1: pix_feratures_temp->my_pix_type = PIX_SKL;
				break;
		case 2: pix_feratures_temp->my_pix_type = PIX_BRANCH;
				break;
		case 3: pix_feratures_temp->my_pix_type = PIX_CROSS;
				break;
	}

	return pix_feratures_temp;
}





/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: extractSkeletonFeatures                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* ----------                --------                   -   -------------------------------------            *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The features of each pixel in the skeleton.                                                               *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::extractSkeletonFeatures()
{
	if (!my_boundaries) {
		computeBoundaries();
	}

	const unsigned char reference_table[] = {
		0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 0, 0, 0, 0, 2, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 0, 0, 0, 0, 2, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 2, 0, 2, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	double *skl_temp = new double[(my_height + 2)*(my_width + 2)];
	memcpy(skl_temp, my_skeleton, (my_height + 2)*(my_width + 2) * sizeof(double));

	/* Look for a ending point inside the skeleton */
	bool *was_visited = new bool[my_height*my_width];
	memset(was_visited, 0, my_height*my_width * sizeof(bool));

	int start_x, start_y, xy = my_width + 2;
	unsigned char resp;
	do {
		xy++;
		start_x = xy % (my_width + 2);
		start_y = (int)((double)xy / ((double)my_width + 2.0));
		resp = sklMask(skl_temp, start_x, start_y);
	} while (reference_table[resp] != (unsigned char)1);

	*(skl_temp + xy) = 2.0;
	int deep_level = 0;

	my_skeleton_features = computeSkeletonGraph(skl_temp, start_x, start_y, &deep_level, reference_table, was_visited);
	my_skeleton_graph_deep = deep_level;


	delete[] was_visited;
	delete[] skl_temp;
}




/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: deleteSkeletonGraph                                                                        *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* graph_root                PIX_PAIR                   I   The realtive root of a PIX_PAIR graph.           *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::deleteSkeletonGraph(PIX_PAIR *graph_root)
{
	if (graph_root->my_n_children > 0) {
		deleteSkeletonGraph(graph_root->my_children[0]);
		delete graph_root->my_children[0];

		if (graph_root->my_n_children > 1) {
			deleteSkeletonGraph(graph_root->my_children[1]);
			delete graph_root->my_children[1];

			if (graph_root->my_n_children > 2) {
				deleteSkeletonGraph(graph_root->my_children[2]);
				delete graph_root->my_children[2];
			}
		}
	}
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: sklMask                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* skl_temp                  double *                   I   An array whith the temporary skeleton            *
* pos_x                     const unsigned int         I   The position in the X-axis.                      *
* pos_y                     const unsigned int         I   The position in the Y-axis.                      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The fit, hit or miss response of the skeleton mask                                                        *
*                                                                                                           *
************************************************************************************************************/
inline unsigned char IMGCONT::sklMask(double * skl_temp, const unsigned int pos_x, const unsigned int pos_y) {
	return   1 * (*(skl_temp + (pos_x - 1) + (pos_y - 1)*my_width) > 0.0) + /* P2 */
		2 * (*(skl_temp + pos_x + (pos_y - 1)*(my_width + 2)) > 0.0) + /* P3 */
		4 * (*(skl_temp + (pos_x + 1) + (pos_y - 1)*(my_width + 2)) > 0.0) + /* P4 */
		8 * (*(skl_temp + (pos_x + 1) + pos_y  *(my_width + 2)) > 0.0) + /* P5 */
		16 * (*(skl_temp + (pos_x + 1) + (pos_y + 1)*(my_width + 2)) > 0.0) + /* P6 */
		32 * (*(skl_temp + pos_x + (pos_y + 1)*(my_width + 2)) > 0.0) + /* P7*/
		64 * (*(skl_temp + (pos_x - 1) + (pos_y + 1)*(my_width + 2)) > 0.0) + /* P8 */
		128 * (*(skl_temp + (pos_x - 1) + pos_y  *(my_width + 2)) > 0.0);  /* P9 */
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: computeSkeleton                                                                            *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Extracts the morphologic skeleton of the image and stores it in the same image data (a black and white    *
* image is required).                                                                                       *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::computeSkeleton()
{
	if (!my_skeleton) {
		my_skeleton = new double[(my_height + 2)*(my_width + 2)];
	}

	double * skl_temp = new double[(my_height + 2)*(my_width + 2)];
	double * swap_temp;

	const unsigned char reference_table[] = {
		0, 0, 0, 1, 0, 0, 1, 3, 0, 0, 3, 1, 1, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 3, 0, 3, 3,
		0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 2, 2,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0, 2, 0,
		0, 0, 3, 1, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
		3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		2, 3, 1, 3, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		2, 3, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 1, 0, 0, 0, 0, 2, 2, 0, 0, 2, 0, 0, 0 };
	int n_erased;

	do {
		n_erased = 0;
		for (unsigned int y = 1; y <= (my_height + 1); y++) {
			for (unsigned int x = 1; x <= (my_width + 1); x++) {
				if (*(my_skeleton + x + y*(my_width + 2)) > 0.0) {
					const unsigned char resp = *(reference_table + sklMask(my_skeleton, x, y));
					if ((resp == 1) || (resp == 3)) {
						*(skl_temp + x + y * (my_width + 2)) = 0.0;
						n_erased++;
					}
				}
			}
		}

		swap_temp = skl_temp;
		skl_temp = my_skeleton;
		my_skeleton = swap_temp;

		// Segundo paso:
		for (unsigned int y = 1; y <= (my_height + 1); y++) {
			for (unsigned int x = 1; x <= (my_width + 1); x++) {
				if (*(my_skeleton + x + y*(my_width + 2)) > 0.0) {
					unsigned char resp = *(reference_table + sklMask(my_skeleton, x, y));
					if (resp == 2 || resp == 3) {
						*(skl_temp + x + y*(my_width + 2)) = 0.0;
						n_erased++;
					}
				}
			}
		}

		swap_temp = skl_temp;
		skl_temp = my_skeleton;
		my_skeleton = swap_temp;

	} while (n_erased > 0);

	delete[] skl_temp;

}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: getSkeleton                                                                                *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Returns a pointer to the skeleton of the image data.                                                      *
*                                                                                                           *
************************************************************************************************************/
double * IMGCONT::getSkeleton()
{
	if (!my_skeleton) {
		computeSkeleton();
	}

	return my_skeleton;
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: getSkeletonFeatures                                                                        *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Returns a pointer to the skeleton features of each pixel in a graph structure.                            *
*                                                                                                           *
************************************************************************************************************/
IMGCONT::PIX_PAIR * IMGCONT::getSkeletonFeatures()
{
	if (!my_skeleton_features) {
		extractSkeletonFeatures();
	}

	return my_skeleton_features;
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: getSkeletonFeatures                                                                        *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Returns a pointer to the skeleton features of each pixel in a graph structure.                            *
*                                                                                                           *
************************************************************************************************************/
int IMGCONT::getSkeletonFeaturesDeep()
{
	if (!my_skeleton_features) {
		extractSkeletonFeatures();
	}

	return my_skeleton_graph_deep;
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: computeMask                                                                                *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Defines the mask from the image data.                                                                     *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::computeMask()
{

	my_FOV_mask = (double*)malloc(my_height * my_width * sizeof(double));

	computeMaskFOV();
	fillMask();

}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: getMask                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ------------               -   ------------------------------------------       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* A pointer to the FOV mask of the image data.                                                              *
*                                                                                                           *
************************************************************************************************************/
double * IMGCONT::getMask()
{
	if (!my_FOV_mask) {
		computeMask();
	}

	return my_FOV_mask;
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: threshold_by_Otsu                                                                          *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* min_intensity             const double               I   The threshold algorithm                          *
* max_intensity             const double               I   A defined threshold value                        *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The threshold value estimated by the Otsu's thresholding method.                                          *
*                                                                                                           *
************************************************************************************************************/
double IMGCONT::threshold_by_Otsu(const double min, const double max)
{
    const int n_classes = 256;
    double * freq_histogram = new double[n_classes];
	memset(freq_histogram, 0, n_classes*sizeof(double));

    /* Calculate the frequencies histogram */
    double intensities_sum = 0.0;
    const double fraction = 1.0 / (double)(my_height * my_width);

    for( int xy = 0; xy < my_height * my_width; xy ++){
        const int class_i = (int)((n_classes-1) * (*(my_img_data + xy) - min)/(max - min + 1e-12));
        freq_histogram[class_i] += fraction;
        intensities_sum += (double)(class_i +1);
    }

    intensities_sum *= fraction;

    double background_intensities_sum = 0;
    double background_weight = 0.0;
    double between_var, max_between_var = -1.0;
    double threshold_value;

    for( int k = 0; k < n_classes; k++){
        background_weight += (double)freq_histogram[k];
		background_intensities_sum += (double)(freq_histogram[k]*(k+1));

        // Calcular la varianza entre el fore pos_y back ground:
        between_var = (intensities_sum*background_weight - background_intensities_sum);
        between_var *= between_var;
        between_var /= (background_weight*(1.0 - background_weight));

        if( between_var > max_between_var ){
            max_between_var = between_var;
            threshold_value = ((double)k / (double)(n_classes-1));
        }
    }

    delete[] freq_histogram;

    return threshold_value;
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: threshold_by_Ridler_and_Calvard                                                            *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* min_intensity             const double               I   The threshold algorithm                          *
* max_intensity             const double               I   A defined threshold value                        *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The threshold value estimated by the Ridler and Calvard's method.                                         *
*                                                                                                           *
************************************************************************************************************/
double IMGCONT::threshold_by_Ridler_and_Calvard(const double min_intensity, const double max_intensity)
{
	double new_threshold_value = 0.0;
	const double range = 1.0 / (max_intensity - min_intensity);
	const double fraction = 1.0 / (double)(my_height * my_width);

	for (int xy = 0; xy < my_height * my_width; xy++) {
		new_threshold_value += (*(my_img_data + xy) - min_intensity) * range;
	}
	new_threshold_value *= fraction;

	double threshold_value;
	while (1) {
		const double previous_threshold_value = new_threshold_value;
		threshold_value = previous_threshold_value;

		double upper_mean = 0.0;
		double lower_mean = 0.0;
		int upper_count = 0;
		int lower_count = 0;

		for (int xy = 0; xy < my_height * my_width; xy++) {

			if ((*(my_img_data + xy) - min_intensity) * range <= threshold_value) {
				upper_mean += (*(my_img_data + xy) - min_intensity) * range;
				upper_count++;
			}
			else {
				lower_mean += (*(my_img_data + xy) - min_intensity) * range;
				lower_count++;
			}
		}

		upper_mean /= (double)upper_count;
		lower_mean /= (double)lower_count;

		new_threshold_value = (upper_mean + upper_mean) / 2.0;

		if (fabs(new_threshold_value - previous_threshold_value) < 1e-5) {
			break;
		}
	}

	return threshold_value;
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: threshold                                                                                  *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* my_threshold_alg          const THRESHOLD_ALG        I   The threshold algorithm                          *
* threshold_value           const double               I   A defined threshold value                        *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The trhesholded image inside the same image data array.                                                   *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::threshold(const THRESHOLD_ALG my_threshold_alg, const double threshold_value)
{
	double my_max_intensity = -MY_INF;
	double my_min_intensity = MY_INF;

	for (int xy = 0; xy < my_height * my_width; xy++) {
		if (*(my_img_data + xy) < my_min_intensity) {
			my_min_intensity = *(my_img_data + xy);
		}
		if (*(my_img_data + xy) > my_max_intensity) {
			my_max_intensity = *(my_img_data + xy);
		}
	}

	double threshold_estimated_value = 0.0;

	switch (my_threshold_alg) {
	case THRESH_LEVEL:
		threshold_estimated_value = threshold_value;
		break;
	case THRESH_OTSU:
		threshold_estimated_value = threshold_by_Otsu(my_min_intensity, my_max_intensity);
		break;
	case THRESH_RIDLER_CALVARD:
		threshold_estimated_value = threshold_by_Ridler_and_Calvard(my_min_intensity, my_max_intensity);
		break;
	}

	// Se umbraliza la imagen con el valor optimo encontrado:
	for (int xy = 0; xy < my_height * my_width; xy++) {
		*(my_img_data + xy) = (((*(my_img_data + xy) - my_min_intensity) / (my_max_intensity - my_min_intensity)) 
			>= threshold_estimated_value)
			? 1.0 : 0.0;
	}
}
