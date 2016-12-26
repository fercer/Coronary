/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* FILE NAME: IMGCONT.cpp                                                                                    *
*                                                                                                           *
* PURPOSE: Implementation of the IMGCONT class for image loading and processing.                            *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date           Author        Change Id    Release    Description Of Change                                *
* 25/Dic/2016    Fernando C.   0            1.0        Creation                                             *
*                                                                                                           *
************************************************************************************************************/


#include "IMGCONT.h"




/************************************************************************************************************
*                                                                                                           *
* VOID CONSTRUCTOR                                                                                          *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ----                       -   ----                                             *
*                                                                                                           *
************************************************************************************************************/

IMGCONT::IMGCONT()
{
	my_height = 0;
	my_width = 0;
		
	my_img_data = NULL;
}




/************************************************************************************************************
*                                                                                                           *
* STRANDARD CONSTRUCTOR                                                                                     *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_height                const unsigned int         I   Height of the new image                          *
* new_width                 const unsigned int         I   Width of the new image                           *
* init_val                  const double               I   The initial value in the whole image             *
*                                                                                                           *
************************************************************************************************************/
IMGCONT::IMGCONT(const unsigned int new_height, const unsigned int new_width, const double init_val = 0.0)
{
	my_height = new_height;
	my_width = new_width;
	
	my_img_data = (double *)malloc((my_height * my_width) * sizeof(double));

	for (unsigned int xy = 0; xy < (my_height * my_width); xy++) {
		*(my_img_data + xy) = init_val;
	}
}




/************************************************************************************************************
*                                                                                                           *
* STRANDARD CONSTRUCTOR                                                                                     *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_height                const unsigned int         I   Height of the new image                          *
* new_width                 const unsigned int         I   Width of the new image                           *
* src_data                  const double *             I   An array to form the image data.                 *
*                                                                                                           *
************************************************************************************************************/
IMGCONT::IMGCONT(const unsigned int new_height, const unsigned int new_width, const double * src_data)
{
	my_height = new_height;
	my_width = new_width;

	my_img_data = (double *)malloc((my_height * my_width) * sizeof(double));
	memcpy(my_img_data, src_data, my_width * my_height * sizeof(double));
}








/************************************************************************************************************
*                                                                                                           *
* COPY CONSTRUCTOR                                                                                          *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* img_src                   const IMGCONT&             I   An image containter to copy to this new image.   *
*                                                                                                           *
************************************************************************************************************/
IMGCONT::IMGCONT(const IMGCONT& img_src) {
	my_height = img_src.my_height;
	my_width = img_src.my_width;

	if (my_img_data) {
		free(my_img_data);
	}

	my_img_data = (double *)malloc((my_height * my_width) * sizeof(double));

	memcpy(my_img_data, img_src.my_img_data, (my_height * my_width) * sizeof(double));
}








/************************************************************************************************************
*                                                                                                           *
* DESTRUCTOR                                                                                                *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* --------                  ----                       -   ----                                             *
*                                                                                                           *
************************************************************************************************************/
IMGCONT::~IMGCONT() {
	if (my_img_data) {
		free(my_img_data);
	}
}







/************************************************************************************************************
*                                                                                                           *
* COPY OPERATOR                                                                                             *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* img_src                   const IMGCONT&             I   An image containter to copy to this new image.   *
*                                                                                                           *
************************************************************************************************************/
IMGCONT & IMGCONT::operator= (const IMGCONT & img_src)
{
	my_height = img_src.my_height;
	my_width = img_src.my_width;

	if (my_img_data) {
		free(my_img_data);
	}

	if (my_img_data) {
		free(my_img_data);
	}
	my_img_data = (double*)malloc(my_height * my_width * sizeof(double));
	memcpy(my_img_data, img_src.my_img_data, my_height * my_width * sizeof(double));

	return *this;
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: writeLog                                                                                   *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* message                   const char *               I   A message to write as log.                       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The message in the starndard output stream.                                                               *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::writeLog( const char *message )
{
    std::cout << message;
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: LoadPNG                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* src_path                  const char *               I   The source path from the image is loaded.        *
* level                     const unsigned int         -   Not used                                         *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The .png image loaded in the image data array.                                                            *
*                                                                                                           *
************************************************************************************************************/
int IMGCONT::LoadPNG(const char *src_path, const unsigned int level)
{
	return -1;
}








/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: LoadPGM                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* src_path                  const char *               I   The source path from the image is loaded.        *
* level                     const unsigned int         -   Not used                                         *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The .pgm image loaded in the image data array.                                                            *
*                                                                                                           *
************************************************************************************************************/
int IMGCONT::LoadPGM(const char *src_path, const unsigned int level)
{
	FILE *img_file = NULL;
#if defined(_WIN32) || defined(_WIN64)
	fopen_s(&img_file, src_path, "r");
#else
	img_file = fopen(ruta_origen, "r");
#endif

	char temp_str[512] = "";

	/* Read the 'Magic number' */
	fgets(temp_str, 512, img_file);

	if ((temp_str[0] != 'P') || (temp_str[1] != '2')) {
		fclose(img_file);
		return -1;
	}

	/* Read commentary, if there is not any commentary, jump to read the image dimensions */
	fpos_t position;
	fgetpos(img_file, &position);

	temp_str[0] = getc(img_file);
	if (temp_str[0] == '#') {
		/* Read the remaining of the commentary */
		fgets(temp_str, 512, img_file);
		fgetpos(img_file, &position);
	}

	double max_intensity;
	unsigned int height;
	unsigned int width;

	fsetpos(img_file, &position);

#if defined(_WIN32) || defined(_WIN64)
	fscanf_s(img_file, "%i", &width);
	fscanf_s(img_file, "%i", &height);
	fscanf_s(img_file, "%lf", &max_intensity);
#else
	fscanf(img_file, "%i", &width);
	fscanf(img_file, "%i", &height);
	fscanf(img_file, "%lf", &max_intensity);
#endif

	IMGCONT new_img(height, width);

	int read_intensity;
	for (int xy = 0; xy < (height* width); xy++) {
#if defined(_WIN32) || defined(_WIN64)
		fscanf_s(img_file, "%i", &read_intensity);
#else
		fscanf(img_file, "%i", &read_intensity);
#endif
		*(new_img.my_img_data + xy) = (double)read_intensity / max_intensity;
	}

	fclose(img_file);

	*this = new_img;

	return 0;
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: LoadDICOM                                                                                  *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* src_path                  const char *               I   The source path from the image is loaded.        *
* level                     const unsigned int         I   Level of the dicom file loaded as image.         *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The DICOOM image loaded in the image data array.                                                          *
*                                                                                                           *
************************************************************************************************************/
int IMGCONT::LoadDICOM(const char *src_path, const unsigned int level)
{
#ifdef BUILD_GDCM_VERSION
	gdcm::ImageReader DICOMreader;
	DICOMreader.SetFileName(src_path);

	DICOMreader.Read();

	gdcm::File &file = DICOMreader.GetFile();
	gdcm::DataSet &ds = file.GetDataSet();

	/* Extract SID */
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x18, 0x1110));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			std::string strm(bv->GetPointer(), bv->GetLength());
			SID = atof(strm.c_str());
		}
	}
	/* Extract SOD */
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x18, 0x1111));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			std::string strm(bv->GetPointer(), bv->GetLength());
			SOD = atof(strm.c_str());
		}
	}

	DDP = SID - SOD;
	const double Magnification = SID / SOD;

	/* Extract pixY\pixX */
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x18, 0x1164));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			std::string strm(bv->GetPointer(), bv->GetLength());
			char *pixXYstr = new char[bv->GetLength()];
			memcpy(pixXYstr, strm.c_str(), bv->GetLength() * sizeof(char));
			char *tmp = strchr(pixXYstr, '\\');
			pixXYstr[tmp - pixXYstr] = '\0';
			pixY = atof(pixXYstr);// / Magnification;
			pixX = atof(tmp + 1);// / Magnification;
			delete[] pixXYstr;
		}
	}
	/* Extract LAO/RAO left(-) to right(+) */
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x18, 0x1510));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			std::string strm(bv->GetPointer(), bv->GetLength());
			LAORAO = atof(strm.c_str());
		}
	}
	/* Extract CRA/CAU cranial(-) to caudal(+) */
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x18, 0x1511));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			std::string strm(bv->GetPointer(), bv->GetLength());
			CRACAU = atof(strm.c_str());
		}
	}

	/* Extract Window Width */
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x18, 0x7030));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			std::string strm(bv->GetPointer(), bv->GetLength());
		}
	}

	/* Extraer Source to Isocenter */
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x21, 0x1017));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			gdcm::Element<gdcm::VR::SL, gdcm::VM::VM1_n> el;
			el.Set(de.GetValue());
			const double SISO = el.GetValue();
			DISO = SID - SISO;
		}
	}
	/* Extract Window Center */
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x28, 0x1050));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			std::string strm(bv->GetPointer(), bv->GetLength());
			WCenter = atof(strm.c_str());
		}
	}
	/* Extract Window Width */
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x28, 0x1051));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			std::string strm(bv->GetPointer(), bv->GetLength());
			WWidth = atof(strm.c_str());
		}
	}
	/* Extraer ECG info */
	int ecg_dim = 0, ecg_np = 0;
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x5000, 0x0005));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			gdcm::Element<gdcm::VR::US, gdcm::VM::VM1_n> el;
			el.Set(de.GetValue());
			ecg_dim = el.GetValue();
		}
	}
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x5000, 0x0010));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			gdcm::Element<gdcm::VR::US, gdcm::VM::VM1_n> el;
			el.Set(de.GetValue());
			ecg_np = el.GetValue();
		}
	}

	int n_levels = 1;
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x8, 0x2143));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			std::string strm(bv->GetPointer(), bv->GetLength());
			n_levels = atof(strm.c_str());
		}
	}
	{
		const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x5000, 0x3000));
		const gdcm::ByteValue *bv = de.GetByteValue();
		if (bv) {
			gdcm::Element<gdcm::VR::OW, gdcm::VM::VM1_n> el;
			el.Set(de.GetValue());
			char filepath_ecg[18];
			strcpy(filepath_ecg, ruta_origen + ruta_l - 7);
			sprintf(filepath_ecg, "%s_%i.dat", filepath_ecg, n_levels);
			FILE *fp = fopen(filepath_ecg, "w");

			for (int i = 0; i < ecg_np; i++) {
				fprintf(fp, "%i %i\n", i, el.GetValue(i));
			}
			fclose(fp);
		}
	}
	/*
	std::string num_points_str = file.GetEntryString(0x5000,0x0010);
	unsigned short num_points;
	convert.clear();
	convert.str(num_points_str);
	convert >> num_points;
	DEB_MSG("Number of Points: " << num_points);

	std::string data_type = file.GetEntryString(0x5000,0x0020);
	DEB_MSG("Type of Data: " << data_type);

	std::string curve_desc = file.GetEntryString(0x5000,0x0022);
	DEB_MSG("Curve Description: " << curve_desc);

	std::string data_rep_str = file.GetEntryString(0x5000,0x0103);
	unsigned short data_rep;
	convert.clear();
	convert.str(data_rep_str);
	convert >> data_rep;

	gdcm::DocEntry *pCurveDataDoc = file.GetDocEntry(0x5000, 0x3000);
	gdcm::DataEntry *pCurveData = dynamic_cast<gdcm::DataEntry *>(pCurveDataDoc);
	uint8_t *curve_data = pCurveData->GetBinArea();

	for(int i = 0; i < num_points; i++){
	DEB_MSG("Pt(" << i <<  ") = " << ((unsigned short*)curve_data)[i]);
	}
	*/

	/* Read images */
	DICOMreader.Read();
	const gdcm::Image &gimage = DICOMreader.GetImage();
	char *buffer = new char[gimage.GetBufferLength()];
	gimage.GetBuffer(buffer);

	const unsigned int* dimension = gimage.GetDimensions();
	unsigned int width = dimension[0];
	unsigned int height = dimension[1];
	const int my_levels = dimension[2];
	const unsigned int my_rows_cols = heigth * width;

	IMGCONT new_img = new IMGCONT(height, width);

	gdcm::PhotometricInterpretation scl_comps = gimage.GetPhotometricInterpretation();
	gdcm::PixelFormat pix_format = gimage.GetPixelFormat();

	switch (scl_comps) {
	case gdcm::PhotometricInterpretation::RGB: {
		if (pix_format == gdcm::PixelFormat::UINT8) {
			for (int y = 0; y < my_height; y++) {
				for (int x = 0; x < my_width; x++) {
					const double pixR = (double)(unsigned char)*(buffer + 3 * x + y*my_width * 3 + level*my_rows_cols * 3) - WCenter + 0.5;
					const double pixG = (double)(unsigned char)*(buffer + 3 * x + 1 + y*my_width * 3 + level*my_rows_cols * 3) - WCenter + 0.5;
					const double pixB = (double)(unsigned char)*(buffer + 3 * x + 2 + y*my_width * 3 + level*my_rows_cols * 3) - WCenter + 0.5;
					double pix = (0.297)*pixR + (0.589)*pixG + (0.114)*pixB;
					if (pix <= -((WWidth - 1) / 2)) {
						pix = 0.0;
					}
					else if (pix > ((WWidth - 1) / 2)) {
						pix = 0.0;
					}
					else {
						pix = pix / (WWidth - 1) + 0.5;
					}
					*(new_img.my_img_data + (my_height - y - 1)*my_height + x) = pix; // 255.0;
				}
			}
		}
		else {
			writeLog(COLOR_RED"<<ERROR AL LEER ARCHIVO DICOM:" COLOR_BACK_YELLOW " Formato de imagen RGB no soportado>>" COLOR_RESET "\n");
			return NULL;
		}
		break;
	}
	case gdcm::PhotometricInterpretation::MONOCHROME1:
	case gdcm::PhotometricInterpretation::MONOCHROME2: {
		if (pix_format == gdcm::PixelFormat::UINT8) {
			for (int y = 0; y < my_height; y++) {
				for (int x = 0; x < my_width; x++) {
					double pix = (double)(unsigned char)*(buffer + level*my_rows_cols + x + y*my_width) - WCenter + 0.5;

					if (pix <= -((WWidth - 1) / 2)) {
						pix = 0.0;
					}
					else if (pix > ((WWidth - 1) / 2)) {
						pix = 1.0;
					}
					else {
						pix = pix / (WWidth - 1) + 0.5;
					}
					*(new_img.my_img_data + (my_height - y - 1)*my_height + x) = pix; // 255.0;
				}
			}

		}
		else if (pix_format == gdcm::PixelFormat::UINT16) {
			unsigned short *buffer16 = (unsigned short*)buffer;
			for (int y = 0; y < my_height; y++) {
				for (int x = 0; x < my_width; x++) {
					const double pixR = (double)((unsigned char)*(buffer16 + 3 * x + y*my_width * 3 + level*my_rows_cols * 3) / 16) - WCenter + 0.5;
					const double pixG = (double)((unsigned char)*(buffer16 + 3 * x + 1 + y*my_width * 3 + level*my_rows_cols * 3) / 16) - WCenter + 0.5;
					const double pixB = (double)((unsigned char)*(buffer16 + 3 * x + 2 + y*my_width * 3 + level*my_rows_cols * 3) / 16) - WCenter + 0.5;
					double pix = (0.297)*pixR + (0.589)*pixG + (0.114)*pixB;
					if (pix <= -((WWidth - 1) / 2)) {
						pix = 0.0;
					}
					else if (pix > ((WWidth - 1) / 2)) {
						pix = 0.0;
					}
					else {
						pix = pix / (WWidth - 1) + 0.5;
					}
					*(new_img.my_img_data + (my_height - y - 1)*my_height + x) = pix; // 255.0;
				}
			}

		}
		else {
			writeLog(COLOR_RED"<<ERROR AL LEER ARCHIVO DICOM:" COLOR_BACK_YELLOW " Formato de imagen RGB no soportado>>" COLOR_RESET "\n");
		}
		break;
	}
	}
	delete[] buffer;

	*this = new_img;

	return 0;
#else
	return -1;
#endif
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: Load                                                                                       *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* src_path                  const char *               I   The source path from the image is loaded.        *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The image loaded in the data array.                                                                       *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::Load(const char *src_path, const unsigned int level)
{
	if (LoadPNG(src_path) == 0) {
		return;
	}
	else if (LoadPGM(src_path) == 0) {
		return;
	}
	else if (LoadDICOM(src_path, level) == 0) {
		return;
	}
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: SavePGM                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* out_path                  const char *               I   The source path from the image is loaded.        *
* my_min                    const double               I   Minimum intensity in the image.                  *
* my_max                    const double               I   Maximum intensity in the image.                  *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::SavePGM(const char * out_path, const double my_min, const double my_max)
{
	FILE *img_file = NULL;

#if defined(_WIN32) || defined(_WIN64)
	fopen_s(&img_file, out_path, "w");
#else
	img_file = fopen(out_path, "w");
#endif
	fprintf(img_file, "P2\n");
	fprintf(img_file, "# by FerCer\n");

	fprintf(img_file, "%i %i\n", my_width, my_height);
	fprintf(img_file, "255\n");

	int intensity;
	for (unsigned int y = 0; y < my_height; y++) {
		for (unsigned int x = 0; x < my_width; x++) {
			intensity = (unsigned char)(255.0 * (*(my_img_data + x + y*my_width) - my_min) / (my_max - my_min));
			fprintf(img_file, "%i\n", intensity);
		}
	}

	fclose(img_file);
}







/************************************************************************************************************
* IMGCONT::PRIVATE                                                                                          *
*                                                                                                           *
* FUNCTION NAME: SavePNG                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* out_path                  const char *               I   The source path from the image is loaded.        *
* my_min                    const double               I   Minimum intensity in the image.                  *
* my_max                    const double               I   Maximum intensity in the image.                  *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::SavePNG(const char * out_path, const double my_min, const double my_max)
{
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: Save                                                                                       *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* out_path                  const char *               I   The path where the image will be saved.          *
* output_type               IMG_TYPE                   I   The file type the image will be saved.           *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::Save(const char * out_path, const IMG_TYPE output_type)
{
	double my_min = MY_INF;
	double my_max = -MY_INF;

	for (int y = 0; y < my_height; y++) {
		for (int x = 0; x < my_width; x++) {
			if (my_min > *(my_img_data + x + y*my_width)) {
				my_min = *(my_img_data + x + y*my__width);
			}
			if (my_max < *(my_img_data + x + y*my_width)) {
				my_max = *(my_img_data + x + y*my_width);
			}
		}
	}

	if (fabs(my_max - my_min) < 1e-12) {
		my_max = 1.0;
		my_min = 0.0;
	}

	switch (output_type) {
	case IMGPGM:
		DEB_MSG("As pgm");
		SavePGM(out_path, my_min, my_max);
		break;
	case IMGPNG:
		DEB_MSG("As png");
		break;
	}
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: getPix                                                                                     *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* row_y                     const unsigned int         I   Row position.                                    *
* col_x                     const unsigned int         I   Column position.                                 *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The pixel value in the defined position.                                                                  *
*                                                                                                           *
************************************************************************************************************/
double IMGCONT::getPix(const unsigned int row_y, const unsigned int col_x)
{
	return *(my_img_data + row_y + col_x);
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: setPix                                                                                     *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* row_y                     const unsigned int         I   Row position.                                    *
* col_x                     const unsigned int         I   Column position.                                 *
* new_val                   const double               I   New value to assign to the current position.     *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void IMGCONT::setPix(const unsigned int row_y, const unsigned int col_x, const double new_val)
{
	*(my_img_data + row_y + col_x) = new_val;
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: getHeight                                                                                  *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* -----                     -----                      -   -------------                                    *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The image height.                                                                                         *
*                                                                                                           *
************************************************************************************************************/
int IMGCONT::getHeight()
{
	return my_height;
}







/************************************************************************************************************
* IMGCONT::PUBLIC                                                                                           *
*                                                                                                           *
* FUNCTION NAME: getWidth                                                                                   *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* -----                     -----                      -   -------------                                    *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The image width.                                                                                          *
*                                                                                                           *
************************************************************************************************************/
int IMGCONT::getWidth()
{
	return my_width;
}