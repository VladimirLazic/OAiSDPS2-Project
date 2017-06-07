#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include <math.h>


void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	uchar *Y_buff = new uchar[xSize * ySize];
	uchar *Y_output = new uchar[newXSize * newYSize];

	char *V_buff = new char[xSize * ySize / 4];
	char *U_buff = new char[xSize * ySize / 4];
	char *V_output = new char[newXSize * newYSize / 4];
	char *U_output = new char[newXSize * newYSize / 4];

	//Convertig from RGB to YUV420
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	//Scaling factors
	const double F_horizontal = (double) newXSize / xSize;
	const double F_vertical = (double) newYSize / ySize;

	//Sample and hold algorithm
	for (int i = 0; i < newYSize; i++) {
		for (int j = 0; j < newXSize; j++) {

			int ii;
			int jj;

			ii = (i - 1) / F_vertical;
			jj = (j - 1) / F_horizontal;

			if (ii < ySize - 1)
				ii = (i - 1) / F_vertical + 1;

			if (jj < xSize - 1)
				jj = (j - 1) / F_horizontal + 1;

			Y_output[i * newXSize + j] = Y_buff[ii * xSize + jj];
		}
	}

	for (int i = 0; i < newYSize / 2; i++) {
		for (int j = 0; j < newXSize / 2; j++) {

			int ii;
			int jj;

			ii = (i - 1) / F_vertical;
			jj = (j - 1) / F_horizontal;

			if (ii < ySize / 2 - 1)
				ii = (i - 1) / F_vertical + 1;

			if (jj < xSize / 2 - 1)
				jj = (j - 1) / F_horizontal + 1;


			U_output[i * newXSize / 2 + j] = U_buff[ii * xSize / 2 + jj];
			V_output[i * newXSize / 2 + j] = V_buff[ii * xSize / 2 + jj];
		}
	}

	YUV420toRGB(Y_output, U_output, V_output, newXSize, newYSize, output);

	delete[] Y_buff;
	delete[] Y_output;
	delete[] U_buff;
	delete[] U_output;
	delete[] V_buff;
	delete[] V_output;
}

void bilinearInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	uchar *Y_buff = new uchar[xSize * ySize];
	uchar *Y_output = new uchar[newXSize * newYSize];

	char *V_buff = new char[xSize * ySize / 4];
	char *U_buff = new char[xSize * ySize / 4];
	char *V_output = new char[newXSize * newYSize / 4];
	char *U_output = new char[newXSize * newYSize / 4];

	//Convertig from RGB to YUV420
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	//Scaling factors
	const double F_horizontal = (double) newXSize / xSize;
	const double F_vertical = (double) newYSize / ySize;

	for (int i = 0; i < newXSize; i++) {				//horizontal
		for (int j = 0; j < newYSize; j++) {			//vertical
			double a = i / F_horizontal - floor(i / F_horizontal);
			double b = j / F_vertical - floor(j / F_vertical);

			int ii = i / F_horizontal;
			int jj = j / F_vertical;

			int iii = ii, jjj = jj;
			if (ii < xSize - 1)
				iii = ii + 1;
			if (jj < ySize - 1)
				jjj = jj + 1;

			Y_output[j * newXSize + i] = (1 - a) * (1 - b) * Y_buff[jj * xSize + ii] + a * (1 - b) *Y_buff[jj * xSize + (iii)]
				+ (1 - a) * b * Y_buff[(jjj) * xSize + ii] + a * b * Y_buff[(jjj) * xSize + (iii)];
		}
	}

	for (int i = 0; i < newXSize / 2; i++) {				//horizontal
		for (int j = 0; j < newYSize / 2; j++) {			//vertical
			double a = i / F_horizontal - floor(i / F_horizontal);
			double b = j / F_vertical - floor(j / F_vertical);

			int ii = i / F_horizontal;
			int jj = j / F_vertical;

			int iii = ii, jjj = jj;

			if (ii < xSize / 2 - 1)
				iii = ii + 1;
			if (jj < ySize / 2 - 1)
				jjj = jj + 1;

			U_output[j * newXSize / 2 + i] = (1 - a)*(1 - b)*U_buff[jj * xSize / 2 + ii] + (a)*(1 - b)*U_buff[jj * xSize / 2 + (iii)]
				+ (1 - a)*(b)*U_buff[(jjj) * xSize / 2 + ii] + (a)*(b)*U_buff[(jjj) * xSize / 2 + (iii)];

			V_output[j * newXSize / 2 + i] = (1 - a)*(1 - b)*V_buff[jjj * xSize / 2 + iii] + (a)*(1 - b)*V_buff[jj * xSize / 2 + (iii)]
				+ (1 - a)*(b)*V_buff[(jjj) * xSize / 2 + ii] + (a)*(b)*V_buff[(jjj) * xSize / 2 + (iii)];
		}
	}


	// Back to RGB
	YUV420toRGB(Y_output, U_output, V_output, newXSize, newYSize, output);

	delete[] Y_buff;
	delete[] Y_output;
	delete[] U_buff;
	delete[] U_output;
	delete[] V_buff;
	delete[] V_output;
}

void bicubicInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */
}

void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
}

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
}