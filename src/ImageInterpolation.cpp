#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include "ImageFilter.h"
#include <math.h>

#define PI 3.141592653589793238462643383279502884L


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


/* A help function for bicubic interpolation. Works for unsigned char */
uchar cubicInterpolate(uchar p[4], double x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

/* A help function for bicubic interpolation. Works for signed chars*/
char cubicInterpolate(char p[4], double x) {
	int retVal = p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));

	if (retVal > 127) {
		retVal = 127;
	} else if (retVal < -127) {
		retVal = -127;
	}

	return (char) retVal;
}

void bicubicInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	uchar *Y_buff = new uchar[xSize * ySize];
	uchar *Y_output = new uchar[newXSize * newYSize];

	char *V_buff = new char[xSize * ySize / 4];
	char *U_buff = new char[xSize * ySize / 4];
	char *V_output = new char[newXSize * newYSize / 4];
	char *U_output = new char[newXSize * newYSize / 4];

	//Buffers for single cubic interpolation
	uchar *Y_horizontal = new uchar[4];
	char *U_horizontal = new char[4];
	char *V_horizontal = new char[4];

	uchar *Y_vertical = new uchar[4];
	char *U_vertical = new char[4];
	char *V_vertical = new char[4];

	//Convertig from RGB to YUV420
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	//Scaling factors
	const double F_horizontal = (double)newXSize / xSize;
	const double F_vertical = (double)newYSize / ySize;

	//Extended buffers for correct pixel calculation
	uchar *Y_extnd = new uchar[(xSize + 4) * (ySize + 4)];
	char *U_extnd = new char[(xSize + 16) * (ySize + 16) / 4];
	char *V_extnd = new char[(xSize + 16) * (ySize + 16) / 4];

	extendBorders(Y_buff, xSize, ySize, Y_extnd, 2);
	extendBorders(U_buff, xSize / 2, ySize / 2, U_extnd, 2);
	extendBorders(V_buff, xSize / 2, ySize / 2, V_extnd, 2);

	for (int i = 0; i < newYSize; i++) {
		for (int j = 0; j < newXSize; j++) {
			double b = i / F_vertical - floor(i / F_vertical);
			double a = j / F_horizontal - floor(j / F_horizontal);

			int new_x = j / F_horizontal;
			int new_y = i / F_vertical;

			int v = 0;
			for (int h = new_y - 1; h < new_y + 3; h++, v++) {
				Y_horizontal[0] = Y_extnd[h * (xSize + 4) + new_x - 1];
				Y_horizontal[1] = Y_extnd[h * (xSize + 4) + new_x];
				Y_horizontal[2] = Y_extnd[h * (xSize + 4) + new_x + 1];
				Y_horizontal[3] = Y_extnd[h * (xSize + 4) + new_x + 2];

				Y_vertical[v] = cubicInterpolate(Y_horizontal, a);
			}

			Y_output[i * newXSize + j] = cubicInterpolate(Y_vertical, b);
		}
	}

	for (int i = 0; i < newYSize / 2; i++) {
		for (int j = 0; j < newXSize / 2; j++) {
			double b = i / F_vertical - floor(i / F_vertical);
			double a = j / F_horizontal - floor(j / F_horizontal);

			int new_x = j / F_horizontal;
			int new_y = i / F_vertical;

			int v = 0;
			for (int h = new_y - 1; h < new_y + 3; h++, v++) {
				U_horizontal[0] = U_extnd[h * (xSize / 2 + 4) + new_x - 1];
				U_horizontal[1] = U_extnd[h * (xSize / 2 + 4) + new_x];
				U_horizontal[2] = U_extnd[h * (xSize / 2 + 4) + new_x + 1];
				U_horizontal[3] = U_extnd[h * (xSize / 2 + 4) + new_x + 2];

				V_horizontal[0] = V_extnd[h * (xSize / 2 + 4) + new_x - 1];
				V_horizontal[1] = V_extnd[h * (xSize / 2 + 4) + new_x];
				V_horizontal[2] = V_extnd[h * (xSize / 2 + 4) + new_x + 1];
				V_horizontal[3] = V_extnd[h * (xSize / 2 + 4) + new_x + 2];

				U_vertical[v] = cubicInterpolate(U_horizontal, a);
				V_vertical[v] = cubicInterpolate(V_horizontal, a);
			}

			U_output[i * newXSize / 2 + j] = cubicInterpolate(U_vertical, b);
			V_output[i * newXSize / 2 + j] = cubicInterpolate(V_vertical, b);
		}
	}

	YUV420toRGB(Y_output, U_output, V_output, newXSize, newYSize, output);

	delete[] Y_extnd;
	delete[] U_extnd;
	delete[] V_extnd;

	delete[] Y_vertical;
	delete[] U_vertical;
	delete[] V_vertical;

	delete[] Y_horizontal;
	delete[] U_horizontal;
	delete[] V_horizontal;

	delete[] Y_output;
	delete[] U_output;
	delete[] V_output;
}

void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	//YUV420 buffers
	uchar *Y_buff = new uchar[xSize * ySize];
	char *U_buff = new char[xSize * ySize / 4];
	char *V_buff = new char[xSize * ySize / 4];

	uchar *Y_output = new uchar[xSize * ySize];
	char *U_output = new char[xSize * ySize / 4];
	char *V_output = new char[xSize * ySize / 4];

	//Convertig from RGB to YUV420
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	double rotation_angle = PI * angle / 180.00;

	for (int i = 0; i < ySize; i++) {
		for (int j = 0; j < xSize; j++) {
			int new_Y = (int)(i * cos(rotation_angle) + j * sin(rotation_angle) 
				- m * sin(rotation_angle) - n * cos(rotation_angle) + n);
			int new_X = (int)(j * cos(rotation_angle) - i * sin(rotation_angle)
				- m * cos(rotation_angle) + n * sin(rotation_angle) + m);

			if (new_X < 0 || new_X > xSize - 1 || new_Y < 0 || new_Y > ySize - 1) {
				Y_output[i * xSize + j] = 0;
			}
			else {
				Y_output[i * xSize + j] = Y_buff[new_Y * xSize + new_X];
			}
		}
	}

	for (int i = 0; i < ySize / 2; i++) {
		for (int j = 0; j < xSize / 2; j++) {
			int new_Y = (int)(i * cos(rotation_angle) + j * sin(rotation_angle)
				- m/2 * sin(rotation_angle) - n/2 * cos(rotation_angle) + n/2);
			int new_X = (int)(j * cos(rotation_angle) - i * sin(rotation_angle)
				- m/2 * cos(rotation_angle) + n/2 * sin(rotation_angle) + m/2);

			if (new_X < 0 || new_X > xSize / 2 - 1 || new_Y < 0 || new_Y > ySize / 2 - 1) {
				U_output[i * xSize / 2 + j] = 0;
				V_output[i * xSize / 2 + j] = 0;
			}
			else {
				U_output[i * xSize / 2 + j] = U_buff[new_Y * xSize / 2 + new_X];
				V_output[i * xSize / 2 + j] = V_buff[new_Y * xSize / 2 + new_X];
			}
		}
	}

	YUV420toRGB(Y_output, U_output, V_output, xSize, ySize, output);
	
	delete[] Y_buff;
	delete[] Y_output;
	delete[] U_buff;
	delete[] U_output;
	delete[] V_buff;
	delete[] V_output;
}

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	//YUV420 buffers
	uchar *Y_buff = new uchar[xSize * ySize];
	char *U_buff = new char[xSize * ySize / 4];
	char *V_buff = new char[xSize * ySize / 4];

	uchar *Y_output = new uchar[xSize * ySize];
	char *U_output = new char[xSize * ySize / 4];
	char *V_output = new char[xSize * ySize / 4];

	uchar *Y_extnd = new uchar[(xSize + 2) * (ySize + 2)];
	char *U_extnd = new char[(xSize + 8) * (ySize + 8) / 4];
	char *V_extnd = new char[(xSize + 8) * (ySize + 8) / 4];

	//Convertig from RGB to YUV420
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	//Extending borders for exgde portection
	extendBorders(Y_buff, xSize, ySize, Y_extnd, 1);
	extendBorders(U_buff, xSize / 2, ySize / 2, U_extnd, 1);
	extendBorders(V_buff, xSize / 2, ySize / 2, V_extnd, 1);

	double rotation_angle = PI * angle / 180.00;

	for (int i = 0; i < ySize; i++) {
		for (int j = 0; j < xSize; j++) {
			double new_Y = (i * cos(rotation_angle) + j * sin(rotation_angle)
				- m * sin(rotation_angle) - n * cos(rotation_angle) + n);
			double new_X = (j * cos(rotation_angle) - i * sin(rotation_angle)
				- m * cos(rotation_angle) + n * sin(rotation_angle) + m);

			int first_X = floor(new_X);
			int first_Y = floor(new_Y);

			int second_X = first_X + 1;
			int second_Y = first_Y + 1;

			double a = new_Y - first_Y;
			double b = new_X - first_X;

			if (first_X < 0 || first_X > xSize - 1 || first_Y < 0 || first_Y > ySize - 1) {
				Y_output[i * xSize + j] = 0;
			}
			else {
				Y_output[i * xSize + j] = (1 - a) * (1 - b) * Y_extnd[first_Y * (xSize + 2) + first_X]
					+ (1 - a) * b * Y_extnd[first_Y * (xSize + 2) + second_X]
					+ a * (1 - b) * Y_extnd[second_Y * (xSize + 2) + first_X]
					+ a * b * Y_extnd[second_Y * (xSize + 2) + second_X];
			}
		}
	}

	for (int i = 0; i < ySize / 2; i++) {
		for (int j = 0; j < xSize / 2; j++) {
			double new_Y = (i * cos(rotation_angle) + j * sin(rotation_angle)
				- m / 2 * sin(rotation_angle) - n / 2 * cos(rotation_angle) + n / 2);
			double new_X = (j * cos(rotation_angle) - i * sin(rotation_angle)
				- m / 2 * cos(rotation_angle) + n / 2 * sin(rotation_angle) + m / 2);

			int first_X = floor(new_X);
			int first_Y = floor(new_Y);

			int second_X = first_X + 1;
			int second_Y = first_Y + 1;

			double a = new_Y - first_Y;
			double b = new_X - first_X;
			if (new_X < 0 || new_X > xSize / 2 - 1 || new_Y < 0 || new_Y > ySize / 2 - 1) {
				U_output[i * xSize / 2 + j] = 0;
				V_output[i * xSize / 2 + j] = 0;
			}
			else {
				U_output[i * xSize / 2 + j] = (1 - a) * (1 - b) * U_extnd[first_Y * (xSize / 2 + 2) + first_X]
					+ (1 - a) * b * U_extnd[first_Y * (xSize / 2 + 2) + second_X]
					+ a * (1 - b) * U_extnd[second_Y * (xSize / 2 + 2) + first_X]
					+ a * b * U_extnd[second_Y * (xSize / 2 + 2) + second_X];

				V_output[i * xSize / 2 + j] = (1 - a) * (1 - b) * V_extnd[first_Y * (xSize / 2 + 2) + first_X]
					+ (1 - a) * b * V_extnd[first_Y * (xSize / 2 + 2) + second_X]
					+ a * (1 - b) * V_extnd[second_Y * (xSize / 2 + 2) + first_X]
					+ a * b * V_extnd[second_Y * (xSize / 2 + 2) + second_X];
			}
		}
	}

	YUV420toRGB(Y_output, U_output, V_output, xSize, ySize, output);

	delete[] Y_buff;
	delete[] Y_output;
	delete[] U_buff;
	delete[] U_output;
	delete[] V_buff;
	delete[] V_output;


}