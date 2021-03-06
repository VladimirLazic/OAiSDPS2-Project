#include "ImageProcessing.h"
#include "ColorSpaces.h"
#include "ImageFilter.h"
#include "NoiseReduction.h"
#include "ImageInterpolation.h"
#include <math.h>

#include <QDebug>

int round_by(int x, int y) {
	return (x + y - 1) & ~(y - 1);
}

void imageProcessingFun(const QString& progName, QImage* const outImgs, const QImage* const inImgs, const QVector<double>& params) 
{
	int X_SIZE = inImgs->width();
	int Y_SIZE = inImgs->height();
		
	if(progName == "Sample and hold") 
	{	
		double vertical_factor = params[0], horizontal_factor = params[1];
		int newX_SIZE = round_by(horizontal_factor * X_SIZE , 4);
		int newY_SIZE = round_by(vertical_factor * Y_SIZE , 4);

		/* Create empty output image */
		*outImgs = *(new QImage(newX_SIZE, newY_SIZE, inImgs->format()));

		sampleAndHold(inImgs->bits() , X_SIZE , Y_SIZE , outImgs->bits() , newX_SIZE , newY_SIZE);
	}
	else if (progName == "Bilinear") 
	{
		double vertical_factor = params[0], horizontal_factor = params[1];
		int newX_SIZE = round_by(horizontal_factor * X_SIZE, 4);
		int newY_SIZE = round_by(vertical_factor * Y_SIZE, 4);

		/* Create empty output image */
		*outImgs = *(new QImage(newX_SIZE, newY_SIZE, inImgs->format()));

		/* Perform Bilinear interpolation  */
		bilinearInterpolate(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), newX_SIZE, newY_SIZE);
	}
	else if (progName == "Bicubic")
	{
		double vertical_factor = params[0], horizontal_factor = params[1];
		int newX_SIZE = round_by(horizontal_factor * X_SIZE, 4);
		int newY_SIZE = round_by(vertical_factor * Y_SIZE, 4);

		/* Create empty output image */
		*outImgs = *(new QImage(newX_SIZE, newY_SIZE, inImgs->format())); 

		/* Perform Bicubic interpolation  */
		bicubicInterpolate(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), newX_SIZE, newY_SIZE);
	}
	else if(progName == "Rotation") 
	{	
		*outImgs = *(new QImage(X_SIZE, Y_SIZE, inImgs->format()));
		imageRotate(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), X_SIZE / 2, Y_SIZE / 2, params[0]);	
	}
	else if (progName == "Rotation Bilinear") 
	{
		*outImgs = *(new QImage(X_SIZE, Y_SIZE, inImgs->format()));
		imageRotateBilinear(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), X_SIZE / 2, Y_SIZE / 2, params[0]);
	}

}

