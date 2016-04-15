#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "control.h"

using std::min;

/*
 * convolve with a box filter
 */
Image* ip_blur_box (Image* src, int size)
{
    double** kernel = new double*[size];
    double kernelValue = 1./(size*size);
    
    for (int i = 0; i < size; ++i) {
        kernel[i] = new double[size];
        for (int j = 0; j < size; ++j) {
            kernel[i][j] = kernelValue;
        }
    }
    
    Image* dest = ip_convolve(src, size, kernel);
    
    for (int i = 0; i < size; ++i) {
        delete[] kernel[i];
    }
    delete[] kernel;
    
    return dest;
}

double gaussian(int i, int j, double sigma)
{
    return exp(-(i*i + j*j)/(2*sigma*sigma));
}

/*
 * convolve with a gaussian filter
 */
Image* ip_blur_gaussian (Image* src, int size, double sigma)
{
    double** kernel = new double*[size];
    int halfEdge = (size - 1) / 2;
    
    double kernelSum = 0;
    
    for (int i = 0; i < size; ++i) {
        kernel[i] = new double[size];
        for (int j = 0; j < size; ++j) {
            kernel[i][j] = gaussian(i-halfEdge, j-halfEdge, sigma);
            kernelSum += kernel[i][j];
        }
    }
    
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            kernel[i][j] /= kernelSum;
        }
    }
    
    Image* dest = ip_convolve(src, size, kernel);
    
    for (int i = 0; i < size; ++i) {
        delete[] kernel[i];
    }
    delete[] kernel;
    
    return dest;
}


/*
 * convolve with a triangle filter
 */
Image* ip_blur_triangle (Image* src, int size)
{
	
    cerr << "This filter has not been implemented.\n";
    return NULL;
}


/*
 * interpolate with a black image
 */
Image* ip_brighten (Image* src, double alpha)
{
    Image blackImg = Image(*src);
    
    for (int i = 0; i < blackImg.getWidth(); ++i) {
        for (int j = 0; j < blackImg.getHeight(); ++j) {
            Pixel p = Pixel(0, 0, 0);
            blackImg.setPixel(i, j, p);
        }
    }

    Image* dest = ip_interpolate (src, &blackImg, alpha);
    return dest;
}

Image* ip_color_shift(Image* src)
{
    Image* dest = new Image(*src);

    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            Pixel p = Pixel(
                src->getPixel(i, j, GREEN),
                src->getPixel(i, j, BLUE),
                src->getPixel(i, j, RED)
                            );
            dest->setPixel(i, j, p);
        }
    }
    return dest;
}

/*
 * interpolate with the average intensity of the src image
 */
Image* ip_contrast (Image* src, double alpha)
{
    Image greyImg = Image(*src);
    
    for (int i = 0; i < greyImg.getWidth(); ++i) {
        for (int j = 0; j < greyImg.getHeight(); ++j) {
            Pixel p = Pixel(.5, .5, .5);
            greyImg.setPixel(i, j, p);
        }
    }
    
    Image* dest = ip_interpolate (src, &greyImg, alpha);
    return dest;
}


/*
 * convolve an image with another image
 */
Image* ip_convolve (Image* src, int size, double** kernel)
{
    Image* dest = new Image(*src);
    const int halfEdge = (size - 1)/2;
    
    const int width = dest->getWidth();
    const int height = dest->getHeight();
    
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            double destRed   = 0;
            double destGreen = 0;
            double destBlue  = 0;
            
            int dummyi = i - halfEdge;
            int dummyj = j - halfEdge;
            
            for (int ii = 0; ii < size; ++ii) {
                for (int jj = 0; jj < size; ++jj) {
                    if (dummyi + ii >= 0 && dummyi + ii < width &&
                        dummyj + jj >= 0 && dummyj + jj < height) {
                        
                        destRed   += kernel[ii][jj] * src->getPixel(dummyi + ii, dummyj + jj, RED);
                        destGreen += kernel[ii][jj] * src->getPixel(dummyi + ii, dummyj + jj, GREEN);
                        destBlue  += kernel[ii][jj] * src->getPixel(dummyi + ii, dummyj + jj, BLUE);
                    }
                }
            }
            
            Pixel p = Pixel(clamp(destRed, 0, 1),
                            clamp(destGreen, 0, 1),
                            clamp(destBlue, 0, 1)
                            );
            dest->setPixel(i, j, p);
        }
    }
    
    return dest;
}


/*
 * use a mask image for a per-pixel alpha value to perform
 * interpolation with a second image
 */
Image* ip_composite(Image* src1, Image* src2, Image* mask)
{
    Image* dest = new Image(*src1);
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            const double maskValue = mask->getPixel(i, j, RED);
            const double destRed   = src1->getPixel(i, j, RED) * maskValue + src2->getPixel(i, j, RED) * (1 - maskValue);
            const double destGreen = src1->getPixel(i, j, GREEN) * maskValue + src2->getPixel(i, j, GREEN) * (1 - maskValue);
            const double destBlue  = src1->getPixel(i, j, BLUE) * maskValue + src2->getPixel(i, j, BLUE) * (1 - maskValue);
            Pixel p = Pixel(clamp(destRed, 0, 1),
                            clamp(destGreen, 0, 1),
                            clamp(destBlue, 0, 1)
                            );
            dest->setPixel(i, j, p);
        }
    }
    
    return dest;
}


/*
 * cut away all but a subset of the image
 */
Image* ip_crop (Image* src, int x0, int y0, int x1, int y1)
{
    Image* dest = new Image(x1-x0, y1-y0);
    
    for (int i = x0; i < x1; ++i) {
        for (int j = y0; j < y1; ++j) {
            Pixel p = src->getPixel(i, j);
            dest->setPixel(i - x0, j - y0, p);
        }
    }
    return dest;
}


/*
 * convolve with an edge detection kernel
 */
Image* ip_edge_detect (Image* src)
{
    double** kernel = new double*[3];
    kernel[0] = new double[3]{ -1, -1, -1 };
    kernel[1] = new double[3]{ -1, 8, -1 };
    kernel[2] = new double[3]{ -1, -1, -1 };
    
    Image* dest = ip_convolve(src, 3, kernel);
    
    for (int i = 0; i < 3; ++i) {
        delete[] kernel[i];
    }
    delete[] kernel;
    
    return dest;
}


/*
 * create a new image with color values from one channel of src
 */
Image* ip_extract (Image* src, int channel)
{
    Image* dest = new Image(*src);
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            Pixel p = Pixel(
                            channel == RED ? src->getPixel(i, j, RED) : 0,
                            channel == GREEN ? src->getPixel(i, j, GREEN) : 0,
                            channel == BLUE ? src->getPixel(i, j, BLUE) : 0
                            );
            dest->setPixel(i, j, p);
        }
    }
    return dest;
}

/*
 * perform do what you like
 * you must query user for parameters here
 */
Image* ip_fun_warp (Image* src, int mode)
{
    Image* dest = new Image(*src);
    
    double sigma = 0;
    double size = 0;
    
    if (mode == I_GAUSSIAN) {
        sigma = getDouble("sigma");
        size = getDouble("size");
    }
    double strength = getDouble("strength");
    double denominator = getPositiveDouble("denominator");
    
    int width = dest->getWidth();
    int height = dest->getHeight();
    
    int cx = width/2;
    int cy = height/2;
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            double r = sqrt(sqr(i-cx) + sqr(j-cy));
            double a = atan((i-cx)/(j == cy ? 1 : j - cy));
            double rn = pow(r, strength)/denominator;
            
            double xp = rn*cos(a) + cx;
            double yp = rn*sin(a) + cy;
            
            if (xp >= src->getWidth() || xp < 0 || yp >= src->getHeight() || yp < 0) {
                continue;
            }
            
            Pixel p;
            switch (mode) {
                case I_NEAREST:
                    p = ip_resample_nearest(src, xp, yp);
                    break;
                case I_BILINEAR:
                    p = ip_resample_bilinear(src, xp, yp);
                    break;
                case I_GAUSSIAN:
                    p = ip_resample_gaussian(src, xp, yp, size, sigma);
                    break;
                default:
                    throw std::logic_error("This shouldn't happen");
            }
            
            dest->setPixel(i, j, p);
            
        }
    }
    
    return dest;
}


/*
 * create a new image with values equal to the psychosomatic intensities
 * of the source image
 */
Image* ip_grey (Image* src)
{
    Image* dest = new Image(*src);
    
    const double redWeight = .2126;
    const double greenWeight = .7152;
    const double blueWeight = .0722;
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            const double totalColor =
                src->getPixel(i, j, RED) * redWeight +
                src->getPixel(i, j, GREEN) * greenWeight +
                src->getPixel(i, j, BLUE) * blueWeight;
            Pixel p = Pixel(totalColor, totalColor, totalColor);
            dest->setPixel(i, j, p);
        }
    }
    return dest;
}

/*
 * shift image by dx, dy
 *
*/
Image* ip_image_shift(Image* src, int dx, int dy)
{
    Image* dest = new Image(*src);
    
    int width = dest->getWidth();
    int height = dest->getHeight();
    
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            Pixel p = src->getPixel(i, j);
            dest->setPixel((i + dx) % width, (j + dy) % height, p);
        }
    }

    return dest;
}

/*
 * interpolate an image with another image
 */
Image* ip_interpolate (Image* src1, Image* src2, double alpha)
{
    Image* dest = new Image(*src1);

    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            const double destRed   = clamp(alpha * src1->getPixel(i, j, RED)   + (1 - alpha) * src2->getPixel(i, j, RED), 0 , 1);
            const double destGreen = clamp(alpha * src1->getPixel(i, j, GREEN) + (1 - alpha) * src2->getPixel(i, j, GREEN), 0 , 1);
            const double destBlue  = clamp(alpha * src1->getPixel(i, j, BLUE)  + (1 - alpha) * src2->getPixel(i, j, BLUE), 0 , 1);

            Pixel p = Pixel(destRed, destGreen, destBlue);
            dest->setPixel(i, j, p);
        }
    }
    return dest;
}

/*
 * subtract the image from a white image
 */
Image* ip_invert (Image* src)
{
    Image greyImg = Image(*src);
    
    for (int i = 0; i < greyImg.getWidth(); ++i) {
        for (int j = 0; j < greyImg.getHeight(); ++j) {
            Pixel p = Pixel(.5, .5, .5);
            greyImg.setPixel(i, j, p);
        }
    }
    
    Image* dest = ip_interpolate (src, &greyImg, -1);
    return dest;
}

/*
 * median filter
 */
Image* ip_median (Image* src, int n)
{
    int size = n*n;
    double** neighborhood = new double*[3]{new double[size], new double[size], new double[size]};
    
    Image* dest = new Image(*src);
    
    const int width = dest->getWidth();
    const int height = dest->getHeight();
    
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            for (int k = 0; k < 3; ++k) {
                std::fill(&neighborhood[k][0], &neighborhood[k][size], std::numeric_limits<double>::max());
            }
            
            int dummyi = i - n;
            int dummyj = j - n;
            
            int count = 0;
            for (int ii = 0; ii < n; ++ii) {
                for (int jj = 0; jj < n; ++jj) {
                    if (dummyi + ii >= 0 && dummyi + ii < width &&
                        dummyj + jj >= 0 && dummyj + jj < height) {
                        
                        neighborhood[RED][ii*n + jj]   = src->getPixel(dummyi + ii, dummyj + jj, RED);
                        neighborhood[GREEN][ii*n + jj] = src->getPixel(dummyi + ii, dummyj + jj, GREEN);
                        neighborhood[BLUE][ii*n + jj]  = src->getPixel(dummyi + ii, dummyj + jj, BLUE);
                        ++count;
                    }
                }
            }
            
            std::sort(&neighborhood[RED][0], &neighborhood[RED][size]);
            std::sort(&neighborhood[GREEN][0], &neighborhood[GREEN][size]);
            std::sort(&neighborhood[BLUE][0], &neighborhood[BLUE][size]);
            
            Pixel p;
            if (count % 2 == 0) {
                double red = (neighborhood[RED][count/2] + neighborhood[RED][count/2+1])/2;
                double green = (neighborhood[GREEN][count/2] + neighborhood[GREEN][count/2+1])/2;
                double blue = (neighborhood[BLUE][count/2] + neighborhood[BLUE][count/2+1])/2;
                
                p = Pixel(clamp(red, 0, 1),
                          clamp(green, 0, 1),
                          clamp(blue, 0, 1)
                          );
            } else {                
                p = Pixel(clamp(neighborhood[RED][count/2+1], 0, 1),
                          clamp(neighborhood[GREEN][count/2+1], 0, 1),
                          clamp(neighborhood[BLUE][count/2+1], 0, 1)
                          );
            }
            
            dest->setPixel(i, j, p);
        }
    }

    
    for (int i = 0; i < 3; ++i) {
        delete[] neighborhood[i];
    }
    
    delete[] neighborhood;
    
    return dest;
}

/* misc
*/
Image* ip_misc(Image* src) 
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}

/*
 * round each pixel to the nearest value in the new number of bits
 */
Image* ip_quantize_simple (Image* src, int bitsPerChannel)
{
    Image* dest = new Image(src->getWidth(), src->getHeight(), bitsPerChannel);
    
    for (int i = 0; i < dest->getHeight(); ++i) {
        for (int j = 0; j < dest->getWidth(); ++j) {
            Pixel p = src->getPixel(j, i);
            dest->setPixel(j, i, p);
        }
    }
    
    return dest;
}


/*
 * dither each pixel to the nearest value in the new number of bits
 * using a static 2x2 matrix
 */
Image* ip_quantize_ordered (Image* src, int bitsPerChannel)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}


/*
 * dither each pixel to the nearest value in the new number of bits
 * using error diffusion
 */
Image* ip_quantize_fs (Image* src, int bitsPerChannel)
{
    int size = src->getWidth() * src->getHeight();
    double** imageError = new double*[3]{new double[size], new double[size], new double[size]};
    for (int i = 0; i < 3; ++i) {
        std::fill(&imageError[i][0], &imageError[i][size], 0);
    }
    
    Image* dest = new Image(src->getWidth(), src->getHeight(), bitsPerChannel);
    
    double alpha = 7./16;
    double beta = 13./16;
    double gamma = 5./16;
    double delta = 1./16;
    
    int width = dest->getWidth();
    int height = dest->getHeight();
    
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            Pixel oldPixel = src->getPixel(i, j);
            Pixel newPixel(
                           clamp(oldPixel.getColor(RED) + imageError[RED][i*width + j], 0, 1),
                           clamp(oldPixel.getColor(GREEN) + imageError[GREEN][i*width + j], 0, 1),
                           clamp(oldPixel.getColor(BLUE) + imageError[BLUE][i*width + j], 0, 1)
            );
            dest->setPixel(i, j, newPixel);
            Pixel bitAdjustedNewPixel = dest->getPixel(i, j);
            
            int redError   = bitAdjustedNewPixel.getColor(RED) - oldPixel.getColor(RED);
            int greenError = bitAdjustedNewPixel.getColor(GREEN) - oldPixel.getColor(GREEN);
            int blueError  = bitAdjustedNewPixel.getColor(BLUE) - oldPixel.getColor(BLUE);
            
            if (i+1 < width) {
                imageError[RED][(i+1)*width + j]   += alpha * redError;
                imageError[GREEN][(i+1)*width + j] += alpha * greenError;
                imageError[BLUE][(i+1)*width + j]  += alpha * blueError;
            }
            
            if (i - 1 > 0 && j + 1 < height) {
                imageError[RED][(i-1)*width + j + 1]   += beta * redError;
                imageError[GREEN][(i-1)*width + j + 1] += beta * greenError;
                imageError[BLUE][(i-1)*width + j + 1]  += beta * blueError;
            }
            
            if (j + 1 < height) {
                imageError[RED][i*width + j + 1]   += gamma * redError;
                imageError[GREEN][i*width + j + 1] += gamma * greenError;
                imageError[BLUE][i*width + j + 1]  += gamma * blueError;
            }
            
            if (i + 1 < width && j + 1 < height) {
                imageError[RED][(i+1)*width + j + 1]   += delta * redError;
                imageError[GREEN][(i+1)*width + j + 1] += delta * greenError;
                imageError[BLUE][(i+1)*width + j + 1]  += delta * blueError;
            }
        }
    }
    
    for (int i = 0; i < 3; ++i) {
        delete[] imageError[i];
    }
    
    delete[] imageError;
    
    return dest;
}

/* helper functions you may find useful for resampling */

/*
 * nearest neighbor sample
 */
Pixel ip_resample_nearest(Image* src, double x, double y)
{
    int ix = round(x);
    int iy = round(y);
    
    if (ix >= src->getWidth() || ix < 0 || iy >= src->getHeight() || iy < 0) {
        return Pixel(0, 0, 0);
    }
    return src->getPixel(round(x), round(y));
}

/*
 * bilinear resample
 */

Pixel ip_resample_bilinear(Image* src, double x, double y)
{
    int xfloor = floor(x);
    int yfloor = floor(y);
    
    double alpha = x - xfloor;
    double beta = y - yfloor;
    
    if (xfloor >= src->getWidth() || xfloor < 0 || xfloor+1 >= src->getWidth() || xfloor+1 < 0 ||
        yfloor >= src->getHeight() || yfloor < 0 || yfloor+1 >= src->getHeight() || yfloor+1 < 0) {
        return Pixel(0, 0, 0);
    }
    
    Pixel p1 = src->getPixel(xfloor, yfloor);
    Pixel p2 = src->getPixel(xfloor+1, yfloor);
    Pixel p3 = src->getPixel(xfloor, yfloor+1);
    Pixel p4 = src->getPixel(xfloor+1, yfloor+1);
    
    Pixel p5 = Pixel(
                     clamp(p1.getColor(RED) * (1-alpha) + p2.getColor(RED) * alpha, 0, 1),
                     clamp(p1.getColor(GREEN) * (1-alpha) + p2.getColor(GREEN) * alpha, 0, 1),
                     clamp(p1.getColor(BLUE) * (1-alpha) + p2.getColor(BLUE) * alpha, 0, 1)
                     );
    Pixel p6 = Pixel(
                     clamp(p3.getColor(RED) * (1-alpha) + p4.getColor(RED) * alpha, 0, 1),
                     clamp(p3.getColor(GREEN) * (1-alpha) + p4.getColor(GREEN) * alpha, 0, 1),
                     clamp(p3.getColor(BLUE) * (1-alpha) + p4.getColor(BLUE) * alpha, 0, 1)
                     );
    
    Pixel newPixel = Pixel(
                           clamp(p5.getColor(RED) * (1-beta) + p6.getColor(RED) * beta, 0, 1),
                           clamp(p5.getColor(GREEN) * (1-beta) + p6.getColor(GREEN) * beta, 0, 1),
                           clamp(p5.getColor(BLUE) * (1-beta) + p6.getColor(BLUE) * beta, 0, 1)
                           );
    
    return newPixel;
}

/*
 * gausian sample
 */
Pixel ip_resample_gaussian(Image* src, double x, double y, int size, double sigma) 
{
    int ix = round(x);
    int iy = round(y);
    
    if (ix >= src->getWidth() || ix < 0 || iy >= src->getHeight() || iy < 0) {
        return Pixel(0, 0, 0);
    }
    
    double** kernel = new double*[size];
    int halfEdge = size / 2;
    
    double kernelSum = 0;
    
    for (int i = 0; i < size; ++i) {
        kernel[i] = new double[size];
        for (int j = 0; j < size; ++j) {
            kernel[i][j] = gaussian(i-halfEdge, j-halfEdge, sigma);
            kernelSum += kernel[i][j];
        }
    }
    
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            kernel[i][j] /= kernelSum;
        }
    }

    int width = src->getWidth();
    int height = src->getHeight();
    
    double destRed   = 0;
    double destGreen = 0;
    double destBlue  = 0;
    
    int dummyi = ix - halfEdge;
    int dummyj = iy - halfEdge;
    
    for (int ii = 0; ii < size; ++ii) {
        for (int jj = 0; jj < size; ++jj) {
            if (dummyi + ii >= 0 && dummyi + ii < width &&
                dummyj + jj >= 0 && dummyj + jj < height) {
                
                destRed   += kernel[ii][jj] * src->getPixel(dummyi + ii, dummyj + jj, RED);
                destGreen += kernel[ii][jj] * src->getPixel(dummyi + ii, dummyj + jj, GREEN);
                destBlue  += kernel[ii][jj] * src->getPixel(dummyi + ii, dummyj + jj, BLUE);
            }
        }
    }
    
    Pixel p = Pixel(clamp(destRed, 0, 1),
                    clamp(destGreen, 0, 1),
                    clamp(destBlue, 0, 1)
                    );
    
    for (int i = 0; i < 3; ++i) {
        delete[] kernel[i];
    }
    delete[] kernel;

    return p;
}


/*
 * rotate image using one of three sampling techniques
 */
Image* ip_rotate (Image* src, double theta, int x, int y, int mode, 
                  int size, double sigma)
{
    Image* dest = new Image(src->getWidth(), src->getHeight());
    
    double ctheta = cos(-theta*M_PI/180.);
    double stheta = sin(-theta*M_PI/180.);
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            Pixel p = Pixel(0, 0, 0);
            dest->setPixel(i, j, p);
        }
    }
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            double xp = x + (i-x)*ctheta - (j-y)*stheta;
            double yp = y + (i-x)*stheta + (j-y)*ctheta;

            if (xp >= src->getWidth() || xp < 0 || yp >= src->getHeight() || yp < 0) {
                continue;
            }
            
            Pixel p;
            switch (mode) {
                case I_NEAREST:
                    p = ip_resample_nearest(src, xp, yp);
                    break;
                case I_BILINEAR:
                    p = ip_resample_bilinear(src, xp, yp);
                    break;
                case I_GAUSSIAN:
                    p = ip_resample_gaussian(src, xp, yp, size, sigma);
                    break;
                default:
                    throw std::logic_error("This shouldn't happen");
            }
            
            dest->setPixel(i, j, p);
        }
    }
    
    return dest;
}

/*
 * interpolate with the greyscale version of the image
 */
Image* ip_saturate (Image* src, double alpha)
{
    Image* greyImg = ip_grey(src);
    Image* dest = ip_interpolate (src, greyImg, alpha);
    delete greyImg;
    return dest;
}


/*
 * scale image using one of three sampling techniques
 */
Image* ip_scale (Image* src, double xFac, double yFac, int mode, 
                 int size, double sigma)
{
    Image* dest = new Image(src->getWidth()*xFac, src->getHeight()*yFac);
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            double xp = i*1./xFac;
            double yp = j*1./yFac;
            
            if (xp >= src->getWidth() || xp < 0 || yp >= src->getHeight() || yp < 0) {
                continue;
            }
            
            Pixel p;
            switch (mode) {
                case I_NEAREST:
                    p = ip_resample_nearest(src, xp, yp);
                    break;
                case I_BILINEAR:
                    p = ip_resample_bilinear(src, xp, yp);
                    break;
                case I_GAUSSIAN:
                    p = ip_resample_gaussian(src, xp, yp, size, sigma);
                    break;
                default:
                    throw std::logic_error("This shouldn't happen");
            }
            
            dest->setPixel(i, j, p);
        }
    }
    
    return dest;
}

/*
 * create a new sepia tones image
 */
Image* ip_sepia (Image* src)
{
    Image* dest = new Image(*src);
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            const double srcRed = src->getPixel(i, j, RED);
            const double srcGreen = src->getPixel(i, j, GREEN);
            const double srcBlue = src->getPixel(i, j, BLUE);
            
            double destRed = srcRed*.393 + srcGreen*.769 + srcBlue*.189;
            double destGreen = srcRed*.349 + srcGreen*.686 + srcBlue*.168;
            double destBlue = srcRed*.272 + srcGreen*.534 + srcBlue*.131;
            
            Pixel p = Pixel(
                            clamp(destRed, 0, 1),
                            clamp(destGreen, 0, 1),
                            clamp(destBlue, 0, 1)
                            );
            dest->setPixel(i, j, p);
        }
    }
    return dest;
}


/*
 * create a new one bit/channel image with the intensity 0 or 1
 * depending on whether the input value is above or below the 
 * threshold
 */
Image* ip_threshold (Image* src, double cutoff)
{
    Image* dest = new Image(*src);
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            Pixel p = Pixel(
                            src->getPixel(i, j, RED) > cutoff ? 1 : 0,
                            src->getPixel(i, j, GREEN) > cutoff ? 1 : 0,
                            src->getPixel(i, j, BLUE) > cutoff ? 1 : 0
                            );
            dest->setPixel(i, j, p);
        }
    }
    return dest;
}




