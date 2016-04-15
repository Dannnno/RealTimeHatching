#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>



/*
 * convolve with a box filter
 */
Image* ip_blur_box (Image* src, int size)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}

/*
 * convolve with a gaussian filter
 */
Image* ip_blur_gaussian (Image* src, int size, double sigma)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
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
    cerr << "This filter has not been implemented.\n";
    return NULL;
}

Image* ip_color_shift(Image* src)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}

/*
 * interpolate with the average intensity of the src image
 */
Image* ip_contrast (Image* src, double alpha)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}


/*
 * convolve an image with another image
 */
Image* ip_convolve (Image* src, int size, double** kernel)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}


/*
 * use a mask image for a per-pixel alpha value to perform
 * interpolation with a second image
 */
Image* ip_composite(Image* src1, Image* src2, Image* mask)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}


/*
 * cut away all but a subset of the image
 */
Image* ip_crop (Image* src, int x0, int y0, int x1, int y1)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}


/*
 * convolve with an edge detection kernel
 */
Image* ip_edge_detect (Image* src)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}


/*
 * create a new image with color values from one channel of src
 */
Image* ip_extract (Image* src, int channel)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}

/*
 * perform do what you like
 * you must query user for parameters here
 */
Image* ip_fun_warp (Image* src,int samplingMethod)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}


/*
 * create a new image with values equal to the psychosomatic intensities
 * of the source image
 */
Image* ip_grey (Image* src)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}

/*
 * shift image by dx, dy
 *
*/
Image* ip_image_shift(Image* src, int dx, int dy)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}

/*
 * interpolate an image with another image
 */
Image* ip_interpolate (Image* src1, Image* src2, double alpha)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}

/*
 * subtract the image from a white image
 */
Image* ip_invert (Image* src)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}

/*
 * median filter
 */
Image* ip_median (Image* src)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
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
    cerr << "This filter has not been implemented.\n";
    return NULL;
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
    cerr << "This filter has not been implemented.\n";
    return NULL;
}

/* helper functions you may find useful for resampling */

/*
 * nearest neighbor sample
 */

Pixel ip_resample_nearest(Image* src, double x, double y)
{
    cerr << "This filter has not been implemented.\n";
    return Pixel(0,0,0);
}

/*
 * bilinear resample
 */

Pixel ip_resample_bilinear(Image* src, double x, double y)
{

    cerr << "This filter has not been implemented.\n";
    return Pixel(0,0,0);
}

/*
 * gausian sample
 */
Pixel ip_resample_gaussian(Image* src, double x, double y, int size, double sigma) 
{
    cerr << "This filter has not been implemented.\n";
    return Pixel(0,0,0);
}


/*
 * rotate image using one of three sampling techniques
 */
Image* ip_rotate (Image* src, double theta, int x, int y, int mode, 
                  int size, double sigma)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}

/*
 * interpolate with the greyscale version of the image
 */
Image* ip_saturate (Image* src, double alpha)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}


/*
 * scale image using one of three sampling techniques
 */
Image* ip_scale (Image* src, double xFac, double yFac, int mode, 
                 int size, double sigma)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}

/*
 * create a new sepia tones image
 */
Image* ip_sepia (Image* src)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}


/*
 * create a new one bit/channel image with the intensity 0 or 1
 * depending on whether the input value is above or below the 
 * threshold
 */
Image* ip_threshold (Image* src, double cutoff)
{
    cerr << "This filter has not been implemented.\n";
    return NULL;
}




