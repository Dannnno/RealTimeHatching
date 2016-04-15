#include "image.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <setjmp.h>
#include "SOIL.h"



static const int maxvals[9] = { 0, 1, 3, 7, 15, 31, 63, 127, 255 };

/*  val2bits												*/
/*  converts a double v into an integer in range [0,255].	*/
/*	rounding is done so as to ensure that only				*/
/*  only 2^b - 1 choices in the range [0,255] are possible.	*/

static inline unsigned char val2bits (double v, int b)
{

	return (unsigned char) (floor(v * maxvals[b] + 0.5) / maxvals[b] * 255);
}

Pixel::Pixel()
{
	col[0]=0;
	col[1]=0;
	col[2]=0;
}
Pixel::Pixel(double r, double g, double b)
{
	col[0]=r;
	col[1]=g;
	col[2]=b;
}



Pixel::~Pixel()
{
}

// copy constructor
Pixel::Pixel (const Pixel& toCopy)
{
	col[0] = toCopy.col[0];
	col[1] = toCopy.col[1];
	col[2] = toCopy.col[2];
}

// assignment operator
Pixel& Pixel::operator = (const Pixel& toCopy)
{

	if (this != &toCopy) {
		col[0] = toCopy.col[0];
		col[1] = toCopy.col[1];
		col[2] = toCopy.col[2];
	}
	return *this;
}


double Pixel::getColor(int chan)
{
	assert (chan>=0 && chan <=2);
	return col[chan];
}
void Pixel::setColor(int chan, double val)
{
	assert (chan>=0 && chan <=2);
	val=max(val,0);
	val=min(val,1);
	col[chan]=val;
}

// write
ostream &operator<<(ostream &out_file, Pixel& thePixel)
{
	out_file << thePixel.getColor(0) << ", " << thePixel.getColor(1) << ", " << thePixel.getColor(2);
	return out_file;
}


Image::Image ()
{
	width    = 0;
	height   = 0;
	bits     = 0;
	maxValue = 0;
	pixels   = NULL;
}


Image::Image (int width_, int height_)
{
	width    = width_;
	height   = height_;

	bits     = 8;
	maxValue = 255;

	// note we don't handle alpha channels
	assert(width > 0 && height > 0 &&
		bits > 0 && bits < 9
		);

	pixels   = new unsigned char[width*height*3];;
	memset(pixels, 0, width*height*3);
}


Image::Image (int width_, int height_, int bits_)
{
	width    = width_;
	height   = height_;
	bits     = bits_;
	maxValue = 255;

	assert(width > 0 && height > 0 &&
		bits > 0 && bits < 9
		);


	pixels   = new unsigned char[width*height*3];

	memset(pixels, 0, width*height*3);
}




Image::Image (const char* filename)
{
	width    = 0;
	height   = 0;
	bits     = 0;
	maxValue = 0;
	pixels   = NULL;

	readBMP(filename);
}


Image::~Image ()
{
	if (pixels) delete[] pixels;
}


Image::Image (const Image& image)
{
	width    = image.width;
	height   = image.height;
	bits     = image.bits;
	maxValue = image.maxValue;
	pixels   = new unsigned char[width*height*3];

	for (int i = 0; i < width*height*3; ++i)
		pixels[i] = image.pixels[i];
}


Image& Image::operator= (const Image& image)
{
	if (&image == this) return *this;

	if (pixels) delete[] pixels;

	width    = image.width;
	height   = image.height;
	bits     = image.bits;
	maxValue = image.maxValue;
	pixels   = new unsigned char[width*height*3];

	for (int i = 0; i < width*height*3; ++i)
		pixels[i] = image.pixels[i];

	return *this;
}


bool Image::good ()
{
	return (width > 0 && height > 0 &&
		bits > 0 && bits < 9 && pixels);
}


bool Image::bad ()
{
	return !good();
}


void Image::clear ()
{
	memset(pixels, 0, width*height*3);
}


int Image::index (int x, int y, int c)
{
	return (((height - y - 1) * width + x) * 3 + c);
}


double Image::getPixel (int x, int y, int channel)
{
	assert(good());
	assert((x >= 0)       &&
		(x <  width)   &&
		(y >= 0)       &&
		(y <  height)  &&
		(channel >= 0) &&
		(channel < 3));

	return pixels[index(x,y,channel)] / 255.0;
}


double Image::getPixel_ (int x, int y, int channel)
{
	if (!good()        ||
		(x <  0)       ||
		(x >= width)   ||
		(y <  0)       ||
		(y >= height)  ||
		(channel < 0)  ||
		(channel >= 3))
		return 0.0;

	return getPixel(x,y,channel);
}


Pixel Image::getPixel (int x, int y)
{
	assert(good());
	assert((x >= 0)       &&
		(x <  width)   &&
		(y >= 0)       &&
		(y <  height));

	Pixel pixel;


	pixel.setColor(BLUE, pixels[index(x,y,BLUE)]  / 255.0);
	pixel.setColor(GREEN, pixels[index(x,y,GREEN)] / 255.0);
	pixel.setColor(RED,pixels[index(x,y,RED)]   / 255.0);

	return pixel;
}


Pixel Image::getPixel_ (int x, int y)
{
	if (!good()        ||
		(x <  0)       ||
		(x >= width)   ||
		(y <  0)       ||
		(y >= height))
	{
		Pixel pixel;
		memset(&pixel, 0, sizeof(Pixel));
		return pixel;
	}

	return getPixel(x,y);
}


Pixel& Image::getPixel (int x, int y, Pixel& pixel)
{
	assert(good());
	assert((x >= 0)       &&
		(x <  width)   &&
		(y >= 0)       &&
		(y <  height));


	pixel.setColor(BLUE,pixels[index(x,y,BLUE)]  / 255.0);
	pixel.setColor(GREEN,pixels[index(x,y,GREEN)] / 255.0);
	pixel.setColor(RED, pixels[index(x,y,RED)]   / 255.0);

	return pixel;
}


Pixel& Image::getPixel_ (int x, int y, Pixel& pixel)
{
	if (!good()        ||
		(x <  0)       ||
		(x >= width)   ||
		(y <  0)       ||
		(y >= height))
	{
		return pixel;
	}

	return getPixel(x,y,pixel);
}


/* set pixel value.  throw error if any problem with parameters */
/* including value outside of [0,1] */
void Image::setPixel (int x, int y, int channel, double value)
{
	assert(good());
	assert((x >= 0)       &&
		(x <  width)   &&
		(y >= 0)       &&
		(y <  height)  &&
		(channel >= 0) &&
		(channel < 3));
	assert((value >= 0.0) && 
		(value <= 1.0));  

	pixels[index(x,y,channel)] = val2bits(value, bits);
}


/*  set pixel clamps value clamping to [0,1]	*/
/*  if any problem with parameters it RETURNS	*/
void Image::setPixel_ (int x, int y, int channel, double value)
{
	if (!good()        ||
		(x <  0)       ||
		(x >= width)   ||
		(y <  0)       ||
		(y >= height)  ||
		(channel < 0)  ||
		(channel >= 3))
		return;
	value=clamp(value, 0, 1);
	setPixel(x,y,channel,value);
}



/* set pixel;  throws error if problems with parameters */
/* including pixel value outside of [0,1]			    */
void Image::setPixel (int x, int y, Pixel& pixel)
{
	assert(good());
	assert((x >= 0)       &&
		(x <  width)   &&
		(y >= 0)       &&
		(y <  height));


	assert((pixel.getColor(BLUE) >= 0.0) && 
		(pixel.getColor(BLUE)<= 1.0));  
	assert((pixel.getColor(GREEN) >= 0.0) && 
		(pixel.getColor(GREEN) <= 1.0));  
	pixels[index(x,y,BLUE)]  = val2bits(pixel.getColor(BLUE), bits);
	pixels[index(x,y,GREEN)] = val2bits(pixel.getColor(GREEN), bits);

	assert((pixel.getColor(RED) >= 0.0) && 
		(pixel.getColor(RED) <= 1.0));  
	pixels[index(x,y,RED)]   = val2bits(pixel.getColor(RED), bits);


}


/*  set pixel value, clamping to [0,1]    */
/*  if any problem with parameters RETURN */
void Image::setPixel_ (int x, int y, Pixel& pixel)
{
	if (!good()        ||
		(x <  0)       ||
		(x >= width)   ||
		(y <  0)       ||
		(y >= height))
		return;
	pixel.setColor(RED,clamp(pixel.getColor(RED),0,1));
	pixel.setColor(GREEN,clamp(pixel.getColor(GREEN),0,1));
	pixel.setColor(BLUE,clamp(pixel.getColor(BLUE),0,1));
	setPixel(x,y,pixel);
}


void Image::glReadPixelsWrapper ()
{
	assert(good());
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

}


void Image::glDrawPixelsWrapper ()
{
	assert(good());

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);


	glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);


}


void Image::glTexImage2DWrapper ()
{
	assert(good());

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);


	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0,
		GL_RGB, GL_UNSIGNED_BYTE, pixels);

}


int Image::readBMP (const char* filename)
{
    // check for png
    unsigned long len = strlen(filename);
    const char* ext = &(filename[len-4]);
    
    // to change to read png as well replace the if condition with the following:
    //if (strncmp(ext, ".bmp", 4) == 0 || strncmp(ext, ".png", 4) == 0 ) {
    if (strncmp(ext, ".bmp", 4) == 0) {
        int w, h, c;
        unsigned char* buffer = SOIL_load_image(filename, &w, &h, &c, 3);
        if (w>0 && h>0 && c==3 && buffer) {
            if (pixels) {
                delete[] pixels;
            }
            pixels = new unsigned char[w*h*c];
            
            // flip
            for (int j=0; j<h/2; j++) {
                for (int i=0; i<w; i++) {
                    for (int ch=0; ch<c; ch++) {
                        char tmp = buffer[c*(i + j*w) + ch];
                        buffer[c*(i + j*w) + ch] = buffer[c*(i + (h-1-j)*w) + ch];
                        buffer[c*(i + (h-1-j)*w) + ch] = tmp;
                    }
                }
            }
            pixels = buffer;
            width = w;
            height = h;
            bits = 8;
            pixels = buffer;
            return 1;
        }
        
    }
    return -1;
}


int Image::writeBMP (const char* filename)
{
	unsigned long len = strlen(filename);
	const char* ext = &(filename[len-4]);

	if (strncmp(ext, ".bmp", 4) == 0)
	{
		cerr << "Writing BMP " << filename << endl;
        unsigned char* buffer= new unsigned char[width*height*bits];
        for (int j=0; j<(height+1)/2; j++) {
            for (int i=0; i<width; i++) {
                for (int c=0; c<3; c++) {
                    buffer[3*(i + j*width) + c] = pixels[3*(i + (height-1-j)*width) + c];
                    buffer[3*(i + (height-1-j)*width) + c] = pixels[3*(i+j*width)+c];
                }
            }
        }

        SOIL_save_image (filename, SOIL_SAVE_TYPE_BMP, width, height, 3, buffer);
        delete[] buffer;
        
		//return writeBMP(filename);
        return 1;
	}
	else
		return -1;
}


//#if !defined(WIN32) || defined(_CONSOLE)

typedef unsigned char BYTE;		/* 8 bits */
typedef unsigned short int WORD;	/* 16-bit unsigned integer. */
typedef unsigned int DWORD;		/* 32-bit unsigned integer */
typedef int LONG;			/* 32-bit signed integer */



typedef struct tagBITMAPFILEHEADER {
	WORD bfType;
	DWORD bfSize;
	WORD bfReserved1;
	WORD bfReserved2;
	DWORD bfOffBits;
} BITMAPFILEHEADER;



typedef struct tagBITMAPINFOHEADER {
	DWORD biSize;
	LONG biWidth;
	LONG biHeight;
	WORD biPlanes;
	WORD biBitCount;
	DWORD biCompression;
	DWORD biSizeImage;
	LONG biXPelsPerMeter;
	LONG biYPelsPerMeter;
	DWORD biClrUsed;
	DWORD biClrImportant;
} BITMAPINFOHEADER;



/* constants for the biCompression field */
#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L



typedef struct tagRGBTRIPLE {
	BYTE rgbtBlue;
	BYTE rgbtGreen;
	BYTE rgbtRed;
} RGBTRIPLE;



typedef struct /*tagRGBQUAD*/ {
	BYTE rgbBlue;
	BYTE rgbGreen;
	BYTE rgbRed;
	BYTE rgbReserved;
} RGBQUAD;


//#endif // !defined(WIN32) || defined(_CONSOLE)


/* Some magic numbers */

#define BMP_BF_TYPE 0x4D42
/* word BM */

#define BMP_BF_OFF_BITS 54
/* 14 for file header + 40 for info header (not sizeof(), but packed size) */

#define BMP_BI_SIZE 40
/* packed size of info header */


/* Reads a WORD from a file in little endian format */
static WORD WordReadLE(ifstream& fp)
{
	WORD lsb, msb;

	lsb = fp.get();
	msb = fp.get();
	return (msb << 8) | lsb;
}



/* Writes a WORD to a file in little endian format */
static void WordWriteLE(WORD x, ofstream& fp)
{
	BYTE lsb, msb;

	lsb = (BYTE) (x & 0x00FF);
	msb = (BYTE) (x >> 8);
	fp.put(lsb);
	fp.put(msb);
}



/* Reads a DWORD word from a file in little endian format */
static DWORD DWordReadLE(ifstream& fp)
{
	DWORD b1, b2, b3, b4;

	b1 = fp.get();
	b2 = fp.get();
	b3 = fp.get();
	b4 = fp.get();
	return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



/* Writes a DWORD to a file in little endian format */
static void DWordWriteLE(DWORD x, ofstream& fp)
{
	unsigned char b1, b2, b3, b4;

	b1 = (x & 0x000000FF);
	b2 = ((x >> 8) & 0x000000FF);
	b3 = ((x >> 16) & 0x000000FF);
	b4 = ((x >> 24) & 0x000000FF);
	fp.put(b1);
	fp.put(b2);
	fp.put(b3);
	fp.put(b4);
}



/* Reads a LONG word from a file in little endian format */
static LONG LongReadLE(ifstream& fp)
{
	LONG b1, b2, b3, b4;

	b1 = fp.get();
	b2 = fp.get();
	b3 = fp.get();
	b4 = fp.get();
	return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



/* Writes a LONG to a file in little endian format */
static void LongWriteLE(LONG x, ofstream& fp)
{
	char b1, b2, b3, b4;

	b1 = (x & 0x000000FF);
	b2 = ((x >> 8) & 0x000000FF);
	b3 = ((x >> 16) & 0x000000FF);
	b4 = ((x >> 24) & 0x000000FF);
	fp.put(b1);
	fp.put(b2);
	fp.put(b3);
	fp.put(b4);
}


int bitcount (DWORD w)
{
	w = (0x55555555 & w) + (0x55555555 & (w>> 1));
	w = (0x33333333 & w) + (0x33333333 & (w>> 2));
	w = (0x0f0f0f0f & w) + (0x0f0f0f0f & (w>> 4));
	w = (0x00ff00ff & w) + (0x00ff00ff & (w>> 8));
	w = (0x0000ffff & w) + (0x0000ffff & (w>>16));
	return w;
}



