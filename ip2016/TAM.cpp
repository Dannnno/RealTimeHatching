//
//  TAM.cpp
//  ip2016
//
//  Created by swarren on 4/15/16.
//  Copyright © 2016 Z Sweedyk. All rights reserved.
//

#include "TAM.hpp"

TAM::TAM(int numTones, int numRes, Image* stroke)
: images(numTones, std::vector<Image*>(numRes)) {
    std::vector<Image*> candidates(numRes);
    
    
    RandomStroke bestStroke;
    int score;
    
    int size = 512;
    
    Image* square = new Image (size, size);
    drawStroke(stroke, getRandomStroke(), square);
    images[0][0] = square;
    
}



TAM::RandomStroke TAM::getRandomStroke() const{
    static std::uniform_real_distribution<double> pos_dist(0, 1);
    static std::uniform_real_distribution<double> length_dist(.3, 1);
    static std::uniform_int_distribution<int> rot_dist(-2, 2);
    
    RandomStroke stroke;
    
    stroke.length = length_dist(rng);
    stroke.x = pos_dist(rng);
    stroke.y = pos_dist(rng);
    stroke.rot = deg2rad(rot_dist(rng));
    
    return stroke;
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
 * rotate image using one of three sampling techniques
 */
Image* ip_rotate (Image* src, double theta)
{
    Image* dest = new Image(src->getWidth(), src->getHeight());
    double x = src->getWidth()/2;
    double y = src->getHeight()/2;
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
            
            Pixel p = ip_resample_bilinear(src, xp, yp);
            
            dest->setPixel(i, j, p);
        }
    }
    
    return dest;
}

/*
 * scale image using one of three sampling techniques
 */
Image* ip_scale (Image* src, double xFac, double yFac)
{
    Image* dest = new Image(src->getWidth()*xFac, src->getHeight()*yFac);
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            double xp = i*1./xFac;
            double yp = j*1./yFac;
            
            if (xp >= src->getWidth() || xp < 0 || yp >= src->getHeight() || yp < 0) {
                continue;
            }
            
            Pixel p = ip_resample_bilinear(src, xp, yp);
            
            dest->setPixel(i, j, p);
        }
    }
    
    return dest;
}

/*
 * use a mask image for a per-pixel alpha value to perform
 * interpolation with a second image
 */
void ip_composite(Image* dest, Image* strokeImage, double x, double y)
{
    static const double whiteThreshold = 0.9;
    for (int i = 0; i < strokeImage->getWidth(); ++i) {
        for (int j = 0; j < strokeImage->getHeight(); ++j) {
            Pixel strokePixel = strokeImage->getPixel(i, j);
            const bool isWhite = strokePixel.getColor(RED) > whiteThreshold;
            const double destRed   = isWhite ? dest->getPixel(i, j, RED) : strokeImage->getPixel(i, j, RED);
            const double destGreen = isWhite ? dest->getPixel(i, j, RED) : strokeImage->getPixel(i, j, RED);
            const double destBlue  = isWhite ? dest->getPixel(i, j, RED) : strokeImage->getPixel(i, j, RED);
            Pixel p = Pixel(destRed, destGreen, destBlue);
            dest->setPixel(x + i, y + j, p);
        }
    }
}

void TAM::drawStroke(Image* strokeImage, RandomStroke strokeModifiers, Image* dest){
    Image* scaledImage = ip_scale(strokeImage, strokeModifiers.length, 1);
    Image* rotImage = ip_rotate(scaledImage, strokeModifiers.rot);
    ip_composite(dest, strokeImage, strokeModifiers.x*dest->getWidth(), strokeModifiers.y*dest->getHeight());
    delete scaledImage;
    delete rotImage;
}