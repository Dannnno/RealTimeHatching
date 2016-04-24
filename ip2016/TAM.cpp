//
//  TAM.cpp
//  ip2016
//
//  Created by swarren on 4/15/16.
//  Copyright © 2016 Z Sweedyk. All rights reserved.
//

#include "TAM.hpp"
std::mt19937 rng;

double getActualStrokeLength(double actualPosition, double actualLength, double imageSize) {
    double normalizedStrokeLength;
    if(actualPosition < 0){
        normalizedStrokeLength = max(actualPosition + actualLength, 0)/2;
    }
    else if (actualPosition + actualLength >= imageSize) {
        normalizedStrokeLength = imageSize - actualPosition;
    } else {
        normalizedStrokeLength = actualLength;
    }
    return normalizedStrokeLength / imageSize;
}

TAM::TAM(int numTones, int numRes, Image* stroke)
: images(numTones, std::vector<Image*>(numRes)) {
    std::vector<Image> candidates(numRes);
    
    for (int tone = 0; tone < numTones; ++tone) {
        int imageSize = 32;
        for (int resolution = 0; resolution < numRes; ++resolution) {
            images[tone][resolution] = new Image(imageSize, imageSize);
            static Pixel white = Pixel(1,1,1);
            for (int i=0; i < imageSize; i++) {
                for (int j=0; j < imageSize; j++) {
                    images[tone][resolution]->setPixel(i, j, white);
                }
            }
            if (tone == 0) {
                candidates[resolution] = Image(imageSize, imageSize);
            }
            imageSize *= 2;
        }
    }
    
    double darkestTone =.9;
    double lightestTone = .15;
    double toneInterval;
    if (numTones == 1) {
        toneInterval = darkestTone - lightestTone;
    } else {
        toneInterval = (darkestTone - lightestTone) / (numTones - 1);
    }
    for (int toneLevel = 0; toneLevel < numTones; ++toneLevel) {
        double maxTone = lightestTone + toneInterval * toneLevel;
        while (fabs(maxTone - images[toneLevel][numRes-1]->getTone()) > (.1/(numRes))) {
            RandomStroke bestStroke;
            double bestTone = -10000;
            for (int i = 0; i < 25; ++i) {
                double toneSum = 0;
                bool horizontalStroke = maxTone <= .7;
                RandomStroke currentStroke = getRandomStroke(horizontalStroke);
                for (int resolution = 0; resolution < numRes; ++resolution) {
                    if (fabs(maxTone - images[toneLevel][resolution]->getTone()) > (.1/(resolution+1.))) {
                        drawStroke(stroke, currentStroke, &candidates[resolution]);
                        // Get the effective length of the stroke, given that it might "run off" the edge
                        
                        double withStroke = fabs(maxTone - candidates[resolution].getTone());
                        double withoutStroke = fabs(maxTone - images[toneLevel][resolution]->getTone());
                        if (withoutStroke < withStroke) {
                            continue;
                        }
                        double toneContribution = candidates[resolution].getTone() - images[toneLevel][resolution]->getTone();
                        // All images are squares, side lengths are equal
                        const double imageSize = images[toneLevel][resolution]->getWidth();
                        
                        double actualStrokeLength = currentStroke.length * stroke->getWidth();
                        double normalizedStrokeLength;
                        double actualStrokeX = currentStroke.x*imageSize;
                        double actualStrokeY = currentStroke.y*imageSize;
                        if(horizontalStroke){
                            normalizedStrokeLength = getActualStrokeLength(actualStrokeX, actualStrokeLength, imageSize);
                        }
                        else {
                            normalizedStrokeLength = getActualStrokeLength(actualStrokeY, actualStrokeLength, imageSize);
                        }
                        toneContribution /= normalizedStrokeLength;
                        toneSum += toneContribution;
                        candidates[resolution] = *images[toneLevel][resolution];
                    }
                }
                
                if (toneSum > bestTone) {
                    bestTone = toneSum;
                    bestStroke = currentStroke;
                }
            }
            
            for (int candidate = 0; candidate < numRes; ++candidate) {
                if (fabs(maxTone - images[toneLevel][candidate]->getTone()) > (.1/(candidate + 1.))) {
                    drawStroke(stroke, bestStroke, images[toneLevel][candidate]);
                    candidates[candidate] = *images[toneLevel][candidate];
                }
            }
        }
        
        if (toneLevel != numTones - 1) {
            for (int resolution = 0; resolution < numRes; ++resolution) {
                *images[toneLevel+1][resolution] = *images[toneLevel][resolution];
                candidates[resolution] = *images[toneLevel][resolution];
            }
        }
    }
}



TAM::RandomStroke TAM::getRandomStroke(bool isHorizontal) const{
    static std::uniform_real_distribution<double> positive_dist(0, 1);
    static std::uniform_real_distribution<double> total_dist(-.2, 1);

    static std::uniform_real_distribution<double> length_dist(.3, 1);
    static std::uniform_int_distribution<int> rot_dist(-2, 2);
    
    RandomStroke stroke;
    
    stroke.length = length_dist(rng);
    if(isHorizontal){
        stroke.x= total_dist(rng);
        stroke.y = positive_dist(rng);
    } else {
        stroke.x = positive_dist(rng);
        stroke.y = total_dist(rng);
    }
    
    stroke.rot =deg2rad(rot_dist(rng));
    stroke.horizontal = isHorizontal;
    
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
        return Pixel(1, 1, 1);
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

Image* rotate_90(Image*src)
{
    Image* dest = new Image(src->getHeight(), src->getWidth());
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            Pixel p = src->getPixel(j, i);
            dest->setPixel(i, j, p);
        }
    }

    return dest;
}

/*
 * rotate image
 */
Image* ip_rotate (Image* src, double theta)
{
    Image* dest = new Image(src->getWidth(), src->getHeight());
    double x = src->getWidth()/2;
    double y = src->getHeight()/2;
    double ctheta = cos(-theta);
    double stheta = sin(-theta);
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            Pixel p = Pixel(1, 1, 1);
            dest->setPixel(i, j, p);
        }
    }
    
    for (int i = 0; i < dest->getWidth(); ++i) {
        for (int j = 0; j < dest->getHeight(); ++j) {
            double xp = x + (i-x)*ctheta - (j-y)*stheta;
            double yp = y + (i-x)*stheta + (j-y)*ctheta;
            
//            if (xp >= src->getWidth() || xp < 0 || yp >= src->getHeight() || yp < 0) {
//                continue;
//            }
            
            Pixel p = ip_resample_bilinear(src, xp, yp);
            
            dest->setPixel(i, j, p);
        }
    }
    
    return dest;
}

/*
 * scale image
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
void ip_composite(Image* dest, Image* strokeImage, int x, int y)
{
    static const double whiteThreshold = 1;
    for (int i = 1; i < strokeImage->getWidth()-1; ++i) {
        for (int j = 1; j < strokeImage->getHeight()-1; ++j) {
            if((x+i)<dest->getWidth() && (y+j)<dest->getHeight() && (x + i) >= 0 && (y+j)>=0) {
                Pixel strokePixel = strokeImage->getPixel(i, j);
                Pixel imagePixel = dest->getPixel(x+i, y+j);
                if(strokePixel.getColor(0)<imagePixel.getColor(0)){
                    dest->setPixel(x+i, y+j, strokePixel);
                }
//                const bool isWhite = strokePixel.getColor(RED) >= whiteThreshold;
//                const double destRed   = isWhite ? dest->getPixel(x+i, y+j, RED) : strokeImage->getPixel(i, j, RED);
//                const double destGreen = isWhite ? dest->getPixel(x+i, y+j, GREEN) : strokeImage->getPixel(i, j, GREEN);
//                const double destBlue  = isWhite ? dest->getPixel(x+i, y+j, BLUE) : strokeImage->getPixel(i, j, BLUE);
                //Pixel p = Pixel(destRed, destGreen, destBlue);
                //dest->setPixel(x + i, y + j, p);
            }
        }
    }
}

void TAM::drawStroke(Image* strokeImage, RandomStroke strokeModifiers, Image* dest){
    Image* scaledImage = ip_scale(strokeImage, strokeModifiers.length, 1);
    Image* rotImage = ip_rotate(scaledImage, strokeModifiers.rot);
    
    if (!strokeModifiers.horizontal) {
        Image* old = rotImage;
        rotImage = rotate_90(rotImage);
        delete old;
    }
    
    ip_composite(dest, rotImage, strokeModifiers.x*dest->getWidth(), strokeModifiers.y*dest->getHeight());
    delete scaledImage;
    delete rotImage;
}