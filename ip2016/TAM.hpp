//
//  TAM.hpp
//  ip2016
//
//  Created by swarren on 4/15/16.
//  Copyright © 2016 Z Sweedyk. All rights reserved.
//

#ifndef TAM_hpp
#define TAM_hpp

#include <stdio.h>
#include <random>
#include <vector>
#include "image.h"

void ip_composite(Image* dest, Image* strokeImage, int x, int y);
Image* ip_scale (Image* src, double xFac, double yFac);
class TAM{
struct RandomStroke;
public:
    TAM (int numTones, int numRES, Image* stroke);
    void drawStroke(Image* strokeImage, RandomStroke strokeModifiers, Image* dest);

    
    std::vector<std::vector<Image*>> images;
private:
    struct RandomStroke {
        double length;
        double x, y;
        double rot;
        bool horizontal;
    };
    
    RandomStroke getRandomStroke(bool isHorizontal) const;
    
    

};

#endif /* TAM_hpp */
