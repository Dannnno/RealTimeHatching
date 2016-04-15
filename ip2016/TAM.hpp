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

std::mt19937 rng;

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
    };
    
    RandomStroke getRandomStroke() const;
    
    

};

#endif /* TAM_hpp */
