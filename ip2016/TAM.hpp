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
public:
    TAM (int numTones, int numRES, Image* stroke);
    
private:
    struct RandomStroke {
        double length;
        double x, y;
        int rot;
    };
    
    RandomStroke getRandomStroke() {
        static std::uniform_real_distribution<double> pos_dist(0, 1);
        static std::uniform_real_distribution<double> length_dist(.3, 1);
        static std::uniform_int_distribution<int> rot_dist(-2, 2);
        
        RandomStroke stroke;
        
        stroke.length = length_dist(rng);
        stroke.x = pos_dist(rng);
        stroke.y = pos_dist(rng);
        stroke.rot = rot_dist(rng);
        
        return stroke;
    }
    
    std::vector<std::vector<Image*>> images;

};

#endif /* TAM_hpp */
