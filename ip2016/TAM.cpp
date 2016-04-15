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
    
    
}

TAM::drawStroke(Image* strokeImage, RandomStroke strokeModifiers, Image* dest){
    
}