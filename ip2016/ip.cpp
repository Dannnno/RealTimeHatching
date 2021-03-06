#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "TAM.hpp"

Image* ip_generate_TAM ()
{

    const char* filepath = "/Users/swarren/Downloads/ip2016skeleton/SampleStrokeBMP.bmp";
    Image stroke = Image(filepath);
    Image* scaledstroke = ip_scale(&stroke, .7, .7);
    TAM* t = new TAM(NUMBER_OF_TONES, NUMBER_OF_RESOLUTIONS, scaledstroke);
    delete scaledstroke;
    
    for(int a=0; a<t->images[0].size();a++){
        for( int b =0; b<t->images.size(); b++) {
            char filename[100];
            sprintf(filename, "/Users/swarren/Downloads/ip2016skeleton/TAM_tone%d_resolution%d.bmp", b, a);
            t->images[b][a]->writeBMP(filename);
        }
    }
    
    int height = 0;
    for( int i=0; i<t->images[0].size(); i++){
        height+= t->images[0][i]->getHeight();
        height+=10;
    }
    
    
    int colWidth = t->images[0][t->images[0].size()-1]->getWidth();
    int width = (colWidth+10)*t->images.size();
    
    Pixel white = Pixel (1,1,1);
    Image* result = new Image(width+1, height+1);
    result->fillImage(white);
    int yStart = 1;
    for(int y=0; y<t->images[0].size(); y++){
        int xStart=1;
        for(int x=0; x<t->images.size(); x++){
            ip_composite(result, t->images[x][y], xStart, yStart);
            xStart += 10 + colWidth;
        }
        yStart+= 10 + t->images[0][y]->getHeight();
    }
    
    delete t;
    
    char filename[100];
    sprintf(filename, "/Users/swarren/Downloads/ip2016skeleton/TAM.bmp");
    result->writeBMP(filename);
    return result;
}






