#include "control.h"
#include "ip.h"
#include "main.h"
#include <stdlib.h>

enum {
    M_QUIT = 0,
    M_HELP = 1,
    
    M_FILE_OPEN =2,
    M_FILE_SAVE =3,
    M_GEN_TAM =4,
    M_SHOW_TAM =5,
    M_LOAD_CUBE =6,
    M_LAST_ENUM
} MENU_ITEMS;


int make_menu ()
{
    int main = glutCreateMenu(menu_func);
    glutAddMenuEntry( "Save...",		M_FILE_SAVE);
    glutAddMenuEntry( "Generate TAM",		M_GEN_TAM);
    glutAddMenuEntry( "Show TAM",		M_SHOW_TAM);
    glutAddMenuEntry( "Load Cube",		M_LOAD_CUBE);
    glutAddMenuEntry( "Help",		M_HELP);
    glutAddMenuEntry( "Quit",		M_QUIT);
    
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    
    return main;
}


static inline void checkStream (const istream& in)
{
    if (in.fail())
    {
        cerr << "Fatal error: stream failed!" << endl;
        exit(-1);
    }
}


void menu_func (int value)
{
    // variables used in the switch statement
    char filename[MAX_LINE];
    
//    if (currentImage) delete currentImage;
    
    switch (value)
    {
        case M_QUIT:  // enum #0
            exit(0);
            break;
        case M_HELP:  // enum #1
            menu_help();
            break;
//        case M_FILE_SAVE:   // enum #3
//            cerr << "Save as (string - no spaces) : ";
//            cin  >> filename;
//            checkStream(cin);
//            image_save(filename);
//            break;
////
//        case M_GEN_TAM:
//            currentImage = ip_generate_TAM();
//            
//            if (currentImage->getWidth()  != window_width    ||
//                currentImage->getHeight() != window_height)
//            {
//                reshape(currentImage->getWidth(), currentImage->getHeight());
//            }
//            
//            if (!quietMode)
//            {
//                cerr << "done!" << endl;
//            }
//            
//            if (!textMode)
//            {
//                glutPostRedisplay();
//            }
//            
//            break;
//            
        default:
            cerr << "Tried to use value " << value << " but not yet implemented" << endl;
    }
    return;
}

void keyboard_func (unsigned char key, int x, int y)
{
    switch (key)
    {
        case 'H':
        case 'h':
            menu_help();
            break;
        case 'Q':
        case 'q':
            exit(0);
            break;
//        case 'z':
//            ++zoomFactor;
//            break;
//        case 'Z':
//            --zoomFactor;
//            break;
    }
    
    glutPostRedisplay();
}


void menu_help ()
{
    cerr << endl
    << "hmc cs155 image processor" << endl
    << "please see the ip manual for usage and algorithm information" << endl
    << "http://www.cs.hmc.edu/courses/2002/fall/cs155/proj1/doc/ip_manual.html"
    << endl << endl;
}


#define MENUOP(num, tag)	cerr << " " << num << ") " << tag << endl;



void image_load (const char* filename)
{
//    if (currentImage)
//        delete currentImage;
//    if (originalImage)
//        delete originalImage;
//    currentImage  = NULL;
//    originalImage = NULL;
//    
//    originalImage = new Image();
//    originalImage->readBMP(filename);
//    
//    if (originalImage->good())
//    {
//        currentImage = new Image(*originalImage);
//        reshape(currentImage->getWidth(), currentImage->getHeight());
//    }
//    else
//    {
//        delete originalImage;
//        originalImage = NULL;
//        cerr << "Couldn't load image " << filename << "!" << endl;
//        return;
//    }
//    
//    if (!textMode)
//        glutPostRedisplay();
//    
//    if (!quietMode)
//        cerr << "done!" << endl;
}


void image_save (const char* filename)
{
//    if (currentImage)
//    {
//        if (currentImage->writeBMP(filename) == 0)
//        {
//            //delete originalImage;
//            //originalImage = new Image(*currentImage);
//        }
//    }
//    else if (originalImage)
//    {
//        originalImage->writeBMP(filename);
//    }
//    else
//    {
//        cerr << "No image!" << endl;
//        return;
//    }
//    
//    if (!quietMode)
//        cerr << "done!" << endl;
}
