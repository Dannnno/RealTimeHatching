#include "main.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "control.h"
#include "SOIL/SOIL.h"

int  window_width  = 300;
int  window_height = 300;

double zoomFactor = 1.;

Image* currentImage  = NULL;
Image* originalImage = NULL;

bool quietMode = false;
bool textMode  = false;



void checkGLError() {
    GLenum errCode = glGetError();
    if (errCode != GL_NO_ERROR)
    {
        const GLubyte* errString = gluErrorString(errCode);
        cerr << "OpenGL error: " << errString << endl;    } else {
        cerr << "No error :)" << endl;
    }
}


GLuint textures[NUMBER_OF_TONES];

void loadTAM(){
    
    int width, height;
    unsigned char* image;
    std::vector<std::vector<std::string>> tamNames(NUMBER_OF_TONES, std::vector<std::string>(NUMBER_OF_RESOLUTIONS));
    
    for (int tone = 0; tone < NUMBER_OF_TONES; ++tone) {
        for (int res = 0; res < NUMBER_OF_RESOLUTIONS; ++res) {
            char filename[100];
            sprintf(filename, "/Users/swarren/Downloads/ip2016skeleton/TAM_tone%d_resolution%d.bmp", tone, res);
            tamNames[tone][res] = std::string(filename);
        }
    }
    
    for (int tone = 0; tone < NUMBER_OF_TONES; ++tone) {
        glActiveTexture(GL_TEXTURE0 + tone);
        textures[tone] = tone;
        glBindTexture(GL_TEXTURE_2D, textures[tone]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, NUMBER_OF_RESOLUTIONS - 1);
        
        for (int resolution = 0; resolution < tamNames[tone].size(); ++resolution) {
            image = SOIL_load_image(tamNames[tone][resolution].c_str(), &width, &height, 0, SOIL_LOAD_RGB);
            glTexImage2D(GL_TEXTURE_2D, resolution, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
            SOIL_free_image_data(image);
        }
    }
    
}


int main (int argc, char** argv)
{
    // initialize parameters
    char* toLoad = init(argc, argv);
    
    if (textMode)
    {
        if (toLoad)
            image_load(toLoad);
        textMenuLoop();
    }
    else
    {
        // set up the window
        glutInit(&argc, &argv[0]);
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
        glutInitWindowPosition(100,100);
        glutInitWindowSize(window_width, window_height);
        glutCreateWindow("real-time hatching");
        
        // register call back functions
        glutDisplayFunc(display);
        glutReshapeFunc(unreshape);
        
        glClearColor(0.0,0.0,0.0,0.0);
        glDisable(GL_DEPTH_TEST);
        
        // setup main menu
        make_menu();
        
        // register keyboard callback function
        glutKeyboardFunc(keyboard_func);
        
        if (toLoad)
            image_load(toLoad);
        
        
        loadTAM();
        
        // wait for something to happen
        glutMainLoop();
    }
    return 0;
}


char* init (int argc, char** argv)
{
    // init random number generator
    //srand48(time(0));
    
    char* toLoad = NULL;
    
    // parse the command line options
    bool noMoreArgs  = false;
    bool noMoreFlags = false;
    if (argc > 1)
    {
        for (int i = 1; i < argc; i++)
        {
            if (noMoreArgs)
                usage();
            
            if (!noMoreFlags && argv[i][0] == '-')
            {
                switch (argv[i][1])
                {
                    case 't':
                        textMode = true;
                        break;
                        
                    case 'q':
                        quietMode = true;
                        break;
                        
                    case '-':
                        if (argv[i][2] == '\0')
                            noMoreFlags = true;
                        else
                            usage();
                        break;
                        
                    default:
                        usage();
                }
            }
            else
            {
                noMoreArgs = true;
                //        image_load(argv[i]);
                toLoad = argv[i];
            }
        }
    }
    
    return toLoad;
}


void usage ()
{
    cerr << "usage: ./ip [ -t ] [ -q ] [ -- ] [ file ]" << endl;
    exit(-1);
}

void draw_cube(int size, int x, int y, int z) {
    glPushMatrix();
    glBegin(GL_QUADS);
    glTexCoord2i(0, 0); glVertex3f(x, y, z);
    glTexCoord2i(0, 1); glVertex3f(x, y+size, z);
    glTexCoord2i(1, 1); glVertex3f(x+size, y+size, z);
    glTexCoord2i(1, 0); glVertex3f(x+size, y, z);
    glEnd();
    glPopMatrix();
}

void draw_test_cubes() {
//    glEnable(GL_TEXTURE_2D);
    glColor3d(1, 1, 1);
    for (int i = 0; i < NUMBER_OF_TONES; ++i) {
//        glBindTexture(GL_TEXTURE_2D, textures[i]);
        draw_cube(2, i, 0, 0);
    }
//    glDisable(GL_TEXTURE_2D);
}


void display ()
{
    if (textMode)
        return;
    
    // check if there have been any openGL problems
    GLenum errCode = glGetError();
    if (errCode != GL_NO_ERROR)
    {
        const GLubyte* errString = gluErrorString(errCode);
        cerr << "OpenGL error: " << errString << endl;
    }
    
    // clear the frame buffer
    glClear(GL_COLOR_BUFFER_BIT);
    
    // draw the image
    if (currentImage) {
        currentImage->glDrawPixelsWrapper();
    } else {
        cerr << "here" << endl;
        gluLookAt(10, 2*zoomFactor, 10, 0, 0, 0, 0, 1, 0);
        draw_test_cubes();
    }
    
    // swap buffers
    glutSwapBuffers();
}


void unreshape (int width, int height)
{
    // don't allow user to manuall resize the window
    reshape(window_width, window_height);
}


void reshape (int width, int height)
{
    // set window height and width
    window_width  = max(width,  64);
    window_height = max(height, 64); 
    
    if (textMode)
        return;
    
    // change the actual window's size
    glutReshapeWindow(window_width, window_height);
    
    // the lower left corner of the viewport is 0,0
    // the upper right corner is width, height
    glViewport(0, 0, (GLint) window_width, (GLint) window_height);  
    
    // setup orthographic projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, window_width, 0.0, window_height);
    
    // default mode should be modelview
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}
