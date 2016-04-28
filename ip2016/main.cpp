#include "common.h"
#include <iostream>
#include <vector>
#include <string>
#include <glut/glut.h>
#include <math.h>
#include "SOIL/SOIL.h"
#include "ip.h"

using std::cerr;
using std::endl;


// function prototypes
void display(void);
void reshape(int width, int height);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void init();


GLuint textures[NUMBER_OF_TONES];
GLfloat lightPosition[]={5, 5, 5,1};

void checkGLError() {
    GLenum errCode = glGetError();
    if (errCode != GL_NO_ERROR)
    {
        const GLubyte* errString = gluErrorString(errCode);
        cerr << "OpenGL error: " << errString << endl;    } else {
            cerr << "No error :)" << endl;
        }
}

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
    
    glGenTextures(NUMBER_OF_TONES, textures);
    for (int tone = 0; tone < NUMBER_OF_TONES; ++tone) {
        glBindTexture(GL_TEXTURE_2D, textures[tone]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 8);
        glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

        for (int resolution = 0; resolution < tamNames[tone].size(); ++resolution) {
            image = SOIL_load_image(tamNames[tone][tamNames[tone].size()-1-resolution].c_str(), &width, &height, 0, SOIL_LOAD_RGB);
            glTexImage2D(GL_TEXTURE_2D, resolution, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
            SOIL_free_image_data(image);
        }
    }
    
}


int windowWidth=600;
int windowHeight=600;
double phi = 0;

using namespace std;

int main(int argc, char **argv)
{
    
    // initialize glut
    glutInit(&argc, argv);
    
    // set window size
    glutInitWindowSize(windowWidth,windowHeight);
    
    // establish glut display parameters
    glutInitDisplayMode(GLUT_DOUBLE   | GLUT_RGB  |GLUT_DEPTH);
    
    // create window
    glutCreateWindow("My Second OpenGL program");
    
    // register callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    
    
    // initalize opengl parameters
    init();
   ip_generate_TAM();
    loadTAM();
    
    // loop until something happens
    glutMainLoop();
    return 0;
}

void init()
{
    // initialize viewing system
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(20.0, 1.0, 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    // enable light0 and lighting
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    //glEnable(GL_COLOR_MATERIAL);
    
    // set color of light0
    GLfloat white[] = {1,1,1,0};			// light color
    glLightfv(GL_LIGHT0, GL_DIFFUSE, white);	// set diffuse light color
    glLightfv(GL_LIGHT0, GL_SPECULAR, white);	// set specular light color
    
    // initialize background color to black
    glClearColor(0,0,0,0);
    
    // enable depth buffering
    glEnable(GL_DEPTH_TEST);
}

void mouse(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
        cout << "Left click with cursor at" << x << " " << y << endl;
}

void motion(int x, int y)
{
    //cout << "Mouse at" << x << " " << y << endl;
    
    
    static int oldX = 0;
    static int oldY = 0;
    
    int deltaX = x - oldX;
    int deltaY = y - oldY;
    
    if (deltaX > 20 || deltaX < -20){
        oldX = x;
        oldY = y;
        return;
    }
    
    phi += (deltaX*1.0/windowWidth)*180;
    
    oldX = x;
    oldY = y;
    //cout << phi << endl;
    glutPostRedisplay();
}

void reshape(int width, int height)
{
    if (width<height)
        glViewport(0,0,width,width);
    else
        glViewport(0,0,height,height);
    
}

//void drawSphere(int slices, int stacks, int radius)

//void drawQuad(float v1 [3], float v2 [3], float v3 [3], float v4 [3]){
//    
//}

void display()
{
    
    GLfloat white[] = {1,1,1,0};	    // white
    GLfloat red[] = {1,0,0,0};              // red
    GLfloat green[] = {0,1,0,0};             // green
    GLfloat purple[] = {1,0,1,0};	    // purple
    
    // clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // initialize modelview matrix
    glMatrixMode(GL_MODELVIEW_MATRIX);
    glLoadIdentity();
    
    // set viewpoint position/orientation
    gluLookAt(20*sin(phi*3.14/180.0),10,20*cos(phi*3.14/180.0),0,0,0,0,1,0);
    // position of light0
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    
    glNormal3f(0,0,1);
    glNormal3f(1,0,0);
    
    
    glEnable(GL_TEXTURE_2D);
    glColor3d(1, 1, 1);
    for (int i = -3; i < 3; i += 1) {
        glActiveTexture(textures[i+3]); //GL_TEXTURE0 + i + 3);
        glBindTexture(GL_TEXTURE_2D, textures[i+3]);
        int x = i, y = 0, z = 0, size = 1;
        glPushMatrix();
        glBegin(GL_QUADS);
        glTexCoord2i(0, 0); glVertex3f(x, y, z);
        glTexCoord2i(1, 0); glVertex3f(x, y+size, z);
        glTexCoord2i(1, 1); glVertex3f(x+size, y+size, z);
        glTexCoord2i(0, 1); glVertex3f(x+size, y, z);
        glEnd();
        glPopMatrix();
    }
        glDisable(GL_TEXTURE_2D);

//    // draw a red triangle
//    glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
//    glMaterialfv(GL_FRONT, GL_SPECULAR, red);
//    glMateriali(GL_FRONT,GL_SHININESS,0);
//    glBegin(GL_TRIANGLES);
//    glVertex3f(-3,-1,-8);
//    glVertex3f(3,-1,-10);
//    glVertex3f(0,3,-9);
//    glEnd();
//    
//    glNormal3f(0,1,0);
//    
//    
//    // draw Triangle strip
//    glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
//    glMaterialfv(GL_FRONT, GL_SPECULAR, green);
//    glMateriali(GL_FRONT,GL_SHININESS,0);			// gree
//    //glColor4fv(green);
//    glBegin(GL_TRIANGLE_STRIP);
//    glVertex3f(-10.0,-1.0,-10.0);
//    glVertex3f(-10.0,-1.0,10.0);
//    glVertex3f(10.0,-1.0,-10.0);
//    glVertex3f(10.0,-1.0,10.0);
//    glEnd();
//    
//    // draw Sphere
//    glMaterialfv(GL_FRONT, GL_DIFFUSE, purple);
//    glMaterialfv(GL_FRONT, GL_SPECULAR, purple);
//    glMateriali(GL_FRONT,GL_SHININESS,50);			// purple
//    //glColor4fv(purple);
//    glutSolidSphere(.75,100,100);
    
    // swap buffers
    glutSwapBuffers();
    
    
    
}