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

#define X .525731112119133606
#define Z .850650808352039932

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



static GLfloat vdata[12][3] = {
    {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},
    {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
    {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}
};

static GLuint tindices[20][3] = {
    {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
    {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},
    {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
    {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11} };
void normalize(GLfloat *a) {
    GLfloat d=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    a[0]/=d; a[1]/=d; a[2]/=d;
}

double dot (GLfloat left[3], GLfloat right[3]) {
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

int getTextureIndex(GLfloat colorValue){
    //colorValue = 1 - colorValue;
    double darkestTone =.9;
    double lightestTone = .15;
    double toneInterval;
    int closestIndex = -1;
    double smallestDiff= 100;
    
    if (NUMBER_OF_TONES == 1) {
        toneInterval = darkestTone - lightestTone;
    } else {
        toneInterval = (darkestTone - lightestTone) / (NUMBER_OF_TONES - 1);
    }
    for (int toneLevel = 0; toneLevel < NUMBER_OF_TONES; ++toneLevel) {
        double maxTone = lightestTone + toneInterval * toneLevel;
        double newDiff = fabs(maxTone - colorValue);
        if (newDiff<smallestDiff){
            smallestDiff= newDiff;
            closestIndex=toneLevel;
        }
    }
    return closestIndex;
}

GLfloat getLightValue(GLfloat* normal, GLfloat* vertexPosition, GLfloat* lightDir, GLfloat* cameraPos) {
    
    GLfloat normalizedLight[3] = {lightDir[0], lightDir[1], lightDir[2]};
    GLfloat normalizedNormal[3] = {normal[0], normal[1], normal[2]};
    normalize(normalizedLight);
    normalize(normalizedNormal);
    double angleFactor = -1*dot(normalizedLight, normalizedNormal);
    if (angleFactor<0)
        return 0;  // light is falling on other side of surface
    else
        return angleFactor;

    //implement specular later (maybe)
//    float lightDirDotNormal = 2*angleFactor;
//    
//    GLfloat reflect[3]={lightDir[0]+normal[0]*lightDirDotNormal, lightDir[1]+normal[1]*lightDirDotNormal,
//        lightDir[2]+normal[2]*lightDirDotNormal};
//        //direction - info.normal*(2*direction.dot(info.normal));
//    float reflectionLength = sqrt(pow(reflect[0],2)+pow(reflect[1],2)+pow(reflect[2],2));
//    reflect[0]/=reflectionLength;
//    reflect[1]/=reflectionLength;
//    reflect[2]/=reflectionLength;
//    double angleFactor = - reflect.dot(info.theRay.getDir());
//    if (angleFactor<0)
//    return Color3d(0,0,0);
//    else {
//    angleFactor = pow(angleFactor,info.material->getKshine());
//    return   color * info.material->getSpecular() * angleFactor;
}

void drawCylinder(int slices, int height, int radius){
    double angleFraction = (M_PI*2)/slices;
    for(int i=0; i<=slices; i++){
        
        // Draw the top and bottom circles
        glBegin(GL_TRIANGLES);
        glNormal3f(0, 1, 0);
        glVertex3f(0, height, 0);
        glVertex3f(radius*cos(i*angleFraction), height, radius*sin(i*angleFraction));
        glVertex3f(radius*cos((i+1)*angleFraction), height, radius*sin((i+1)*angleFraction));
        
        glNormal3f(0, -1, 0);
        glVertex3f(0, 0, 0);
        glVertex3f(radius*cos(i*angleFraction), 0, radius*sin(i*angleFraction));
        glVertex3f(radius*cos((i+1)*angleFraction), 0, radius*sin((i+1)*angleFraction));
        glEnd();
        
        // Draw the side
        GLfloat v0[3] = {
            static_cast<GLfloat>(radius*cos(i*angleFraction)),
            static_cast<GLfloat>(height),
            static_cast<GLfloat>(radius*sin(i*angleFraction))
        };
        GLfloat normal[3] = {
            static_cast<GLfloat>(cos((i+.5)*angleFraction)),
            static_cast<GLfloat>(0),
            static_cast<GLfloat>(sin((i+.5)*angleFraction))
        };
        
        GLfloat lightValue = getLightValue(v0, normal, lightPosition, nullptr);
        int texIndex = getTextureIndex(lightValue);
        glEnable(GL_TEXTURE_2D);
        glActiveTexture(textures[texIndex]);
        cout << lightValue << ": " << texIndex << endl; 
        glBindTexture(GL_TEXTURE_2D, textures[texIndex]);
        glBegin(GL_QUADS);
        glNormal3f(normal[0], normal[1], normal[2]);
        glTexCoord2i(0, 1); glVertex3f(v0[0], v0[1], v0[2]);
        glTexCoord2i(1, 1); glVertex3f(radius*cos((i+1)*angleFraction), height, radius*sin((i+1)*angleFraction));
        glTexCoord2i(1, 0); glVertex3f(radius*cos((i+1)*angleFraction), 0, radius*sin((i+1)*angleFraction));
        glTexCoord2i(0, 0); glVertex3f(radius*cos(i*angleFraction), 0, radius*sin(i*angleFraction));
        glDisable(GL_TEXTURE_2D);
        glEnd();
        
    }
}

GLfloat** sortVertices(GLfloat *a, GLfloat* b, GLfloat* c, int axis) {
    GLfloat* vertices[3] = {a, b, c};
    std::sort(vertices, vertices + 3, [=](GLfloat* l, GLfloat* r) { return l[axis] < r[axis]; });
    return vertices;
}


void drawtri(GLfloat *a, GLfloat *b, GLfloat *c, int div, float r) {
    if (div<=0) {
        glBlendFunc(GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
        glEnable(GL_BLEND);
        glEnable(GL_COLOR_MATERIAL);
        
        glColor4f(1.0f, 1.0f, 1.0f, (.4f));
        GLfloat lightVal[3] = {getLightValue(a, a, lightPosition, NULL), getLightValue(b, b, lightPosition, NULL),getLightValue(c, c, lightPosition, NULL)};
        int texIndex[3] = {getTextureIndex(lightVal[0]), getTextureIndex(lightVal[1]), getTextureIndex(lightVal[2])};
        for (int i=0; i<3; i++) {
            glActiveTexture(textures[texIndex[i]]);
            glBindTexture(GL_TEXTURE_2D, textures[texIndex[i]]);
            glBegin(GL_TRIANGLES);
            // Sort the textures by y coordinate to determine orientation
            GLfloat** sortedByYCoordinates = sortVertices(a, b, c, 1);
            GLfloat* v0 = sortedByYCoordinates[0];
            GLfloat* v1 = sortedByYCoordinates[1];
            GLfloat* v2 = sortedByYCoordinates[2];
            if (fabs(v0[1] - v1[1]) < fabs(v1[1] - v2[1])) {
                if (v0[0] < v1[0]) {
                    // We look like this:
                    //      2
                    //   0     1
                    
                    glNormal3fv(v0); glTexCoord2i(0, 0); glVertex3f(v0[0]*r, v0[1]*r, v0[2]*r);
                    glNormal3fv(v1); glTexCoord2i(1, 0); glVertex3f(v1[0]*r, v1[1]*r, v1[2]*r);
                    glNormal3fv(v2); glTexCoord2i(.5, 1); glVertex3f(v2[0]*r, v2[1]*r, v2[2]*r);
                } else {
                    // We look like this:
                    //      2
                    //   1     0
                    
                    glNormal3fv(v1); glTexCoord2i(0, 0); glVertex3f(v1[0]*r, v1[1]*r, v1[2]*r);
                    glNormal3fv(v0); glTexCoord2i(1, 0); glVertex3f(v0[0]*r, v0[1]*r, v0[2]*r);
                    glNormal3fv(v2); glTexCoord2i(.5, 1); glVertex3f(v2[0]*r, v2[1]*r, v2[2]*r);
                }
            } else {
                if (v2[0] < v1[0]) {
                    // We look like this:
                    //     2     1
                    //        0
                    glNormal3fv(v0); glTexCoord2i(.5, 0); glVertex3f(v0[0]*r, v0[1]*r, v0[2]*r);
                    glNormal3fv(v1); glTexCoord2i(1, 1); glVertex3f(v1[0]*r, v1[1]*r, v1[2]*r);
                    glNormal3fv(v2); glTexCoord2i(0, 1); glVertex3f(v2[0]*r, v2[1]*r, v2[2]*r);
                } else {
                    // We look like this:
                    //     1     2
                    //        0
                    glNormal3fv(v0); glTexCoord2i(.5, 0); glVertex3f(v0[0]*r, v0[1]*r, v0[2]*r);
                    glNormal3fv(v2); glTexCoord2i(1, 1); glVertex3f(v2[0]*r, v2[1]*r, v2[2]*r);
                    glNormal3fv(v1); glTexCoord2i(0, 1); glVertex3f(v1[0]*r, v1[1]*r, v1[2]*r);
                }
            }
            glEnd();
        }
        glDisable(GL_BLEND);
        glDisable(GL_COLOR_MATERIAL);
    } else {
        GLfloat ab[3], ac[3], bc[3];
        for (int i=0;i<3;i++) {
            ab[i]=(a[i]+b[i])/2;
            ac[i]=(a[i]+c[i])/2;
            bc[i]=(b[i]+c[i])/2;
        }
        normalize(ab); normalize(ac); normalize(bc);
        drawtri(a, ab, ac, div-1, r);
        drawtri(b, bc, ab, div-1, r);
        drawtri(c, ac, bc, div-1, r);
        drawtri(ab, bc, ac, div-1, r);  //<--Comment this line and sphere looks really cool!
    }
}
void drawsphere(int ndiv, float radius=1.0) {
    glEnable(GL_TEXTURE_2D);
    
//    glBegin(GL_TRIANGLES);
    for (int i=0;i<20;i++)
        drawtri(vdata[tindices[i][0]], vdata[tindices[i][1]], vdata[tindices[i][2]], ndiv, radius);
    glDisable(GL_TEXTURE_2D);
//    glEnd();
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
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
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
   //ip_generate_TAM();
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
    gluLookAt(10*sin(phi*3.14/180.0),5,10*cos(phi*3.14/180.0),0,0,0,0,1,0);
    // position of light0
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    
    glNormal3f(0,0,1);
    glNormal3f(1,0,0);
    glColor3d(1, 1, 1);

//    drawsphere(1);
    
    drawCylinder(20, 1, 1);
    
    
//    glEnable(GL_TEXTURE_2D);
//    for (int i = -3; i < 3; i += 1) {
//        glActiveTexture(textures[i+3]); //GL_TEXTURE0 + i + 3);
//        glBindTexture(GL_TEXTURE_2D, textures[i+3]);
//        int x = i, y = 0, z = 0, size = 1;
//        glPushMatrix();
//        glBegin(GL_QUADS);
//        glTexCoord2i(0, 0); glVertex3f(x, y, z);
//        glTexCoord2i(1, 0); glVertex3f(x, y+size, z);
//        glTexCoord2i(1, 1); glVertex3f(x+size, y+size, z);
//        glTexCoord2i(0, 1); glVertex3f(x+size, y, z);
//        glEnd();
//        glPopMatrix();
//    }
//        glDisable(GL_TEXTURE_2D);
    

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