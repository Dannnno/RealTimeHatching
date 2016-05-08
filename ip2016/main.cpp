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
        cerr << "OpenGL error: " << errString << endl;
    } else {
        cerr << "No error :)" << endl;
    }
}

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

void setSphereTexCoord(GLfloat coords[3], GLfloat radius) {
    GLfloat s = acosf(coords[2] / radius) / M_PI;
    GLfloat t = (atan2f(coords[1], coords[0]) - M_PI_2) / (2 * M_PI);
    if( s<.01 && s>-.01){
        cout << s << endl;
    }
    glTexCoord2d(t, s);
}

GLfloat getLightValue(GLfloat vertexPosition[3], GLfloat normal[3], GLfloat lightDir[3], GLfloat* cameraPos) {
    
//    if (lightDir[4] == 1) {
//        // Point light
//        return getPointLightValue(vertexPosition...)
//    }
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

void drawCone(double height, float baseRadius, int numStacks, int numSlices) {
    float angleFraction = (M_PI*2)/numSlices;
    double stackHeight = height/numStacks;
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    
    // draw the tip
    float curStackY= height; //((2*radius*i)/numStacks)-radius;
    float nextStackY = stackHeight*(numSlices-2);
    float nextStackRadius= baseRadius - (((numSlices-1)*baseRadius)/(numStacks-1.));
        for(int j=0; j<numSlices; j++) {
            GLfloat v0 [3] = {static_cast<GLfloat>(0, curStackY, 0)};
            GLfloat v1 [3] =  {static_cast<GLfloat>(nextStackRadius*cos(j*angleFraction)), nextStackY, static_cast<GLfloat>(nextStackRadius*sin(j*angleFraction))};
            GLfloat v2 [3] = {static_cast<GLfloat>(nextStackRadius*cos((j+1)*angleFraction)), nextStackY, static_cast<GLfloat>(nextStackRadius*sin((j+1)*angleFraction))};
            GLfloat m[3] = {
                static_cast<GLfloat>(sqrt(pow(v0[0], 2) + pow(v0[2], 2))),
                static_cast<GLfloat>(sqrt(pow(v1[0], 2) + pow(v1[2], 2))),
                static_cast<GLfloat>(sqrt(pow(v2[0], 2) + pow(v2[2], 2)))
            };
            GLfloat n0 [3] = {0, 1, 0};            GLfloat n1 [3] = {
                static_cast<GLfloat>((v1[0] / m[0]) * height / baseRadius),
                static_cast<GLfloat>(baseRadius / height),
                static_cast<GLfloat>((v1[2] / m[0]) * height / baseRadius),
            };
            GLfloat n2 [3] = {
                static_cast<GLfloat>((v2[0] / m[0]) * height / baseRadius),
                static_cast<GLfloat>(baseRadius / height),
                static_cast<GLfloat>((v2[2] / m[0]) * height / baseRadius),
            };
            GLfloat lightVal[3] = {
                getLightValue(v0, n0, lightPosition, NULL),
                getLightValue(v1, n1, lightPosition, NULL),
                getLightValue(v2, n2, lightPosition, NULL)
            };
            int texIndex[3] = {
                getTextureIndex(lightVal[0]),
                getTextureIndex(lightVal[1]),
                getTextureIndex(lightVal[2])
            };
            //
            //            glBlendFunc(GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            //            glEnable(GL_BLEND);
            //            glEnable(GL_COLOR_MATERIAL);
            
            glColor4f(1.0f, 1.0f, 1.0f, .25f);
            
            for (int k = 0; k < 3; ++k) {
                glEnable(GL_TEXTURE_2D);
                glActiveTexture(textures[texIndex[k]]);
                glBindTexture(GL_TEXTURE_2D, textures[texIndex[k]]);
                glBegin(GL_TRIANGLES);
                glNormal3fv(n0); glTexCoord2d((j+1)*angleFraction/(2*M_PI), curStackY/height);glVertex3fv(v0);
                glNormal3fv(n1); glTexCoord2d(j*angleFraction/(2*M_PI), curStackY/height); glVertex3fv(v1);
                glNormal3fv(n2); glTexCoord2d(j*angleFraction/(2*M_PI), nextStackY/height); glVertex3fv(v2);
                glEnd();
                glDisable(GL_TEXTURE_2D);
            }
            
            //            glDisable(GL_BLEND);
            //            glDisable(GL_COLOR_MATERIAL);
            
        }


    for(int i=1; i<numStacks-1; i++){
        float curStackY= stackHeight*i; //((2*radius*i)/numStacks)-radius;
        float curStackRadius = baseRadius - ((i*baseRadius)/(numStacks-1.));
        float nextStackY = stackHeight*(i+1);
        float nextStackRadius= baseRadius - (((i+1)*baseRadius)/(numStacks-1.));
        for(int j=0; j<numSlices; j++) {
            GLfloat v0 [3] = {static_cast<GLfloat>(curStackRadius*cos((j+1)*angleFraction)), curStackY, static_cast<GLfloat>(curStackRadius*sin((j+1)*angleFraction))};
            GLfloat v1 [3] = {static_cast<GLfloat>(curStackRadius*cos(j*angleFraction)), curStackY, static_cast<GLfloat>(curStackRadius*sin(j*angleFraction))};
            GLfloat v2 [3] =  {static_cast<GLfloat>(nextStackRadius*cos(j*angleFraction)), nextStackY, static_cast<GLfloat>(nextStackRadius*sin(j*angleFraction))};
            GLfloat v3 [3] = {static_cast<GLfloat>(nextStackRadius*cos((j+1)*angleFraction)), nextStackY, static_cast<GLfloat>(nextStackRadius*sin((j+1)*angleFraction))};
            GLfloat m[4] = {
                static_cast<GLfloat>(sqrt(pow(v0[0], 2) + pow(v0[2], 2))),
                static_cast<GLfloat>(sqrt(pow(v1[0], 2) + pow(v1[2], 2))),
                static_cast<GLfloat>(sqrt(pow(v2[0], 2) + pow(v2[2], 2))),
                static_cast<GLfloat>(sqrt(pow(v3[0], 2) + pow(v3[2], 2)))
            };
            GLfloat n0 [3] = {
                static_cast<GLfloat>((v0[0] / m[0]) * height / baseRadius),
                static_cast<GLfloat>(baseRadius / height),
                static_cast<GLfloat>((v0[2] / m[0]) * height / baseRadius),
            };
            GLfloat n1 [3] = {
                static_cast<GLfloat>((v1[0] / m[0]) * height / baseRadius),
                static_cast<GLfloat>(baseRadius / height),
                static_cast<GLfloat>((v1[2] / m[0]) * height / baseRadius),
            };
            GLfloat n2 [3] = {
                static_cast<GLfloat>((v2[0] / m[0]) * height / baseRadius),
                static_cast<GLfloat>(baseRadius / height),
                static_cast<GLfloat>((v2[2] / m[0]) * height / baseRadius),
            };
            GLfloat n3 [3] = {
                static_cast<GLfloat>((v3[0] / m[0]) * height / baseRadius),
                static_cast<GLfloat>(baseRadius / height),
                static_cast<GLfloat>((v3[2] / m[0]) * height / baseRadius),
            };
            GLfloat lightVal[4] = {
                getLightValue(v0, n0, lightPosition, NULL),
                getLightValue(v1, n1, lightPosition, NULL),
                getLightValue(v2, n2, lightPosition, NULL),
                getLightValue(v3, n3, lightPosition, NULL)
            };
            int texIndex[4] = {
                getTextureIndex(lightVal[0]),
                getTextureIndex(lightVal[1]),
                getTextureIndex(lightVal[2]),
                getTextureIndex(lightVal[3])
            };
//            
//            glBlendFunc(GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
//            glEnable(GL_BLEND);
//            glEnable(GL_COLOR_MATERIAL);
            
            glColor4f(1.0f, 1.0f, 1.0f, .25f);
            
            for (int k = 0; k < 4; ++k) {
                glEnable(GL_TEXTURE_2D);
                glActiveTexture(textures[texIndex[k]]);
                glBindTexture(GL_TEXTURE_2D, textures[texIndex[k]]);
                glBegin(GL_QUADS);
                glNormal3fv(n0); glTexCoord2d((j+1)*angleFraction/(2*M_PI), curStackY/height);glVertex3fv(v0);
                glNormal3fv(n1); glTexCoord2d(j*angleFraction/(2*M_PI), curStackY/height); glVertex3fv(v1);
                glNormal3fv(n2); glTexCoord2d(j*angleFraction/(2*M_PI), nextStackY/height); glVertex3fv(v2);
                glNormal3fv(n3); glTexCoord2d((j+1)*angleFraction/(2*M_PI), nextStackY/height); glVertex3fv(v3);
                glEnd();
                glDisable(GL_TEXTURE_2D);
            }
            
//            glDisable(GL_BLEND);
//            glDisable(GL_COLOR_MATERIAL);
            
        }
    }
    
    //draw base
    for( int a=0; a<numSlices; a++) {
        glBegin(GL_TRIANGLES);
        glNormal3f(0, -1, 0);
        glVertex3f(0, 0, 0);
        glVertex3f(baseRadius*cos(a*angleFraction), 0, baseRadius*sin(a*angleFraction));
        glVertex3f(baseRadius*cos((a+1)*angleFraction), 0, baseRadius*sin((a+1)*angleFraction));
        glEnd();
    }
    glDisable(GL_CULL_FACE);
}

void drawCubeWithQuads(double sideLength, int numSlices) {
    const double halfLength = sideLength / 2;
    const double quadSideLength = sideLength / numSlices;
    
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    // draw the top and bottom
    for (double x = -halfLength; x < halfLength; x += quadSideLength){
        for (double z = -halfLength; z < halfLength; z += quadSideLength){
            glBegin(GL_QUADS);
            // top
//            GLfloat v0[3] = {static_cast<GLfloat>(x), static_cast<GLfloat>(halfLength), static_cast<GLfloat>(z)};
//            GLfloat v1[3] = {x, halfLength, z + quadSideLength};
//            GLfloat v2[3] = {x+quadSideLength, halfLength, z};
//            GLfloat v3[3] = {x+quadSideLength, halfLength, z + quadSideLength};
//            GLfloat normal[3] = {0, 1, 0};
//            
//            GLfloat lightVal[4] = {
//                getLightValue(v0, normal, lightPosition, NULL),
//                getLightValue(v1, normal, lightPosition, NULL),
//                getLightValue(v2, normal, lightPosition, NULL),
//                getLightValue(v3, normal, lightPosition, NULL)
//            };
//            int texIndex[4] = {
//                getTextureIndex(lightVal[0]),
//                getTextureIndex(lightVal[1]),
//                getTextureIndex(lightVal[2]),
//                getTextureIndex(lightVal[3])
//            };
            glNormal3f(0, 1, 0);
            glTexCoord2d(x/sideLength, z/sideLength);glVertex3f(x, halfLength, z);
            glTexCoord2d(x/sideLength, (z+quadSideLength)/sideLength);glVertex3f(x, halfLength, z + quadSideLength);
            glTexCoord2d((x+quadSideLength)/sideLength, (z+quadSideLength)/sideLength);glVertex3f(x + quadSideLength, halfLength, z + quadSideLength);
            glTexCoord2d((x+quadSideLength)/sideLength, z/sideLength);glVertex3f(x + quadSideLength, halfLength, z);
            
            // bottom
            glNormal3f(0, -1, 0);
            glTexCoord2d(x/sideLength, z/sideLength);glVertex3f(x, -halfLength, z);
            glTexCoord2d(x/sideLength, (z+quadSideLength)/sideLength);glVertex3f(x, -halfLength, z + quadSideLength);
            glTexCoord2d((x+quadSideLength)/sideLength, (z+quadSideLength)/sideLength);glVertex3f(x + quadSideLength, -halfLength, z + quadSideLength);
            glTexCoord2d((x+quadSideLength)/sideLength, z/sideLength);glVertex3f(x + quadSideLength, -halfLength, z);
            glEnd();
        }
    }
    
    // draw the +/- x sides
    for (double y = -halfLength; y < halfLength; y += quadSideLength){
        for (double z = -halfLength; z < halfLength; z += quadSideLength){
            glBegin(GL_QUADS);
            // +x
            glNormal3f(1, 0, 0);
            glVertex3f(halfLength, y, z + quadSideLength);
            glVertex3f(halfLength, y, z);
            glVertex3f(halfLength, y + quadSideLength, z);
            glVertex3f(halfLength, y + quadSideLength, z + quadSideLength);
            
            // -x
            glNormal3f(-1, 0, 0);
            glVertex3f(-halfLength, y, z);
            glVertex3f(-halfLength, y, z + quadSideLength);
            glVertex3f(-halfLength, y + quadSideLength, z + quadSideLength);
            glVertex3f(-halfLength, y + quadSideLength, z);
            glEnd();
        }
    }
    
    // draw the +/- z sides
    for (double y = -halfLength; y < halfLength; y += quadSideLength){
        for (double x = -halfLength; x < halfLength; x += quadSideLength){
            glBegin(GL_QUADS);
            // +z
            glNormal3f(0, 0, 1);
            glVertex3f(x, y, halfLength);
            glVertex3f(x + quadSideLength, y, halfLength);
            glVertex3f(x + quadSideLength, y + quadSideLength, halfLength);
            glVertex3f(x, y + quadSideLength, halfLength);
            
            // -z
            glNormal3f(0, 0, -1);
            glVertex3f(x + quadSideLength, y, -halfLength);
            glVertex3f(x, y, -halfLength);
            glVertex3f(x, y + quadSideLength, -halfLength);
            glVertex3f(x + quadSideLength, y + quadSideLength, -halfLength);
            glEnd();
        }
    }
    
    glDisable(GL_CULL_FACE);
}



void drawSphereWithQuads(float radius, int numSlices, int numStacks){
    float angleFraction = (M_PI*2)/numSlices;
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    //draw sphere with center at origin from bottom up

    for(int i=0; i<numStacks; i++){
        float curStackY= ((2*radius*i)/numStacks)-radius;
        float curStackRadius = sqrt(pow(radius, 2)-pow(curStackY, 2));
        float nextStackY = ((2*radius*(i+1))/numStacks)-radius;
        float nextStackRadius= sqrt(pow(radius, 2)-pow(nextStackY, 2));
        for(int j=0; j<numSlices; j++) {

            GLfloat v0 [3] = {static_cast<GLfloat>(curStackRadius*cos((j+1)*angleFraction)), curStackY, static_cast<GLfloat>(curStackRadius*sin((j+1)*angleFraction))};
            GLfloat v1 [3] = {static_cast<GLfloat>(curStackRadius*cos(j*angleFraction)), curStackY, static_cast<GLfloat>(curStackRadius*sin(j*angleFraction))};
            GLfloat v2 [3] =  {static_cast<GLfloat>(nextStackRadius*cos(j*angleFraction)), nextStackY, static_cast<GLfloat>(nextStackRadius*sin(j*angleFraction))};
            GLfloat v3 [3] = {static_cast<GLfloat>(nextStackRadius*cos((j+1)*angleFraction)), nextStackY, static_cast<GLfloat>(nextStackRadius*sin((j+1)*angleFraction))};
            GLfloat lightVal[4] = {
                getLightValue(v0, v0, lightPosition, NULL),
                getLightValue(v1, v1, lightPosition, NULL),
                getLightValue(v2, v2, lightPosition, NULL),
                getLightValue(v3, v3, lightPosition, NULL)
            };
            int texIndex[4] = {
                getTextureIndex(lightVal[0]),
                getTextureIndex(lightVal[1]),
                getTextureIndex(lightVal[2]),
                getTextureIndex(lightVal[3])
            };
            

            glEnable(GL_BLEND);
            glEnable(GL_TEXTURE_2D);
            glBlendFunc(GL_DST_ALPHA, GL_ONE);
            //glEnable(GL_COLOR_MATERIAL);
            
//            glColor4f(1.0f, 1.0f, 1.0f, .25f);
            
            for (int k = 0; k < 4; ++k) {

                glActiveTexture(textures[texIndex[k]]);
                glBindTexture(GL_TEXTURE_2D, textures[texIndex[k]]);
                glBegin(GL_QUADS);
                glNormal3fv(v0); setSphereTexCoord(v0, radius); glVertex3fv(v0);
                glNormal3fv(v1); setSphereTexCoord(v1, radius); glVertex3fv(v1);
                glNormal3fv(v2); setSphereTexCoord(v2, radius); glVertex3fv(v2);
                glNormal3fv(v3); setSphereTexCoord(v3, radius); glVertex3fv(v3);
                glEnd();
       
            }
            glDisable(GL_TEXTURE_2D);
            glDisable(GL_BLEND);
            glDisable(GL_COLOR_MATERIAL);

        }
    }
    glDisable(GL_CULL_FACE);
}

void drawCylinder(int slices, int height, double radius){
    double angleFraction = (M_PI*2)/slices;
    int numHorizontalTiles = 21;
    double sliceWidth = sqrt(pow(radius*cos((1)*angleFraction)-radius*cos((2)*angleFraction),2)+pow(radius*sin((1)*angleFraction)-radius*sin((2)*angleFraction),2));
    double realCircumference = sliceWidth * slices;
    double realTextureWidth = static_cast<double>(slices) / numHorizontalTiles;
    double realTextureHeight = height * numHorizontalTiles / realCircumference;
    
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    
    double topTriangleTexHeight = radius/sliceWidth;
    
    for(int i=0; i<slices; i++){
        
        // Draw the top and bottom circles
        glBegin(GL_TRIANGLES);
        glNormal3f(0, 1, 0);
        glTexCoord2i(.5,0);glVertex3f(0, height, 0);
        glTexCoord2i(topTriangleTexHeight,0); glVertex3f(radius*cos((i+1)*angleFraction), height, radius*sin((i+1)*angleFraction));
        glTexCoord2i(topTriangleTexHeight,1);glVertex3f(radius*cos(i*angleFraction), height, radius*sin(i*angleFraction));
        
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
        glBindTexture(GL_TEXTURE_2D, textures[texIndex]);
        glBegin(GL_QUADS);
        glNormal3f(normal[0], normal[1], normal[2]);
        glTexCoord2i(i/realTextureWidth, realTextureHeight); glVertex3f(v0[0], v0[1], v0[2]);
        glTexCoord2i((i+1)/realTextureWidth, realTextureHeight); glVertex3f(radius*cos((i+1)*angleFraction), height, radius*sin((i+1)*angleFraction));
        glTexCoord2i((i+1)/realTextureWidth, 0); glVertex3f(radius*cos((i+1)*angleFraction), 0, radius*sin((i+1)*angleFraction));
        glTexCoord2i(i/realTextureWidth, 0); glVertex3f(radius*cos(i*angleFraction), 0, radius*sin(i*angleFraction));
        glDisable(GL_TEXTURE_2D);
        glEnd();
        
    }
}

void loadTAM(){
    
    int width, height;
    unsigned char* image;
    std::vector<std::vector<std::string>> tamNames(NUMBER_OF_TONES, std::vector<std::string>(NUMBER_OF_RESOLUTIONS));
    
    for (int tone = 0; tone < NUMBER_OF_TONES; ++tone) {
        for (int res = 0; res < NUMBER_OF_RESOLUTIONS; ++res) {
            char filename[100];
            //uncomment line below and comment line after to have correct path on Sam's comp, line below
            //is correct file path for github code 
//            sprintf(filename, "/Users/swarren/Downloads/ip2016skeleton/TAM_tone%d_resolution%d.bmp", tone, res);
            sprintf(filename, "/Users/swarren/Downloads/ip2016skeleton/ip2016/TAMimages/TAM_tone%d_resolution%d.bmp", tone, res);

            tamNames[tone][res] = std::string(filename);
        }
    }
    
    glGenTextures(NUMBER_OF_TONES, textures);
    for (int tone = 0; tone < NUMBER_OF_TONES; ++tone) {
        glBindTexture(GL_TEXTURE_2D, textures[tone]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, NUMBER_OF_RESOLUTIONS-1);
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
//    glEnable(GL_LIGHT0);
//    glEnable(GL_LIGHTING);
    //glEnable(GL_COLOR_MATERIAL);
    
    // set color of light0
    GLfloat white[] = {1,1,1,0};			// light color
    glLightfv(GL_LIGHT0, GL_DIFFUSE, white);	// set diffuse light color
    glLightfv(GL_LIGHT0, GL_SPECULAR, white);	// set specular light color
    
    // initialize background color to black
    glClearColor(0,0,0,.25);
    
    // enable depth buffering
    glDisable(GL_DEPTH_TEST);
}

void mouse(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN);
//        cout << "Left click with cursor at" << x << " " << y << endl;
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

void draw_quad()
{
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glEnable(GL_BLEND);
//    glTexEnvf(, <#GLenum pname#>, )
    glBlendFunc(GL_DST_ALPHA, GL_ONE);
    glEnable(GL_TEXTURE_2D);
//    glColor3f(1,1,1);
    for (int k = 0; k < 4; ++k) {
//        if (k==0) {
//          
//            glColor3f(1,0,0);
//        }
//        if (k==1)
//        {
//            glColor3f(0,1,0);
//        }
//        if (k==2) {
//            glColor3f(0, 0, 1);
//        }
    glBindTexture(GL_TEXTURE_2D, textures[0]);
        
//        printf("here\n");
//        checkGLError();
//        printf("there\n");
    glNormal3f(0, 0, 1);
    glBegin(GL_QUADS);
    glTexCoord2i(0, 0); glVertex3f(0, 0, 0);
    glTexCoord2i(1, 0); glVertex3f(1, 0, 0);
    glTexCoord2i(1, 1); glVertex3f(1, 1, 0);
    glTexCoord2i(0, 1); glVertex3f(0, 1, 0);
    glEnd();
//    }
//    glBegin(GL_QUADS);
//    glVertex3f(0, 0, 0);
//    glVertex3f(1, 0, 0);
//    glVertex3f(1, 1, 0);
//    glVertex3f(0, 1, 0);
//    glEnd();
}
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
}


void display()
{
    // clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // initialize modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    // set viewpoint position/orientation
    gluLookAt(10*sin(phi*3.14/180.0),5,10*cos(phi*3.14/180.0),0,0,0,0,1,0);
    // position of light0
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

    
//    drawCylinder(20, 1, 1);
    drawSphereWithQuads(1, 50, 50);
//    drawCone(2, 1, 15, 15);
//    drawCubeWithQuads(2, 50);
//    draw_quad();
    
    // swap buffers
    glutSwapBuffers();
    
    
    
}