#include<GL/gl.h>
#include<GL/glu.h>
#include<GL/glut.h>
#include <bits/stdc++.h>
#include <windows.h>
#include <stdlib.h>

using namespace std;
#define PI 3.1416

class BmpLoader
{
public:
    unsigned char* textureData;
    int iWidth, iHeight;

    BmpLoader(const char*);
    ~BmpLoader();

private:
    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
};

BmpLoader::BmpLoader(const char* filename)
{
    FILE *file=0;
    file=fopen(filename, "rb");
    if(!file)
        cout<<"File not found"<<endl;
    fread(&bfh, sizeof(BITMAPFILEHEADER),1,file);
    if(bfh.bfType != 0x4D42)
        cout<<"Not a valid bitmap"<<endl;
    fread(&bih, sizeof(BITMAPINFOHEADER),1,file);
    if(bih.biSizeImage==0)
        bih.biSizeImage=bih.biHeight*bih.biWidth*3;
    textureData = new unsigned char[bih.biSizeImage];
    fseek(file, bfh.bfOffBits, SEEK_SET);
    fread(textureData, 1, bih.biSizeImage, file);
    unsigned char temp;
    for(int i=0; i<bih.biSizeImage; i+=3)
    {
        temp = textureData[i];
        textureData[i] = textureData[i+2];
        textureData[i+2] = temp;

    }

    iWidth = bih.biWidth;
    iHeight = bih.biHeight;
    fclose(file);
}

BmpLoader::~BmpLoader()
{
    delete [] textureData;
}

void animate_car();


const float height=1500;
const float width=800;
GLfloat eye_x=1,eye_y=12,eye_z=115;
GLfloat look_x=1,look_y=17,look_z=0,theta=0;
float rm_h=1000,rm_w=0.4,rm_l=1000,
      window_rn=0,gate_rn=0,door_rn=0,fan_rn=0,rn,
      new_x,new_z,new_angle=1,new_r=0,
      car_x=0,car_z=80,car_angle=1,car_rn=-90,car_r=0,
     switch_baranda=0,side_check=0,stair_check=0,window_check=0,room_check=0,
    remove_upper_stairblock=1,main_floor_baranda_check=0;
bool start_angle=true,start_car=true,
flagSpotLight1=false,flagSpotLight2=false,flagSpotLight3=false,
flagLight1=true,flagLight2=false,flagLight3=false,
flagCar=false,
flagLeaf1=true,flagLeaf2=true,flagNode1=true,flagNode2=true;
unsigned int n=100;

//simple cube
static GLfloat v[8][3]=
{
    {-1,1,1},{1,1,1},{1,-1,1},{-1,-1,1},
    {-1,1,-1},{1,1,-1},{1,-1,-1},{-1,-1,-1}
};


static GLfloat v1[8][3]=
{
    {-1,1,1},{1,1,0.6},{1,-1,0.6},{-1,-1,1},
    {-1,1,-1},{1,1,-0.6},{1,-1,-0.6},{-1,-1,-1}
};



static GLfloat v2[8][3]=
{
    {-1,1,1},{1,1,1},{1,-1,1},{-0.5,-1,1},
    {-1,1,-1},{1,1,-1},{1,-1,-1},{-0.5,-1,-1}
};

//kuet design
static GLfloat v3[8][3] =
{
    {-1,1,1},{1,1,1},{2,-1,2},{-2,-1,2},
    {-1,1,-1},{1,1,-1},{2,-1,-2},{-2,-1,-2}
};

//Gate Triangle
static GLfloat v4[8][3]=
{
    {0,1,1},{0,1,1},{1,-1,1},{-1,-1,1},
    {0,1,-1},{0,1,-1},{1,-1,-1},{-1,-1,-1}
};

static GLfloat v_car[8][3] =
{
    {-1,1,1},{1,1.4,1},{1,-1,1},{-1,-1,1},
    {-1,1,-1},{1,1.4,-1},{1,-1,-1},{-1,-1,-1}
};

static GLfloat v_car2[8][3] =
{
    {-1,1,1},{1,1,1},{1.4,-1,1},{-1.4,-1,1},
    {-1,1,-1},{1,1,-1},{1.4,-1,-1},{-1.4,-1,-1}
};

static GLubyte index[6][4]=
{
    {0,3,2,1},{4,5,6,7},{0,4,7,3},
    {1,2,6,5},{0,1,5,4},{7,6,2,3}
};

static GLfloat c[18][3]=
{
    {0.404, 0.592, 0.176},// 0-->Floor Background
    {0.34,0.843,0.96},// 1-->Sky Background
    {0.55,0.62,0.71},//2-->Main Gate wall
    {0.97,0.81,0.68},//3-->Right & Left Wall
    {0.62,0.62,0.62},//4-->Floor Line
    {1.0,1.0,1.0},// 5-->White Board Base
    {0.8,0.45,0.21},//6--> Window BackUp
    {0.219, 0.22, 0.196},//7-->Window Grill
    {0.145,0.14,0.133},//8-->Fan
    {0.36,0.255,0.204},//9--> Fan Base
    {0.0,0.0,0.0},//10-->Black
    {0.0,1.0,0.0},//11-->Green
    {0.59,0.19,0.0},//12-->Table Base
    {0.796,0.6,0.0},//13-->Table Leg
    {0.57,0.37,0.26},//14-->Chair Base
    {0.16,0.17,0.19},//15-->Road
    {1.0,0.0,0.0},//16-->Red
    {0.196,0.58,0.192}//17-->Main_field Background

};

static GLuint texture_id[100]= {};


//<----------------------------------------------Curve Start------------------------------------------------->

const int L=4;
const int dgre=3;
int ncpt=L+1;
int clikd=0;
const int nt = 40;				//number of slices along x-direction
const int ntheta = 14;


GLfloat ctrlpoints[L+1][3] =
{

    { 0.2, 0.5, 0.0},
    { 0.1, 0.5, 0.0},
    { 0.0, 0.5, 0.0},
    { -0.1, 0.5, 0.0},
    { -0.2, 0.5, 0.0}
};

float wcsClkDn[3],wcsClkUp[3];
///////////////////////////////
class point1
{
public:
    point1()
    {
        x=0;
        y=0;
    }
    int x;
    int y;
}
clkpt[2];
int flag=0;
GLint viewport[4];
GLdouble modelview[16];
GLdouble projection[16];

//////////////////////////
void scsToWcs(float sx,float sy, float wcsv[3] );
void processMouse(int button, int state, int x, int y);
void matColor(float kdr, float kdg, float kdb,  float shiny, int frnt_Back=0, float ambFactor=1.0, float specFactor=1.0);
///////////////////////////

void scsToWcs(float sx,float sy, float wcsv[3] )
{

    GLfloat winX, winY, winZ;
    GLdouble worldX, worldY, worldZ;

    glGetDoublev( GL_PROJECTION_MATRIX, projection );
    glGetIntegerv( GL_VIEWPORT, viewport );

    winX = sx;
    winY = (float)viewport[3] - (float)sy;
    winZ = 0;
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &worldX, &worldY, &worldZ);
    wcsv[0]=worldX;
    wcsv[1]=worldY;
    wcsv[2]=worldZ;
}

//control points
long long nCr(int n, int r)
{
    if(r > n / 2)
        r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++)
    {
        ans *= n - r + i;
        ans /= i;
    }

    return ans;
}

//polynomial interpretation for N points
void BezierCurve ( double t,  float xy[2])
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L; i++)
    {
        int ncr=nCr(L,i);
        double oneMinusTpow=pow(1-t,double(L-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints[i][0];
        y+=coef*ctrlpoints[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);
}

///////////////////////
void setNormal(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(-Nx,-Ny,-Nz);
}

void wheelBezier()
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurve( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);
        }
        glEnd();
        x = x1;
        r = r1;
    }
}


//<---------------------------------Stop Curve----------------------------------------->
static void getNormal3p(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat x3, GLfloat y3, GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(Nx,Ny,Nz);
}

void cube(GLfloat C[],float a=0,float b=1,float c=0,float d=0,float e=1,float f=0,float g=1,float h=1,bool emp=false)
{
    float x=C[0],y=C[1],z=C[2];
    GLfloat m_amb[]= {x,y,z,1};
    GLfloat m_diff[]= {x,y,z,1};
    GLfloat m_spec[]= {1,1,1,1};
    GLfloat m_sh[]= {30};
    GLfloat m_em[] = {x,y,z,1};
    GLfloat m_no[] = {0, 0, 0, 1.0};

    glMaterialfv(GL_FRONT,GL_AMBIENT,m_amb);
    glMaterialfv(GL_FRONT,GL_DIFFUSE,m_diff);
    glMaterialfv(GL_FRONT,GL_SPECULAR,m_spec);
    glMaterialfv(GL_FRONT,GL_SHININESS,m_sh);

     if(emp) glMaterialfv(GL_FRONT, GL_EMISSION, m_em);
    else glMaterialfv(GL_FRONT, GL_EMISSION, m_no);

    glBegin(GL_QUADS);
    for(GLint i=0; i<6; i++)
    {
        getNormal3p(v[index[i][0]][0], v[index[i][0]][1], v[index[i][0]][2],
                    v[index[i][1]][0], v[index[i][1]][1], v[index[i][1]][2],
                    v[index[i][2]][0], v[index[i][2]][1], v[index[i][2]][2]);

        glTexCoord2f(a,b);
        glVertex3fv(&v[index[i][0]][0]);
        glTexCoord2f(c,d);
        glVertex3fv(&v[index[i][1]][0]);
        glTexCoord2f(e,f);
        glVertex3fv(&v[index[i][2]][0]);
        glTexCoord2f(g,h);
        glVertex3fv(&v[index[i][3]][0]);

    }

    glEnd();
}


void cube1(GLfloat C[],float a=0,float b=1,float c=0,float d=0,float e=1,float f=0,float g=1,float h=1)
{
    float x=C[0],y=C[1],z=C[2];
    GLfloat m_amb[]= {x,y,z,1};
    GLfloat m_diff[]= {x,y,z,1};
    GLfloat m_spec[]= {1,1,1,1};
    GLfloat m_sh[]= {30};

    glMaterialfv(GL_FRONT,GL_AMBIENT,m_amb);
    glMaterialfv(GL_FRONT,GL_DIFFUSE,m_diff);
    glMaterialfv(GL_FRONT,GL_SPECULAR,m_spec);
    glMaterialfv(GL_FRONT,GL_SHININESS,m_sh);

    glBegin(GL_QUADS);
    for(GLint i=0; i<6; i++)
    {
        getNormal3p(v1[index[i][0]][0], v1[index[i][0]][1], v1[index[i][0]][2],
                    v1[index[i][1]][0], v1[index[i][1]][1], v1[index[i][1]][2],
                    v1[index[i][2]][0], v1[index[i][2]][1], v1[index[i][2]][2]);

        glTexCoord2f(a,b);
        glVertex3fv(&v1[index[i][0]][0]);
        glTexCoord2f(c,d);
        glVertex3fv(&v1[index[i][1]][0]);
        glTexCoord2f(e,f);
        glVertex3fv(&v1[index[i][2]][0]);
        glTexCoord2f(g,h);
        glVertex3fv(&v1[index[i][3]][0]);

    }
    glEnd();
}


void cube2(GLfloat C[],float a=0,float b=1,float c=0,float d=0,float e=1,float f=0,float g=1,float h=1)
{
    float x=C[0],y=C[1],z=C[2];
    GLfloat m_amb[]= {x,y,z,1};
    GLfloat m_diff[]= {x,y,z,1};
    GLfloat m_spec[]= {1,1,1,1};
    GLfloat m_sh[]= {30};

    glMaterialfv(GL_FRONT,GL_AMBIENT,m_amb);
    glMaterialfv(GL_FRONT,GL_DIFFUSE,m_diff);
    glMaterialfv(GL_FRONT,GL_SPECULAR,m_spec);
    glMaterialfv(GL_FRONT,GL_SHININESS,m_sh);

    glBegin(GL_QUADS);
    for(GLint i=0; i<6; i++)
    {
        getNormal3p(v2[index[i][0]][0], v2[index[i][0]][1], v2[index[i][0]][2],
                    v2[index[i][1]][0], v2[index[i][1]][1], v2[index[i][1]][2],
                    v2[index[i][2]][0], v2[index[i][2]][1], v2[index[i][2]][2]);

        glTexCoord2f(a,b);
        glVertex3fv(&v2[index[i][0]][0]);
        glTexCoord2f(c,d);
        glVertex3fv(&v2[index[i][1]][0]);
        glTexCoord2f(e,f);
        glVertex3fv(&v2[index[i][2]][0]);
        glTexCoord2f(g,h);
        glVertex3fv(&v2[index[i][3]][0]);

    }
    glEnd();
}


void cube3(GLfloat C[],float a=0,float b=1,float c=0,float d=0,float e=1,float f=0,float g=1,float h=1)
{

    float x=C[0],y=C[1],z=C[2];
    GLfloat m_amb[]= {x,y,z,1};
    GLfloat m_diff[]= {x,y,z,1};
    GLfloat m_spec[]= {1,1,1,1};
    GLfloat m_sh[]= {30};

    glMaterialfv(GL_FRONT,GL_AMBIENT,m_amb);
    glMaterialfv(GL_FRONT,GL_DIFFUSE,m_diff);
    glMaterialfv(GL_FRONT,GL_SPECULAR,m_spec);
    glMaterialfv(GL_FRONT,GL_SHININESS,m_sh);

    glBegin(GL_QUADS);
    for(GLint i=0; i<6; i++)
    {
        getNormal3p(v3[index[i][0]][0], v3[index[i][0]][1], v3[index[i][0]][2],
                    v3[index[i][1]][0], v3[index[i][1]][1], v3[index[i][1]][2],
                    v3[index[i][2]][0], v3[index[i][2]][1], v3[index[i][2]][2]);

        glTexCoord2f(a,b);
        glVertex3fv(&v3[index[i][0]][0]);
        glTexCoord2f(c,d);
        glVertex3fv(&v3[index[i][1]][0]);
        glTexCoord2f(e,f);
        glVertex3fv(&v3[index[i][2]][0]);
        glTexCoord2f(g,h);
        glVertex3fv(&v3[index[i][3]][0]);
    }
    glEnd();
}


void cube4(GLfloat C[],float a=0,float b=1,float c=0,float d=0,float e=1,float f=0,float g=1,float h=1)
{
    float x=C[0],y=C[1],z=C[2];
    GLfloat m_amb[]= {x,y,z,1};
    GLfloat m_diff[]= {x,y,z,1};
    GLfloat m_spec[]= {1,1,1,1};
    GLfloat m_sh[]= {30};

    glMaterialfv(GL_FRONT,GL_AMBIENT,m_amb);
    glMaterialfv(GL_FRONT,GL_DIFFUSE,m_diff);
    glMaterialfv(GL_FRONT,GL_SPECULAR,m_spec);
    glMaterialfv(GL_FRONT,GL_SHININESS,m_sh);

    glBegin(GL_QUADS);
    for(GLint i=0; i<6; i++)
    {
        getNormal3p(v4[index[i][0]][0], v4[index[i][0]][1], v4[index[i][0]][2],
                    v4[index[i][1]][0], v4[index[i][1]][1], v4[index[i][1]][2],
                    v4[index[i][2]][0], v4[index[i][2]][1], v4[index[i][2]][2]);

        glTexCoord2f(a,b);
        glVertex3fv(&v4[index[i][0]][0]);
        glTexCoord2f(c,d);
        glVertex3fv(&v4[index[i][1]][0]);
        glTexCoord2f(e,f);
        glVertex3fv(&v4[index[i][2]][0]);
        glTexCoord2f(g,h);
        glVertex3fv(&v4[index[i][3]][0]);

    }
    glEnd();
}


void cube5(GLfloat C[],float a=0,float b=0,float c=0,float d=1,float e=1,float f=1,float g=1,float h=0)
{
    float x=C[0],y=C[1],z=C[2];
    GLfloat m_amb[]= {x,y,z,1};
    GLfloat m_diff[]= {x,y,z,1};
    GLfloat m_spec[]= {1,1,1,1};
    GLfloat m_sh[]= {30};

    glMaterialfv(GL_FRONT,GL_AMBIENT,m_amb);
    glMaterialfv(GL_FRONT,GL_DIFFUSE,m_diff);
    glMaterialfv(GL_FRONT,GL_SPECULAR,m_spec);
    glMaterialfv(GL_FRONT,GL_SHININESS,m_sh);

    glBegin(GL_QUADS);
    for(GLint i=0; i<6; i++)
    {
        getNormal3p(v_car[index[i][0]][0], v_car[index[i][0]][1], v_car[index[i][0]][2],
                    v_car[index[i][1]][0], v_car[index[i][1]][1], v_car[index[i][1]][2],
                    v_car[index[i][2]][0], v_car[index[i][2]][1], v_car[index[i][2]][2]);

        glVertex3fv(&v_car[index[i][0]][0]);
        glTexCoord2f(a,b);
        glVertex3fv(&v_car[index[i][1]][0]);
        glTexCoord2f(c,d);
        glVertex3fv(&v_car[index[i][2]][0]);
        glTexCoord2f(e,f);
        glVertex3fv(&v_car[index[i][3]][0]);
        glTexCoord2f(g,h);
    }
    glEnd();
}

void cube6(GLfloat C[],float a=0,float b=0,float c=0,float d=1,float e=1,float f=1,float g=1,float h=0)
{
    float x=C[0],y=C[1],z=C[2];
    GLfloat m_amb[]= {x,y,z,1};
    GLfloat m_diff[]= {x,y,z,1};
    GLfloat m_spec[]= {1,1,1,1};
    GLfloat m_sh[]= {30};

    glMaterialfv(GL_FRONT,GL_AMBIENT,m_amb);
    glMaterialfv(GL_FRONT,GL_DIFFUSE,m_diff);
    glMaterialfv(GL_FRONT,GL_SPECULAR,m_spec);
    glMaterialfv(GL_FRONT,GL_SHININESS,m_sh);

    glBegin(GL_QUADS);
    for(GLint i=0; i<6; i++)
    {
        getNormal3p(v_car2[index[i][0]][0], v_car2[index[i][0]][1], v_car2[index[i][0]][2],
                    v_car2[index[i][1]][0], v_car2[index[i][1]][1], v_car2[index[i][1]][2],
                    v_car2[index[i][2]][0], v_car2[index[i][2]][1], v_car2[index[i][2]][2]);

        glVertex3fv(&v_car2[index[i][0]][0]);
        glTexCoord2f(a,b);
        glVertex3fv(&v_car2[index[i][1]][0]);
        glTexCoord2f(c,d);
        glVertex3fv(&v_car2[index[i][2]][0]);
        glTexCoord2f(e,f);
        glVertex3fv(&v_car2[index[i][3]][0]);
        glTexCoord2f(g,h);
    }
    glEnd();
}


void animate()
{
    fan_rn++;
    if(fan_rn>=360)
        fan_rn=1;
    glutPostRedisplay();
}

void myKeyboardFunc(unsigned char key,int x,int y)
{
    GLfloat dx,dz,cx,cz;
    if(start_angle)
    {
        new_angle=atan((look_z-eye_z)/(look_x-eye_x));
        new_r=(look_z-eye_z)*(look_z-eye_z)+(look_x-eye_x)*(look_x-eye_x);
        new_r=sqrt(new_r);
        start_angle=false;
    }
    if(start_car)
    {
        car_angle=atan((look_z-car_z)/(look_x-car_x));
        car_r=(look_z-car_r)*(look_z-car_z)+(look_x-car_x)*(look_x-car_x);
        car_r=sqrt(car_r);
        start_car=false;
    }
    switch(key)
    {
    case 'd':
        new_angle*=57.296;
         if(new_angle>=360)
            new_angle=0;

        new_angle+=1;
        new_angle/=57.296;
        look_x=new_r*cos(new_angle)+eye_x;
        look_z=new_r*sin(new_angle)+eye_z;

        break;
    case 'a':
        new_angle*=57.296;
          if(new_angle<=-360)
            new_angle=0;
        new_angle-=1;
        new_angle/=57.296;
        look_x=new_r*cos(new_angle)+eye_x;
        look_z=new_r*sin(new_angle)+eye_z;

        break;
    case 'w':
        if(new_angle>=360)
            new_angle/=360;

        dx=cos(new_angle);
        dz=sin(new_angle);
        new_angle*=57.296;

        eye_x+=dx;
        eye_z+=dz;
        look_x+=dx;
        look_z+=dz;

        new_angle/=57.296;
        break;
    case 's':
        if(new_angle>=360)
            new_angle/=360;

        dx=cos(new_angle);
        dz=sin(new_angle);
        new_angle*=57.296;

         eye_x-=dx;
        eye_z-=dz;
        look_x-=dx;
        look_z-=dz;

        new_angle/=57.296;
        break;
    case '-':
        look_y-=0.1;
        eye_y-=0.1;
        break;
    case '+':
        look_y+=0.1;
        eye_y+=0.1;
        break;
    case 'i':
        look_y+=0.1;
        break;
    case 'k':
        look_y-=0.1;
        break;
    case 'y':
        rn++;
        break;
    case 't':
        rn--;
        break;
    case '1':
        if(gate_rn<90)
            gate_rn++;
        break;
    case '2':
        if(gate_rn>0)
            gate_rn--;
        break;
    case '3':
        if(door_rn<90)
            door_rn++;
        break;
    case '4':
        if(door_rn>0)
            door_rn--;
        break;
    case '5':
         if(window_rn<90)
            window_rn++;
        break;
    case '6':
        if(window_rn>0)
            window_rn--;
        break;
      case 'b':
        flagLight1=1-flagLight1;
        break;
     case 'c':
        flagLight2=1-flagLight2;
        break;
    case 'v':
        flagLight3=1-flagLight3;
        break;
     case 'n':
        flagSpotLight1=1-flagSpotLight1;
        break;
    case 'W':
        if(car_angle>=360)
            car_angle/=360;

        cx=cos(car_angle);
        cz=sin(car_angle);
        car_angle*=57.296;

        car_x+=cx;
        car_z+=cz;
        eye_x+=cx;
        eye_z+=cz;
        look_x+=cx;
        look_z+=cz;

        car_angle/=57.296;
        break;
    case 'S':
        if(car_angle>=360)
            car_angle/=360;

        dx=cos(car_angle);
        dz=sin(car_angle);
        car_angle*=57.296;


        car_x-=cx;
        car_z-=cz;
        eye_x-=cx;
        eye_z-=cz;
        look_x-=cx;
        look_z-=cz;

        car_angle/=57.296;
        break;
    case 'D':
        car_rn-=90;
        car_angle*=57.296;
        if(car_angle>=360)
            car_angle/=360;

        car_angle+=90;
        car_angle/=57.296;
        look_x=new_r*cos(car_angle)+car_x;
        look_z=new_r*sin(car_angle)+car_z;
        eye_x=car_x-30;
        eye_z=car_z;

        break;
    case 'A':
        car_rn+=90;
        car_angle*=57.296;
        if(car_angle>=360)
            car_angle/=360;

        car_angle-=90;
        car_angle/=57.296;
        look_x=new_r*cos(car_angle)+car_x;
        look_z=new_r*sin(car_angle)+car_z;
        eye_x=car_x;
        eye_z=car_z+30;
        break;

    case 'L':
        if(car_angle<-1)
        {
            car_x++,eye_x++,look_x++;
            car_z--,eye_z--,look_z--;
        }
        else
        {
            car_x++,eye_x++,look_x++;
            car_z++,eye_z++,look_z++;
        }
        break;
    case 'J':
        if(car_angle<-1)
        {
            car_x--,eye_x--,look_x--;
            car_z--,eye_z--,look_z--;
        }
        else
        {
            car_x++,eye_x++,look_x++;
            car_z--,eye_z--,look_z--;
        }
        break;
    case 9:
        exit(1);
    }
    cout<<"EyeX: "<<eye_x<<" EyeY: "<<eye_y<<" EyeZ: "<<eye_z<<endl;
    cout<<"LookX: "<<look_x<<" LookY: "<<look_y<<" LookZ: "<<look_z<<endl;
    cout<<"dx: "<<dx<<"  dz: "<<dz<<endl;
    cout<<"Angle: "<<new_angle*57.296<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    glutPostRedisplay();
}

static void resz(int width, int height)
{
    float rto=1.0*(width/height);
    glViewport(0, 0, width, width/rto);
}

void light1(float x,float y,float z)
{
    glEnable(GL_LIGHTING);

    GLfloat l_no[]= {0.0,0.0,0.0,1.0};
    GLfloat l_amb[] = {0.2,0.2,0.2,1};
    GLfloat l_diff[] = {1,1,1,1};
    GLfloat l_spec[]= {1,1,1,1};
    GLfloat l_pos[]= {x,y,z,1};

    glEnable(GL_LIGHT0);

    if(flagLight1)
    {
        glLightfv(GL_LIGHT0,GL_AMBIENT,l_amb);
        glLightfv(GL_LIGHT0,GL_DIFFUSE,l_diff);
        glLightfv(GL_LIGHT0,GL_SPECULAR,l_spec);
    }
    else if(!flagLight1)
    {
        glLightfv(GL_LIGHT0,GL_AMBIENT,l_no);
        glLightfv(GL_LIGHT0,GL_DIFFUSE,l_no);
        glLightfv(GL_LIGHT0,GL_SPECULAR,l_no);
    }
    glLightfv(GL_LIGHT0,GL_POSITION,l_pos);
}

void light2(float x,float y,float z)
{
    //<----Light-->
    glEnable(GL_LIGHTING);

    GLfloat l_no[]= {0.0,0.0,0.0,1.0};
    GLfloat l_amb[] = {0.2,0.2,0.2,1};
    GLfloat l_diff[] = {1,1,1,1};
    GLfloat l_spec[]= {1,1,1,1};
    GLfloat l_pos[]= {x,y,z,1};

    glEnable(GL_LIGHT1);

    if(flagLight2)
    {
        glLightfv(GL_LIGHT1,GL_AMBIENT,l_amb);
        glLightfv(GL_LIGHT1,GL_DIFFUSE,l_diff);
        glLightfv(GL_LIGHT1,GL_SPECULAR,l_spec);
    }
    else if(!flagLight2)
    {
        glLightfv(GL_LIGHT1,GL_AMBIENT,l_no);
        glLightfv(GL_LIGHT1,GL_DIFFUSE,l_no);
        glLightfv(GL_LIGHT1,GL_SPECULAR,l_no);
    }
    glLightfv(GL_LIGHT1,GL_POSITION,l_pos);
}


void light3(float x,float y,float z)
{
    //<----Light-->
    glEnable(GL_LIGHTING);

    GLfloat l_no[]= {0.0,0.0,0.0,1.0};
    GLfloat l_amb[] = {0.2,0.2,0.2,1};
    GLfloat l_diff[] = {1,1,1,1};
    GLfloat l_spec[]= {1,1,1,1};
    GLfloat l_pos[]= {x,y,z,1};

    glEnable(GL_LIGHT2);

    if(flagLight3)
    {
        glLightfv(GL_LIGHT2,GL_AMBIENT,l_amb);
        glLightfv(GL_LIGHT2,GL_DIFFUSE,l_diff);
        glLightfv(GL_LIGHT2,GL_SPECULAR,l_spec);
    }
    else if(!flagLight3)
    {
        glLightfv(GL_LIGHT2,GL_AMBIENT,l_no);
        glLightfv(GL_LIGHT2,GL_DIFFUSE,l_no);
        glLightfv(GL_LIGHT2,GL_SPECULAR,l_no);
    }
    glLightfv(GL_LIGHT2,GL_POSITION,l_pos);
}


void spot_light1(float x, float y, float z)
{
     //<----Light-->
    glEnable(GL_LIGHTING);

    GLfloat l_no[]= {0.0,0.0,0.0,1.0};
    GLfloat l_amb[] = {0.2,0.2,0.2,1};
    GLfloat l_diff[] = {1,1,1,1};
    GLfloat l_spec[]= {1,1,1,1};
    GLfloat l_pos[]= {x,y,z,1};

    glEnable( GL_LIGHT3);
     if(flagSpotLight1)
    {
        glLightfv(GL_LIGHT3,GL_AMBIENT,l_amb);
        glLightfv(GL_LIGHT3,GL_DIFFUSE,l_diff);
        glLightfv(GL_LIGHT3,GL_SPECULAR,l_spec);
    }
    else if(!flagSpotLight1)
    {
        glLightfv(GL_LIGHT3,GL_AMBIENT,l_no);
        glLightfv(GL_LIGHT3,GL_DIFFUSE,l_no);
        glLightfv(GL_LIGHT3,GL_SPECULAR,l_no);
    }
   glLightfv(GL_LIGHT3,GL_POSITION,l_pos);

    GLfloat direction[]= {0,-1,0,1};
    GLfloat cut_off=30;
    glLightfv(GL_LIGHT3,GL_SPOT_DIRECTION,direction);
    glLightf(GL_LIGHT3,GL_SPOT_CUTOFF,cut_off);

}

///<-----------------------------------------------MAIN BACKGROUND------------------------------------------------------------------>
void main_background(float h,float w,float l)
{
    //<------------------------Floor--------------------------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[0]);

    glPushMatrix();
    glScalef(h+w,w,l);
    cube(c[5],0,500,0,0,500,0,500,500);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);

    //<------------------------Front Side--------------------------------->

    glPushMatrix();
    glTranslatef(0,(h+w),-(h-w));
    glScalef(h,l,w);
    cube(c[1]);
    glPopMatrix();

    //<------------------------Ceil--------------------------------->
    glPushMatrix();
    glTranslatef(0,2*(h+w),0);
    glScalef(h+w,w,l);
    cube(c[1]);
    glPopMatrix();

    //<------------------------Left Side--------------------------------->
    glPushMatrix();
    glTranslatef(-h,(h+w),0);
    glScalef(w,h,l);
    cube(c[1]);
    glPopMatrix();

    //<------------------------Right Side--------------------------------->
    glPushMatrix();
    glTranslatef(h,(h+w),0);
    glScalef(w,h,l);
    cube(c[1]);
    glPopMatrix();

    //<------------------------Back Side--------------------------------->
    glPushMatrix();
    glTranslatef(0,(h+w),(h-w));
    glScalef(h,l,w);
    cube(c[1]);
    glPopMatrix();
}


///<-----------------------------------------------MAIN ROAD---------------------------------------------------------------------------->

void road_side_lamp()
{
    float light_h=1,light_w=0.2,light_l=1;
    //<---------TubeLight------------->
    glPushMatrix();
    glScalef(light_h,light_w,light_l);
    cube(c[5],true);
    glPopMatrix();
}

void main_road()
{
    float x,y,z,road_h=10,road_w=0.3,road_l=100;


    //<--------------------------straight1------------------------------>
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[2]);

    y=road_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(road_h,road_w,road_l);
    cube(c[5],0,0,1,0,1,10,0,10);
    glPopMatrix();

    //<--------------------------Turn Right1------------------------------>

    x=road_l-road_h;
    z=(road_h+road_l);
    glPushMatrix();
    glTranslatef(x,y,-z);
    glScalef(road_l,road_w,road_h);
    cube(c[5],1,0,1,10,0,10,0,0);
    glPopMatrix();

    //<-------------------------straight2-------------------------------------->
    x=2*road_l;
    z=2*road_l;
    for(int i=0; i<3; i++)
    {
        glPushMatrix();
        glTranslatef(x,y,-z);
        glScalef(road_h,road_w,road_l);
        cube(c[5],0,0,1,0,1,10,0,10);
        glPopMatrix();
        z+=(2*road_l);
    }

    glDisable(GL_TEXTURE_2D);

}


///<-----------------------------------------------KUET DESIGN---------------------------------------------------------------------------->

//<--------------'K' Main  Design----------------->
void k_create()
{
    float main_h=0.5,main_w=3,main_l=1,
          side_h=2,side_w=0.5,side_l=1;

    //<--------------------Straight----------------->
    glPushMatrix();
    glTranslatef(0,main_w,0);
    glScalef(main_h,main_w,main_l);
    cube(c[5]);
    glPopMatrix();

    //<--------------Upper ----------------->
    glPushMatrix();
    glTranslatef(side_h,main_w+1.2,0);
    glRotatef(40,0,0,1);
    glScalef(side_h,side_w,side_l);
    cube(c[5]);
    glPopMatrix();

    //<--------------Lower ----------------->
    glPushMatrix();
    glTranslatef(side_h,main_w-1.2,0);
    glRotatef(-45,0,0,1);
    glScalef(side_h,side_w,side_l);
    cube(c[5]);
    glPopMatrix();

}

//<--------------'U' Main  Design----------------->
void u_create()
{
    float x,y,z, main_h=0.5,main_w=3,main_l=1,
                 side_h=1.5,side_w=0.5,side_l=1;

    //<---------------------Left & Right Side------------->
    x=main_h+side_h;
    y=main_w;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,y,0);
        glScalef(main_h,main_w,main_l);
        cube(c[5]);
        glPopMatrix();
        x*=(-1);
    }

    //<---------------------Lower------------->
    y=side_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(side_h,side_w,side_l);
    cube(c[5]);
    glPopMatrix();

}

//<--------------'E' Main  Design----------------->
void e_create()
{
    float x,y,z, main_h=0.5,main_w=3,main_l=1,
                 side_h=1.5,side_w=0.5,side_l=1;

    //<---------------------Left  Side------------->
    x=main_h+side_h;
    y=main_w;
    glPushMatrix();
    glTranslatef(-x,y,0);
    glScalef(main_h,main_w,main_l);
    cube(c[5]);
    glPopMatrix();

//<---------------------Upper & Lower & Middle------------->
    y=side_w;
    for(int i=0; i<3; i++)
    {
        glPushMatrix();
        glTranslatef(0,y,0);
        glScalef(side_h,side_w,side_l);
        cube(c[5]);
        glPopMatrix();
        y+=2.5;
    }

}

 //<--------------'T' Main  Design----------------->
void t_create()
{
    float x,y,z, main_h=0.5,main_w=3,main_l=1,
                 side_h=2.5,side_w=0.5,side_l=1;

    //<---------------------Middle Main  ------------->
    y=main_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(main_h,main_w,main_l);
    cube(c[5]);
    glPopMatrix();

//<---------------------Upper ------------->
    y=2*main_w-side_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(side_h,side_w,side_l);
    cube(c[5]);
    glPopMatrix();
    y+=2.5;
}


void KUET(float base_h,float base_w,float base_l)
{
    float x,y,z;
    //<--------------K--------------------->
    x=base_h*0.9;
    y=2*base_w;
    glPushMatrix();
    glTranslatef(-x,y,0);
    k_create();
    glPopMatrix();

    //<--------------U--------------------->
    x*=(0.5);
    glPushMatrix();
    glTranslatef(-x,y,0);
    u_create();
    glPopMatrix();

    //<--------------E--------------------->
    glPushMatrix();
    glTranslatef(0,y,0);
    e_create();
    glPopMatrix();

    //<--------------T--------------------->
    x*=(0.8);
    glPushMatrix();
    glTranslatef(x,y,0);
    t_create();
    glPopMatrix();
}

void kuet_design()
{
    float x,y,z, base_h=15,base_w=3,base_l=5;


    //<-------------Base---------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[6]);

    y=base_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(base_h,base_w,base_l);
    cube3(c[5],0,3,0,0,3,0,3,3);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);

    //<-------------------------KUET-------------------------->
    glPushMatrix();
    glTranslatef(3,0,0);
    KUET(base_h,base_w,base_l);
    glPopMatrix();

}



///<-----------------------------------------------MAIN GATE------------------------------------------------------------------>

void triangle(float h,float w,float l)
{
    float gate_h=h,gate_w=w,gate_l=l;
    glPushMatrix();
    glScalef(gate_h,gate_w,gate_l);
    cube4(c[5]);
    glPopMatrix();
}


void triangle_with_logo()
{
    glPushMatrix();
    glTranslatef(12,0,29);
    //<-------------------Triangle--------------->
     glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[43]);
    glPushMatrix();
    glTranslatef(0,12,0);
    glRotatef(90,0,0,1);
    triangle(12,5,1);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

     //<-------------------kuet logo--------------->
     glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[41]);
    glPushMatrix();
    glTranslatef(2,14,1);
   glScalef(1,1,0.01);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

     //<-------------------kuet title--------------->
     glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[42]);
    glPushMatrix();
    glTranslatef(1,11,1);
   glScalef(4,1.5,0.01);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glPopMatrix();
}

void gate()
{
    float x,y,z,left_side_h=0.1,left_side_w=30,left_side_l=1,
    right_side_h=6.5,right_side_w=0.1,right_side_l=1;
    //<-----------------Triangle Logo Part--------------------->
    glPushMatrix();
    triangle_with_logo();
    glPopMatrix();

    //<----------------Left Side--------------------------->
    glPushMatrix();
    glTranslatef(-20,0,29);
     glRotatef(-45,0,0,1);
    glTranslatef(0,30,0);
    glScalef(left_side_h,left_side_w,left_side_l);
    cube(c[5]);
    glPopMatrix();

    //<------------------Right Middle------------------------->

    glPushMatrix();
    glTranslatef(10.5,24,29);
    glScalef(right_side_h,right_side_w,right_side_l);
    cube(c[5]);
    glPopMatrix();

    //<------------------Right Lower------------------------->
    right_side_h=6.15;
    glPushMatrix();
    glTranslatef(4,24,29);
    glRotatef(-75.5,0,0,1);
    glTranslatef(right_side_h,0,0);
    glScalef(right_side_h,right_side_w,right_side_l);
    cube(c[5]);
    glPopMatrix();

    //<------------------Right Upper------------------------->
    right_side_h=5.5;
    glPushMatrix();
    glTranslatef(17,24,29);
    glRotatef(-78,0,0,1);
    glTranslatef(-right_side_h,0,0);
    glScalef(right_side_h,right_side_w,right_side_l);
    cube(c[5]);
    glPopMatrix();



}


void orginal_gate_function(float gate_h,float gate_w,float gate_l)
{
    float x=4.65,y,z,
          gate_backup_h=4.65,gate_backup_w=0.05,gate_backup_l=0.03;
    //<--------------Left & Right Side--------->
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,gate_w,0);
        glScalef(gate_h,gate_w,gate_l);
        cube(c[5]);
        glPopMatrix();
        x*=(-1);
    }

    //<--------------Middle BackUp--------->
    y=0;
    while(y<2*gate_w)
    {
        glPushMatrix();
        glTranslatef(0,y,0);
        glScalef(gate_backup_h,gate_backup_w,gate_backup_l);
        cube(c[5]);
        glPopMatrix();
        y+=0.2;
    }

}

void main_gate()
{
    float x,y,z, gate_h=0.2,gate_w=4,gate_l=0.2,
                 side_road_h=10,side_road_w=0.5,side_road_l=34,
                 road_seat_h=7,road_seat_w=1,road_seat_l=30,
                 road_grass_h=5,road_grass_w=0.05,road_grass_l=25,
                 side_wall_h=50,side_wall_w=4.5,side_wall_l=1,
                 logo_h=4,logo_w=1,logo_l=0.1;

    //<------------Gate Side Wall------------------>
    x=side_wall_h+10;
    y=rm_w+side_wall_w;
    z=side_wall_l+30;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(-x,y,z);
        glScalef(side_wall_h,side_wall_w,side_wall_l);
        cube(c[2]);
        glPopMatrix();
        x*=(-1);
    }

    //<------------Gate Side Road------------------>
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[4]);

    glPushMatrix();
    x=side_road_h+10;
    y=rm_w+side_road_w;
    z=2*side_wall_l+side_road_l+30;
    glTranslatef(-x,y,z);
    glScalef(side_road_h,side_road_w,side_road_l);
    cube(c[5],0,0,5,0,5,5,0,5);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);

    //<------------Gate Side Road Seat------------------>
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[3]);

    glPushMatrix();
    y=rm_w+2*side_road_w+road_seat_w;
    glTranslatef(-x,y,z);
    glScalef(road_seat_h,road_seat_w,road_seat_l);
    cube(c[5],0,0,10,0,10,10,0,10);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);

    //<------------Gate Side Road Grass------------------>
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[6]);

    glPushMatrix();
    y=rm_w+2*(side_road_w+road_seat_w)+road_grass_w;
    glTranslatef(-x,y,z);
    glScalef(road_grass_h,road_grass_w,road_grass_l);
    cube(c[5],0,0,1,0,1,5,0,5);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<---------------------Gate Side KUET Logo-------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[5]);

    glPushMatrix();
    z=2*side_wall_l+side_road_l+5;
    y=rm_w+2*(side_road_w+road_seat_w+road_grass_w)+logo_w;
    glTranslatef(-x,y,z);
    glRotatef(20,0,1,0);
    glScalef(logo_h,logo_w,logo_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<-----------------------Main Gate Left Part---------------------------->
    x=9.5;
    y=rm_w+0.6;
    z=side_wall_l+30;

    glPushMatrix();
    glTranslatef(-x,y,z);
    glRotatef(gate_rn,0,1,0);
    glTranslatef(4.65,0,0);
    orginal_gate_function(gate_h,gate_w,gate_l);
    glPopMatrix();

    //<-----------------------Main Gate Right Part---------------------------->
    glPushMatrix();
    glTranslatef(x,y,z);
    glRotatef(-gate_rn,0,1,0);
    glTranslatef(-4.65,0,0);
    orginal_gate_function(gate_h,gate_w,gate_l);
    glPopMatrix();

    //<-----------------------Gate Design----------------------------->
    glPushMatrix();
    gate();
    glPopMatrix();

}



///<-----------------------------------------------CENTRAL SHAHID MINAR---------------------------------------------------------------------------->

void minar_circle(float h,float w,float l)
{
    glPushMatrix();
    for(float i=0; i<=360; i+=0.05)
    {
        glPushMatrix();
        glRotatef(i,0,0,1);
        glScalef(h,l,w);
        cube(c[16]);
        glPopMatrix();
    }
    glPopMatrix();
}

void minar_part()
{
    float x,y,z,main_base_h=18,main_base_w=1,main_base_l=8,
    minar_h,minar_w,minar_l,
    upper_h,upper_w,upper_l,
    circle_h,circle_w,circle_l;

     //<--------------Main Base of Shahid Miner------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[4]);

    y=main_base_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(main_base_h,main_base_w,main_base_l);
    cube(c[5],0,0,5,0,5,3,0,3);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);

    //<------------------Main Part of Shahid Minar -> Middle----------------------------->
    minar_h=0.3;
    minar_w=10;
    minar_l=0.3;
    //<---------------(1)Minar Part------------------>
    x=-2;
    y=2*main_base_w+minar_w;
    for(int i=0; i<3; i++)
    {
        glPushMatrix();
        glTranslatef(x,y,0);
        glScalef(minar_h,minar_w,minar_l);
        cube(c[5]);
        x+=2;
        glPopMatrix();
    }

    upper_h=minar_h*1.5+2;
    upper_w=0.3;
    upper_l=0.3;
    //<---------------(2)Upper Part------------------>
    y=2*(main_base_w+minar_w);
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(upper_h,upper_w,upper_l);
    cube(c[5]);
    glPopMatrix();

    circle_h=5;
    circle_w=0.05;
    circle_l=0.05;
    //<-----------------(3)Circle Part------------------>
    y=2*main_base_w+minar_w;
    z=minar_l+circle_l;
    glPushMatrix();
    glTranslatef(0,y,-z);
    minar_circle(circle_h,circle_w,circle_l);
    glPopMatrix();

    //<------------------Main Part Right Side1----------------------------->
    minar_h=0.3;
    minar_w=8;
    minar_l=0.3;
    //<---------------(1)Minar Part------------------>
    x=4;
    y=2*main_base_w+minar_w;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,y,0);
        glScalef(minar_h,minar_w,minar_l);
        cube(c[5]);
        x+=2;
        glPopMatrix();
    }

    upper_h=minar_h+1;
    upper_w=0.3;
    upper_l=0.3;
    //<---------------(2)Upper Part------------------>
    x=5;
    y= 2*(main_base_w+minar_w);
    glPushMatrix();
    glTranslatef(x,y,0);
    glScalef(upper_h,upper_w,upper_l);
    cube(c[5]);
    glPopMatrix();

    //<------------------Main Part Left Side1----------------------------->
    minar_h=0.3;
    minar_w=8;
    minar_l=0.3;
    //<---------------(1)Minar Part------------------>
    x=-4;
    y=2*main_base_w+minar_w;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,y,0);
        glScalef(minar_h,minar_w,minar_l);
        cube(c[5]);
        x-=2;
        glPopMatrix();
    }

    upper_h=minar_h+1;
    upper_w=0.3;
    upper_l=0.3;
    //<---------------(2)Upper Part------------------>
    x=-5;
    y=2*(main_base_w+minar_w);
    glPushMatrix();
    glTranslatef(x,y,0);
    glScalef(upper_h,upper_w,upper_l);
    cube(c[5]);
    glPopMatrix();


    //<------------------Main Part Right Side2----------------------------->
    minar_h=0.3;
    minar_w=5;
    minar_l=0.3;
    //<---------------(1)Minar Part------------------>
    x=8;
    y=2*main_base_w+minar_w;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,y,0);
        glScalef(minar_h,minar_w,minar_l);
        cube(c[5]);
        x+=2;
        glPopMatrix();
    }

    upper_h=minar_h+1;
    upper_w=0.3;
    upper_l=0.3;
    //<---------------(2)Upper Part------------------>
    x=9;
    y=2*(main_base_w+minar_w);
    glPushMatrix();
    glTranslatef(x,y,0);
    glScalef(upper_h,upper_w,upper_l);
    cube(c[5]);
    glPopMatrix();

    //<------------------Main Part Left Side2----------------------------->
    minar_h=0.3;
    minar_w=5;
    minar_l=0.3;
    //<---------------(1)Minar Part------------------>
    x=-8;
    y=2*main_base_w+minar_w;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,y,0);
        glScalef(minar_h,minar_w,minar_l);
        cube(c[5]);
        x-=2;
        glPopMatrix();
    }

    upper_h=minar_h+1;
    upper_w=0.3;
    upper_l=0.3;
    //<---------------(2)Upper Part------------------>
    x=-9;
    y=2*(main_base_w+minar_w);
    glPushMatrix();
    glTranslatef(x,y,0);
    glScalef(upper_h,upper_w,upper_l);
    cube(c[5]);
    glPopMatrix();

}

void central_shahid_minar()
{
    float x,y,z,main_stair_h=20,main_stair_w=0.3,main_stair_l=50,
    middle_stair_h=20,middle_stair_w=0.2,middle_stair_l=35,
    lower_stair_h=20,lower_stair_w=0.2,lower_stair_l=15,check;


    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[3]);
    //<------------------Main Stair--------------------->
    y=main_stair_w;
    z=main_stair_l-lower_stair_l;
    glPushMatrix();
    glTranslatef(0,y,z);
    glScalef(main_stair_h,main_stair_w,main_stair_l);
    cube(c[5],0,0,10,0,10,20,0,20);
    glPopMatrix();


    //<------------------Middle Stair--------------------->
    y=2*main_stair_w+middle_stair_w;
    z=middle_stair_l-lower_stair_l;
    check=middle_stair_l;
    while(middle_stair_l>=33)
    {
        glPushMatrix();
        glTranslatef(0,y,z);
        glScalef(middle_stair_h,middle_stair_w,middle_stair_l);
        cube(c[5],0,0,10,0,10,20,0,20);
        glPopMatrix();
        middle_stair_l-=1;
        z-=1;
        y+=(2*middle_stair_w);
    }
    middle_stair_l=check;


    //<------------------Lower Stair--------------------->
    y=2*main_stair_w+6*middle_stair_w+lower_stair_w;
    z=0;
    check=lower_stair_l;
    while(lower_stair_l>=10)
    {
        glPushMatrix();
        glTranslatef(0,y,z);
        glScalef(lower_stair_h,lower_stair_w,lower_stair_l);
        cube(c[5],0,0,10,0,10,8,0,8);
        glPopMatrix();
        lower_stair_l-=1;
        z-=1;
        y+=(2*lower_stair_w);
    }
    lower_stair_l=check;
    glDisable(GL_TEXTURE_2D);


   glPushMatrix();
    y=2*main_stair_w+6*middle_stair_w+12*lower_stair_w;
    z=lower_stair_l-8;
    glTranslatef(0,y,-z);
    minar_part();
    glPopMatrix();
}




///<-----------------------------------------------KHAN JAHAN ALI HALL START--------------------------------------->


//<-----------------------------------------------SINGLE FAN--------------------------------------->
void room_fan()
{
    float x,y,z,stand_h=0.02,stand_w=0.4,stand_l=0.02,
                base_h=0.2,base_w=0.04,base_l=0.2,
                leg_h=0.6,leg_w=0.015,leg_l=0.15;

    //<-------------------Stand---------------------->
    y=stand_w+base_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(stand_h,stand_w,stand_l);
    cube(c[7]);
    glPopMatrix();

    //<-------------------Base---------------------->
    glPushMatrix();
    glScalef(base_h,base_w,base_l);
    cube(c[9]);
    glPopMatrix();

    //<-------------------Leg---------------------->
    float rotatn=10;
    x=leg_h;
    for(int i=0; i<3; i++)
    {
        glPushMatrix();
        glRotatef(rotatn,0,1,0);
        glTranslatef(-x,0,0);
        glScalef(leg_h,leg_w,leg_l);
        cube1(c[8]);
        glPopMatrix();
        rotatn+=120;
    }
}


//<-----------------------------------------------SINGLE WINDOW--------------------------------------->
void room_window()
{
    float x,y,z,window_h=1.3,window_w=2.5,window_l=0.02,
                support_h=2.6,support_w=0.1,support_l=0.12,
                middle_h=2.5,middle_w=0.08,middle_l=0.1,
                helper_h=0.05,helper_w=0.2,helper_l=0.1,
                l_side_h=0.65,l_side_w=2.6,l_side_l=0.1,
                u_side_h=4,u_side_w=1.2,u_side_l=0.1;

    //<---------------Left & Right Support------------------->
    x=support_h;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,0,0);
        glScalef(support_w,support_h,support_l);
        cube(c[6]);
        glPopMatrix();
        x*=(-1);
    }

    //<---------------Lower & Upper Support------------------->
    y=support_h-support_w;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(0,y,0);
        glScalef(support_h,support_w,support_l);
        cube(c[6]);
        glPopMatrix();
        y*=(-1);
    }


    //<----------------------Middle Part1------------------->
    float var1=-(middle_h-2*helper_w-middle_w),var2=abs(var1);
    while(var1<=var2)
    {
        glPushMatrix();
        glTranslatef(0,var1,0);

        glPushMatrix();
        glScalef(middle_h,middle_w,middle_l);
        cube(c[7]);
        glPopMatrix();

        //<----------------------Middle Upper Part-------------->
        y=helper_w;
        glPushMatrix();
        glTranslatef(0,y,0);
        glScalef(helper_h,helper_w,helper_l);
        cube(c[7]);
        glPopMatrix();

        //<----------------------Middle Lower Part-------------->
        y=helper_w;
        x=-1.5;
        for(int i=0; i<2; i++)
        {
            glPushMatrix();
            glTranslatef(x,-y,0);
            glScalef(helper_h,helper_w,helper_l);
            cube(c[7]);
            glPopMatrix();
            x*=(-1);
        }

        //<-------------------Middle Part2------------------>
        y=2*helper_w+middle_w;
        glPushMatrix();
        glTranslatef(0,y,0);
        glScalef(middle_h,middle_w,middle_l);
        cube(c[7]);
        glPopMatrix();

        glPopMatrix();
        var1+=(4*helper_w+2*middle_w);

    }


    if(window_check==1)
    {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D,texture_id[21]);
        //<-----------------Main Left window --------------------->
        x=support_h;
        z=support_l;
        glPushMatrix();
        glTranslatef(-x,0,z);
        glRotatef(-window_rn,0,1,0);
        glTranslatef(window_h,0,0);
        glScalef(window_h,window_w,window_l);
        cube(c[5]);
        glPopMatrix();

        //<-----------------Main Right window  --------------------->
        glPushMatrix();
        glTranslatef(x,0,z);
        glRotatef(window_rn,0,1,0);
        glTranslatef(-window_h,0,0);
        glScalef(window_h,window_w,window_l);
        cube(c[5]);
        glPopMatrix();

        glDisable(GL_TEXTURE_2D);
    }


    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[20]);
    //<-----------------Left & Right of the Window--------------->
    x=support_h+support_w+l_side_h;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,0,0);
        glScalef(l_side_h,l_side_w,l_side_l);
        cube(c[5]);
        glPopMatrix();
        x*=(-1);
    }

    //<-----------------Upper & Lower of the Window--------------->
    y=support_h+u_side_w;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(0,y,0);
        glScalef(u_side_h,u_side_w,u_side_l);
        cube(c[5]);
        glPopMatrix();
        y*=(-1);
    }
    glDisable(GL_TEXTURE_2D);

}


//<-------------------------------------------------------SINGLE DOOR---------------------------------------------->
void room_door()
{
    float x,y,z,door_h=2,door_w=4,door_l=0.1,
                support_h=2,support_w=4,support_l=0.12,
                upper_h=2.2,upper_w=1,upper_l=0.1;

    //<----------Upper & Lower Support-------------->
    y=support_w;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(0,y,0);
        glScalef(support_h+support_l,support_l,support_l);
        cube(c[4]);
        glPopMatrix();
        y*=(-1);
    }
    //<----------Left & Right Support-------------->
    x=support_h;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,0,0);
        glScalef(support_l,support_w,support_l);
        cube(c[4]);
        glPopMatrix();
        x*=(-1);
    }

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[22]);
    //<-----------------Main Door--------------------->
    x=support_h;
    z=support_l;
    glPushMatrix();
    glTranslatef(x,0,-z);
    glRotatef(-door_rn,0,1,0);
    glTranslatef(-door_h,0,0);
    glScalef(door_h,door_w,door_l);
    cube(c[5]);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[20]);
    //<---------------------Upper Side of the Door------------->

    y=door_w+upper_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(upper_h,upper_w,upper_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

}


//<-----------------------------------------------------ROOM BACKGROUND-------------------------------->
void room_background()
{
    float x,y,z,room_h=10,room_w=0.1,room_l=12,room_y=5;

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[19]);
    //<-----------------Floor----------------->
    y=room_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(room_h+room_w,room_w,room_l);
    cube(c[5],0,0,10,0,10,10,0,10);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[23]);
    //<-----------------Ceil----------------->
    y=2*room_y;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(room_h+room_w,room_w,room_l);
    cube(c[5]);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[20]);
    //<-----------------Left Side---------------->
    x=room_h;
    y=room_y;
    glPushMatrix();
    glTranslatef(-x,y,0);
    glScalef(room_w,room_y,room_l);
    cube(c[5]);
    glPopMatrix();

    //<-----------------Right Side---------------->
    glPushMatrix();
    glTranslatef(x,y,0);
    glScalef(room_w,room_y,room_l);
    cube(c[5]);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);

    //<-------------------Front Part Left & Right-------------------->
    x=room_h-4;
    z=room_l;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(-x,y,z);
        room_window();
        glPopMatrix();
        x*=(-1);
    }

    //<--------------------Front Middle Part------------------------->

    y--;
    glPushMatrix();
    glTranslatef(0,y,z);
    room_door();
    glPopMatrix();


    //<-------------------Back  Left & Right & Middle Part-------------------->
    x=room_h-4;
    y++;
    z=room_l;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(-x,y,-z);
        glRotatef(180,0,1,0);
        room_window();
        glPopMatrix();
        x*=(-1);
    }

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[20]);
    glPushMatrix();
    glTranslatef(0,y,-z);
    glScalef(2,5,0.1);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}



//<-------------------------------------------------------ROOM CHAIR------------------------------------->
void room_chair()
{
    float x,y,z,
          chair_h=0.6,chair_w=0.05,chair_l=0.6,
          leg_h=0.05,leg_w=1.2,leg_l=0.05,
          backup_h=0.5,backup_w=0.05,backup_l=0.05,
          side_h=0.05,side_w=0.05,side_l=0.5;

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[15]);

    glPushMatrix();
    glTranslatef(0,leg_w,0);

    ///<----------------Base of Chair------------------->
    glPushMatrix();
    glScalef(chair_h,chair_w,chair_l);
    cube(c[5]);
    glPopMatrix();

    ///<-------------------Back Leg--------------------->
    x=chair_h-leg_h;
    z=chair_l-leg_l;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,0,-z);
        glScalef(leg_h,leg_w,leg_l);
        cube(c[5]);
        glPopMatrix();
        x*=(-1);
    }
    ///<-------------------Front Leg--------------------->
    leg_w/=2;
    x=chair_h-leg_h;
    y=chair_w+leg_w;
    z=chair_l-leg_l;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,-y,z);
        glScalef(leg_h,leg_w,leg_l);
        cube(c[5]);
        glPopMatrix();
        x*=(-1);
    }

    ///<----------------------------Upper BackUp--------------->
    y=0.25;
    z=chair_l-leg_l;
    for(int i=0; i<4; i++)
    {
        glPushMatrix();
        glTranslatef(0,y,-z);
        glScalef(backup_h,backup_w,backup_l);
        cube(c[5]);
        glPopMatrix();
        y+=0.25;
    }

    ///<----------------------------Lower Side--------------->
    x=chair_h-side_h;
    y=0.75;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,-y,0);
        glScalef(side_h,side_w,side_l);
        cube(c[5]);
        glPopMatrix();
        x*=(-1);
    }

    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

}


//<-------------------------------------------------------ROOM TABLE------------------------------------->
void room_laptop()
{
    float x,y,z,base_h=0.5,base_w=0.01,base_l=0.3,
                front_h=0.5,front_w=0.3,front_l=0.005;

    //<----------------Base--------------------->
    y=base_w;
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[36]);
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(base_h,base_w,base_l);
    cube(c[5],0,0,1,0,1,1,0,1);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<---------------Front-------------------->
    y=0;
    z=base_l+front_l;
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[37]);
    glPushMatrix();
    glTranslatef(0,0,-z);
    glRotatef(rn,1,0,0);
    glTranslatef(0,front_w,0);
    glScalef(front_h,front_w,front_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);


}
void room_table()
{
    float   x,y,z,table_h=1.5,table_w=0.1,table_l=1.5,
                  leg_h=0.05,leg_w=0.8,leg_l=0.05,
                  leg_support_h=table_h,leg_support_w=0.05,leg_support_l=0.05,
                  upper_h=table_h,upper_w=table_w/2,upper_l=table_l/2,
                  back_h=table_h,back_w=0.85,back_l=0.05,
                  side_h=0.05,side_w=0.85,side_l=table_l/2,
                  book_h=0.05,book_w=0.22,book_l=table_l/2;


    glPushMatrix();
    glTranslatef(0,leg_w,0);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[14]);
    //<------------Table Base----------->
    y=table_w+leg_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(table_h,table_w,table_l);
    cube(c[5]);
    glPopMatrix();

    //<-------------Leg Side---------------->
    x=table_h-leg_h;
    z=table_l-leg_l;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<2; j++)
        {
            glPushMatrix();
            glTranslatef(x,0,z);
            glScalef(leg_h,leg_w,leg_l);
            cube(c[5]);
            glPopMatrix();
            x*=(-1);
        }
        z*=(-1);
    }

    //<-------------Leg Support Left & Right Side---------------->


    //<--------------Upper Part---------------->
    y=leg_w+2*table_w+upper_w+0.8;
    z=upper_l;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(0,y,-z);
        glScalef(upper_h,upper_w,upper_l);
        cube(c[5]);
        glPopMatrix();
        y+=(0.8+upper_w);
    }

    //<--------------Upper Back backup---------------->
    y=leg_w+2*table_w+back_w;
    z=table_l-back_l;

    glPushMatrix();
    glTranslatef(0,y,-z);
    glScalef(back_h,back_w,back_l);
    cube(c[5]);
    glPopMatrix();

    //<--------------Upper Side backup---------------->
    x=table_h-side_h;
    y=leg_w+2*table_w+side_w;
    z=table_l-side_l;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(-x,y,-z);
        glScalef(side_h,side_w,side_l);
        cube(c[5]);
        glPopMatrix();
        x*=(-1);
    }
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    y=2*(table_w+leg_w);
    glPushMatrix();
    glTranslatef(0,y,0);
    room_laptop();
    glPopMatrix();
}


//<-------------------------------------------------------ROOM BED------------------------------------->
void room_bed()
{
    float   x,y,z,bed_h=1.5,bed_w=0.1,bed_l=3,
                  leg_h=0.1,leg_w=0.5,leg_l=0.1,
                  backup_h=bed_h,backup_w=bed_w,backup_l=0.1,
                  backup_leg_h=0.1,backup_leg_w=0.25,backup_leg_l=0.1,
                  balish_h=0.75,balish_w=0.12,balish_l=0.5;

    glPushMatrix();
    glTranslatef(0,leg_w,0);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[17]);
    ///<------------Bed Base----------->
    y=leg_w+bed_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(bed_h,bed_w,bed_l);
    cube(c[5]);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[16]);
    ///<-------------Leg Side---------------->
    x=bed_h-leg_h;
    z=bed_l-leg_l;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<2; j++)
        {
            glPushMatrix();
            glTranslatef(x,0,z);
            glScalef(leg_h,leg_w,leg_l);
            cube(c[5]);
            glPopMatrix();
            x*=(-1);
        }
        z*=(-1);
    }

    ///<-------------BackUp  Side---------------->
    y=leg_w+2*bed_w+backup_w+0.5;
    z=bed_l-backup_l;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(0,y,-z);
        glScalef(backup_h,backup_w,backup_l);
        cube(c[5]);
        glPopMatrix();
        z*=(-1);
    }
    ///<-------------BackUp Leg---------------->
    x=bed_h-backup_leg_h;
    y=leg_w+2*bed_w+backup_leg_w;
    z=bed_l-backup_l;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<2; j++)
        {
            glPushMatrix();
            glTranslatef(x,y,-z);
            glScalef(backup_leg_h,backup_leg_w,backup_leg_l);
            cube(c[5]);
            glPopMatrix();
            x*=(-1);
        }
        z*=(-1);
    }
    glDisable(GL_TEXTURE_2D);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[18]);
    ///<------------------------Balish----------------------->
    y= y=leg_w+2*bed_w+balish_w;
    z=bed_l-balish_l-bed_w;
    glPushMatrix();
    glTranslatef(0,y,-z);
    glScalef(balish_h,balish_w,balish_l);
    cube(c[5]);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);


    glPopMatrix();
}


//<-------------------------------------------------COMBINE CHAIR & TABLE & BED & FAN-------------------------------->
void room_elements()
{
    float x,y,z;
    glPushMatrix();
    x=-4.3;
    for(int i=0; i<2; i++)
    {
        z=-9;
        for(int j=0; j<3; j++)
        {
            glPushMatrix();
            glTranslatef(x,rm_w+0.2,z);
            //<---------------------------Room Chair--------------------------->
            glPushMatrix();
            glTranslatef(0,0,1.5);
            glRotatef(180,0,1,0);
            room_chair();
            glPopMatrix();

            //<---------------------------Room Table--------------------------->
            glPushMatrix();
            room_table();
            glPopMatrix();
            z+=8;
            glPopMatrix();
        }
        x*=(-1);
    }


    //<---------------------------Room All Beds --------------------------->

    x=-7.8;
    for(int j=0; j<2; j++)
    {
        z=-8;
        for(int i=0; i<3; i++)
        {
            glPushMatrix();
            glTranslatef(x,rm_w+0.2,z);
            room_bed();
            glPopMatrix();
            z+=8;
        }
        x*=(-1);
    }
    glPopMatrix();

    //<---------------------------Room Fans--------------------------->
    x=-4;
    y=9.5;
    z=-5;
    glPushMatrix();
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<2; j++)
        {
            glPushMatrix();
            glTranslatef(x,y,z);
            glRotatef(fan_rn,0,1,0);
            room_fan();
            glPopMatrix();
            x*=(-1);
        }
        z*=(-1);
    }
    glPopMatrix();
}

//<-------------------------------------------------------SINGLE ROOM----------------------------------------->
void single_room()
{
    float x,y,z;


    //<---------------------------Room BackGround--------------------------->
    glPushMatrix();
    glTranslatef(0,rm_w,0);
    room_background();
    glPopMatrix();

    //<---------------------------Room Elements--------------------------->
    glPushMatrix();
    room_elements();
    glPopMatrix();

}


//<-------------------------------------------------------ROOM CARROM BOARD------------------------------------->
void room_carrom()
{
    float   x,y,z,base_h=1.5,base_w=0.05,base_l=1.5,
                  leg_h=0.05,leg_w=0.8,leg_l=0.05,
                  carrom_h=1.5,carrom_w=0.01,carrom_l=1.5;

    glPushMatrix();
    glTranslatef(-3.5,leg_w+0.2,-8);

    ///<----------------------Carrom Base--------------------->
    y=base_w+leg_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(base_h,base_w,base_l);
    cube(c[10]);
    glPopMatrix();

    ///<-------------------Main Carrom Board---------------------->

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[26]);
    y=2*base_w+leg_w+carrom_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(carrom_h,carrom_w,carrom_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);


    ///<---------------------------Leg Side----------------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[14]);
    x=base_h-leg_h-0.1;
    z=base_l-leg_l-0.1;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<2; j++)
        {
            glPushMatrix();
            glTranslatef(x,0,z);
            glScalef(leg_h,leg_w,leg_l);
            cube(c[5]);
            glPopMatrix();
            x*=(-1);
        }
        z*=(-1);
    }
    glDisable(GL_TEXTURE_2D);

    glPopMatrix();

}


//<--------------------------------------------------------ROOM TABLE TENNIS--------------------------------------->
void room_table_tennis()
{
    float   x,y,z,base_h=1.5,base_w=0.05,base_l=2.5,
                  leg_h=0.05,leg_w=0.8,leg_l=0.05,
                  tt_h=1.5,tt_w=0.01,tt_l=2.5,
                  net_h=1.5,net_w=0.2,net_l=0.005;

    glPushMatrix();
    glTranslatef(-6.5,leg_w+0.2,0);

    ///<----------------------TT Base--------------------->
    y=base_w+leg_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(base_h,base_w,base_l);
    cube(c[10]);
    glPopMatrix();

    ///<-------------------Main TT Board---------------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[27]);
    y=2*base_w+leg_w+tt_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(tt_h,tt_w,tt_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    ///<-------------------TT Net---------------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[28]);
    y=2*base_w+leg_w+2*tt_w+net_w;
    z=net_l;
    glPushMatrix();
    glTranslatef(0,y,-z);
    glScalef(net_h,net_w,net_l);
    cube(c[5],5,1,0,1,0,0,5,0);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[28]);
    y=2*base_w+leg_w+2*tt_w+net_w;
    glPushMatrix();
    glTranslatef(0,y,z);
    glScalef(net_h,net_w,net_l);
    cube(c[5],0,1,0,0,5,0,5,1);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    ///<---------------------------Leg Side----------------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[14]);
    x=base_h-leg_h-0.1;
    z=base_l-leg_l-0.1;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<2; j++)
        {
            glPushMatrix();
            glTranslatef(x,0,z);
            glScalef(leg_h,leg_w,leg_l);
            cube(c[5]);
            glPopMatrix();
            x*=(-1);
        }
        z*=(-1);
    }
    glDisable(GL_TEXTURE_2D);

    glPopMatrix();

}


//<-------------------------------------------------------- COMBINE TV ROOM TV & MAT & CHAIRS-------------------------------------------->
void tv_mat_and_chairs()
{
    float x,y,z,tv_h=1.5,tv_w=1,tv_l=0.1,
                mat_h=5,mat_w=0.02,mat_l=1.5;

    glPushMatrix();
    //<------------------------TV------------------------------->
    y=4*tv_w;
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[29]);
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(tv_h,tv_w,tv_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<------------------------Mat------------------------------->
    y=mat_w;
    z=mat_l;
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[30]);
    glPushMatrix();
    glTranslatef(0,y,z);
    glScalef(mat_h,mat_w,mat_l);
    cube(c[5],0,0,1,0,1,1,0,1);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<---------------------Multiple Chairs---------------------->
    z=3*mat_l+1.2;
    for(int i=0; i <3; i++)
    {
        x=0;
        for(int j=0; j<4; j++)
        {
            glPushMatrix();
            glTranslatef(x,0,z);
            glRotatef(180,0,1,0);
            room_chair();
            glPopMatrix();

            glPushMatrix();
            glTranslatef(-x,0,z);
            glRotatef(180,0,1,0);
            room_chair();
            glPopMatrix();
            x+=2;
        }
        z+=1.5;
    }
    glPopMatrix();
}

//<-------------------------------------------------------TV ROOM  ALL ELEMENTS----------------------------------------->
void tv_room_elements()
{
    float x,y,z;

    //<----------------------------TV and MAT----------------------------->
    glPushMatrix();
    glTranslatef(9.5,0.2,0);
    glRotatef(-90,0,1,0);
    tv_mat_and_chairs();
    glPopMatrix();

    //<---------------------------TT--------------------------->
    glPushMatrix();
    room_table_tennis();
    glPopMatrix();

    //<---------------------------CARROM--------------------------->
    glPushMatrix();
    room_carrom();
    glPopMatrix();

    //<---------------------------Room Fans--------------------------->
    x=-4;
    y=9.3;
    z=-5;
    glPushMatrix();
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<2; j++)
        {
            glPushMatrix();
            glTranslatef(x,y,z);
            glRotatef(fan_rn,0,1,0);
            room_fan();
            glPopMatrix();
            x*=(-1);
        }
        z*=(-1);
    }
    glPopMatrix();
}

//<-------------------------------------------------TV ROOM------------------------------------------>
void tv_room()
{
    //<------------------------Room Background---------------------------->
    glPushMatrix();
    glTranslatef(0,rm_w,0);
    room_background();
    glPopMatrix();
    //<----------------------------TV Room Elements------------------------------>
    window_check=1;
    glPushMatrix();
    glTranslatef(0,rm_w,0);
    tv_room_elements();
    glPopMatrix();


}


//<----------------------------------------------------SOFA SET----------------------------------------------------->
void room_sofaset()
{
    float x,y,z,sofa_h=1.8,sofa_w=0.2,sofa_l=1,
                base_h=2,base_w=0.1,base_l=1,
                leg_h=0.05,leg_w=0.05,leg_l=0.05,
                side_h=0.1,side_w=0.2,side_l=1,
                upper_side_h=0.2,upper_side_w=0.2,upper_side_l=1,
                back_h=1.6,back_w=0.5,back_l=0.1,
                balish_h=1.6,balish_w=0.5,balish_l=0.01;

    glPushMatrix();
    glTranslatef(0,0.2,0);

    //<-------------------Leg of Sofa---------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[32]);
    x=base_h-leg_h;
    y=leg_w;
    z=base_l-leg_l;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<2; j++)
        {
            glPushMatrix();
            glTranslatef(x,y,z);
            glScalef(leg_h,leg_w,leg_l);
            cube(c[5],0,2,0,0,2,0,2,2);
            glPopMatrix();
            x*=(-1);
        }
        z*=(-1);
    }

    //<------------------Sofa Base----------------->
    y=2*leg_w+base_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(base_h,base_w,base_l);
    cube(c[5],0,2,0,0,2,0,2,2);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<----------------Main Sofa------------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[31]);
    y=2*leg_w+2*base_w+sofa_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(sofa_h,sofa_w,sofa_l);
    cube(c[5],0,5,0,0,5,0,5,5);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<--------------Side of Sofa------------------------>
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[32]);
    x=sofa_h+side_h;
    y=2*leg_w+2*base_w+side_w;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,y,0);
        glScalef(side_h,side_w,side_l);
        cube(c[5],0,2,0,0,2,0,2,2);
        glPopMatrix();
        x*=(-1);
    }

    //<--------------Upper Side of Sofa-------------------->
    x=base_h-upper_side_h;
    y=2*leg_w+2*base_w+2*sofa_w+upper_side_w;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,y,0);
        glScalef(upper_side_h,upper_side_w,upper_side_l);
        cube(c[5],0,2,0,0,2,0,2,2);
        glPopMatrix();
        x*=(-1);
    }

    glDisable(GL_TEXTURE_2D);
    //<--------------Back of Sofa-------------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[31]);
    y=2*leg_w+2*base_w+2*sofa_w+back_w;
    z=base_l-back_l;
    glPushMatrix();
    glTranslatef(0,y,-z);
    glScalef(back_h,back_w,back_l);
    cube(c[5],0,5,0,0,5,0,5,5);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<---------------Balish of Sofa------------------>
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[31]);
    y=2*leg_w+2*base_w+2*sofa_w+balish_w;
    z=base_l-(2*back_l+balish_l);
    glPushMatrix();
    glTranslatef(0,y,-z);
    glScalef(balish_h,balish_w,balish_l);
    cube(c[5],0,5,0,0,5,0,5,5);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glPopMatrix();
}

//<--------------------------------------------------ROOM ALMARI----------------------------------------->
void room_almari()
{
    float x,y,z,almari_front_h=1.5,almari_front_w=3,almari_front_l=0.01,
                almari_side_h=0.01,almari_side_w=3,almari_side_l=0.5,
                almari_top_h=1.5,almari_top_w=0.01,almari_top_l=0.5;

    glPushMatrix();
    glTranslatef(0,0.2,0);
    //<-----------------------Front & Back----------------------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[34]);
    y=almari_front_w;
    z=almari_side_l;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(0,y,z);
        glScalef(almari_front_h,almari_front_w,almari_front_l);
        cube(c[5]);
        glPopMatrix();
        z*=(-1);
    }
    glDisable(GL_TEXTURE_2D);

    //<-----------------------Left & Right Side----------------------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[35]);
    x=almari_front_h;
    y=almari_side_w;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,y,0);
        glScalef(almari_side_h,almari_side_w,almari_side_l);
        cube(c[5]);
        glPopMatrix();
        x*=(-1);
    }
    //<---------------------------Top----------------------------------->
    y=2*almari_front_w+almari_top_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(almari_top_h,almari_top_w,almari_top_l);
    cube(c[5]);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);

    glPopMatrix();
}

//<--------------------------------------------COMBINE SOFA & ALMARI & FANS---------------------------------------->
void guest_room_elements()
{
    float x,y,z,rotatn;

    //<-------------------------Room Sofaset--------------------------->
    glPushMatrix();
    glTranslatef(0,0,-4);
    room_sofaset();
    glPopMatrix();

    glPushMatrix();
    x=4;
    rotatn=-90;
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(x,0,0);
        glRotatef(rotatn,0,1,0);
        room_sofaset();
        glPopMatrix();
        x*=(-1);
        rotatn*=(-1);
    }
    glPopMatrix();

    //<--------------------------Room Almari--------------------------->
    glPushMatrix();
    glTranslatef(-7,0,-9);
    glRotatef(15,0,1,0);
    room_almari();
    glPopMatrix();

    //<---------------------------Room Fans--------------------------->
    x=-4;
    y=9.2;
    z=-5;
    glPushMatrix();
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<2; j++)
        {
            glPushMatrix();
            glTranslatef(x,y,z);
            glRotatef(fan_rn,0,1,0);
            room_fan();
            glPopMatrix();
            x*=(-1);
        }
        z*=(-1);
    }
    glPopMatrix();
}

//<-------------------------------------------------------GUEST ROOM----------------------------------------->
void guest_room()
{
    float x,y,z;

    //<--------------------------- Guest Room BackGround--------------------------->
    glPushMatrix();
    room_background();
    glPopMatrix();

    //<----------------------------Guest Room Elements-------------------->
    glPushMatrix();
    guest_room_elements();
    glPopMatrix();

}



//<-------------------------------------------------------SINGLE BLOCK----------------------------------------->
void single_block()
{
    float x,y,z;
    glPushMatrix();
    //<---------------------------Single Room--------------------------->
    glPushMatrix();
    single_room();
    glPopMatrix();
    x=20.2;
    //<---------------------------Block --------------------------->
    for(int i=0; i<2; i++)
    {
        glPushMatrix();
        glTranslatef(-x,0,0);
        if(i==1 && room_check==0)
        {
            window_check=0;
            tv_room();
        }
        else
        {
            single_room();
        }
        glPopMatrix();

        glPushMatrix();
        glTranslatef(x,0,0);
        if(i==1 && room_check==0)
        {
            window_check=0;
            glPushMatrix();
            glTranslatef(0,rm_w,0);
            guest_room();
            glPopMatrix();
        }
        else
        {
            single_room();
        }
        glPopMatrix();
        x+=20.2;

    }
    room_check=1;
    window_check=1;
    glPopMatrix();

}

///<--------------------------------------------------ROOM BARANDA--------------------------------->


void room_baranda_type1()
{
    float x,y,z,baranda_h=50.5,baranda_w=0.1,baranda_l=16,
                front_h=50.5,front_w=1.25,front_l=0.2,
                side_h=0.2,side_w=1.25,side_l=4;

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[24]);

    glPushMatrix();
    glTranslatef(0,0,4);
    //<--------------------------Baranda----------------------->
    glPushMatrix();
    glTranslatef(0,baranda_w,0);
    glScalef(baranda_h,baranda_w,baranda_l);
    cube(c[5]);
    glPopMatrix();

    //<--------------------------Front----------------------->
    y=front_w;
    z=baranda_l-front_l;
    glPushMatrix();
    glTranslatef(0,y,z);
    glScalef(front_h,front_w,front_l);
    cube(c[5]);
    glPopMatrix();

    //<--------------------------Side----------------------->
    x=baranda_h-side_h;
    y=side_w;
    z=baranda_l-side_l;
    glPushMatrix();
    if(side_check==0)
        glTranslatef(-x,y,z);
    else
        glTranslatef(x,y,z);
    glScalef(side_h,side_w,side_l);
    cube(c[5]);
    glPopMatrix();

    glPopMatrix();

    glDisable(GL_TEXTURE_2D);
}


///<------------------------------------GROUND FLOOR BARANDA--------------------------------------->

void grill_helper(float h,float w,float l)
{
    float x=0,n,increment;
    if(main_floor_baranda_check==0)
        n=51,increment=0.91;
    else
        n=57,increment=0.9;
    for(int k=0; k<n; k++)
    {
        glPushMatrix();
        glTranslatef(x,0,0);
        glScalef(h,w,l);
        cube(c[5]);
        glPopMatrix();

        glPushMatrix();
        glTranslatef(-x,0,0);
        glScalef(h,w,l);
        cube(c[5]);
        glPopMatrix();
        x+=increment;
    }
}

void room_baranda_grill(float h,float w,float l)
{
    float x,y=0,z,support_h=0.1,support_w=0.35,support_l=0.1;
    glPushMatrix();


    for(int i=0; i<7; i++)
    {
        glPushMatrix();
        glTranslatef(0,y,0);
        glScalef(h,w,l);
        cube(c[5]);
        glPopMatrix();

        glPushMatrix();
        glTranslatef(0,y-support_w,0);
        grill_helper(support_h,support_w,support_l);
        glPopMatrix();

        glPushMatrix();
        glTranslatef(0,-y,0);
        glScalef(h,w,l);
        cube(c[5]);
        glPopMatrix();

        glPushMatrix();
        glTranslatef(0,-(y-support_w),0);
        grill_helper(support_h,support_w,support_l);
        glPopMatrix();

        y+=0.7;
    }
    glPopMatrix();
}

void room_baranda_ground_floor()
{
    float x,y,z,baranda_h=50.5,baranda_w=0.1,baranda_l=16,
                front_h=50.5,front_w=1.25,front_l=0.2,
                helper_h=50.5,helper_w=0.1,helper_l=0.1,
                side_h=0.2,side_w=1.25,side_l=4;

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[24]);

    glPushMatrix();
    glTranslatef(0,0,4);
    //<--------------------------Baranda----------------------->
    glPushMatrix();
    glTranslatef(0,baranda_w,0);
    glScalef(baranda_h,baranda_w,baranda_l);
    cube(c[5]);
    glPopMatrix();

    //<--------------------------Front----------------------->
    x=0;
    if(main_floor_baranda_check==0)
    {
        front_h=helper_h=45.5;
        x=baranda_h-front_h;
    }
    y=front_w;
    z=baranda_l-front_l;
    glPushMatrix();
    glTranslatef(-x,y,z);
    glScalef(front_h,front_w,front_l);
    cube(c[5]);
    glPopMatrix();

    //<--------------------------Side----------------------->
    x=baranda_h-side_h;
    y=side_w;
    z=baranda_l-side_l;
    glPushMatrix();
    if(side_check==0)
        glTranslatef(-x,y,z);
    else
        glTranslatef(x,y,z);
    glScalef(side_h,side_w,side_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<--------------------------Grill-------------------------------------->
    if(main_floor_baranda_check==0)
        x=baranda_h-front_h;
    else
        x=0;
    y=front_w+5.8;
    z=baranda_l-front_l;
    glPushMatrix();
    glTranslatef(-x,y,z);
    room_baranda_grill(helper_h,helper_w,helper_l);
    glPopMatrix();

    glPopMatrix();

}


///<------------------------------------SINGLE BLOCK & BARANDA------------------------------------------>

void single_block_baranda()
{
    glPushMatrix();
    glTranslatef(0,rm_w/2,0);
    single_block();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,rm_w,0);
    if(switch_baranda==0)
        room_baranda_ground_floor();
    else
        room_baranda_type1();
    glPopMatrix();
}


///<---------------------------------------FULL ONE PART OF THE HALL------------------------------------------------------->
void full_one_part_of_hall()
{
    float y=0;
    glPushMatrix();
    for(int i=0; i<3; i++)
    {
        //<---------------------------Full Block--------------------------->
        glPushMatrix();
        glTranslatef(0,y,0);
        single_block_baranda();
        glPopMatrix();
        y+=10.2;
        switch_baranda++;
    }
    glPushMatrix();
    glTranslatef(0,y+0.4,4);
    glScalef(51,0.5,16);
    cube(c[5]);
    glPopMatrix();
    switch_baranda=0;
    glPopMatrix();
}



///<-------------------------------------------------STAIRCASE------------------------------------------------------------------------>

//<------------------------------------------------HALL STAIRCASE-------------------------------------->

void staircase_background()
{
    float x,y,z,background_h=5,background_w=0.1,background_l=12;

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[20]);
    //<----------------Right---------------->
    y=background_h+2*background_w;
    x=background_h-background_w;
    glPushMatrix();
    glTranslatef(x,y,0);
    glScalef(background_w,background_h,background_l);
    cube(c[5]);
    glPopMatrix();

    //<----------------Back---------------->
    y=background_h+2*background_w;
    z=background_l-background_w;
    glPushMatrix();
    glTranslatef(0,y,-z);
    glScalef(background_h,background_h,background_w);
    cube(c[5]);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);
}

void staircase()
{
    float x,y,z,stair_h=1,stair_w=0.12,stair_l=1.5,stair_z,
                front_h=0.01,front_w=0.12,front_l=1.5,
                stand_h=0.05,stand_w=0.4,stand_l=0.05,
                helper_h=3,helper_w=0.07,helper_l=0.07;


    //<-------------------Lower Stair-------------------->
    x=0,y=stair_w;
    for(int i=0; i<10; i++)
    {
        glPushMatrix();
        glTranslatef(x,y,0);

        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D,texture_id[14]);
        glPushMatrix();
        glScalef(stair_h,stair_w,stair_l);
        cube2(c[5]);
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);

        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D,texture_id[25]);
        glPushMatrix();
        glTranslatef(stair_h+front_h,0,0);
        glScalef(front_h,front_w,front_l);
        cube(c[7]);
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);

        glPushMatrix();
        glTranslatef(stair_h-stand_h,stand_w,stair_l-stand_l);
        glScalef(stand_h,stand_w,stand_l);
        cube(c[7]);
        glPopMatrix();
        y+=2*stair_w;
        x-=0.5;
        glPopMatrix();
    }
    //<----------------------helper-------------------------->
    y=2*stand_w-helper_w;
    z=stair_l-helper_l;
    glPushMatrix();
    glTranslatef(helper_h/2,y,z);
    glRotatef(-25,0,0,1);
    glTranslatef(-helper_h,0,0);
    glScalef(helper_h,helper_w,helper_l);
    cube(c[7]);
    glPopMatrix();

}

//<--------------------------------------------HALL STAIRCASE BLOCK-------------------------------->

void block_staircase()
{

    glPushMatrix();
    //<---------------------------Main StairCase--------------------------->
    glPushMatrix();
    glTranslatef(0.9,rm_w,0);
    glRotatef(-90,0,1,0);
    staircase();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-1.6,rm_w+2.5,-7);
    glScalef(4,0.1,3);
    cube(c[7]);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-4.1,rm_w+2.5,-3);
    glRotatef(90,0,1,0);
    staircase();
    glPopMatrix();

    glPopMatrix();
}

///<----------------------------------------------STAIRCASE FRONT------------------------------------------------->

void staircase_front()
{
    float x,y,z,front_h=4.8,front_w=1.25,front_l=0.2;

    if(stair_check>0)
    {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D,texture_id[24]);
        glPushMatrix();
        glTranslatef(-1.2,rm_w+0.1,3.4);
        glScalef(4.8,0.1,10.5);
        cube(c[5]);
        glPopMatrix();

        glPushMatrix();
        glTranslatef(-1.3,rm_w+front_w,13.7);
        glScalef(front_h,front_w,front_l);
        cube(c[5]);
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);
    }
}

void combine_staircase_front()
{
    //<---------------------------Staircase Front Floor--------------------------------->
    glPushMatrix();
    glTranslatef(56.5,0,11.1);
    staircase_front();
    glPopMatrix();
}
///<--------------------------------------------COMBINE ALL FRONT OF THE STAIRCASE-------------------------------->
void full_front_staircase()
{
    float y=0;
    stair_check=0;
    for(int i=0; i<3; i++)
    {
        glPushMatrix();
        glTranslatef(0,y,0);
        combine_staircase_front();
        glPopMatrix();
        y+=10.2;
        stair_check++;
    }
}

///<---------------------------------------------STAIRCASE COMBINE---------------------------------->

void combine_background_and_staircase()
{

    glPushMatrix();
    glTranslatef(56.2,0,11.1);

    //<----------------------------Background--------------------------------->
    glPushMatrix();
    glTranslatef(-1.3,rm_w,-6);
    staircase_background();
    glPopMatrix();


    //<---------------------------Staircase--------------------------------->
    if(remove_upper_stairblock<=2)
    {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D,texture_id[24]);
        glPushMatrix();
        glTranslatef(0,0,-8);
        glScalef(1,2.01,1);
        block_staircase();
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);
    }
    glPopMatrix();
}


void full_loop_of_staircase()
{
    float y=0;
    remove_upper_stairblock=1;
    for(int i=0; i<3; i++)
    {
        glPushMatrix();
        glTranslatef(0,y,0);
        combine_background_and_staircase();
        glPopMatrix();
        y+=10;
        remove_upper_stairblock++;
    }
}

void complete_staircase()
{
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[24]);
    glPushMatrix();
    glTranslatef(55.2,rm_w+0.05,9);
    glScalef(4.7,0.2,16);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glPushMatrix();
    glTranslatef(55.2,rm_w+30.6,9);
    glScalef(4.7,0.5,16);
    cube(c[5]);
    glPopMatrix();

    glPushMatrix();
    full_loop_of_staircase();
    glPopMatrix();

    glPushMatrix();
    full_front_staircase();
    glPopMatrix();
}


void main_block_of_hall()
{
    //<---------------------------Combine Hall Staircase------------------------->
    glPushMatrix();
    full_one_part_of_hall();
    glPopMatrix();

    //<---------------------------Combine Hall Staircase------------------------->
    glPushMatrix();
    glTranslatef(0,0,-5);
    complete_staircase();
    glPopMatrix();
}

///<------------------------------------------------------------------------------ALL THREE BLOCK COMBINATOR----------------------------------------->


void block_combinator()
{
    float x,y,z,h=10,w=0.1,l=63,
                side_h=10,side_w=1.25,side_l=0.2,
                left_h=0.2,left_w=1.25,left_l=21.5;

    x=70;
    y=0.1;
    z=63;

    for(int i=0; i<3; i++)
    {
        //<---------------------Base----------------------------------->
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D,texture_id[24]);
        glPushMatrix();
        glTranslatef(70,y,0);
        glScalef(h,w,l);
        cube(c[5]);
        glPopMatrix();

        //<---------------------Front & Back ----------------------------------->
        glPushMatrix();
        glTranslatef(x,y+side_w,z);
        glScalef(side_h,side_w,side_l);
        cube(c[5]);
        glPopMatrix();

        glPushMatrix();
        glTranslatef(x,y+side_w,-z);
        glScalef(side_h,side_w,side_l);
        cube(c[5]);
        glPopMatrix();

        //<--------------------- Front left & Right ----------------------------------->
        glPushMatrix();
        glTranslatef(60,y+left_w,41.5);
        glScalef(left_h,left_w,left_l);
        cube(c[5]);
        glPopMatrix();


        left_l=31;
        glPushMatrix();
        glTranslatef(80,y+left_w,0);
        glScalef(left_h,left_w,left_l);
        cube(c[5]);
        glPopMatrix();

        //<--------------------------Back Left Part Only----------------------------------->
        left_l=25.5;
        glPushMatrix();
        glTranslatef(60,y+left_w,-37.5);
        glScalef(left_h,left_w,left_l);
        cube(c[5]);
        glPopMatrix();


        glDisable(GL_TEXTURE_2D);

        y+=10.1;
        left_l=21.5;
    }


}

void khaja_gate()
{
    float x,y,z,
    left_h=0.5,left_w=5.2,left_l=0.5,
    upper_h=10,upper_w=1.25,upper_l=0.5;


    //<-------------Left & Right    ------------------>
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[20]);
    x=9.5;
    y=left_w+rm_w;
    for(int i=0;i<2;i++)
    {
        glPushMatrix();
        glTranslatef(x,y,0);
        glScalef(left_h,left_w,left_l);
        cube(c[5]);
        glPopMatrix();
        x*=(-1);
    }
    glDisable(GL_TEXTURE_2D);
    //<----------------Upper------------------------------>
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[38]);
    y=2*left_w+rm_w+upper_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(upper_h,upper_w,upper_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}


///<----------------------------------KHAJA FRONT ROAD & CUBE--------------------------------->

void khaja_front_cube()
{
    float x,y,z,base_h=4,base_w=0.5,base_l=5,
    back_h=3,back_w=5,back_l=0.5,
    img_h=3,img_w=5,img_l=0.1;

    //<------------------base part----------------->
     glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[4]);
    y=base_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(base_h,base_w,base_l);
    cube(c[5]);
    glPopMatrix();

     y=2*base_w+back_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(back_h,back_w,back_l);
    cube(c[5]);
    glPopMatrix();

     glDisable(GL_TEXTURE_2D);

     glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[13]);
    z=back_l+img_l;
    glPushMatrix();
    glTranslatef(0,y,z);
    glScalef(img_h,img_w,img_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}


void khaja_front_road()
{
    float x,y,z,front_h=10,front_w=0.3,front_l=5,
    road_h=10,road_w=0.2,road_l=45;

    glPushMatrix();
    //<------------------front------------------------>
      glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[39]);
    glPushMatrix();
    glTranslatef(49.8,front_w,25);
    glScalef(front_h,front_w,front_l);
      cube(c[5],0,0,10,0,10,3,0,3);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<------------------road------------------------>
     glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[2]);
    glPushMatrix();
    glTranslatef(49.8,front_w,75);
    glScalef(road_h,road_w,road_l);
      cube(c[5],0,0,1,0,1,10,0,10);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glPopMatrix();

}

void khaja_front_combine_road_and_cube()
{
    //<---------------------------Image in Picture--------------------------->
     glPushMatrix();
     glTranslatef(65,rm_w,100);
     glRotatef(-75,0,1,0);
     khaja_front_cube();
     glPopMatrix();

     //<-----------------------Front Road--------------------------------------->
     glPushMatrix();
     glTranslatef(0,rm_w,0);
     khaja_front_road();
     glPopMatrix();
}

///<----------------------------------KHAJA FRONT ROAD & CUBE--------------------------------->


///<----------------------------------KHAN JAHAN ALI HALL--------------------------------->
void khan_jahan_ali_hall()
{
    //<------------------------Hall All Block---------------------------->
    side_check=0,main_floor_baranda_check=0;
    glPushMatrix();
    main_block_of_hall();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,rm_w,0);
    block_combinator();
    glPopMatrix();

    side_check=1,main_floor_baranda_check=1;
    glPushMatrix();
    glTranslatef(130.5,0,-51);
    full_one_part_of_hall();
    glPopMatrix();

    side_check=0;
    glPushMatrix();
    glTranslatef(130.5,0,51);
    glRotatef(180,0,1,0);
    full_one_part_of_hall();
    glPopMatrix();

     //<------------------------Block of Staircase---------------------------->
     glPushMatrix();
    glTranslatef(0,-rm_w/2,0);
    glScalef(1,2,1);
    block_staircase();
    glPopMatrix();

    //<---------------------------Khaja Gate------------------------------->
    glPushMatrix();
    glTranslatef(49.5,0,21);
    khaja_gate();
    glPopMatrix();

    //<-------------Front Road and Cube------------------->
    glPushMatrix();
    khaja_front_combine_road_and_cube();
    glPopMatrix();

}




///<----------------------------------------ALL TYPES OF TREE------------------------------------------------------------>
void single_leaf_creation()
{
    float x,y,z,leaf_h=0.04,leaf_w=0.4,leaf_l=0.04;

    float rotatn=45;
    y=leaf_w;
    //<----------------Rotate with z axis----------------->
    if(flagLeaf1)
    {
        for(int i=0; i<2; i++)
        {
            glEnable(GL_TEXTURE_2D);
            glBindTexture(GL_TEXTURE_2D,texture_id[11]);
            glPushMatrix();
            glRotatef(-rotatn,0,0,1);
            glTranslatef(0,y,0);
            glScalef(leaf_h,leaf_w,leaf_l);
            cube(c[5]);
            glPopMatrix();
            rotatn*=(-1);
            glDisable(GL_TEXTURE_2D);
        }
    }
    //<----------------Rotate with x axis----------------->
    if(flagLeaf2)
    {
        for(int i=0; i<2; i++)
        {
            glEnable(GL_TEXTURE_2D);
            glBindTexture(GL_TEXTURE_2D,texture_id[11]);
            glPushMatrix();
            glRotatef(-rotatn,1,0,0);
            glTranslatef(0,y,0);
            glScalef(leaf_h,leaf_w,leaf_l);
            cube(c[5]);
            glPopMatrix();
            rotatn*=(-1);
            glDisable(GL_TEXTURE_2D);
        }
    }
}

void leaf_leg_combine()
{
    float x,y,z,leg_h=0.05,leg_w=1.5,leg_l=0.05;

    //<------------------leg Part----------------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[10]);
    y=leg_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(leg_h,leg_w,leg_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<--------------------leaf Part----------------------->
    y=0.6;

    while(y<=2*leg_w)
    {
        glPushMatrix();
        glTranslatef(0,y,0);
        single_leaf_creation();
        glPopMatrix();
        y+=0.3;
    }

}

void single_node_creation()
{
    float x,y,z,rotatn=45;
    //<----------------Rotate with z axis----------------->
    if(flagNode1)
    {
        for(int i=0; i<2; i++)
        {
            glPushMatrix();
            glRotatef(-rotatn,0,0,1);
            leaf_leg_combine();
            glPopMatrix();
            rotatn*=(-1);
        }
    }
    //<----------------Rotate with x axis----------------->
    if(flagNode2)
    {
        for(int i=0; i<2; i++)
        {
            glPushMatrix();
            glRotatef(-rotatn,1,0,0);
            leaf_leg_combine();
            glPopMatrix();
            rotatn*=(-1);
        }
    }
}

void complex_node_creation()
{
    float x,y,z,node_h=0.15,node_w=3,node_l=0.15;

    //<------------------Node Part----------------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[10]);
    y=node_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(node_h,node_w,node_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<--------------------Leg Part----------------------->
    y=2;

    while(y<=2*node_w)
    {
        glPushMatrix();
        glTranslatef(0,y,0);
        single_node_creation();
        glPopMatrix();
        y+=1.8;
    }

}

void tree_type1()
{
    flagLeaf1=flagNode1=true;
    flagLeaf2=flagNode2=true;
    glPushMatrix();
    complex_node_creation();
    glPopMatrix();
}

void tree_type2()
{
    float x,y,z,base_h=0.5,base_w=2.5,base_l=0.5,rotatn=45;

    //<-------------------------Base---------------------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[10]);
    y=base_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(base_h,base_w,base_l);
    cube(c[5],0,2,0,0,1,0,1,2);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    y=1.2*base_w;
    for(int j=0; j<2; j++)
    {
        flagLeaf1=flagNode1=false;
        flagLeaf2=flagNode2=true;
        for(int i=0; i<2; i++)
        {
            glPushMatrix();
            glTranslatef(0,y,0);
            glRotatef(-rotatn,0,0,1);
            complex_node_creation();
            glPopMatrix();
            rotatn*=(-1);
        }
        flagLeaf2=flagNode2=false;
        flagLeaf1=flagNode1=true;
        for(int i=0; i<2; i++)
        {
            glPushMatrix();
            glTranslatef(0,y,0);
            glRotatef(-rotatn,1,0,0);
            complex_node_creation();
            glPopMatrix();
            rotatn*=(-1);
        }
        y=1.7*base_w;
        rotatn=35;
    }


}

void loop_of_tree()
{
     float z=0,rotatn=30;
    for(int i=0;i<2;i++)
    {
         glPushMatrix();
        glTranslatef(0,rm_w,z);
        glRotatef(rotatn,0,1,0);
        tree_type2();
        glPopMatrix();
        z-=10;
        glPushMatrix();
        glTranslatef(0,rm_w,z);
        tree_type1();
        glPopMatrix();
        z-=10;
        rotatn+=5;
    }
}


void trees()
{

}


///<---------------------------------AUDITORIUM--------------------------------------->


//<-----------------------------------------------------AUDI BACKGROUND-------------------------------->
void audi_door()
{
    float x,y,z,door_h=3,door_w=6.7,door_l=0.2,
    side_h=13.5,side_w=7,side_l=0.2;

    //<------------------Side---------------->
    x=30-side_h;
    y=side_w;
    z=14;
    for(int i=0;i<2;i++)
    {
         glEnable(GL_TEXTURE_2D);
     glBindTexture(GL_TEXTURE_2D,texture_id[46]);
        glPushMatrix();
        glTranslatef(x,y,z);
        glScalef(side_h,side_w,side_l);
       cube(c[5],0,10,0,0,10,0,10,10);
        glPopMatrix();
        x*=(-1);
        glDisable(GL_TEXTURE_2D);
    }

    //<-----------------Door-------------------->
    x=door_h;
    y=door_w+2*door_l;
     glEnable(GL_TEXTURE_2D);
     glBindTexture(GL_TEXTURE_2D,texture_id[22]);
    glPushMatrix();
    glTranslatef(x,y,z);
    glRotatef(rn,0,1,0);
    glTranslatef(-door_h,0,0);
    glScalef(door_h,door_w,door_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

}
void audi_background()
{
    float x,y,z,room_h=30,room_w=0.2,room_l=14,room_y=7;

    //<-----------------Floor----------------->
     glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[19]);
    y=room_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(room_h+room_w,room_w,room_l);
    cube(c[5],0,0,10,0,10,10,0,10);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);


    //<-----------------Ceil----------------->
     glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[23]);
    y=2*room_y;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(room_h+room_w,room_w,room_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);


    //<-----------------Left Side---------------->
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[46]);
    x=room_h;
    y=room_y;
    glPushMatrix();
    glTranslatef(-x,y,0);
    glScalef(room_w,room_y,room_l);
    cube(c[5],0,10,0,0,10,0,10,10);
    glPopMatrix();

    //<-----------------Right Side---------------->
    glPushMatrix();
    glTranslatef(x,y,0);
    glScalef(room_w,room_y,room_l);
    cube(c[5],0,10,0,0,10,0,10,10);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);


    //<-------------------Back   Part-------------------->
    z=room_l-room_w;
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[46]);
    glPushMatrix();
    glTranslatef(0,y,-z);
    glScalef(room_h,room_y,room_w);
    cube(c[5],0,10,0,0,10,0,10,10);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<----------------------------Front Part----------------->
    glPushMatrix();
    audi_door();
    glPopMatrix();
}

//<---------------------------------------------AUDI ELEMENTS------------------------>
void audi_elements()
{
    float x,y,z,bannar_h=6,bannar_w=1.5,bannar_l=0.02,
    base_h=12,base_w=2,base_l=3,
    stair_h=2,stair_w=0.4,stair_l=2;

    glPushMatrix();
    //<------------------------Front Bannar------------------------------->
    y=6*bannar_w;
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[45]);
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(bannar_h,bannar_w,bannar_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

     //<------------------------Base------------------------------->
    y=base_w;
    z=base_l;
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[19]);
    glPushMatrix();
    glTranslatef(0,y,z);
    glScalef(base_h,base_w,base_l);
    cube(c[5],0,0,10,0,10,10,0,10);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    //<---------------------- Left Stair----------------------------->
     glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[25]);
    x=base_h-stair_h;
    y=stair_w;
    z=2*base_l+stair_l;
    glPushMatrix();
    glTranslatef(-x,y,z);
    float ty=0,tz=0,sz=stair_l;
     for(int i=0;i<5;i++)
    {

        glPushMatrix();
        glTranslatef(0,ty,tz);
        glScalef(stair_h,stair_w,sz);
        cube(c[5]);
        glPopMatrix();
        sz-=0.4;
        tz-=0.4;
        ty+=(2*stair_w);
    }
    glPopMatrix();

    //<---------------------- Right Stair----------------------------->
    glPushMatrix();
    glTranslatef(x,y,z);
     ty=0,tz=0,sz=stair_l;
     for(int i=0;i<5;i++)
    {
        glPushMatrix();
        glTranslatef(0,ty,tz);
        glScalef(stair_h,stair_w,sz);
        cube(c[5]);
        glPopMatrix();
        sz-=0.4;
        tz-=0.4;
        ty+=(2*stair_w);
    }
    glPopMatrix();
   glDisable(GL_TEXTURE_2D);
    //<--------------------FrontLine Sofa Set----------------------->
    x=0;
    z=5*base_l;
    for(int i=0;i<3;i++)
    {
        glPushMatrix();
        glTranslatef(x,0,z);
        glRotatef(180,0,1,0);
        room_sofaset();
        glPopMatrix();

        glPushMatrix();
        glTranslatef(-x,0,z);
        glRotatef(180,0,1,0);
        room_sofaset();
        glPopMatrix();
        x+=4.1;
    }


    //<---------------------Multiple Chairs---------------------->
    z=5*base_l+stair_l;
     for(int i=0;i <14;i++)
     {
         x=0;
         for(int j=0;j<5;j++)
         {
             glPushMatrix();
             glTranslatef(x,0,z);
             glRotatef(180,0,1,0);
             room_chair();
             glPopMatrix();

             glPushMatrix();
             glTranslatef(-x,0,z);
             glRotatef(180,0,1,0);
             room_chair();
             glPopMatrix();
             x+=2;
         }
         z+=3;
     }
     glPopMatrix();
}

void audi_tubelight()
{
    float light_h=3,light_w=0.2,light_l=0.1;
    //<---------TubeLight------------->
    glPushMatrix();
    glScalef(light_h,light_w,light_l);
    cube(c[5],true);
    glPopMatrix();
}
//<----------------------------------------------AUDITORIUM MAIN FUNCTION-------------------------------------->

void auditorium()
{
    float x,y,z,upper_h=30,upper_w=2,upper_l=15,
    board_h=5,board_w=2,board_l=0.05;

    //<-------------------------Upper Side of Audi----------------------->
    y=14.2+upper_w;
    glPushMatrix();
    glTranslatef(0,y,1);
    glScalef(upper_h,upper_w,upper_l);
    cube(c[7]);
    glPopMatrix();

    //<-------------------------Upper Side Board of Audi----------------------->
    y=14.2+upper_w;
    z=upper_l+2*board_l+1;
    glEnable(GL_TEXTURE_2D);
     glBindTexture(GL_TEXTURE_2D,texture_id[44]);
    glPushMatrix();
    glTranslatef(0,y,z);
    glScalef(board_h,board_w,board_l);
    cube(c[5]);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
     //<-------------------------Audi BackGround----------------------->
    glPushMatrix();
    audi_background();
    glPopMatrix();

    //<-------------------------Audi Elements----------------------->
    glPushMatrix();
    glTranslatef(29.7,0.2,-1);
    glRotatef(-90,0,1,0);
    audi_elements();
    glPopMatrix();

    //<-----------------------------Audi TubeLight------------------------>
    glPushMatrix();
    glTranslatef(6,12,-13.5);
    audi_tubelight();
    light2(0,0,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,12,13.5);
    audi_tubelight();
    light3(0,0,-1);
    glPopMatrix();

}



///<--------------------------------------------CAR--------------------------------------------------->

//<-------------------------------------------------------CAR DESIGN---------------------------------------->

void circle_of_wheel(float wheel_h,float wheel_w,float wheel_l)
{
    glPushMatrix();
    for(float i=0; i<360; i+=0.04)
    {
        glPushMatrix();
        glRotatef(i,0,0,1);
        glTranslatef(wheel_h,0,0);
        glScalef(wheel_h,wheel_w,wheel_l);
        cube(c[10]);
        glPopMatrix();
    }
    glPopMatrix();
}

void wheel()
{
     glPushMatrix();
     glRotatef(90,0,1,0);
    glScalef(1,1,1);
    wheelBezier();
    glPopMatrix();
}

void car()
{
    float x,y,z,
          middle_h=2.5,middle_w=0.6,middle_l=1.5,
          upper_h=1.8,upper_w=0.5,upper_l=1.5,
          front_h=1,front_w=0.5,front_l=1.5,
          wheel_h=0.3,wheel_w=0.01,wheel_l=0.1,
          window_h=1.6,window_w=0.5,window_l=0.01;

    glPushMatrix();
    glTranslatef(0,wheel_h,0);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id[7]);
    //<------------Middle Part--------------->
    y=middle_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(middle_h,middle_w,middle_l);
    cube(c[5],0,0,5,0,5,5,0,5);
    glPopMatrix();

    //<------------Front Part--------------->
    x=middle_h+front_h;
    y=front_w;
    glPushMatrix();
    glTranslatef(-x,y,0);
    glScalef(front_h,front_w,front_l);
    cube5(c[5],0,0,5,0,5,5,0,5);
    glPopMatrix();

    //<------------Back Part--------------->
    x=middle_h+front_h;
    y=front_w;
    glPushMatrix();
    glTranslatef(x,y,0);
    glRotatef(-180,0,1,0);
    glScalef(front_h,front_w,front_l);
    cube5(c[5],0,0,5,0,5,5,0,5);
    glPopMatrix();

    //<------------Upper Part--------------->
    y=2*middle_w+upper_w;
    glPushMatrix();
    glTranslatef(0,y,0);
    glScalef(upper_h,upper_w,upper_l);
    cube6(c[5],0,0,5,0,5,5,0,5);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

     //<------------------Left & Rigth Side Window------------------------->
     y=2*middle_w+upper_w;
     z=upper_l+window_l;
     for(int i=0;i<2;i++)
     {
         glPushMatrix();
         glTranslatef(0,y,z);
         glScalef(window_h,window_w,window_l);
         cube6(c[5]);
         glPopMatrix();
         z*=(-1);
     }

     //<------------------Front & Back Side Window------------------------->
     window_h=0.01;
     window_w=0.35;
     window_l=1.2;
     y=2*middle_w+upper_w;
     x=middle_h;
     float rotatn=35;
     for(int i=0;i<2;i++)
     {
         glPushMatrix();
         glTranslatef(x,y,0);
         glRotatef(rotatn,0,0,1);
         glTranslatef(0,window_w,0);
         glScalef(window_h,window_w,window_l);
         cube5(c[5]);
         glPopMatrix();
         x*=(-1);
         rotatn*=(-1);
     }



    //<--------------------Wheel------------------------->
    y=2*wheel_h;
    z=upper_l+wheel_l;
    for(int i=0; i<2; i++)
    {
        x=3;
        for(int j=0; j<2; j++)
        {
            glPushMatrix();
            glTranslatef(x,y,z);
            wheel();
            glPopMatrix();
            x*=(-1);
        }
        z*=(-1);
    }


}


///<---------------------------MAIN PROJECT--------------------------->

void main_project()
{

    //<---------------------------Main BackGround--------------------------->
    glPushMatrix();
    main_background(rm_h,rm_w,rm_l);
    glPopMatrix();

    //<---------------------------Main Road--------------------------->
    glPushMatrix();
    glTranslatef(0,rm_w,0);
    main_road();
    glPopMatrix();

      //<---------------------------KUET Design--------------------------->
    glPushMatrix();
    glTranslatef(0,rm_w,-150);
    kuet_design();
    glPopMatrix();

    //<---------------------------Main Gate--------------------------->
    glPushMatrix();
    main_gate();
    glPopMatrix();

    //<---------------------------Central Shahid Minar--------------------------->
    glPushMatrix();
    glTranslatef(300,rm_w,-170);
    glRotatef(-90,0,1,0);
    central_shahid_minar();
    glPopMatrix();

     //<----------------------------Main Hall_Block------------------------------>
    room_check=0;
    glPushMatrix();
     glTranslatef(70,0,-300);
    glRotatef(90,0,1,0);
   khan_jahan_ali_hall();
    glPopMatrix();

    //<----------------------------Trees-------------------------->
   glPushMatrix();
   glTranslatef(180,0,-250);
    loop_of_tree();
    glPopMatrix();

     glPushMatrix();
   glTranslatef(160,0,-250);
    loop_of_tree();
    glPopMatrix();


     glPushMatrix();
   glTranslatef(185,0,-400);
    loop_of_tree();
    glPopMatrix();

     glPushMatrix();
   glTranslatef(165,0,-400);
    loop_of_tree();
    glPopMatrix();

    //<------------------------Auditorium------------------------------>
    glPushMatrix();
    glTranslatef(250,rm_w,-450);
    glRotatef(-90,0,1,0);
    auditorium();
    glPopMatrix();

    //<--------------------Road Side Lamp----------------------------->
    glPushMatrix();
    //glTranslatef(0,10,0);
    //road_side_lamp();
   //spot_light1(0,9,0);
    glPopMatrix();

  /*  //<-----------------------Car---------------------------------->
    glPushMatrix();
    glTranslatef(car_x,rm_w+1,car_z);
    glRotatef(car_rn,0,1,0);
    glRotatef(rn,0,1,0);
   // car();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,5,0);
    glRotatef(rn,0,1,0);
     //wheel();
     glPopMatrix();*/

}


void display()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glFrustum(-1,1,-1,1,2,5000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    gluLookAt(eye_x,eye_y,eye_z,look_x,look_y,look_z, 0,1,0);
    light1(0,1000,0);


    glPushMatrix();
    main_project();
    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}

void LoadTexture(const char*filename,GLuint texture_index)
{
    glGenTextures(1, &texture_id[texture_index]);
    glBindTexture(GL_TEXTURE_2D,texture_id[texture_index]);
    glPixelStorei(GL_UNPACK_ALIGNMENT, texture_id[texture_index]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    BmpLoader bl(filename);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, bl.iWidth, bl.iHeight, GL_RGB, GL_UNSIGNED_BYTE, bl.textureData );
}


int main (int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    glutInitWindowPosition(10,10);
    glutInitWindowSize(height, width);
    glutCreateWindow("Niloy Banik_1607038");

    glShadeModel( GL_SMOOTH );
    glEnable( GL_DEPTH_TEST );
    glEnable(GL_NORMALIZE);
    glEnable(GL_BLEND);

    LoadTexture("G:\\z_bitmap\\grass.bmp",0); /// Floor BackGround
    LoadTexture("G:\\z_bitmap\\sky.bmp",1);/// Sky BackGround
    LoadTexture("G:\\z_bitmap\\road.bmp",2);/// Road
    LoadTexture("G:\\z_bitmap\\minar_tiles.bmp",3);/// Shahid Minar Red Tiles
    LoadTexture("G:\\z_bitmap\\minar_black_tiles.bmp",4);/// Shahid Minar Black Tiles
    LoadTexture("G:\\z_bitmap\\logo.bmp",5);/// Kuet Logo
    LoadTexture("G:\\z_bitmap\\grs.bmp",6);/// Car Wheel
    LoadTexture("G:\\z_bitmap\\car.bmp",7);/// Car Black Background
    LoadTexture("G:\\z_bitmap\\wheel.bmp",8);/// Car Wheel
    LoadTexture("G:\\z_bitmap\\building.bmp",9);/// Building
    LoadTexture("G:\\z_bitmap\\tree.bmp",10);/// Tree Base
    LoadTexture("G:\\z_bitmap\\leaf.bmp",11);/// Tree Leaf
    LoadTexture("G:\\z_bitmap\\white.bmp",12);/// White
    LoadTexture("G:\\z_bitmap\\img.bmp",13);/// Img in Cube
    LoadTexture("G:\\z_bitmap\\wood.bmp",14);/// Table Wood
    LoadTexture("G:\\z_bitmap\\chair_wood.bmp",15);/// Chair Wood
    LoadTexture("G:\\z_bitmap\\bed_wood.bmp",16);/// Bed Wood
    LoadTexture("G:\\z_bitmap\\bed_sheet1.bmp",17);/// Bed Sheet1
    LoadTexture("G:\\z_bitmap\\balish1.bmp",18);/// Balish1
    LoadTexture("G:\\z_bitmap\\floor.bmp",19);/// Room Floor
    LoadTexture("G:\\z_bitmap\\wall.bmp",20);/// Room Wall
    LoadTexture("G:\\z_bitmap\\window.bmp",21);/// Window
    LoadTexture("G:\\z_bitmap\\door.bmp",22);/// Door
    LoadTexture("G:\\z_bitmap\\room_ceil.bmp",23);/// Room Ceil
    LoadTexture("G:\\z_bitmap\\baranda.bmp",24);/// Baranda
    LoadTexture("G:\\z_bitmap\\lower_stair.bmp",25);/// Lower Stair
    LoadTexture("G:\\z_bitmap\\carrom.bmp",26);/// Carrom Board
    LoadTexture("G:\\z_bitmap\\tt_board.bmp",27);/// TT Board
    LoadTexture("G:\\z_bitmap\\net.bmp",28);/// TT Net
    LoadTexture("G:\\z_bitmap\\tv.bmp",29);/// TV
    LoadTexture("G:\\z_bitmap\\mat.bmp",30);/// Mat
    LoadTexture("G:\\z_bitmap\\sofa_main.bmp",31);/// Sofa_Main
    LoadTexture("G:\\z_bitmap\\sofa.bmp",32);/// Sofa
    LoadTexture("G:\\z_bitmap\\sofa_balish.bmp",33);/// Sofa_balish
    LoadTexture("G:\\z_bitmap\\almari_front.bmp",34);/// Almari Front
    LoadTexture("G:\\z_bitmap\\almari_side.bmp",35);/// Almari Side
    LoadTexture("G:\\z_bitmap\\laptop_keyboard.bmp",36);/// Laptop KeyBoard
    LoadTexture("G:\\z_bitmap\\laptop_screen.bmp",37);///Laptop Screen
    LoadTexture("G:\\z_bitmap\\khaja.bmp",38);/// Laptop Front
    LoadTexture("G:\\z_bitmap\\tiles.bmp",39);/// Khaja Hall Road Tiles
    LoadTexture("G:\\z_bitmap\\road1.bmp",40);/// Khaja Hall Road
    LoadTexture("G:\\z_bitmap\\kuet_logo.bmp",41);/// KUET logo
    LoadTexture("G:\\z_bitmap\\kuet_title.bmp",42);/// KUET title
   LoadTexture("G:\\z_bitmap\\cover.bmp",43);/// KUET Triangle Cover
    LoadTexture("G:\\z_bitmap\\audi_title.bmp",44);/// Auditorium Title
    LoadTexture("G:\\z_bitmap\\bannar.bmp",45);/// Auditorium Bannar
    LoadTexture("G:\\z_bitmap\\audi_tex.bmp",46);/// Auditorium Texture
    LoadTexture("G:\\z_bitmap\\audi_texture1.bmp",47);/// Auditorium Texture Others



    cout<<"1607038:- Final Project"<<endl;

    cout<<"Press 'd' for Move Right\n"<<endl;
    cout<<"Press 'a' for Move Left\n"<<endl;
    cout<<"Press 'w' for Go Forward\n"<<endl;
    cout<<"Press 's' for Go Back\n"<<endl;
    cout<<"Press '+' Go Up\n"<<endl;
    cout<<"Press '-' Go Down\n"<<endl;

    cout<<"Press 'D' for CAR Move Right\n"<<endl;
    cout<<"Press 'A' for CAR Move Left\n"<<endl;
    cout<<"Press 'W' for CAR Forward\n"<<endl;
    cout<<"Press 'S' for CAR Back\n"<<endl;

    cout<<"Press 'L' for CAR Translate Right\n"<<endl;
    cout<<"Press 'J' for CAR Translate Left\n"<<endl;

    cout<<"Press 'i' for Look at Up\n"<<endl;
    cout<<"Press 'k' for Look at Down\n"<<endl;

    cout<<"Press '1' for  Open \n"<<endl;
    cout<<"Press '2' for Close\n"<<endl;

    cout<<"Press 'b' for Control Light 1\n"<<endl;
    cout<<"Press 'c' for Control Light 2\n"<<endl;
    cout<<"Press 'v' for Control Light 3\n"<<endl;


    glutDisplayFunc(display);
    glutKeyboardFunc(myKeyboardFunc);
    glutReshapeFunc(resz);
     glutIdleFunc(animate);
    glutMainLoop();

    return 0;
}

