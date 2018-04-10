#define GL_GLEXT_PROTOTYPES
#ifdef __APPLE__
#include<OpenGl/glu.h>
#include<OpenGl/gl.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

 #include <stdio.h>
#include <stdarg.h>
#include <math.h>



#define PI 3.14159265

int convert_2d_to_3d = 0;
int view = 0;

class point
{
 public:

 float x; /*!< x co-ordinate*/
 float y; /*!< y co-ordinate*/
 float z; /*!< z co-ordinate*/
 /**
  Constructor for class point 
  @param a x co-ordinate
  @param b y co-ordinate
  @param c z co-ordinate
 */
 point()
 {
  x=0; 
  y=0;
  z=0;
 }
 /**
 prints the point co-ordinates

*/
 void print()
 {
  printf("\n point co-ordinates are (%f,%f,%f)\n",x,y,z );
 }
};

/**
 returns midpoint if given 2 points
 @param a point 
 @param b point 
 @return  mid_point

*/
point midpoint(point a, point b)
{
  float x = 1/2*(a.x+b.x);
  float y = 1/2*(a.y+b.y);
  float z = 1/2*(a.z+b.z);
  point c;
  c.x = x; c.y = y; c.z = z;
  return c;
}
/**
 rotates the point about x-axis
 @param alpha angle 
 @param p point
 @return new rotated point

*/

point rotation_x(float alpha,point p)
{
  float cos_x,sin_x,z,z1,y,y1;
  cos_x = cos ( alpha * PI / 180.0 );
  sin_x = sin ( alpha * PI / 180.0 );
  z = p.z;
  y = p.y;
  y1 = y*cos_x - z*sin_x;
  z1 = y*sin_x + z*cos_x;
  point news;
  news.x = p.x;
  news.y = y1;
  news.z = z1;
  return news;

}
/**
 rotates the point about y-axis
 @param alpha angle 
 @param p point
 @return news rotated point

*/
point rotation_y(float alpha,point p)
{
  float cos_x,sin_x,z,z1,x,x1;
  cos_x = cos ( alpha * PI / 180.0 );
  sin_x = sin ( alpha * PI / 180.0 );
  z = p.z;
  x = p.x;
  z1 = z*cos_x - x*sin_x;
  x1 = z*sin_x + x*cos_x;
  point news;
  news.x = x1;
  news.y = p.y;
  news.z = z1;
  return news;

}
/**
 rotates the point about z-axis
 @param alpha angle 
 @param p point
 @return news rotated point

*/
point rotation_z(float alpha,point p)
{
  float cos_x,sin_x,x,x1,y,y1;
  cos_x = cos ( alpha * PI / 180.0 );
  sin_x = sin ( alpha * PI / 180.0 );
  x = p.x;
  y = p.y;
  x1 = x*cos_x - y*sin_x;
  y1 = x*sin_x + y*cos_x;
  point news;
  news.x = x1;
  news.y = y1;
  news.z = p.z;
  return news;

}
/**
 produces new point ,given rotation angles about x,y,z-axis respectively
 @param alpha x-angle 
 @param beta y-angle 
 @param gamma_angle z-angle 
 @param p point
 @return news rotated point

*/
point rotation(float alpha,float beta,float gamma_angle ,point p)
{
  point r1,r2,r3;
  r1 = rotation_x(alpha,p);
  r2 = rotation_y(beta,r1);
  r3 = rotation_z(gamma_angle,r2);
  return r3;
}



class equation
{
 public:

 float l; /*!< x coefficient*/
 float m; /*!< y coefficient*/
 float n; /*!< z coefficient*/
 float c; /*!< constant*/
 /**
  Constructor for class point
 */
 equation()
 {
  l=0; 
  m=0;
  n=0;
  c=0;
 }
  /**
 prints the plane equation

*/

 void print()
 {
  printf("\n plane_equation is %d x+ %d y + %d z + %d = 0\n",l,m,n,c );
 }

};


class plane
{
public:
  int n; /*!< number of points*/
  point set_points[100]; /*!< pointer for array of points*/

   /**
  returns centroid of the given face
  @return centroid centre of the given face
 */ 
  point centroid()
  {
      float x,y,z;
      x = 0;
      y = 0;
      z = 0;
      int i;
      for (i = 0; i < n; i++)
      {
        x = x+set_points[i].x;
        y = y+set_points[i].y;
        z = z+set_points[i].z;
      }
      x = x/n;
      y = y/n;
      z = z/n;
      point centroid;
      centroid.x = x;
      centroid.y = y;
      centroid.z = z;
      return centroid;

  }
     /**
  Constructor for class plane
  @param a number of points
  @param A Array of points
 */ 
  plane(int a,point A[100])
  {
    n=a;
    //set_points = point* malloc(a);
    int i = 0;
    for(i=0;i<a;i++)
    {
      set_points[i]=A[i];
    }
  }
    
 /**
  function for equation of a plane object 
  @return A plane equation
 */ 

  equation plane_equation()
  {
    
    /*plane can be represented as lx+my+nz+o=0 ..now l,m,n,o are calculated and put in array A*/
    //ax+by+cz+d=0
    equation A;
    point a = set_points[0];
    point b = set_points[1];
    point c = set_points[2];
    //float *A = new float[4];
    float x1 = a.x;
    float y1 = a.y;
    float z1 = a.z;
    float x2 = b.x;
    float y2 = b.y;
    float z2 = b.z;
    float x3 = c.x;
    float y3 = c.y;
    float z3 = c.z; 
    float l,m,n,o;
    l = (y1-y2)*(z2-z3) - (y2-y3)*(z1-z2);
    m = (z1-z2)*(x2-x3) - (z2-z3)*(x1-x2);
    n = (x1-x2)*(y2-y3) - (x2-x3)*(y1-y2);
    o = -(l*x1+m*y1+n*z1);
    if(l<=0 && m<=0 && n<=0 && o<=0)
    {
      l = -l;
      m = -m;
      n = -n;
      o = -o;
    }
    A.l=l;
    A.m=m;
    A.n=n;
    A.c=o;


    return A;
  } 
  /**
 prints the plane equation

*/

  void print()
  {
        equation A;
    point a = set_points[0];
    point b = set_points[1];
    point c = set_points[2];
    //float *A = new float[4];
    float x1 = a.x;
    float y1 = a.y;
    float z1 = a.z;
    float x2 = b.x;
    float y2 = b.y;
    float z2 = b.z;
    float x3 = c.x;
    float y3 = c.y;
    float z3 = c.z; 
    float l,m,n,o;
    l = (y1-y2)*(z2-z3) - (y2-y3)*(z1-z2);
    m = (z1-z2)*(x2-x3) - (z2-z3)*(x1-x2);
    n = (x1-x2)*(y2-y3) - (x2-x3)*(y1-y2);
    o = -(l*x1+m*y1+n*z1);
    if(l<=0 && m<=0 && n<=0 && o<=0)
    {
      l = -l;
      m = -m;
      n = -n;
      o = -o;
    }
    A.l=l;
    A.m=m;
    A.n=n;
    A.c=o;
    A.print();

  }

};