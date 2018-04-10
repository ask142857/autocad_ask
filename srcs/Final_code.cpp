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
/**
 calculates distance of a point from a plane
 @param a point 
 @param b plane 
 @return distance

*/
float distance_point_plane(point a,plane b)
{
  float d;
  equation A = b.plane_equation();

  float l,m,n,o,x,y,z;
  l = A.l;
  m = A.m;
  n = A.n;
  o = A.c;
  x = a.x;
  y = a.y;
  z = a.z;
  float sq  = l*l+m*m+n*n;
  d = (l*x+m*y+n*z+o)/pow(sq,0.5);
    /**/
  return d;
}

/**
 returns  the projected point
 @param a point input
 @param b plane on to which image is projected
 @output val projected point
 

*/
point projection(point a , plane b)
{
  float x = a.x;
  float y = a.y;
  float z = a.z;
  equation A = b.plane_equation();
  float l = A.l;
  float m = A.m;
  float n = A.n;
  float o = A.c;
  float d = distance_point_plane(a,b);
  float sq = l*l+m*m+n*n;
  d= d/pow(sq,0.5);
  float x0 = x - l*d;
  float y0 = y - m*d;
  float z0 = z - n*d;
  point val;
  val.x = x0;
  val.y = y0;
  val.z = z0;

  
  return val;
}
/**
 calculates distance between line bc and point a
 @param a point input
 @param b starting point
 @param c ending point

 @output k distance
 

*/

float above_the_line (point a , point b , point c)
{
  // x-x2/x-x3 = y-y2/y-y3
  float k = (a.x-b.x)*(a.y-c.y) - (a.x - c.x)*(a.y - b.y);
  return k;
}

/**
 returns 1 if point is above the plane else 0 
 @param a point input
 @param b plane input
 

*/


int above_the_plane(point a , plane b)
{
  point c = projection(a,b);
  if(a.z-c.z>=-0.0001) return 1;
  return 0;
}

/**
 returns difference between z co-ordinates of the original point and projected one
 @param a point input
 @param b plane input
 

*/
float ask_plane(point a , plane b)
{
	point c = projection(a,b);
	return (a.z - c.z);
}
/**
 checks  the inside point if it is hidden according to a plane or not
 @param a point to be checken
 @param b checking plane 
 @param project projection plane
 @return check gives 1 if it is visble , 0 if hidden
 ;
 

*/
int hidden_point_checker(point a,plane b, plane project)
{
  if(above_the_plane(a,b)==1)
    return 1;
  else
  {  

    point centroid = b.centroid();
    //point proj = projection(a,b);

   // point set_points[100] ;
   // set_points= b.set_points;
    int n = b.n;
    int i = 0;
    for(i=0 ; i<n ; i++)
    {
      int k = (i+1)%n;
      point start = b.set_points[i];
      point end = b.set_points[k];

      float a1 = above_the_line(centroid,start,end);
      float b1 = above_the_line(a,start,end);
      if(a1*b1<=0)
      {
        return 1;
      }


    }
    float v = ask_plane(a,b);

      
 //printf("\n case 2 \n");
 printf("\n z difference is %f",v);
 //a.print();  
 //printf(" is below the plane ");   
 //b.print();
 return 0;


  }
    }




// ----------------------------------------------------------
// Function Prototypes
// ----------------------------------------------------------
 
// ----------------------------------------------------------
// Global Variables
// ----------------------------------------------------------
double rotate_z=0;
double rotate_y=0; 
double rotate_x=0;
int n_of_vertices;//number of vertices
int n_of_faces;//number of faces 

int n_of_vertices_face[100];
float vertices [100][3];//let limit for number of vertices be 100
int facets [100][20];// number of faces limit 100 , vertices per face limit 20
//void read_file (char **argv);
int n_of_vertices_2;
int n_of_faces_2;
float vertices_xy[100][2];
float vertices_yz[100][2];
float vertices_zx[100][2];
int n_of_vertices_face_2_xy[100];
int n_of_vertices_face_2_yz[100];
int n_of_vertices_face_2_xz[100];
int facets_xy[100][20];
int facets_yz[100][20];
int facets_xz[100][20];


 
// ----------------------------------------------------------
// display() Callback function
// ----------------------------------------------------------

//plane be lx+my+n+c = 0;
int l_out,m_out,n_out,c_out;
float alpha,beta,gamma_angle;
/**
 checks if both arrays are equal in a way 
 @param size array size
 @param A array 1
 @param B array 2

 

*/
int equality(int size,float A[],float B[])
{
  int i =0;
  int j = 0;
  for(i=0;i<size;i++)
  {
    int offset = i; int count = 1;
    for(j=0;j<size;j++)
    {
      if(A[j] != B[(j+offset)%size])
        count = 0;
    }
    if(count ==1)
    return 1;
  }
  return 0 ;
}
/**
takes 2d input and gives 3d output

 

*/
void display_2d_to_3d()
{
  
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  glRotatef( rotate_x, 1.0, 0.0, 0.0 );
  glRotatef( rotate_y, 0.0, 1.0, 0.0 );
  glRotatef( rotate_z, 0.0, 0.0, 1.0 );

        glColor3f(0.0,0,1.0);
        glLineWidth(0.0001f);
        glBegin(GL_LINES);
        glVertex3f(0,0,0);
        glVertex3f(100000,0,0);
        glEnd();
        glColor3f(0.0,1,0.0);
        glLineWidth(0.0001f);
        glBegin(GL_LINES);
        glVertex3f(0,0,0);
        glVertex3f(0,100000,00000);
        glEnd();
        glColor3f(1.0,0,0.0);
        glLineWidth(0.0001f);
        glBegin(GL_LINES);
        glVertex3f(0,0,0);
        glVertex3f(00000,00000,100000);
        glEnd();
  // actual vertice 1  face 1 xy ..
  int i,j,k,l;
  
  for(i=0;i<n_of_faces_2;i++)
  {
    // facets_xy
   
      // face 1 contains x1 y1 x2 y2 ....
      float y1[n_of_vertices_face_2_xy[i]];
      for(j = 0; j<n_of_vertices_face_2_xy[i] ; j++)
      {
        y1[j] = vertices_xy[facets_xy[i][j]][1];
      }
      
      for(k=0;k<n_of_faces_2;k++)
      {
        if(n_of_vertices_face_2_xy[i] == n_of_vertices_face_2_yz[k])

         {
           float y2[n_of_vertices_face_2_yz[k]];
           for(l = 0; l<n_of_vertices_face_2_yz[k] ; l++)
           {
             y2[l] = vertices_yz[facets_yz[k][l]][0];
           }
           if(equality(n_of_vertices_face_2_xy[k],y1,y2)==1)
           {
              
              
              //glColor3f();
              int i1=0;
              float x[n_of_vertices_face_2_xy[i]];
              float y[n_of_vertices_face_2_xy[i]];
              float z[n_of_vertices_face_2_xy[i]];
              for(i1=0;i1<n_of_vertices_face_2_xy[i];i1++)
              { 

               x[i1] = vertices_xy[facets_xy[i][i1]][0];
               y[i1] = vertices_xy[facets_xy[i][i1]][1];
               z[i1] = vertices_yz[facets_yz[k][i1]][1];
                //glVertex3fv(x1,y1,z1);
              }
              glBegin(GL_POLYGON);
              glColor3f(x[0],y[0],z[0]);
              for(i1=0;i1<n_of_vertices_face_2_xy[i];i1++)
              { 
                glVertex3f(x[i1],y[i1],z[i1]);
              }              

              glEnd();

           }
         }    

       


      }
      
    
  }

  glFlush();
  glutSwapBuffers();
 


}

/**
takes 3d input and gives orthographic views

 

*/

void checker_2()
{

if(convert_2d_to_3d==0)
{
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
 // glRotatef( rotate_x, 1.0, 0.0, 0.0 );
  //glRotatef( rotate_y, 0.0, 1.0, 0.0 );
  //glRotatef( rotate_z, 0.0, 0.0, 1.0 );


          glColor3f(0.0,0,1.0);
        glLineWidth(0.0001f);
        glBegin(GL_LINES);
        glVertex3f(0,0,0);
        glVertex3f(100000,0,0);
        glEnd();
        glColor3f(0.0,1,0.0);
        glLineWidth(0.0001f);
        glBegin(GL_LINES);
        glVertex3f(0,0,0);
        glVertex3f(0,100000,00000);
        glEnd();
        glColor3f(1.0,0,0.0);
        glLineWidth(0.0001f);
        glBegin(GL_LINES);
        glVertex3f(0,0,0);
        glVertex3f(00000,00000,100000);
        glEnd();

int mid_hid[n_of_faces][40];
//point mid_points[100][100];
//int hidmid[100][100];
 for(int ia=0; ia< n_of_faces ; ia++)
 {
  for(int ka=0;ka<40;ka++)
  {
    //hidden[ia][ka]=1;
    mid_hid[ia][ka]=1;
  }
 }


point proj[100];
 point bc;
 bc.x = -2;
 proj[0] = bc;

  point bd;
  bd.y = 3;

 proj[1]=bd;

 point be;
 be.x = 1;
 proj[2]=be;

 plane ortho(3,proj);//xy plane

 int count = 0;
  point vertix[n_of_vertices];
  int ask;
  int check_hidden_points[n_of_vertices];

// points initialized
  for(ask=0;ask<n_of_vertices;ask++)
  {
    point temp;
    temp.x = vertices[ask][0];
    temp.y = vertices[ask][1];
    temp.z = vertices[ask][2];
    if(view==1)
    vertix[ask] = rotation(alpha,beta,gamma_angle,temp);
    else if(view==2)
    vertix[ask] = rotation(alpha,beta+90,gamma_angle,temp);	
    else if(view==3)
    vertix[ask] = rotation(alpha+90,beta,gamma_angle,temp);		
    check_hidden_points[ask] = 1;

  }
int i,j,k,l,m,n,o;
i=0;
j=0;
k=0;
l=0;
m=0;
n=0;
o=0;

// mid point calculation

point mid_points[100][100];
for(i=0;i<n_of_faces;i++)
{
   for(j=0;j<n_of_vertices_face[i];j++)
   {
    k = (j+1)%n_of_vertices_face[i];
    point start = vertix[facets[i][j]];
    point end = vertix[facets[i][k]];
  point mid;  
   mid.x = 0.5*(start.x+end.x);
  mid.y = 0.5*(start.y+end.y);
  mid.z = 0.5*(start.z+end.z);
  
  
    //point mid = midpoint(start,end);
    mid_points[i][j] = mid;
    int u1 = (j+1)%(n_of_vertices_face[i]);

        int start_number = facets[i][j];
        int end_number = facets[i][u1];
       // printf("\nvertex %d is ",start_number);
       // start.print();
       // printf("\nvertex %d is ",end_number);
       // end.print();
                //  printf("\nmidpoint of edge %d___%d is ",start_number,end_number );
         // mid.print();
        // printf("\n");


   }
}

// vertices hidden or not

for(i=0;i<n_of_vertices;i++)
{
  for(j=0;j<n_of_faces;j++)
  {
    point face_vertices[100];
    for(k=0;k<n_of_vertices_face[j];k++)
    {
      face_vertices[k] = vertix[facets[j][k]];
    }

    plane face_checker(n_of_vertices_face[j],face_vertices);

    check_hidden_points[i] = check_hidden_points[i] * hidden_point_checker(vertix[i],face_checker,ortho);
  }// j closes
}// i closes


// midpoints hidden calculation

for(i=0;i<n_of_faces;i++)
{
  for(j=0;j<n_of_vertices_face[i];j++)
  {
    for(k=0;k<n_of_faces;k++)
    {
      point face_vertices[100];
      for(l = 0 ;l<n_of_vertices_face[k];l++)
      {
        face_vertices[l] = vertix[facets[k][l]];
      }
      plane face_checker(n_of_vertices_face[k],face_vertices);
       //printf("\n case 2 \n");
 
 //printf(" is below the plane ");   
 //b.print();

      mid_hid[i][j] = mid_hid[i][j] * hidden_point_checker(mid_points[i][j],face_checker,ortho);
       
       if(hidden_point_checker(mid_points[i][j],face_checker,ortho)==0)
      {
       printf("\n midpoint");
       mid_points[i][j].print(); 
       printf("start vertex %d",facets[i][j]);
       vertix[facets[i][j]].print();
       printf("end vertex %d",facets[i][(j+1)%(n_of_vertices_face[i])]);
       vertix[facets[i][(j+1)%(n_of_vertices_face[i])]].print();       

       printf(" is below the plane ");
       face_checker.print();
      }
    }
  }
}






  //int i =0;int j = 0; int k =0;
 for(i=0;i<n_of_faces;i++)
  {
    //glColor3f(x_color,y_color,z_color);
    for(k=0;k<n_of_vertices_face[i];k++)
    {
       int u1 = (k+1)%(n_of_vertices_face[i]);
       int x = check_hidden_points[facets[i][k]] * check_hidden_points[facets[i][u1]]* mid_hid[i][k];

        int start_number = facets[i][k];
        int end_number = facets[i][u1];
        int a1 = check_hidden_points[start_number];
        int a2 = check_hidden_points[end_number];
        int a12 = mid_hid[i][k];

        if(a1 == 0)
        {
          printf("point %d is hidden\n",start_number );
        }
        if(a2==0)
        {
          printf("point %d is hidden\n",end_number);
        }
        if(a1!=0 && a2!=0 && a12==0)
        {
          printf("midpoint of edge %d___%d is hidden but not corners\n",start_number,end_number );
          mid_points[i][k].print();
        }

        point start = projection(vertix[facets[i][k]] ,ortho);
        point end = projection(vertix[facets[i][u1]] , ortho);
        //float xo1,xo2,yo1,yo2,zo1,zo2;
        float xo1 = vertix[facets[i][k]].x; 
        float yo1 = vertix[facets[i][k]].y ;
        float zo1 = vertix[facets[i][k]].z ;
        float xo2 = vertix[facets[i][u1]].x ;
        float yo2 = vertix[facets[i][u1]].y ;
        float zo2 = vertix[facets[i][u1]].z ;
        float x1 = start.x;
        float y1 = start.y;
        float z1 = start.z;
        float x2 = end.x;
        float y2 = end.y;
        float z2 = end.z;
      if(x==1)
      { 
        //printf("chotu");
        glColor3f(1.0,0,0);
        glLineWidth(5.0f);
        glBegin(GL_LINES);
        glVertex3f(x1,y1,z1);
        glVertex3f(x2,y2,z2);
        glEnd();
      }
      else
      {
        //printf("edge %d -- %d is hidden\n",start_number,end_number );
        glColor3f(1.0,1.00,0.0);
        glLineWidth(0.1f);
        glBegin(GL_LINES);
        glVertex3f(x1,y1,z1);
        glVertex3f(x2,y2,z2);
        glEnd();
        //printf("bhosdi");

      }
    }
  }
   glFlush();
  glutSwapBuffers();
}

else if(convert_2d_to_3d==1)
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  glRotatef( rotate_x, 1.0, 0.0, 0.0 );
  glRotatef( rotate_y, 0.0, 1.0, 0.0 );
  glRotatef( rotate_z, 0.0, 0.0, 1.0 );
  
     glColor3f(0.0,0,1.0);
        glLineWidth(0.0001f);
        glBegin(GL_LINES);
        glVertex3f(0,0,0);
        glVertex3f(100000,0,0);
        glEnd();
        glColor3f(0.0,1,0.0);
        glLineWidth(0.0001f);
        glBegin(GL_LINES);
        glVertex3f(0,0,0);
        glVertex3f(0,100000,00000);
        glEnd();
        glColor3f(1.0,0,0.0);
        glLineWidth(0.0001f);
        glBegin(GL_LINES);
        glVertex3f(0,0,0);
        glVertex3f(00000,00000,100000);
        glEnd();




  // Line

 // polygon
  int i=0;



  for(i=0;i<n_of_vertices;i++)
  {
  	    point p,q;
    	p.x = vertices[i][0];
    	p.y = vertices[i][1];
    	p.z = vertices[i][2];
    	q = rotation(alpha,beta,gamma_angle,p);
    	vertices[i][0] = q.x;
    	vertices[i][1] = q.y;
    	vertices[i][2] = q.z;
    
  }
  for (i = 0; i < n_of_faces; ++i)
  {
    float x_color,y_color,z_color;
    x_color = vertices[i][0];
    y_color = vertices[i][1];
    z_color = vertices[i][2];
    glBegin(GL_POLYGON);
    glColor3f(x_color,y_color,z_color);
    int i1=0;
    
    for(i1=0;i1<n_of_vertices_face[i];i1++)
    {
      glVertex3fv(vertices[facets[i][i1]]);
    }
    glEnd();
  }
  //glLineWidth(2.0f);
  glFlush();
  glutSwapBuffers();
 

}

else if(convert_2d_to_3d == 2)
{
	display_2d_to_3d();
}

} 

//global declaration



 




// ----------------------------------------------------------
// specialKeys() Callback Function
// ----------------------------------------------------------
/**
takes keyboard input and rotates the axes
@param key key
@param x  int
@param y  int
 

*/

void specialKeys( int key, int x, int y ) {
 
  //  Right arrow - increase rotation by 5 degree
  if (key == GLUT_KEY_RIGHT)
    rotate_y += 5;
 
  //  Left arrow - decrease rotation by 5 degree
  else if (key == GLUT_KEY_LEFT)
    rotate_y -= 5;
 
  else if (key == GLUT_KEY_UP)
    rotate_x += 5;
 
  else if (key == GLUT_KEY_DOWN)
    rotate_x -= 5;
  else if (key == GLUT_KEY_F1)
    rotate_z += 5;
  else if (key == GLUT_KEY_F2)
    rotate_z -= -5;
 
  //  Request display update
  glutPostRedisplay();
 
}

// ----------------------------------------------------------
// processKeys() Callback Function
// ----------------------------------------------------------
/**
takes keyboard input and does the changes on screen
@param key key
@param x  int
@param y  int
 

*/
void processKeys( unsigned char key, int x, int y ) {
 
  //  Right arrow - increase rotation by 5 degree
  if (key == 'y')
    
    {
    	rotate_y += 5;
    	beta +=5;
    }
 
  //  Left arrow - decrease rotation by 5 degree
  else if (key == 'b')
   {
    rotate_y -= 5;
    beta -=5;
   }
 
  else if (key == 'x')
    {
    	rotate_x += 5;
    	alpha += 5;
 }
  else if (key == 'a')
    {
     rotate_x -= 5;
     alpha -= 5;
    }
  else if (key == 'z')
    {
    	rotate_z += 5;
    	gamma_angle +=5;
    }
  else if (key == 'c')
    {
    	rotate_z -= 5;
    	gamma_angle -= 5;
    }
  else if (key == 'o')
  {
  	convert_2d_to_3d = 1;

  }  
  else if (key == 'r' )
  {
   convert_2d_to_3d = 0;	
  }
  else if(key == 'q')
  {
  	convert_2d_to_3d = 2;
  }

  else if(key == 'f')
  {
  	view = 1;
  }
  else if(key == 't')
  {
  	view = 2;
  }
  else if(key == 's')
  {
  	view = 3;
  }
 
  //  Request display update
  glutPostRedisplay();
 
}

/**
read file for 2d  to 3d
*/


void read_file_2 (char **argv)
{
  FILE *Map2;
  int loopcounter;
  Map2 = fopen("Map2.txt","r");
  if(Map2 == NULL)
  {
    printf("\ncan't open input file\n");
  }
 
  fscanf(Map2,"%d", &n_of_vertices_2);
  fscanf(Map2,"%d", &n_of_faces_2);

  for (loopcounter = 0 ; loopcounter < n_of_vertices_2 ; loopcounter++ )
  {
    fscanf(Map2,"%f", &vertices_xy[loopcounter][0]);
          //printf("%f ",vertices[loopcounter][0]);
          
    fscanf(Map2,"%f", &vertices_xy[loopcounter][1]);
    //printf("%f ",vertices[loopcounter][1]);
    
    //fscanf(Map,"%f", &vertices[loopcounter][2]);
   // printf("%f\n",vertices[loopcounter][2]);
  }


  for (loopcounter = 0 ; loopcounter < n_of_faces_2 ; loopcounter++ )
  {
    
    fscanf(Map2,"%d", &n_of_vertices_face_2_xy[loopcounter]);
    //printf("%d ",n_of_vertices_face[loopcounter]);
    int count = 0;
    for(count=0;count<n_of_vertices_face_2_xy[loopcounter];count++)
    {
      fscanf(Map2,"%d", &facets_xy[loopcounter][count]);
      //printf("%d ",facets[loopcounter][count]);
    }
  }



  for (loopcounter = 0 ; loopcounter < n_of_vertices_2 ; loopcounter++ )
  {
    fscanf(Map2,"%f", &vertices_yz[loopcounter][0]);
          //printf("%f ",vertices[loopcounter][0]);
          
    fscanf(Map2,"%f", &vertices_yz[loopcounter][1]);
    //printf("%f ",vertices[loopcounter][1]);
    
    //fscanf(Map,"%f", &vertices[loopcounter][2]);
   // printf("%f\n",vertices[loopcounter][2]);
  }
    for (loopcounter = 0 ; loopcounter < n_of_faces_2 ; loopcounter++ )
  {
    
    fscanf(Map2,"%d", &n_of_vertices_face_2_yz[loopcounter]);
    //printf("%d ",n_of_vertices_face[loopcounter]);
    int count = 0;
    for(count=0;count<n_of_vertices_face_2_yz[loopcounter];count++)
    {
      fscanf(Map2,"%d", &facets_yz[loopcounter][count]);
      //printf("%d ",facets[loopcounter][count]);
    }
  }
  for (loopcounter = 0 ; loopcounter < n_of_vertices_2 ; loopcounter++ )
  {
    fscanf(Map2,"%f", &vertices_zx[loopcounter][0]);
          //printf("%f ",vertices[loopcounter][0]);
          
    fscanf(Map2,"%f", &vertices_zx[loopcounter][1]);
    //printf("%f ",vertices[loopcounter][1]);
    
    //fscanf(Map,"%f", &vertices[loopcounter][2]);
   // printf("%f\n",vertices[loopcounter][2]);
  }

  for (loopcounter = 0 ; loopcounter < n_of_faces_2 ; loopcounter++ )
  {
    
    fscanf(Map2,"%d", &n_of_vertices_face_2_xz[loopcounter]);
    //printf("%d ",n_of_vertices_face[loopcounter]);
    int count = 0;
    for(count=0;count<n_of_vertices_face_2_xz[loopcounter];count++)
    {
      fscanf(Map2,"%d", &facets_xz[loopcounter][count]);
      //printf("%d ",facets[loopcounter][count]);
    }
  }


  
  

}
/**
read file for 3d  to 2d
*/
void read_file (char **argv)
{
  FILE *Map;
  int loopcounter;
  Map = fopen("Map.txt", "r");
  
  if (Map == NULL)
  {
    printf("\ncan't open input file\n");
    //exit (1);
  }
    fscanf(Map,"%f", &alpha);
    //printf("%f ",alpha);
    fscanf(Map,"%f", &beta);
    //printf("%f ",beta);
    fscanf(Map,"%f", &gamma_angle);
    //printf("%f ",gamma_angle);
    
    fscanf(Map,"%d", &n_of_vertices);
    //printf("%d ",n_of_vertices);

    fscanf(Map,"%d", &n_of_faces);
    //printf("%d ",n_of_faces);
    //n_of_vertices_face = new int[n_of_faces];


  for (loopcounter = 0 ; loopcounter < n_of_vertices ; loopcounter++ )
  {
    fscanf(Map,"%f", &vertices[loopcounter][0]);
          //printf("%f ",vertices[loopcounter][0]);
          
    fscanf(Map,"%f", &vertices[loopcounter][1]);
    //printf("%f ",vertices[loopcounter][1]);
    
    fscanf(Map,"%f", &vertices[loopcounter][2]);
   // printf("%f\n",vertices[loopcounter][2]);
  }
  
  for (loopcounter = 0 ; loopcounter < n_of_faces ; loopcounter++ )
  {
    
    fscanf(Map,"%d", &n_of_vertices_face[loopcounter]);
    //printf("%d ",n_of_vertices_face[loopcounter]);
    int count = 0;
    for(count=0;count<n_of_vertices_face[loopcounter];count++)
    {
      fscanf(Map,"%d", &facets[loopcounter][count]);
      //printf("%d ",facets[loopcounter][count]);
    }
  }
}
 
// ----------------------------------------------------------
// main() function
// ----------------------------------------------------------
int main(int argc, char* argv[]){
 
  read_file(argv);
  read_file_2(argv);
  //  Initialize GLUT and process user parameters
  glutInit(&argc,argv);
 
  //  Request double buffered true color window with Z-buffer
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
 
  // Create window
  glutCreateWindow("AutoCad_ask_rohan");
 
  //  Enable Z-buffer depth test
  glEnable(GL_DEPTH_TEST);
 
  // Callback functions
  //checker_2();
  //checker();
  //if(convert_2d_to_3d == 1)
  //glutDisplayFunc(display);
  //else	
  glutDisplayFunc(checker_2);
  glutSpecialFunc(specialKeys);
  glutKeyboardFunc(processKeys);
 
  //  Pass control to GLUT for events
  glutMainLoop();
 
  //  Return to OS
  return 0;
 
}