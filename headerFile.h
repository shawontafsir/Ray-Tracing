#include<bits/stdc++.h>
#include <windows.h>
#include <glut.h>
#include "bitmap_image.hpp"
using namespace std;
#define pi (2*acos(0.0))
#define EPSILON 0.0000001

extern int r_level;
bitmap_image texture("floor.bmp");

class Point{
public:
	double x,y,z;
	Point(){}
	Point(double ax, double ay, double az){
        x = ax; y = ay; z = az;
	}

	Point operator + (Point const &other){
        return Point(x+other.x, y+other.y, z+other.z);
	}

	Point operator - (Point const &other){
        return Point(x-other.x, y-other.y, z-other.z);
	}

	Point operator * (double const &other){
        return Point(x*other, y*other, z*other);
	}
};

class Ray{
public:
    Point start, dir;
    Ray(){}
    Ray(Point s, Point d){
        start = s;
        double a = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
        dir.x=d.x/a; dir.y=d.y/a; dir.z=d.z/a;
    }
};

class Object{
public:
    Point reference_point;
    double height, width, length;
    int shine;
    double color[3];
    double co_efficients[4];

    Object(){}

    virtual void draw(){}

    void setColor(double ax, double ay, double az){
        color[0]=ax; color[1]=ay; color[2]=az;
    }

    void setShine(int s){shine=s;}

    void setCoEfficients(double aw,double ax, double ay, double az){
        co_efficients[0]=aw; co_efficients[1]=ax;
        co_efficients[2]=ay; co_efficients[3]=az;
    }

    virtual double intersect(Ray *r, double &colx, double &coly, double &colz, int level){
        return -1;
    }

    double dot(Point p1, Point p2){
        return p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
    }

    Point cross(Point p1, Point p2){
        return Point(p1.y*p2.z-p1.z*p2.y, p1.z*p2.x-p1.x*p2.z, p1.x*p2.y-p1.y*p2.x);
    }

    Point normalize(Point p){
        double d = sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
        return Point(p.x/d, p.y/d, p.z/d);
    }
};

vector<Object*> objects;
vector<Point> lights;

class Sphere:public Object{
public:
    Sphere(){}
    Sphere(Point Center, double Radius){
        reference_point=Center;
        length=Radius;
    }

    void draw(){
        glColor3f(color[0], color[1], color[2]);
        glPushMatrix();
        {
            glTranslatef(reference_point.x, reference_point.y, reference_point.z);
            glutSolidSphere(length,60,20);
        }
        glPopMatrix();
    }

    Point getNormal(Point ip) {return ip-reference_point;}

    void illumination(Point intersectionPoint, Point start, Point re, double &colx, double &coly, double &colz){
        for(unsigned i=0;i<lights.size();i++){
            Point d = lights[i]-intersectionPoint;
            d = normalize(d);
            Point s = intersectionPoint + (d*1);
            Ray *L = new Ray(s, d);

            int nearest = getColorAt(L);

            if(nearest==-1){
                Point m = intersectionPoint-reference_point;
                m = normalize(m);
                Point v = intersectionPoint-start;

                double lambert = max(dot(d,m), 0.0);
                double phong = dot(normalize(v), re);

                colx+= lambert*co_efficients[1]*color[0];
                coly+= lambert*co_efficients[1]*color[1];
                colz+= lambert*co_efficients[1]*color[2];

                colx+= pow(phong, shine)*co_efficients[2]*color[0];
                coly+= pow(phong, shine)*co_efficients[2]*color[1];
                colz+= pow(phong, shine)*co_efficients[2]*color[2];
            }
        }
    }

    Point getReflection(Point dir, Point ip){
        Point normal = getNormal(ip);
        Point rf;
        rf.x = dir.x - 2*dot(dir, normal)*normal.x;
        rf.y = dir.y - 2*dot(dir, normal)*normal.y;
        rf.z = dir.z - 2*dot(dir, normal)*normal.z;

        return rf;
    }

    Point getRefrection(Point dir, Point ip, double eta){
        Point normal = getNormal(ip);
        double N_dot_I = dot(normal, dir);
        double k = 1.f - eta * eta * (1.f - N_dot_I * N_dot_I);
        if (k < 0.f) return Point(0, 0, 0);
        else return (dir*eta)-(normal*(eta * N_dot_I + sqrt(k)));
    }

    int getColorAt(Ray* L){
        int nearest=-1;
        double t_min = 50000;
        for(int k=0;k<objects.size();k++)
        {
            double dcx, dcy, dcz;
            double t = objects[k]->intersect(L, dcx, dcy, dcz, 0);

            if(t<=0) continue;
            if(t<t_min){
                nearest = k;
                t_min = t;
            }
        }
        return nearest;
    }

    double intersect(Ray *r, double &colx, double &coly, double &colz, int level){
        Point start = r->start, dir = r->dir;
        double a=1;
        double b=2*((start.x-reference_point.x)*dir.x+(start.y-reference_point.y)*dir.y+(start.z-reference_point.z)*dir.z);
        double c=(start.x-reference_point.x)*(start.x-reference_point.x)+(start.y-reference_point.y)*(start.y-reference_point.y)+(start.z-reference_point.z)*(start.z-reference_point.z) - length*length;
        double d=b*b-4*a*c;

        double t;
        if(d<0) return -1;
        else{
            double t1,t2;
            d=sqrt(d);
            t1 = (-b+d)/(2*a);
            t2 = (-b-d)/(2*a);

            t = t1<t2? t1:t2;
        }

        if(t<=0) return -1;
        if(level==0)return t;

        Point intersectionPoint = start + (dir*t);
        colx = color[0]; coly = color[1]; colz = color[2];

        colx*=co_efficients[0];
        coly*=co_efficients[0];
        colz*=co_efficients[0];

        if(length>15){
            Point reflection = getReflection(dir, intersectionPoint);
            reflection = normalize(reflection);

            illumination(intersectionPoint, start, reflection, colx, coly, colz);

            if(level<r_level){
                Point s = intersectionPoint + (reflection*1);
                Ray *refRay = new Ray(s, reflection);

                int nearest = getColorAt(refRay);

                if(nearest!=-1){
                    double rfcolx, rfcoly, rfcolz;
                    objects[nearest]->intersect(refRay, rfcolx, rfcoly, rfcolz, level+1);
                    colx += rfcolx*co_efficients[3];
                    coly += rfcoly*co_efficients[3];
                    colz += rfcolz*co_efficients[3];
                }

            }
        }
        else{
            Point refrection = getRefrection(dir, intersectionPoint, 1.3);
            refrection = normalize(refrection);

            illumination(intersectionPoint, start, refrection, colx, coly, colz);

            if(level<r_level){
                Point s = intersectionPoint + (refrection*1);
                Ray *refRay = new Ray(s, refrection);

                int nearest = getColorAt(refRay);

                if(nearest!=-1){
                    double rfcolx, rfcoly, rfcolz;
                    objects[nearest]->intersect(refRay, rfcolx, rfcoly, rfcolz, level+1);
                    colx += rfcolx*co_efficients[3];
                    coly += rfcoly*co_efficients[3];
                    colz += rfcolz*co_efficients[3];
                }
            }
        }


        if(colx>1) colx=1; else if(colx<0) colx=0;
        if(coly>1) coly=1; else if(coly<0) coly=0;
        if(colz>1) colz=1; else if(colz<0) colz=0;


        return t;
    }
};

class Floor:public Object{
public:
    Floor(double floorWidth, double tileWidth){
        reference_point=Point(floorWidth/5, floorWidth/5, 0);
        length=tileWidth;
    }
    void draw(){
        glPushMatrix();
        {
            glTranslatef(-reference_point.x/2, -reference_point.y/2, 0);
            glBegin(GL_QUADS);
            for (unsigned int x=0;x<reference_point.x;x+=length){
                for (unsigned int y=0;y<reference_point.y;y+=length)
                {
                    if (((y/(int)length)+(x/(int)length)) & 1)  glColor3f(1.0f,1.0f,1.0f);
                    else glColor3f(0.0f,0.0f,0.0f);

                    glVertex3f(x, y, 0);
                    glVertex3f((x+length), y, 0);
                    glVertex3f((x+length), (y+length), 0);
                    glVertex3f(x, (y+length), 0);
                }
            }
            glEnd();
        }
        glPopMatrix();
    }

    Point getNormal(){
        return Point(0,0,1);
    }

    bool withinFloor(Point ip){
        if(ip.x>=-reference_point.x/2&&ip.x<=reference_point.x-reference_point.x/2 && ip.y>=-reference_point.y/2&&ip.y<=reference_point.y-reference_point.y/2)
            return true;
        else return false;
    }

    Point getReflection(Point dir){
        Point normal = getNormal();
        Point rf;
        rf.x = dir.x - 2*dot(dir, normal)*normal.x;
        rf.y = dir.y - 2*dot(dir, normal)*normal.y;
        rf.z = dir.z - 2*dot(dir, normal)*normal.z;

        return rf;
    }

    Point getRefrection(Point dir, double eta){
        Point normal = getNormal();
        double N_dot_I = dot(normal, dir);
        double k = 1.f - eta * eta * (1.f - N_dot_I * N_dot_I);
        if (k < 0.f) return Point(0, 0, 0);
        else return (dir*eta)-(normal*(eta * N_dot_I + sqrt(k)));
    }


    void getColor(Point ip, double &cx, double &cy, double &cz){
        int x = (ip.x+reference_point.x/2);
        int y = (ip.y+reference_point.y/2);
        if (((y/(int)length)+(x/(int)length)) & 1) {cx=1;cy=1;cz=1;}
        else {cx=0;cy=0;cz=0;}
    }

    int getColorAt(Ray* L){
        int nearest=-1;
        double t_min = 5000000;
        for(int k=0;k<objects.size();k++)
        {
            double dcx, dcy, dcz;
            double t = objects[k]->intersect(L, dcx, dcy, dcz, 0);

            if(t<=0) continue;
            if(t<t_min){
                nearest = k;
                t_min = t;
            }
        }
        return nearest;
    }

    void illumination(Point intersectionPoint, Point start, Point re, double &colx, double &coly, double &colz){
        for(unsigned i=0;i<lights.size();i++){
            Point d = lights[i]-intersectionPoint;
            d = normalize(d);
            Point s = intersectionPoint + (d*1);
            Ray *L = new Ray(s, d);

            int nearest = getColorAt(L);

            if(nearest==-1){
                Point m = intersectionPoint-reference_point;
                m = normalize(m);
                Point v = intersectionPoint-start;

                double lambert = max(dot(d,m), 0.0);
                double phong = dot(normalize(v), re);

                colx+= lambert*co_efficients[1]*color[0];
                coly+= lambert*co_efficients[1]*color[1];
                colz+= lambert*co_efficients[1]*color[2];

                colx+= pow(phong, shine)*co_efficients[2]*color[0];
                coly+= pow(phong, shine)*co_efficients[2]*color[1];
                colz+= pow(phong, shine)*co_efficients[2]*color[2];
            }
        }
    }


    double intersect(Ray *r, double &colx, double &coly, double &colz, int level){
        Point start = r->start, dir = r->dir;
        double t=(reference_point.z-start.z)/dir.z;
        Point intersectionPoint = start + (dir*t);

        if(!withinFloor(intersectionPoint)) return -1;
        if(level==0) return t;

        getColor(intersectionPoint, colx, coly, colz);

        Point reflection = getReflection(dir);
        reflection = normalize(reflection);

        illumination(intersectionPoint, start, reflection, colx, coly, colz);

        if(level<r_level){
            Point s = intersectionPoint + (reflection*1);
            Ray *refRay = new Ray(s, reflection);

            int nearest = getColorAt(refRay);

            if(nearest!=-1){
                double rfcolx, rfcoly, rfcolz;
                objects[nearest]->intersect(refRay, rfcolx, rfcoly, rfcolz, level+1);
                colx += rfcolx*co_efficients[3];
                coly += rfcoly*co_efficients[3];
                colz += rfcolz*co_efficients[3];
            }
        }

        if(colx>1) colx=1; else if(colx<0) colx=0;
        if(coly>1) coly=1; else if(coly<0) coly=0;
        if(colz>1) colz=1; else if(colz<0) colz=0;


        return t;
    }

};

class Texture:public Object{
public:
    Texture(double floorWidth, double tileWidth){
        reference_point=Point(floorWidth, floorWidth, 0);
        length=tileWidth;
    }
    void draw(){
        glPushMatrix();
        {
            glTranslatef(-reference_point.x/2, -reference_point.y/2, 0);
            glBegin(GL_QUADS);
            for (unsigned int x=0;x<reference_point.x;x+=length){
                for (unsigned int y=0;y<reference_point.y;y+=length)
                {
                    if (((y/(int)length)+(x/(int)length)) & 1)  glColor3f(1.0f,1.0f,1.0f);
                    else glColor3f(0.0f,0.0f,0.0f);

                    glVertex3f(x, y, 0);
                    glVertex3f((x+length), y, 0);
                    glVertex3f((x+length), (y+length), 0);
                    glVertex3f(x, (y+length), 0);
                }
            }
            glEnd();
        }
        glPopMatrix();
    }

    Point getNormal(){
        return Point(0,0,1);
    }

    bool withinFloor(Point ip){
        if(ip.x>=-reference_point.x/2&&ip.x<=reference_point.x-reference_point.x/2 && ip.y>=-reference_point.y/2&&ip.y<=reference_point.y-reference_point.y/2)
            return true;
        else return false;
    }

    Point getReflection(Point dir){
        Point normal = getNormal();
        Point rf;
        rf.x = dir.x - 2*dot(dir, normal)*normal.x;
        rf.y = dir.y - 2*dot(dir, normal)*normal.y;
        rf.z = dir.z - 2*dot(dir, normal)*normal.z;

        return rf;
    }

    Point getRefrection(Point dir, double eta){
        Point normal = getNormal();
        double N_dot_I = dot(normal, dir);
        double k = 1.f - eta * eta * (1.f - N_dot_I * N_dot_I);
        if (k < 0.f) return Point(0, 0, 0);
        else return (dir*eta)-(normal*(eta * N_dot_I + sqrt(k)));
    }


    void getColor(Point ip, double &cx, double &cy, double &cz){
        unsigned i = (ip.x-reference_point.x)*texture.width()/(-2*reference_point.x);
        unsigned j = (ip.y-reference_point.y)*texture.height()/(-2*reference_point.y);
        unsigned char r, g, b;
        texture.get_pixel(i, j, r, g, b);
        cx = r/255.0; cy = g/255.0; cz = b/255.0;
        //cout<<cx<<" "<<cy<<" "<<cz<<endl;
    }

    int getColorAt(Ray* L){
        int nearest=-1;
        double t_min = 5000000;
        for(int k=0;k<objects.size();k++)
        {
            double dcx, dcy, dcz;
            double t = objects[k]->intersect(L, dcx, dcy, dcz, 0);

            if(t<=0) continue;
            if(t<t_min){
                nearest = k;
                t_min = t;
            }
        }
        return nearest;
    }

    void illumination(Point intersectionPoint, Point start, Point re, double &colx, double &coly, double &colz){
        for(unsigned i=0;i<lights.size();i++){
            Point d = lights[i]-intersectionPoint;
            d = normalize(d);
            Point s = intersectionPoint + (d*1);
            Ray *L = new Ray(s, d);

            int nearest = getColorAt(L);

            if(nearest==-1){
                Point m = intersectionPoint-reference_point;
                m = normalize(m);
                Point v = intersectionPoint-start;

                double lambert = max(dot(d,m), 0.0);
                double phong = dot(normalize(v), re);

                colx+= lambert*co_efficients[1]*color[0];
                coly+= lambert*co_efficients[1]*color[1];
                colz+= lambert*co_efficients[1]*color[2];

                colx+= pow(phong, shine)*co_efficients[2]*color[0];
                coly+= pow(phong, shine)*co_efficients[2]*color[1];
                colz+= pow(phong, shine)*co_efficients[2]*color[2];
            }
        }
    }


    double intersect(Ray *r, double &colx, double &coly, double &colz, int level){
        Point start = r->start, dir = r->dir;
        double t=(reference_point.z-start.z)/dir.z;
        Point intersectionPoint = start + (dir*t);

        if(!withinFloor(intersectionPoint)) return -1;
        if(level==0) return t;

        getColor(intersectionPoint, colx, coly, colz);

        Point reflection = getReflection(dir);
        reflection = normalize(reflection);

        illumination(intersectionPoint, start, reflection, colx, coly, colz);

        if(level<r_level){
            Point s = intersectionPoint + (reflection*1);
            Ray *refRay = new Ray(s, reflection);

            int nearest = getColorAt(refRay);

            if(nearest!=-1){
                double rfcolx, rfcoly, rfcolz;
                objects[nearest]->intersect(refRay, rfcolx, rfcoly, rfcolz, level+1);
                colx += rfcolx*co_efficients[3];
                coly += rfcoly*co_efficients[3];
                colz += rfcolz*co_efficients[3];
            }
        }

        if(colx>1) colx=1; else if(colx<0) colx=0;
        if(coly>1) coly=1; else if(coly<0) coly=0;
        if(colz>1) colz=1; else if(colz<0) colz=0;


        return t;
    }

};


class Triangle:public Object{
public:
    Point a, b, c;

    Triangle(){}
    Triangle(Point a, Point b, Point c){this->a=a; this->b=b; this->c=c;}

    void draw(){
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x, a.y, a.z);
            glVertex3f(b.x, b.y, b.z);
            glVertex3f(c.x, c.y, c.z);
        }
        glEnd();
    }

    Point getNormal(){
        Point ab = b-a, ac = c-a;
        return cross(ab, ac);
    }

    double tValue(Ray *r){
        double a,f,u,v;
        Point edge1 = b - this->a, edge2 = c - this->a, h, s, q;

        h = cross(r->dir, edge2);
        a = dot(edge1, h);
        f = 1/a;
        if (a > -EPSILON && a < EPSILON) return -1;

        f = 1/a;
        s = r->start - this->a;
        u = f * dot(s, h);
        if (u < 0.0 || u > 1.0) return -1;

        q = cross(s, edge1);
        v = f * dot(r->dir, q);
        if (v < 0.0 || u + v > 1.0) return -1;

        return f * dot(edge2, q);
    }

    Point getReflection(Point dir){
        Point normal = getNormal();
        Point rf;
        rf.x = dir.x - 2*dot(dir, normal)*normal.x;
        rf.y = dir.y - 2*dot(dir, normal)*normal.y;
        rf.z = dir.z - 2*dot(dir, normal)*normal.z;

        return rf;
    }

    Point getRefrection(Point dir, double eta){
        Point normal = getNormal();
        double N_dot_I = dot(normal, dir);
        double k = 1.f - eta * eta * (1.f - N_dot_I * N_dot_I);
        if (k < 0.f) return Point(0, 0, 0);
        else return (dir*eta)-(normal*(eta * N_dot_I + sqrt(k)));
    }



    int getColorAt(Ray* L){
        int nearest=-1;
        double t_min = 50000;
        for(int k=0;k<objects.size();k++)
        {
            double dcx, dcy, dcz;
            double t = objects[k]->intersect(L, dcx, dcy, dcz, 0);

            if(t<=0) continue;
            if(t<t_min){
                nearest = k;
                t_min = t;
            }
        }
        return nearest;
    }

    void illumination(Point intersectionPoint, Point start, Point re, double &colx, double &coly, double &colz){
        for(unsigned i=0;i<lights.size();i++){
            Point d = lights[i]-intersectionPoint;
            d = normalize(d);
            Point s = intersectionPoint + (d*1);
            Ray *L = new Ray(s, d);

            int nearest = getColorAt(L);

            if(nearest==-1){
                Point m = intersectionPoint-reference_point;
                m = normalize(m);
                Point v = intersectionPoint-start;

                double lambert = max(dot(d,m), 0.0);
                double phong = dot(normalize(v), re);

                colx+= lambert*co_efficients[1]*color[0];
                coly+= lambert*co_efficients[1]*color[1];
                colz+= lambert*co_efficients[1]*color[2];

                colx+= pow(phong, shine)*co_efficients[2]*color[0];
                coly+= pow(phong, shine)*co_efficients[2]*color[1];
                colz+= pow(phong, shine)*co_efficients[2]*color[2];
            }
        }
    }


    double intersect(Ray *r, double &colx, double &coly, double &colz, int level){
        Point start = r->start, dir = r->dir;
        double t = tValue(r);
        if(t<= EPSILON || t==-1) return -1;
        if(level==0) return t;

        Point intersectionPoint = start + (dir * t);
        colx = color[0]; coly = color[1]; colz = color[2];

        colx*=co_efficients[0];
        coly*=co_efficients[0];
        colz*=co_efficients[0];

        Point reflection = getReflection(dir);
        reflection = normalize(reflection);

        illumination(intersectionPoint, start, reflection, colx, coly, colz);

        if(level<r_level){
            Point s = intersectionPoint + (reflection*1);
            Ray *refRay = new Ray(s, reflection);

            int nearest = getColorAt(refRay);

            if(nearest!=-1){
                double rfcolx, rfcoly, rfcolz;
                objects[nearest]->intersect(refRay, rfcolx, rfcoly, rfcolz, level+1);
                colx += rfcolx*co_efficients[3];
                coly += rfcoly*co_efficients[3];
                colz += rfcolz*co_efficients[3];
            }

            Point refrection = getRefrection(dir, 1.3);
            refrection = normalize(refrection);
        }

        if(colx>1) colx=1; else if(colx<0) colx=0;
        if(coly>1) coly=1; else if(coly<0) coly=0;
        if(colz>1) colz=1; else if(colz<0) colz=0;

        return t;
    }
};

class GenQuad:public Object{
public:
    double a, b, c, d, e, f, g, h, i, j;
    GenQuad(){}
    GenQuad(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9, double p10, double p11, double p12, double p13, double p14, double p15, double p16){
        a=p1; b=p2, c=p3, d=p4, e=p5, f=p6, g=p7, h=p8, i=p9, j=p10;
        reference_point.x = p11; reference_point.x = p12; reference_point.x = p13;
        length = p14; width = p15; height = p16;
    }

    Point getNormal(Point ip){
        double x = 2*a*ip.x+d*ip.y+f*ip.z+g;
        double y = 2*b*ip.y+d*ip.x+e*ip.z+h;
        double z = 2*c*ip.z+e*ip.y+f*ip.x+i;

        return Point(x, y, z);
    }

    Point getReflection(Point dir, Point ip){
        Point normal = getNormal(ip);
        Point rf;
        rf.x = dir.x - 2*dot(dir, normal)*normal.x;
        rf.y = dir.y - 2*dot(dir, normal)*normal.y;
        rf.z = dir.z - 2*dot(dir, normal)*normal.z;

        return rf;
    }

    Point getRefrection(Point dir, Point ip, double eta){
        Point normal = getNormal(ip);
        double N_dot_I = dot(normal, dir);
        double k = 1.f - eta * eta * (1.f - N_dot_I * N_dot_I);
        if (k < 0.f) return Point(0, 0, 0);
        else return (dir*eta)-(normal*(eta * N_dot_I + sqrt(k)));
    }

    int intersectingT(Ray *r){
        Point p = r->start, v = r->dir;
        double x = a*v.x*v.x+b*v.y*v.y+c*v.z*v.z+d*v.x*v.y+e*v.y*v.z+f*v.z*v.x;
        double y = 2*a*p.x*v.x+2*b*p.y*v.y+2*c*p.z*v.z+d*p.x*v.y+d*p.y*v.x+e*p.y*v.z+e*p.z*v.y+f*p.z*v.x+f*p.x*v.z+g*v.x+h*v.y+i*v.z;
        double z = a*p.x*p.x+b*p.y*p.y+c*p.z*p.z+d*p.x*p.y+e*p.y*p.z+f*p.z*p.x+g*p.x+h*p.y+i*p.z+j;

        double dit = y*y-4*x*z;
        if(dit<0) return -1;
        dit=sqrt(dit);
        double t1 = (-y+dit)/(2*x);
        double t2 = (-y-dit)/(2*x);

        Point ip1 = p + v*t1;
        Point ip2 = p + v*t2;

        bool check1 = (length>0 && (reference_point.x > ip1.x || ip1.x - reference_point.x > length) ||
                      width>0 && (reference_point.y > ip1.y || ip1.y - reference_point.y > width) ||
                      height>0 && (reference_point.z > ip1.z || ip1.z - reference_point.z > height));

        bool check2 = (length>0 && (reference_point.x > ip2.x || ip2.x- reference_point.x > length) ||
                      width>0 && (reference_point.y > ip2.y || ip2.y - reference_point.y > width) ||
                      height>0 && (reference_point.z > ip2.z || ip2.z - reference_point.z > height));

        if (check1 && check2) return -1;
        else if (check1) return t2;
        else if (check2) return t1;
        else return t1<t2? t1:t2;
    }

    int getColorAt(Ray* L){
        int nearest=-1;
        double t_min = 5000000;
        for(int k=0;k<objects.size();k++)
        {
            double dcx, dcy, dcz;
            double t = objects[k]->intersect(L, dcx, dcy, dcz, 0);

            if(t<=0) continue;
            if(t<t_min){
                nearest = k;
                t_min = t;
            }
        }
        return nearest;
    }

    void illumination(Point intersectionPoint, Point start, Point re, double &colx, double &coly, double &colz){
        for(unsigned i=0;i<lights.size();i++){
            Point d = lights[i]-intersectionPoint;
            d = normalize(d);
            Point s = intersectionPoint + (d*1);
            Ray *L = new Ray(s, d);

            int nearest = getColorAt(L);

            if(nearest==-1){
                Point m = intersectionPoint-reference_point;
                m = normalize(m);
                Point v = intersectionPoint-start;

                double lambert = max(dot(d,m), 0.0);
                double phong = dot(normalize(v), re);

                colx+= lambert*co_efficients[1]*color[0];
                coly+= lambert*co_efficients[1]*color[1];
                colz+= lambert*co_efficients[1]*color[2];

                colx+= pow(phong, shine)*co_efficients[2]*color[0];
                coly+= pow(phong, shine)*co_efficients[2]*color[1];
                colz+= pow(phong, shine)*co_efficients[2]*color[2];
            }
        }
    }


    double intersect(Ray *r, double &colx, double &coly, double &colz, int level){
        Point start = r->start, dir = r->dir;
        double t = intersectingT(r);

        if(t<=0) return -1;
        if(level==0) return t;

        Point intersectionPoint = start + (dir*t);
        colx = color[0]; coly = color[1]; colz = color[2];

        colx*=co_efficients[0];
        coly*=co_efficients[0];
        colz*=co_efficients[0];

        Point reflection = getReflection(dir, intersectionPoint);
        reflection = normalize(reflection);

        illumination(intersectionPoint, start, reflection, colx, coly, colz);

        if(level<r_level){
            Point s = intersectionPoint + (reflection*1);
            Ray *refRay = new Ray(s, reflection);

            int nearest = getColorAt(refRay);

            if(nearest!=-1){
                double rfcolx, rfcoly, rfcolz;
                objects[nearest]->intersect(refRay, rfcolx, rfcoly, rfcolz, level+1);
                colx += rfcolx*co_efficients[3];
                coly += rfcoly*co_efficients[3];
                colz += rfcolz*co_efficients[3];
            }
        }


        if(colx>1) colx=1; else if(colx<0) colx=0;
        if(coly>1) coly=1; else if(coly<0) coly=0;
        if(colz>1) colz=1; else if(colz<0) colz=0;


        return t;
    }

};
