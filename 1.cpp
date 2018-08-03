#include "headerFile.h"

double cameraX, cameraY, cameraZ;
double cameraAngle, lx,ly,lz;
double ux, uy, uz, rx, ry, rz;
int drawaxes;
int image_width=768;
double window_height=500, window_width=500;
double view_angle=80;
double radius=30, a=50;
int r_level=5;
bitmap_image image(image_width, image_width);


void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}

void loadActualData(){
    ifstream in("scene.txt");
    Object *temp;
    string str;
    int N;
    in>>r_level>>image_width>>N;
    cout<<N<<endl;
    while(--N>=0){
        in>>str;
        if(str == "sphere"){
            double x, y, z, r;
            int shine;
            in>>x>>y>>z>>r;

            temp = new Sphere(Point(x, y, z), r);

            in>>x>>y>>z;
            temp->setColor(x, y, z);

            in>>x>>y>>z>>r;
            temp->setCoEfficients(x, y, z, r);

            in>>shine;
            temp->setShine(shine);

            objects.push_back(temp);
        }
        else if(str == "triangle"){
            double x, y, z, r;
            int shine;

            in>>x>>y>>z;
            Point p1(x, y, z);

            in>>x>>y>>z;
            Point p2(x, y, z);

            in>>x>>y>>z;
            Point p3(x, y, z);

            temp = new Triangle(p1, p2, p3);

            in>>x>>y>>z;
            temp->setColor(x, y, z);

            in>>x>>y>>z>>r;
            temp->setCoEfficients(x, y, z, r);

            in>>shine;
            temp->setShine(shine);

            objects.push_back(temp);
        }
        else if(str == "general"){
            double a, b, c, d, e, f, g, h, i, j, rx, ry, rz, length, width, height;
            int shine;

            in>>a>>b>>c>>d>>e>>f>>g>>h>>i>>j>>rx>>ry>>rz>>length>>width>>height;

            temp = new GenQuad(a, b, c, d, e, f, g, h, i, j, rx, ry, rz, length, width, height);

            in>>a>>b>>c;
            temp->setColor(a, b, c);

            in>>a>>b>>c>>d;
            temp->setCoEfficients(a, b, c, d);

            in>>shine;
            temp->setShine(shine);

            objects.push_back(temp);
        }
        else{
            int lsize;
            double x, y, z;
            in>>lsize;
            for(int i=0;i<lsize;i++){
                in>>x>>y>>z;
                lights.push_back(Point(x, y, z));
            }
        }
    }

    temp=new Floor(1000, 20);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(10);
    //objects.push_back(temp);

    temp=new Texture(1000, 20);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(10);
    objects.push_back(temp);
}

void capture(){
    double plane_distance = (window_height/2)/tan(view_angle/2);
    Point topleft;
    topleft.x = cameraX - lx*plane_distance-rx*window_width/2+ux*window_height/2;
    topleft.y = cameraY - ly*plane_distance-ry*window_width/2+uy*window_height/2;
    topleft.z = cameraZ - lz*plane_distance-rz*window_width/2+uz*window_height/2;
    double du=window_width/image_width;
    double dv=du;

    for(int i=0;i<image_width;i++){
        for(int j=0;j<image_width;j++){
            Point corner;
            corner.x = topleft.x+i*rx*du-j*ux*dv;
            corner.y = topleft.y+i*ry*du-j*uy*dv;
            corner.z = topleft.z+i*rz*du-j*uz*dv;

            Ray* ray = new Ray(Point(cameraX,cameraY,cameraZ), Point(corner.x-cameraX, corner.y-cameraY, corner.z-cameraZ));
            int nearest=-1;
            double t_min = 50000;
            for(int k=0;k<objects.size();k++)
            {
                double dcolx, dcoly, dcolz;
                double t =objects[k]->intersect(ray, dcolx, dcoly, dcolz, 0);

                if(t<=0) continue;
                if(t<t_min){
                    nearest = k;
                    t_min = t;
                }
            }

            if(nearest!=-1){
                double colx, coly, colz;
                objects[nearest]->intersect(ray, colx, coly, colz, 1);
                image.set_pixel(i,j,(int)(255*colx),(int)(255*coly),(int)(255*colz));
            }
        }
    }
    cout<<"ok"<<endl;
    image.save_image("test.bmp");
}


Point f3(Point vect, Point perp, int dir)
{
    double c = cos(pi/180);
    double s = dir * sin(pi/180);
    Point point;
    point.x = c * vect.x + s * perp.x;
    point.y = c * vect.y + s * perp.y;
    point.z = c * vect.z + s * perp.z;
    c = sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
    point.x /= c;
    point.y /= c;
    point.z /= c;
    return point;
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
	    case '0':
            capture();
            break;
		case '1': {
            Point l1 = f3(Point(lx,ly,lz), Point(rx,ry,rz), -1);
            Point r = f3(Point(rx,ry,rz), Point(lx,ly,lz), 1); rx=r.x;ry=r.y;rz=r.z;
            lx = l1.x;ly=l1.y;lz=l1.z;
            break;
		}
		case '2': {
            Point l1 = f3(Point(lx,ly,lz), Point(rx,ry,rz), 1);
            Point r = f3(Point(rx,ry,rz), Point(lx,ly,lz), -1); rx=r.x;ry=r.y;rz=r.z;
            lx = l1.x;ly=l1.y;lz=l1.z;
            break;
            break;
		}

		case '3': {
            Point u1 = f3(Point(ux,uy,uz), Point(lx,ly,lz), -1);
            Point l = f3(Point(lx,ly,lz), Point(ux,uy,uz), 1);
            lx = l.x; ly = l.y; lz = l.z;
            ux = u1.x;uy=u1.y;uz=u1.z;
            break;
        }
        case '4': {
            Point u1 = f3(Point(ux,uy,uz), Point(lx,ly,lz), 1);
            Point l = f3(Point(lx,ly,lz), Point(ux,uy,uz), -1);
            lx = l.x; ly = l.y; lz = l.z;
            ux = u1.x;uy=u1.y;uz=u1.z;
            break;
        }
        case '5': {
            Point r1 = f3(Point(rx,ry,rz), Point(ux,uy,uz), -1);
            Point u = f3(Point(ux,uy,uz), Point(rx,ry,rz), 1);
            ux = u.x; uy = u.y; uz = u.z;
            rx = r1.x;ry=r1.y;rz=r1.z;
            break;
        }
        case '6':{
            Point r1 = f3(Point(rx,ry,rz), Point(ux,uy,uz), 1);
            Point u = f3(Point(ux,uy,uz), Point(rx,ry,rz), -1);
            ux = u.x; uy = u.y; uz = u.z;
            rx = r1.x;ry=r1.y;rz=r1.z;
            break;
        }
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
    double x1,x2;
	switch(key){
	    case GLUT_KEY_UP:		// up arrow key
			cameraX+=lx;
			cameraY+=ly;
			cameraZ+=lz;
			break;
		case GLUT_KEY_DOWN:		//down arrow key
			cameraX-=lx;
			cameraY-=ly;
			cameraZ-=lz;
			break;

		case GLUT_KEY_RIGHT:
			cameraX+=rx;
			cameraY+=ry;
			cameraZ+=rz;
			break;
		case GLUT_KEY_LEFT:
		    cameraX-=rx;
			cameraY-=ry;
			cameraZ-=rz;
			break;

		case GLUT_KEY_PAGE_UP:
		    cameraX+=ux;
			cameraY+=uy;
			cameraZ+=uz;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    cameraX-=ux;
			cameraY-=uy;
			cameraZ-=uz;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
		    if (radius<a) radius++;
			break;
		case GLUT_KEY_END:
		    if (radius>0) radius--;
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}

void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();

	gluLookAt(cameraX, cameraY, cameraZ,	cameraX+lx, cameraY+ly, cameraZ+lz,		ux,uy,uz);

	glMatrixMode(GL_MODELVIEW);

	drawAxes();
	//drawSCC();

	for(unsigned i=0;i<objects.size();i++) objects[i]->draw();

	glColor3f(0,0,1);
	glPointSize(10);
	glBegin(GL_POINTS);
	{
        for(unsigned i=0;i<lights.size();i++)
            glVertex3f(lights[i].x, lights[i].y, lights[i].z);
	}
	glEnd();

	glutSwapBuffers();
}

void animate(){
    glutPostRedisplay();
}

void init(){
	drawaxes=0;
	//cameraX = 150, cameraY = 150, cameraZ = 50;
	cameraX = 0, cameraY = -200, cameraZ = 10;
	//lx=-1.0/sqrt(2),ly=-1.0/sqrt(2),lz=0;
	ux=0, uy=0, uz=1;
	//rx= -1/sqrt(2), ry=1/sqrt(2), rz=0;
	lx=0, ly=1, lz=0;
	rx=1, ry=0, rz=0;

	glClearColor(0,0,0,0);

	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(view_angle,	1,	1,	1000.0);
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(window_height, window_width);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();
	//loadTestData();
	loadActualData();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
