#include <windows.h>		// Header File For Windows
#include <fstream>
#include <iostream>
#include <gl.h>			// Header File For The OpenGL32 Library
#include <glu.h>			// Header File For The GLu32 Library
#include <glut.h>
#include <glaux.h>		// Header File For The Glaux Library
#include<texture.h>
#include "math.h"												    // NEW: Needed For Sqrtf
#include "Model_3DS.h"	

using namespace std;
bool        isClicked = false;										// NEW: Clicking The Mouse?
bool        isRClicked = false;
bool        isDragging = false;
int mouseX = 0, mouseY = 0;

#pragma comment(lib,"opengl32.lib")
#pragma comment(lib,"glut32.lib")
#pragma comment(lib,"glu32.lib")
#pragma comment(lib,"GLAUX.LIB")


HDC			hDC = NULL;		// Private GDI Device Context
HGLRC		hRC = NULL;		// Permanent Rendering Context
HWND		hWnd = NULL;		// Holds Our Window Handle
HINSTANCE	hInstance;		// Holds The Instance Of The Application

bool	keys[256];			// Array Used For The Keyboard Routine
bool	active = TRUE;		// Window Active Flag Set To TRUE By Default
bool	fullscreen = TRUE;	// Fullscreen Flag Set To Fullscreen Mode By Default

LRESULT	CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);	// Declaration For WndProc


GLvoid ReSizeGLScene(GLsizei width, GLsizei height)		// Resize And Initialize The GL Window
{
	if (height == 0)										// Prevent A Divide By Zero By
	{
		height = 1;										// Making Height Equal One
	}

	glViewport(0, 0, width, height);						// Reset The Current Viewport

	glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
	glLoadIdentity();									// Reset The Projection Matrix

	// Calculate The Aspect Ratio Of The Window
	gluPerspective(45.0f, (GLfloat)width / (GLfloat)height, 0.1f, 1000.0f);

	glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
	glLoadIdentity();									// Reset The Modelview Matrix
}
float eyeX = 0, eyeY = 0, eyeZ = -10, angel = -1.5;
double xS = 10.0, yS = 0,zS=-15.4;
GLfloat k = 0;
void mouse(int mouseX, int mouseY, bool isClicked, bool isRClicked)
{
	if (isClicked) angel -= 0.0005;;
	if (isRClicked) angel += 0.0005;;
}


double frictionFactor = 0.4;
double frameTime = 1.0 / 60.0;
double width_table= 24.0;
double high_table= 12.0;
double startTableX = -12;
double startTableY = -6;
class Vector
{
public:
	double x, y, z=-15.4;
	Vector(){
		this->z = -15.4;
	}
	Vector(double x, double y){
		this->x = x;
		this->y = y;
		z = -15.4;
	}
	Vector(double x , double y , double z ){
		this->x = x;
		this->y = y;
		this->z = z;
	}
	Vector mul(double d){
		Vector c;
		c.x = x * d;
		c.y = y * d;
		c.z = z ;
		return c;
	}


	Vector sum(Vector v){
		Vector c;
		c.x = x + v.x;
		c.y = y + v.y;
		c.z = z;
		return c;
	}

	Vector sub(Vector v){

		Vector c;
		c.x = x - v.x;
		c.y = y - v.y;
		c.z = z;
		return c;
	}

	double dot(Vector v ){
		return (x*v.x)+(y*v.y);
	}
	Vector cross(Vector vect)
	{
		Vector cross_P;
		cross_P.x = y * vect.z - z * vect.y;
		cross_P.y = z * vect.x- x* vect.z;
		cross_P.z = x * vect.y - y * vect.x;
		return cross_P;
	}
	double value(){
		return sqrt(x*x + y*y );
	}
};
Vector gools[6];
double radiusGool = 1.0;
void goolPosition(){

	gools[0].x = startTableX + radiusGool;
	gools[0].y = startTableY + radiusGool;
	gools[0].z = -16.0;

	gools[1].x = startTableX + width_table/2;
	gools[1].y = startTableY + radiusGool;
	gools[1].z = -16.0;

	gools[2].x = startTableX + width_table-radiusGool;
	gools[2].y = startTableY + radiusGool;
	gools[2].z = -16.0;

	gools[3].x = startTableX + radiusGool;
	gools[3].y = startTableY + high_table-radiusGool;
	gools[3].z = -16.0;

	gools[4].x = startTableX + width_table / 2;
	gools[4].y = startTableY + high_table - radiusGool;
	gools[4].z = -16.0;

	gools[5].x = startTableX + width_table - radiusGool;
	gools[5].y = startTableY + high_table - radiusGool;
	gools[5].z = -16.0;
}
class Ball
{
public:
	Vector position, speed, accelerate, inplusForce, frictionForce,alpha,omiga, theta,inplus_arm,friction_arm;
	double mass, radius,tita;
	bool is_running=true ,is_hit=false,rotate=false;
	Ball(){
		mass = 0.2;
		radius = 0.5;
		speed.x = 0;
		speed.y = 0;
		speed.z = -15.4;

		accelerate.x = 0;
		accelerate.y = 0;
		accelerate.z = -15.4;
	}
	Ball(double x, double y){

		mass = 0.2;
		radius = 0.5;
		position.x = x;
		position.y = y;
		position.z = -15.4;
	}
	Ball(double x, double y, double z,double m,double r){
		mass = m;
		radius = r;
		position.x = x;
		position.y = y;
		position.z = z;
	}
	void forces(Vector diraction, double value){
		Vector norm;
		double f = mass*9.8*frictionFactor;
		double dist = diraction.value();
		norm = diraction.mul(1.0 / dist);
		inplusForce = diraction.mul(value);
		frictionForce = norm.mul(-1.0*f);
	}

	void start(Vector v, double d, Vector eff_point){
		if (mass*9.8*frictionFactor < d){
			is_hit = true;
			forces(v, d);
			Vector segma = inplusForce.sum(frictionForce);
			accelerate = segma.mul(1.0 / mass);
			speed = accelerate.mul(frameTime);
			position = position.sum(speed.mul(frameTime));

			cout << "frictionForce:" << frictionForce.x << "      ,      " << frictionForce.y << "      ,      " << frictionForce.z << endl;
			cout << "accelerate:" << accelerate.x << "      ,      " << accelerate.y << "      ,      " << accelerate.z << endl;
			cout << "speed:" << speed.x << "      ,      " << speed.y << "      ,      " << speed.z << endl;
			cout << position.x << "      ,      " << position.y << "    ,     " << position.z << endl;

			if (is_rotate(v))
			{
				inplus_arm = position.sub(eff_point);
				inplus_arm = inplus_arm.mul(radius / inplus_arm.value());
				Vector fric_point(position.x,position.y,-16);
				friction_arm = position.sub(fric_point);
				alpha = ((frictionForce.cross(friction_arm)).sum(inplusForce.cross(inplus_arm))).mul(5 / (2 * mass*pow(radius, 2)));
				omiga = alpha.mul(frameTime);
				theta = omiga.mul(frameTime);
				tita = theta.value();
			}
		}
	}

	void update(){
		if (is_running && is_hit){		
			Vector norm;
			double f = mass*9.8*frictionFactor;
			double dist = speed.value();
			norm = speed.mul(1.0 / dist);
			frictionForce = norm.mul(-1.0*f);
			accelerate = frictionForce.mul(1.0 / mass);
			speed = speed.sum(accelerate.mul(frameTime));

			shock_table();

			cout << "accelerate:" << accelerate.x << "      ,      " << accelerate.y << "      ,      " << accelerate.z << endl;
			cout << "speed:" << speed.x << "      ,      " << speed.y << "      ,      " << speed.z << endl;
			cout << position.x << "      ,      " << position.y << "    ,     " << position.z << endl;
		    
			if (rotate)
			{
				Vector fric_point(position.x, position.y, -16);
				friction_arm = position.sub(fric_point);
				alpha = (frictionForce.cross(friction_arm)).mul(5/(2*mass*pow(radius,2)));
				omiga = omiga.sum(alpha.mul(frameTime));
				theta =theta.sum( omiga.mul(frameTime));

			}

			if (isInGool())
			{
				position.z = -50.0;
				is_running = false;
			}
			if (stop()){
			omiga.x = 0.0;
			omiga.y = 0.0;
			rotate = false;
			}
		}
	}
	bool is_rotate(Vector v)
	{
		if (v.z == position.z && v.y == position.y)
		{
			rotate = false;
		}
		else
		{
			rotate = true;
		}
		return rotate;
	}
	void shock_table()
	{
		
		if (position.x <= startTableX +radius )
		{
			position.x = startTableX + radius;
			speed.x *= -1;
		}
		if (position.x >= startTableX +width_table- radius)
		{
			position.x = startTableX + width_table - radius;
			speed.x *= -1;
		}
		if (position.y <= startTableY  + radius)
		{
			position.y = startTableY  + radius;
			speed.y *= -1;
			
		}
		if (position.y >= startTableY + high_table - radius)
		{
			position.y = startTableY + high_table - radius;
			speed.y *= -1;
		}
	}
	bool stop(){
		if (abs(speed.x) <= abs(accelerate.mul(frameTime).x) && abs(speed.y) <= abs(accelerate.mul(frameTime).y))
		{
			return true;
		}
		return false;
	}
	bool isInGool(){
		for (int i = 0; i < 6; i++)
		{
			goolPosition();
			double dist = sqrt((pow(position.x - gools[i].x, 2)) + (pow(position.y - gools[i].y, 2)));
			if (dist <= radiusGool )
			{
				return true;
			}
		}
		return false;
	}
};

bool is_shock(Ball ball1, Ball ball2){
	double dist = ((ball1.position.x - ball2.position.x)*(ball1.position.x - ball2.position.x)) + ((ball1.position.y - ball2.position.y)*(ball1.position.y - ball2.position.y));
	if (ball1.radius + ball2.radius >= sqrt(dist) /*distance(ball1.position, ball2.position)*/)
	{
		cout << "sssssssssssssssssssssssssssssssssssssssssssssssssss" << endl;
		return true;
	}
	else{
		return false;
	}
}

void move_after_shock(Ball *ball1, Ball *ball2)
{
	if (ball1->position.x == ball2->position.x && ball1->position.y != ball2->position.y || ball1->position.y == ball2->position.y && ball1->position.x != ball2->position.x){
		if ( ball1->is_hit == true && ball1->is_running==true &&  ball2->is_hit == false/* ball2->speed.value()<1*/){
			Vector temp(ball1->speed.x, ball1->speed.y);
			ball1->is_hit = false;
			ball1->speed.y = 0.000001;
			ball1->speed.x = 0.000001;
			ball2->speed.x = temp.x;
			ball2->speed.y = temp.y;
			ball2->is_hit = true;
			ball2->is_running = true;
			
		}
		else if (ball2->is_hit == true && ball2->is_running == true &&  ball1->is_hit == false /*ball1->speed.value()<1*/){
			Vector temp(ball2->speed.x, ball2->speed.y);
			ball2->is_hit = false;
			ball2->speed.y = 0.000001;
			ball2->speed.x = 0.000001;
			ball1->speed.x = temp.x;
			ball1->speed.y = temp.y;
			ball1->is_hit = true;
			ball1->is_running = true; 
		
		}
	}
	else{
		double d = ((ball1->position.x - ball2->position.x)*(ball1->position.x - ball2->position.x)) + ((ball1->position.y - ball2->position.y)*(ball1->position.y - ball2->position.y));
		double dist = ball1->radius + ball2->radius - sqrt(d);
		if (dist != 0)
		{
			Vector unit = ball1->position.sub(ball2->position);
			unit = unit.mul(1.0 / unit.value());

			ball1->position.x += unit.mul(dist / 2.0).x;
			ball1->position.y += unit.mul(dist / 2.0).y;

			ball2->position.x -= unit.mul(dist / 2.0).x;
			ball2->position.y -= unit.mul(dist / 2.0).y;
		}
		
		ball1->is_hit = true;
		ball2->is_hit = true;

		ball1->rotate = true;
		ball2->rotate = true;
		Vector n(ball1->position.x - ball2->position.x, ball1->position.y - ball2->position.y);
		Vector un = n.mul(1.0 / n.value());
		Vector ut(-1 * un.y, un.x);
		double v1n = un.dot(ball1->speed);
		double v1t = ut.dot(ball1->speed);
		double v2n = un.dot(ball2->speed);
		double v2t = ut.dot(ball2->speed);
		double new_v1n_val = (v1n*(ball1->mass - ball2->mass) + v2n * 2 * ball2->mass) / (ball1->mass + ball2->mass);
		double new_v2n_val = (v2n*(ball2->mass - ball1->mass) + v1n * 2 * ball1->mass) / (ball1->mass + ball2->mass);
		Vector new_v1n_vec = un.mul(new_v1n_val);
		Vector new_v1t_vec = ut.mul(v1t);
		Vector new_v2n_vec = un.mul(new_v2n_val);
		Vector new_v2t_vec = ut.mul(v2t);
		cout << "//////////////////////////////////////";
		cout << ball1->speed.x << "\n" << ball1->speed.y << "\n" << ball1->speed.z;
		cout << "///////////";
		cout << ball2->speed.x << "\n" << ball2->speed.y << "\n" << ball2->speed.z;
		cout << "//////////////////////////////////";
		ball1->speed = new_v1n_vec.sum(new_v1t_vec);
		ball2->speed = new_v2n_vec.sum(new_v2t_vec);

		Vector eff_point(un.mul(ball1->radius).x, un.mul(ball1->radius).y, un.mul(ball1->radius).z);
		ball1->inplus_arm = ball1->position.sub(eff_point);
		ball1->inplus_arm = ball1->inplus_arm.mul(ball1->radius / ball1->inplus_arm.value());
		Vector fric_point(ball1->position.x, ball1->position.y, -16);
		ball1->friction_arm = ball1->position.sub(fric_point);
		ball1->alpha = ((ball1->frictionForce.cross(ball1->friction_arm)).sum(ball1->speed.cross(ball1->inplus_arm))).mul(5 / (2 * ball1->mass*pow(ball1->radius, 2)));
		ball1->omiga = ball1->alpha.mul(frameTime);
		ball1->theta = ball1->omiga.mul(frameTime);
		ball1->tita = ball1->theta.value();

		Vector eff_pointt(un.mul(ball2->radius).x, un.mul(ball2->radius).y, un.mul(ball2->radius).z);
		ball2->inplus_arm = ball2->position.sub(eff_pointt);
		ball2->inplus_arm = ball2->inplus_arm.mul(ball2->radius / ball2->inplus_arm.value());
		Vector fric_pointt(ball2->position.x, ball2->position.y, -16);
		ball2->friction_arm = ball2->position.sub(fric_pointt);
		ball2->alpha = ((ball2->frictionForce.cross(ball2->friction_arm)).sum(ball2->speed.cross(ball2->inplus_arm))).mul(5 / (2 * ball2->mass*pow(ball2->radius, 2)));
		ball2->omiga = ball2->alpha.mul(frameTime);
		ball2->theta = ball2->omiga.mul(frameTime);
		ball2->tita = ball2->theta.value();
	}
}

int b[8], t;

GLUquadric *quadric = gluNewQuadric();

Ball balls[8];
int numberOfBalls;

Vector v;

double d;

int InitGL(GLvoid)										// All Setup For OpenGL Goes Here
{
	glShadeModel(GL_SMOOTH);							// Enable Smooth Shading
	glClearColor(0.0f, 0.0f, 0.0f, 0.50f);				// Black Background
	glClearDepth(1.0f);									// Depth Buffer Setup
	glEnable(GL_DEPTH_TEST);							// Enables Depth Testing
	glDepthFunc(GL_LEQUAL);								// The Type Of Depth Testing To Do
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really Nice Perspective Calculations

	glEnable(GL_TEXTURE_2D);
	b[0] = LoadTexture("data/ball/1.bmp", 255); 
	b[1] = LoadTexture("data/ball/2.bmp", 255); 
	b[2] = LoadTexture("data/ball/3.bmp", 255); 
	b[3] = LoadTexture("data/ball/4.bmp", 255); 
	b[4] = LoadTexture("data/ball/5.bmp", 255);
	t = LoadTexture("data/table.bmp", 255);
	b[5] = LoadTexture("data/ball/6.bmp", 255);
	b[6] = LoadTexture("data/ball/7.bmp", 255);
	b[7] = LoadTexture("data/ball/9.bmp", 255);
	 
	balls[0].position.x = 8 ;
	balls[0].position.y = 0;
	balls[0].position.z = -15.4;

 /*صفة الكرات الفردية*/
double t = startTableX + (width_table / 2);
	for (int i = 1; i < numberOfBalls; i++){
		if (i % 2 != 0){
			balls[i].position.x = t;
			balls[1].position.y = startTableY + (high_table / 2);
			balls[1].position.z = -15.4;
			t-=1.3;
		}
	}
	/*صفة الكرات الزوجية*/
	double s=(startTableX + (width_table / 2)) - 1.3;
	double r = (startTableY + (high_table / 2));
	double q = 1.2;
	int l=2;
	while (l <= numberOfBalls){
			balls[l].position.x = s;
			balls[l].position.y = r+q;
			balls[l].position.z = -15.4;
			l += 2;
			balls[l].position.x = s;
			balls[l].position.y = r - q;
			balls[l].position.z = -15.4;
			l += 2;
			s --;
			q+=0.8;
	}

	return TRUE;										// Initialization Went OK
}

void Key(bool* keys)
{
	float eyeX1 = eyeX;
	float eyeY1 = eyeY;
	float eyeZ1 = eyeZ;

	if (keys[VK_UP])
	{

		eyeX1 += cos(angel)*0.1;
		eyeZ1 += sin(angel)*0.1;
		eyeX = eyeX1;
		eyeZ = eyeZ1;

	}

	if (keys[VK_DOWN])
	{
		eyeX1 -= cos(angel)*0.1;
		eyeZ1 -= sin(angel)*0.1;
		eyeX = eyeX1;
		eyeZ = eyeZ1;
	}

	if (keys[VK_LEFT])
		angel -= 0.0005;

	if (keys[VK_RIGHT])
		angel += 0.0005;
	if (keys['O'])
	{
		eyeY1 += 0.5;
		eyeY = eyeY1;
	}
	if (keys['L'])
	{
		eyeY1 -= 0.5;
		eyeY = eyeY1;

	}
	if (keys['8'])
	{
		yS += 0.02;
	}
	if (keys['2'])
	{
		yS -= 0.02;
	}

	if (keys['4'])
	{
		xS += 0.02;
	}

	if (keys['6'])
	{
		xS -= 0.02;
	}
	if (keys['1'])
	{
		zS -= 0.02;
	}
	if (keys['3'])
	{
		zS += 0.02;
	}
	if (keys['5'])
	{
		Vector eff_point;
		eff_point.x = xS;
		eff_point.y = yS;
		eff_point.z = -15.4;
		v.x = xS - 12;
		v.y = yS;
		v.z = -15.4;
		balls[0].start(v, d, eff_point);
	}

}

int   DrawGLScene(GLvoid)									// Here's Where We Do All The Drawing
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer
	glLoadIdentity(); // Reset The Current Modelview Matrix
	Key(keys);
	//gluLookAt(eyeX, eyeY, eyeZ, eyeX + cos(angel), eyeY, eyeZ + sin(angel), 0, 1, 0);
	glEnable(GL_TEXTURE_2D);
	mouse(mouseX, mouseY, isClicked, isRClicked);
	
	//table
	glBindTexture(GL_TEXTURE_2D, t);
		glBegin(GL_QUADS);
		
		glVertex3f(startTableX, startTableY, -16);
		glTexCoord2f(1, 0);
		glVertex3f(startTableX, startTableY + high_table, -16);
		glTexCoord2f(0, 0);
		glVertex3f(startTableX +width_table , startTableY + high_table, -16);
		glTexCoord2f(0, 1);
		glVertex3f(startTableX + width_table, startTableY, -16);
		glTexCoord2f(1, 1);
		glEnd();
		//stick
		glPushMatrix();

		glBegin(GL_LINES);
		glColor3d(1, 1, 1);
		glVertex3d(12, 0, -15.4);
		glColor3d(1, 1, 1);
		glVertex3d(xS, yS, zS);
		glEnd();
	
		glPopMatrix();
	//balls
	for (int i = 0; i < numberOfBalls; i++)
	{
		glPushMatrix();
		
		glBindTexture(GL_TEXTURE_2D, b[i]);
		glTranslatef(balls[i].position.x, balls[i].position.y, balls[i].position.z);
		balls[i].position.x += balls[i].speed.x*frameTime;
		balls[i].position.y += balls[i].speed.y*frameTime;

	
			glRotatef(balls[i].tita, balls[i].position.x, balls[i].position.y, balls[i].position.z);
			balls[i].tita = balls[i].theta.value();
		
		
		for (int j = i+1; j < numberOfBalls; j++)
		{	
			if (is_shock(balls[i], balls[j]))
			{
				move_after_shock(&balls[i], &balls[j]);
			}
		}
		balls[i].update();
		gluSphere(quadric, 0.5, 10, 10);
		
		glPopMatrix();
		
	}

	return TRUE;
}

GLvoid KillGLWindow(GLvoid)								// Properly Kill The Window
{
	if (fullscreen)										// Are We In Fullscreen Mode?
	{
		ChangeDisplaySettings(NULL, 0);					// If So Switch Back To The Desktop
		ShowCursor(TRUE);								// Show Mouse Pointer
	}

	if (hRC)											// Do We Have A Rendering Context?
	{
		if (!wglMakeCurrent(NULL, NULL))					// Are We Able To Release The DC And RC Contexts?
		{
			MessageBox(NULL, "Release Of DC And RC Failed.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
		}

		if (!wglDeleteContext(hRC))						// Are We Able To Delete The RC?
		{
			MessageBox(NULL, "Release Rendering Context Failed.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
		}
		hRC = NULL;										// Set RC To NULL
	}

	if (hDC && !ReleaseDC(hWnd, hDC))					// Are We Able To Release The DC
	{
		MessageBox(NULL, "Release Device Context Failed.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
		hDC = NULL;										// Set DC To NULL
	}

	if (hWnd && !DestroyWindow(hWnd))					// Are We Able To Destroy The Window?
	{
		MessageBox(NULL, "Could Not Release hWnd.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
		hWnd = NULL;										// Set hWnd To NULL
	}

	if (!UnregisterClass("OpenGL", hInstance))			// Are We Able To Unregister Class
	{
		MessageBox(NULL, "Could Not Unregister Class.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
		hInstance = NULL;									// Set hInstance To NULL
	}
}

BOOL CreateGLWindow(char* title, int width, int height, int bits, bool fullscreenflag)
{
	GLuint		PixelFormat;			// Holds The Results After Searching For A Match
	WNDCLASS	wc;						// Windows Class Structure
	DWORD		dwExStyle;				// Window Extended Style
	DWORD		dwStyle;				// Window Style
	RECT		WindowRect;				// Grabs Rectangle Upper Left / Lower Right Values
	WindowRect.left = (long)0;			// Set Left Value To 0
	WindowRect.right = (long)width;		// Set Right Value To Requested Width
	WindowRect.top = (long)0;				// Set Top Value To 0
	WindowRect.bottom = (long)height;		// Set Bottom Value To Requested Height

	fullscreen = fullscreenflag;			// Set The Global Fullscreen Flag

	hInstance = GetModuleHandle(NULL);				// Grab An Instance For Our Window
	wc.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;	// Redraw On Size, And Own DC For Window.
	wc.lpfnWndProc = (WNDPROC)WndProc;					// WndProc Handles Messages
	wc.cbClsExtra = 0;									// No Extra Window Data
	wc.cbWndExtra = 0;									// No Extra Window Data
	wc.hInstance = hInstance;							// Set The Instance
	wc.hIcon = LoadIcon(NULL, IDI_WINLOGO);			// Load The Default Icon
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);			// Load The Arrow Pointer
	wc.hbrBackground = NULL;									// No Background Required For GL
	wc.lpszMenuName = NULL;									// We Don't Want A Menu
	wc.lpszClassName = "OpenGL";								// Set The Class Name

	if (!RegisterClass(&wc))									// Attempt To Register The Window Class
	{
		MessageBox(NULL, "Failed To Register The Window Class.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;											// Return FALSE
	}

	if (fullscreen)												// Attempt Fullscreen Mode?
	{
		DEVMODE dmScreenSettings;								// Device Mode
		memset(&dmScreenSettings, 0, sizeof(dmScreenSettings));	// Makes Sure Memory's Cleared
		dmScreenSettings.dmSize = sizeof(dmScreenSettings);		// Size Of The Devmode Structure
		dmScreenSettings.dmPelsWidth = width;				// Selected Screen Width
		dmScreenSettings.dmPelsHeight = height;				// Selected Screen Height
		dmScreenSettings.dmBitsPerPel = bits;					// Selected Bits Per Pixel
		dmScreenSettings.dmFields = DM_BITSPERPEL | DM_PELSWIDTH | DM_PELSHEIGHT;

		// Try To Set Selected Mode And Get Results.  NOTE: CDS_FULLSCREEN Gets Rid Of Start Bar.
		if (ChangeDisplaySettings(&dmScreenSettings, CDS_FULLSCREEN) != DISP_CHANGE_SUCCESSFUL)
		{
			// If The Mode Fails, Offer Two Options.  Quit Or Use Windowed Mode.
			if (MessageBox(NULL, "The Requested Fullscreen Mode Is Not Supported By\nYour Video Card. Use Windowed Mode Instead?", "GL template", MB_YESNO | MB_ICONEXCLAMATION) == IDYES)
			{
				fullscreen = FALSE;		// Windowed Mode Selected.  Fullscreen = FALSE
			}
			else
			{
				// Pop Up A Message Box Letting User Know The Program Is Closing.
				MessageBox(NULL, "Program Will Now Close.", "ERROR", MB_OK | MB_ICONSTOP);
				return FALSE;									// Return FALSE
			}
		}
	}

	if (fullscreen)												// Are We Still In Fullscreen Mode?
	{
		dwExStyle = WS_EX_APPWINDOW;								// Window Extended Style
		dwStyle = WS_POPUP;										// Windows Style
		ShowCursor(FALSE);										// Hide Mouse Pointer
	}
	else
	{
		dwExStyle = WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;			// Window Extended Style
		dwStyle = WS_OVERLAPPEDWINDOW;							// Windows Style
	}

	AdjustWindowRectEx(&WindowRect, dwStyle, FALSE, dwExStyle);		// Adjust Window To True Requested Size

	// Create The Window
	if (!(hWnd = CreateWindowEx(dwExStyle,							// Extended Style For The Window
		"OpenGL",							// Class Name
		title,								// Window Title
		dwStyle |							// Defined Window Style
		WS_CLIPSIBLINGS |					// Required Window Style
		WS_CLIPCHILDREN,					// Required Window Style
		0, 0,								// Window Position
		WindowRect.right - WindowRect.left,	// Calculate Window Width
		WindowRect.bottom - WindowRect.top,	// Calculate Window Height
		NULL,								// No Parent Window
		NULL,								// No Menu
		hInstance,							// Instance
		NULL)))								// Dont Pass Anything To WM_CREATE
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL, "Window Creation Error.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	static	PIXELFORMATDESCRIPTOR pfd =				// pfd Tells Windows How We Want Things To Be
	{
		sizeof(PIXELFORMATDESCRIPTOR),				// Size Of This Pixel Format Descriptor
		1,											// Version Number
		PFD_DRAW_TO_WINDOW |						// Format Must Support Window
		PFD_SUPPORT_OPENGL |						// Format Must Support OpenGL
		PFD_DOUBLEBUFFER,							// Must Support Double Buffering
		PFD_TYPE_RGBA,								// Request An RGBA Format
		bits,										// Select Our Color Depth
		0, 0, 0, 0, 0, 0,							// Color Bits Ignored
		0,											// No Alpha Buffer
		0,											// Shift Bit Ignored
		0,											// No Accumulation Buffer
		0, 0, 0, 0,									// Accumulation Bits Ignored
		16,											// 16Bit Z-Buffer (Depth Buffer)  
		0,											// No Stencil Buffer
		0,											// No Auxiliary Buffer
		PFD_MAIN_PLANE,								// Main Drawing Layer
		0,											// Reserved
		0, 0, 0										// Layer Masks Ignored
	};

	if (!(hDC = GetDC(hWnd)))							// Did We Get A Device Context?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL, "Can't Create A GL Device Context.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if (!(PixelFormat = ChoosePixelFormat(hDC, &pfd)))	// Did Windows Find A Matching Pixel Format?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL, "Can't Find A Suitable PixelFormat.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if (!SetPixelFormat(hDC, PixelFormat, &pfd))		// Are We Able To Set The Pixel Format?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL, "Can't Set The PixelFormat.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if (!(hRC = wglCreateContext(hDC)))				// Are We Able To Get A Rendering Context?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL, "Can't Create A GL Rendering Context.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if (!wglMakeCurrent(hDC, hRC))					// Try To Activate The Rendering Context
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL, "Can't Activate The GL Rendering Context.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	ShowWindow(hWnd, SW_SHOW);						// Show The Window
	SetForegroundWindow(hWnd);						// Slightly Higher Priority
	SetFocus(hWnd);									// Sets Keyboard Focus To The Window
	ReSizeGLScene(width, height);					// Set Up Our Perspective GL Screen

	if (!InitGL())									// Initialize Our Newly Created GL Window
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL, "Initialization Failed.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	return TRUE;									// Success
}

LRESULT CALLBACK WndProc(HWND	hWnd,			// Handle For This Window
	UINT	uMsg,			// Message For This Window
	WPARAM	wParam,			// Additional Message Information
	LPARAM	lParam)			// Additional Message Information
{
	static PAINTSTRUCT ps;
	switch (uMsg)									// Check For Windows Messages
	{
	case WM_ACTIVATE:							// Watch For Window Activate Message
	{
													if (!HIWORD(wParam))					// Check Minimization State
													{
														active = TRUE;						// Program Is Active
													}
													else
													{
														active = FALSE;						// Program Is No Longer Active
													}

													return 0;								// Return To The Message Loop
	}

	case WM_SYSCOMMAND:							// Intercept System Commands
	{
													switch (wParam)							// Check System Calls
													{
													case SC_SCREENSAVE:					// Screensaver Trying To Start?
													case SC_MONITORPOWER:				// Monitor Trying To Enter Powersave?
														return 0;							// Prevent From Happening
													}
													break;									// Exit
	}

	case WM_CLOSE:								// Did We Receive A Close Message?
	{
													PostQuitMessage(0);						// Send A Quit Message
													return 0;								// Jump Back
	}

	case WM_KEYDOWN:							// Is A Key Being Held Down?
	{
													keys[wParam] = TRUE;					// If So, Mark It As TRUE
													return 0;								// Jump Back
	}

	case WM_KEYUP:								// Has A Key Been Released?
	{
													keys[wParam] = FALSE;					// If So, Mark It As FALSE
													return 0;								// Jump Back
	}

	case WM_SIZE:								// Resize The OpenGL Window
	{
													ReSizeGLScene(LOWORD(lParam), HIWORD(lParam));  // LoWord=Width, HiWord=Height
													return 0;								// Jump Back
	}
	case WM_MOUSEMOVE:
	{
						 mouseX = (int)LOWORD(lParam);
						 mouseY = (int)HIWORD(lParam);
						 isClicked = (LOWORD(wParam) & MK_LBUTTON) ? true : false;
						 isRClicked = (LOWORD(wParam) & MK_RBUTTON) ? true : false;
						 break;
	}
	case WM_LBUTTONUP:
		isClicked = false; break;
	case WM_RBUTTONUP:
		isRClicked = false; break;
	case WM_LBUTTONDOWN:
		isClicked = true; break;
	case WM_RBUTTONDOWN:
		isRClicked = true; break;
	case WM_PAINT:
		DrawGLScene();
		BeginPaint(hWnd, &ps);

		EndPaint(hWnd, &ps);
		return 0;
	}

	// Pass All Unhandled Messages To DefWindowProc
	return DefWindowProc(hWnd, uMsg, wParam, lParam);
}

int WINAPI WinMain(HINSTANCE    hInstance,HINSTANCE	hPrevInstance,LPSTR	lpCmdLine,int nCmdShow)			// Window Show State
{
	MSG		msg;									// Windows Message Structure
	BOOL	done = FALSE;								// Bool Variable To Exit Loop

	// Ask The User Which Screen Mode They Prefer
	//if (MessageBox(NULL,"Would You Like To Run In Fullscreen Mode?", "Start FullScreen?",MB_YESNO|MB_ICONQUESTION)==IDNO)
	{
		fullscreen = FALSE;							// Windowed Mode
	}

	// Create Our OpenGL Window
	if (!CreateGLWindow("OpenGL template", 640, 480, 16, fullscreen))
	{
		return 0;									// Quit If Window Was Not Created
	}

	while (!done)									// Loop That Runs While done=FALSE
	{
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))	// Is There A Message Waiting?
		{
			if (msg.message == WM_QUIT)				// Have We Received A Quit Message?
			{
				done = TRUE;							// If So done=TRUE
			}
			else									// If Not, Deal With Window Messages
			{
				TranslateMessage(&msg);				// Translate The Message
				DispatchMessage(&msg);				// Dispatch The Message
			}
		}
		else										// If There Are No Messages
		{
			// Draw The Scene.  Watch For ESC Key And Quit Messages From DrawGLScene()
			if (active)								// Program Active?
			{
				if (keys[VK_ESCAPE])				// Was ESC Pressed?
				{
					done = TRUE;						// ESC Signalled A Quit
				}
				else								// Not Time To Quit, Update Screen
				{
					DrawGLScene();					// Draw The Scene
					SwapBuffers(hDC);				// Swap Buffers (Double Buffering)
				}
			}

			if (keys[VK_F1])						// Is F1 Being Pressed?
			{
				keys[VK_F1] = FALSE;					// If So Make Key FALSE
				KillGLWindow();						// Kill Our Current Window
				fullscreen = !fullscreen;				// Toggle Fullscreen / Windowed Mode
				// Recreate Our OpenGL Window
				if (!CreateGLWindow("OpenGL template", 640, 480, 16, fullscreen))
				{
					return 0;						// Quit If Window Was Not Created
				}
			}
		}
	}

	// Shutdown
	KillGLWindow();									// Kill The Window
	return (msg.wParam);							// Exit The Program
}
int main(HINSTANCE hinstance, HINSTANCE hPrevlnstance, LPSTR iPCmdLine, int nCmdShow)
{	
	
	cout << "Enter number of balls (between 1 and 8): " << endl;
	cin >> numberOfBalls;

	while (numberOfBalls > 8 || numberOfBalls < 1)
	{
		cout << "Enter number of balls (between 1 and 8): " << endl;
		cin >> numberOfBalls;
	}

	cout << "Enter value of inplus force: " << endl;
	cin >> d;
	
	return WinMain(hinstance, hPrevlnstance, iPCmdLine, nCmdShow);
}
