#ifndef INCLUDE_H_
#define INCLUDE_H_

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#include "omp.h"
#endif
//#include <OpenEXR/ImathVec.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <omp.h>
#include <fenv.h>

//#include <armadillo>
#include "CHcudaMath.h"
using namespace std;

#define PI 3.141592653589793

real cell_size = .2;
real half_cell = cell_size * .5;
real2 sim_size = R2(5, 2);
real2 grid_size = sim_size / cell_size;
int num_nodes = grid_size.x * grid_size.y;

real density = 1000;
real E = 2e7;
real nu = .333;

real volume = cell_size * cell_size * .25;
real mass = volume * density;
real Vc = sqrt(E / density);
real dt = 1e-5;

int num_particles = 0;

vector<real2> p_pos, p_vel, p_vel_grad, p_acc, p_lpos;
vector<real> p_mass, p_E, p_rho, p_vol;
vector<realtensor2D> p_s, p_F;

vector<int4> p_node;
vector<bool> p_fixed;

vector<real> g_mass(num_nodes, 0);
vector<real2> g_mom(num_nodes, R2(0)), g_fint(num_nodes, R2(0)), g_fext(num_nodes, R2(0)), g_vel(num_nodes, R2(0)), g_acc(num_nodes, R2(0));
vector<bool> g_fixed(num_nodes, 0);
int frame = 0;

real2 max_fint = R2(0, 0);
real2 max_vel = R2(0, 0);

int draw_mode = 1;
int draw_type = 1;
bool paused = 0;
bool deform_grid = 0;

int stress_model = 0;



//OpenGL code
real3 GetColour(double v, double vmin, double vmax) {
	real3 c = { 1.0, 1.0, 1.0 }; // white
	double dv;

	if (v < vmin)
		v = vmin;
	if (v > vmax)
		v = vmax;
	dv = vmax - vmin;

	if (v < (vmin + 0.25 * dv)) {
		c.x = 0;
		c.y = 4 * (v - vmin) / dv;
	} else if (v < (vmin + 0.5 * dv)) {
		c.x = 0;
		c.z = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
	} else if (v < (vmin + 0.75 * dv)) {
		c.x = 4 * (v - vmin - 0.5 * dv) / dv;
		c.z = 0;
	} else {
		c.y = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
		c.z = 0;
	}
	return (c);
}

void changeSize(int w, int h) {
	if (h == 0) {
		h = 1;
	}
	float ratio = 1.0 * w / h;

	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, grid_size.x * cell_size, 0, grid_size.y * cell_size, 0, 1);
	glMatrixMode (GL_MODELVIEW);
}

void initScene() {
	GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glShadeModel (GL_SMOOTH);
	glEnable (GL_COLOR_MATERIAL);

	glEnable (GL_POINT_SMOOTH);
	//glEnable( GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_POINT_SMOOTH_HINT, GL_DONT_CARE);
	//   glEnable(GL_DEPTH_TEST);
	//   glDepthFunc(GL_LESS);
	//glFrontFace(GL_CCW);
	//glCullFace(GL_BACK);
	//glEnable(GL_CULL_FACE);
	//glDepthFunc( GL_LEQUAL);
	//glClearDepth(1.0);
	glPointSize(2);

	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glEnable (GL_LIGHTING);
	glEnable (GL_LIGHT0);
}

void processNormalKeys(unsigned char key, int x, int y) {
	switch (key) {

	case '1':
		draw_mode = 1;
		break;
	case '2':
		draw_mode = 2;
		break;
	case '3':
		draw_mode = 3;
		break;
	case 'p':
		draw_type = 1;
		break;
	case 'e':
		draw_type = 2;
		break;
	case ' ':
		paused = !paused;
		break;
	case 'd':
		deform_grid = !deform_grid;
		break;
	}
}

void drawGrid() {

	glColor3f(.9, .9, .9);
	glBegin (GL_LINES);
	for (int i = 0; i < grid_size.x; i++) {
		glVertex2f(i * cell_size, 0); // origin of the line
		glVertex2f(i * cell_size, grid_size.y * cell_size); // origin of the line
	}
	glEnd();
	glBegin(GL_LINES);
	for (int i = 0; i < grid_size.y; i++) {
		glVertex2f(0, i * cell_size); // origin of the line
		glVertex2f(grid_size.x * cell_size, i * cell_size); // origin of the line
	}
	glEnd();

}

void drawParticles() {
	if (draw_type != 1) {
		return;
	}
	real3 col = R3(0, 0, 0);
	for (int i = 0; i < num_particles; i++) {
		real2 pos = p_pos[i];
		if (draw_mode == 2) {
			real len = sqrt(p_vel[i].x * p_vel[i].x + p_vel[i].y * p_vel[i].y);
			col = GetColour(len, 0, sqrt(max_vel.x * max_vel.x + max_vel.y * max_vel.y));
		} else if (draw_mode == 3) {
			real len = sqrt(p_s[i].xx * p_s[i].xx + p_s[i].yy * p_s[i].yy + p_s[i].xy * p_s[i].yx);
			col = GetColour(len, 0, sqrt(max_fint.x * max_fint.x + max_fint.y * max_fint.y));
		}
		glColor3f(col.x, col.y, col.z);
		glBegin (GL_POINTS);
		glVertex2f(pos.x, pos.y);
		glEnd();
	}
}

void drawElements() {
	if (draw_type != 2) {
		return;
	}

	real3 col = R3(0, 0, 0);
	if (deform_grid == false) {
		glBegin (GL_QUADS);
		for (int i = 0; i < grid_size.x; i++) {
			for (int j = 0; j < grid_size.y; j++) {
				int node = i + grid_size.x * j;
				if (draw_mode == 2) {
					real len = sqrt(g_vel[node].x * g_vel[node].x + g_vel[node].y * g_vel[node].y);
					col = GetColour(len, 0, sqrt(max_vel.x * max_vel.x + max_vel.y * max_vel.y));
				} else if (draw_mode == 3) {
					real len = sqrt(g_fint[node].x * g_fint[node].x + g_fint[node].y * g_fint[node].y);
					col = GetColour(len, 0, sqrt(max_fint.x * max_fint.x + max_fint.y * max_fint.y));
				}
				glColor3f(col.x, col.y, col.z);
				glVertex2f(i * cell_size, j * cell_size);
				glVertex2f(i * cell_size + cell_size, j * cell_size);
				glVertex2f(i * cell_size + cell_size, j * cell_size + cell_size);
				glVertex2f(i * cell_size, j * cell_size + cell_size);

			}

		}
		glEnd();
	}
	if (deform_grid == true) {
		glBegin (GL_POINTS);
		for (int i = 0; i < grid_size.x; i++) {
			for (int j = 0; j < grid_size.y; j++) {
				int node = i + grid_size.x * j;
				//if (draw_mode == 2) {
				real len = sqrt(g_vel[node].x * g_vel[node].x + g_vel[node].y * g_vel[node].y);
				col = GetColour(len, 0, sqrt(max_vel.x * max_vel.x + max_vel.y * max_vel.y));
				//} else if (draw_mode == 3) {
				//	real len = sqrt(g_fint[node].x * g_fint[node].x + g_fint[node].y * g_fint[node].y);
				//	col = GetColour(len, 0, max_stress / 50.);
				//}
				glColor3f(col.x, col.y, col.z);

				real2 p1 = R2(i * cell_size, j * cell_size) + g_vel[i + grid_size.x * j] * dt * 1000;
				real2 p2 = R2(i * cell_size + cell_size, j * cell_size) + g_vel[(i + 1) + grid_size.x * j] * dt * 1000;
				real2 p3 = R2(i * cell_size + cell_size, j * cell_size + cell_size) + g_vel[(i + 1) + grid_size.x * (j + 1)] * dt * 1000;
				real2 p4 = R2(i * cell_size, j * cell_size + cell_size) + g_vel[(i) + grid_size.x * (j + 1)] * dt * 1000;

				glVertex2f(p1.x, p1.y);
				glVertex2f(p2.x, p2.y);
				glVertex2f(p3.x, p3.y);
				glVertex2f(p4.x, p4.y);

			}

		}
		glEnd();

	}
}
void drawVectors() {

	real3 col = R3(.5, .5, .5);
	glColor3f(col.x, col.y, col.z);
	glLineWidth(4);
	glBegin (GL_LINES);

	for (int i = 0; i < grid_size.x; i++) {
		for (int j = 0; j < grid_size.y; j++) {
			int node = i + grid_size.x * j;
			real2 p1 = R2(i * cell_size, j * cell_size);
			real2 vel = g_vel[node];

			real len = sqrt(vel.x * vel.x + vel.y * vel.y);

			real2 p2 = p1 + g_vel[node] / len * .1;
			glVertex2f(p1.x, p1.y);
			glVertex2f(p2.x, p2.y);\
		}
	}

	glEnd();
	glLineWidth(1.1);
}
#endif
