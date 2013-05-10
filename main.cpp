#include "include.h"
#include "model.h"
void computeNodes() {
	//	[D - C]
	//	[A - B]
	for (int i = 0; i < num_particles; i++) {
		real2 p1 = p_pos[i];
		real2 p = p1 / cell_size;
		int x = floor(p.x);
		int y = floor(p.y);
		int4 node;

		//get the node numbers
		node.x = x + grid_size.x * y;
		node.y = node.x + 1;
		node.z = node.y + grid_size.x;
		node.w = node.x + grid_size.x;

		p_node[i] = node;

		//compute local position
		real cx = x * cell_size + cell_size * .5;
		real cy = y * cell_size + cell_size * .5;

		real xi = (p1.x - cx) / half_cell;
		real eta = (p1.y - cy) / half_cell;

		p_lpos[i] = R2(xi, eta);
	}
}

void shapeFunction(real2 lpos, real4& N, real4& dNdx, real4& dNdy) {
	real xi = lpos.x;
	real eta = lpos.y;

	N.x = (1 - xi) * (1 - eta) / 4.0;
	N.y = (1 + xi) * (1 - eta) / 4.0;
	N.z = (1 + xi) * (1 + eta) / 4.0;
	N.w = (1 - xi) * (1 + eta) / 4.0;

	dNdx.x = -1.0 / 4.0 * (1 - eta) / half_cell;
	dNdx.y = 1.0 / 4.0 * (1 - eta) / half_cell;
	dNdx.z = 1.0 / 4.0 * (1 + eta) / half_cell;
	dNdx.w = -1.0 / 4.0 * (1 + eta) / half_cell;

	dNdy.x = -1.0 / 4.0 * (1 - xi) / half_cell;
	dNdy.y = -1.0 / 4.0 * (1 + xi) / half_cell;
	dNdy.z = 1.0 / 4.0 * (1 + xi) / half_cell;
	dNdy.w = 1.0 / 4.0 * (1 - xi) / half_cell;

}

void updateStress(int i, int4 node, real4 dNdx, real4 dNdy) {
	if (stress_model == 0) {
		linearStress(i, node, dNdx, dNdy);
	} else if (stress_model == 1) {
		StVK(i, node, dNdx, dNdy);
	} else if (stress_model == 2) {
		CoRotated(i, node, dNdx, dNdy);
	}
	//NeoHookean(i, node, dNdx, dNdy);

}


void projectToGrid() {
	for (int i = 0; i < num_particles; i++) {

		int4 node = p_node[i];
		real2 lpos = p_lpos[i];
		real4 N, dNdx, dNdy;
		shapeFunction(lpos, N, dNdx, dNdy);

		updateStress(i, node, dNdx, dNdy);

		real mass = p_mass[i];
		real vol = p_vol[i];

		g_mass[node.x] += mass * N.x;
		g_mass[node.y] += mass * N.y;
		g_mass[node.z] += mass * N.z;
		g_mass[node.w] += mass * N.w;

		g_mom[node.x] += p_mass[i] * p_vel[i] * N.x;
		g_mom[node.y] += p_mass[i] * p_vel[i] * N.y;
		g_mom[node.z] += p_mass[i] * p_vel[i] * N.z;
		g_mom[node.w] += p_mass[i] * p_vel[i] * N.w;

		real damping = 1;

		g_fint[node.x].x -= vol * (p_s[i].xx * dNdx.x + p_s[i].yx * dNdy.x)*damping;
		g_fint[node.x].y -= vol * (p_s[i].xy * dNdx.x + p_s[i].yy * dNdy.x)*damping;

		g_fint[node.y].x -= vol * (p_s[i].xx * dNdx.y + p_s[i].yx * dNdy.y)*damping;
		g_fint[node.y].y -= vol * (p_s[i].xy * dNdx.y + p_s[i].yy * dNdy.y)*damping;

		g_fint[node.z].x -= vol * (p_s[i].xx * dNdx.z + p_s[i].yx * dNdy.z)*damping;
		g_fint[node.z].y -= vol * (p_s[i].xy * dNdx.z + p_s[i].yy * dNdy.z)*damping;

		g_fint[node.w].x -= vol * (p_s[i].xx * dNdx.w + p_s[i].yx * dNdy.w)*damping;
		g_fint[node.w].y -= vol * (p_s[i].xy * dNdx.w + p_s[i].yy * dNdy.w)*damping;

	}
}

void solveGrid() {
	for (int i = 0; i < num_nodes; i++) {

		if (g_fixed[i] == false) {
			if (g_mass[i] != 0) {
				g_vel[i] = g_mom[i] / g_mass[i];
				g_acc[i] = (g_fint[i]) / g_mass[i];
				g_vel[i] += g_acc[i] * dt;

				max_fint.x = max(max_fint.x, g_fint[i].x);
				max_fint.y = max(max_fint.y, g_fint[i].y);

			}
		}

	}
}

void integrate() {
	for (int i = 0; i < num_particles; i++) {

		int4 node = p_node[i];
		real2 lpos = p_lpos[i];
		real4 N, dNdx, dNdy;
		shapeFunction(lpos, N, dNdx, dNdy);

		real2 accel = R2(0);
		accel += g_acc[node.x] * N.x;
		accel += g_acc[node.y] * N.y;
		accel += g_acc[node.z] * N.z;
		accel += g_acc[node.w] * N.w;

		real2 velocity = R2(0);
		velocity += g_vel[node.x] * N.x;
		velocity += g_vel[node.y] * N.y;
		velocity += g_vel[node.z] * N.z;
		velocity += g_vel[node.w] * N.w;

		max_vel.x = max(max_vel.x, velocity.x);
		max_vel.y = max(max_vel.y, velocity.y);

		realtensor2D F = GetDefGrad(node, dNdx, dNdy);
		accel*=R2(1);
		p_pos[i] = p_pos[i] + velocity * dt;
		p_vel[i] = p_vel[i] + accel * dt;

		//p_vel_grad[i] = velocity_grad;

//		realtensor2D Fp;
//
//		Fp.xx = F.xx * p_F[i].xx;
//		Fp.yy = F.yy * p_F[i].yy;
//
//		real J = Fp.xx * Fp.yy - Fp.xy * Fp.yx;
//		p_F[i] = Fp;
//
//		if (J <= 0) {
//			cout << "J error Vol " << J << endl;
//			exit(1);
//		}
//		real V0 = cell_size * cell_size * .25;
//		p_vol[i] = V0 * J;

		updateStress(i, node, dNdx, dNdy);

	}

}

void simulate() {
	fill(g_mass.begin(), g_mass.end(), 0);
	fill(g_mom.begin(), g_mom.end(), R2(0));
	fill(g_fint.begin(), g_fint.end(), R2(0));
	fill(g_fext.begin(), g_fext.end(), R2(0));
	fill(g_vel.begin(), g_vel.end(), R2(0));
	fill(g_acc.begin(), g_acc.end(), R2(0));

	computeNodes();

	projectToGrid();
//	fill(p_s.begin(), p_s.end(), realtensor2D(0,0,0,0));
//projectMassVelToGrid();
//projectForceToGrid();

	solveGrid();
	integrate();

}

void createCube(real2 pos, real2 size, real density, real2 vel, real E, real rho) {
	for (int i = 0; i < size.x * density; i++) {
		for (int j = 0; j < size.y * density; j++) {
			p_pos.push_back(R2(i / density + pos.x, j / density + pos.y));
			p_vel.push_back(R2(vel.x, vel.y));
			p_acc.push_back(R2(0, 0));
			p_mass.push_back(mass);
			realtensor2D temp;
			temp.xx = temp.xy = temp.yx = temp.yy = 0;
			p_s.push_back(temp);
			temp.xx = 1;
			temp.xy = 0;
			temp.yx = 0;
			temp.yy = 1;
			p_F.push_back(temp);
			p_node.push_back(I4(0, 0, 0, 0));
			p_lpos.push_back(R2(0, 0));
			p_E.push_back(E);
			p_rho.push_back(rho);
			p_vol.push_back(volume);
			p_fixed.push_back(false);
			p_vel_grad.push_back(R2(0));
			num_particles++;
		}
	}
}

void fixGrid() {
	for (int i = 0; i < grid_size.x; i++) {
		for (int j = 0; j < grid_size.y; j++) {
			int node = i + grid_size.x * j;
			if (i * cell_size >= sim_size.x - cell_size + cell_size || i * cell_size <= -sim_size.x + cell_size + cell_size) {
				g_fixed[node] = true;
			}
			if (j * cell_size >= sim_size.y - cell_size + cell_size || j * cell_size <= -sim_size.y + cell_size + cell_size) {
				g_fixed[node] = true;
			}
		}
	}
}

void initParticles() {
	createCube(R2(3, .5), R2(1, 1), 25, R2(-5, 0), E, density);
	createCube(R2(1, .5), R2(1, 1), 25, R2(5, 0), E, density);

	fixGrid();
}

void renderScene() {
	max_fint.x=max_fint.y=0;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
// glEnable(GL_DEPTH_TEST);
	glDisable (GL_DEPTH_TEST);
	glDepthFunc (GL_LEQUAL);
	glLoadIdentity();
	glPointSize(8);
	glLineWidth(1.1);
	glColor3f(0, 0, 0);
	if (paused == false) {
		simulate();
		cout << "Frame " << frame << " time:" << frame * dt<<" "<<max_fint.x<<"N "<<max_fint.y<<"N" << endl;
	}

	if (frame % 100 == 0) {

		drawElements();
		if (deform_grid == false) {
			drawGrid();
		}
		drawVectors();
		drawParticles();

		glutSwapBuffers();

	}
	frame++;

}
int main(int argc, char* argv[]) {
	stress_model=atoi(argv[1]);
	E=atof(argv[2]);
	density = atof(argv[2]);


	cout << grid_size.x << " " << grid_size.y << endl;
	initParticles();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(500 * 2, 200 * 2);
	glutCreateWindow("MPM");
	glutDisplayFunc(renderScene);
	glutIdleFunc(renderScene);
	glutReshapeFunc(changeSize);
	glutIgnoreKeyRepeat(0);
	glutKeyboardFunc(processNormalKeys);
//glutMouseFunc (mouseButton);
//glutMotionFunc (mouseMove);
	initScene();
	glutMainLoop();

}
