#include "include.h"

realtensor2D GetDefGrad(int4 node, real4 dNdx, real4 dNdy) {

	real4 grid_velx = R4(g_vel[node.w].x, g_vel[node.x].x, g_vel[node.y].x, g_vel[node.z].x);
	real4 grid_vely = R4(g_vel[node.w].y, g_vel[node.x].y, g_vel[node.y].y, g_vel[node.z].y);

	realtensor2D F;
	F.xx = dot(grid_velx * dt, dNdx) + 1;
	F.xy = dot(grid_velx * dt, dNdy);
	F.yx = dot(grid_vely * dt, dNdx);
	F.yy = dot(grid_vely * dt, dNdy) + 1;
	return F;

}

void PolarDecomp(const realtensor2D & F, realtensor2D & R, realtensor2D & U) {
	real J = F.xx * F.yy - F.xy * F.yx;
	if (J <= 0) {
		cout << "J error" << endl;
		exit(1);
	}
	real S = sqrt(0.3e1) * sqrt(0.1e1 / (F.xx * F.xx + F.yx * F.yx + F.xy * F.xy + F.yy * F.yy));
	realtensor2D E, A, X;

	E.xx = F.xx * F.xx + F.yx * F.yx;
	E.xy = F.xx * F.xy + F.yx * F.yy;
	E.yx = F.xx * F.xy + F.yx * F.yy;
	E.yy = F.xy * F.xy + F.yy * F.yy;

	E.xx = S * E.xx / 0.2e1 - 0.1e1 / 0.2e1;
	E.xy = S * E.xy / 0.2e1;
	E.yx = S * E.yx / 0.2e1;
	E.yy = S * E.yy / 0.2e1 - 0.1e1 / 0.2e1;

	S = sqrt(S);
	real ERRZ = E.xx * E.xx + E.yy * E.yy + 2 * E.xy * E.xy;

	A.xx = F.xx * S;
	A.xy = F.xy * S;
	A.yx = F.yx * S;
	A.yy = F.yy * S;

	bool converged = false;
	if (ERRZ + 1.0 == 1.0) {
		converged = true;
	}

	int num_iters = 0;
	while (converged == false) {

		X.xx = A.xx * (1 - E.xx) - A.xy * E.yx;
		X.xy = -A.xx * E.xy + A.xy * (1 - E.yy);
		X.yx = A.yx * (1 - E.xx) - A.yy * E.yx;
		X.yy = -A.yx * E.xy + A.yy * (1 - E.yy);
		A = X;

		E.xx = A.xx * A.xx / 0.2e1 + A.yx * A.yx / 0.2e1 - 0.1e1 / 0.2e1;
		E.xy = A.xx * A.xy / 0.2e1 + A.yx * A.yy / 0.2e1;
		E.yx = A.xx * A.xy / 0.2e1 + A.yx * A.yy / 0.2e1;
		E.yy = A.xy * A.xy / 0.2e1 + A.yy * A.yy / 0.2e1 - 0.1e1 / 0.2e1;

		real ERR = E.xx * E.xx + E.yy * E.yy + 2 * E.xy * E.xy;
		if (ERR >= ERRZ || ERR + 1.0 == 1.0) {
			converged = true;
		}
		real old_ERRZ = ERRZ;
		ERRZ = ERR;
		if (num_iters == 200) {
			cout << "not converged" << endl;
			exit(1);
		}
		num_iters++;
	}

	R = A;
	U.xx = R.xx * F.xx + R.yx * F.yx;
	U.xy = R.xx * F.xy + R.yx * F.yy;
	U.yx = R.xy * F.xx + R.yy * F.yx;
	U.yy = R.xy * F.xy + R.yy * F.yy;

}
void computeInfinitesimalStrain(realtensor2D & strain_rate, int4 node, real4 dNdx, real4 dNdy) {
//	real4 grid_velx = R4(g_vel[node.w].x, g_vel[node.x].x, g_vel[node.y].x, g_vel[node.z].x);
//	real4 grid_vely = R4(g_vel[node.w].y, g_vel[node.x].y, g_vel[node.y].y, g_vel[node.z].y);
//
//	strain_rate.xx = dot(grid_velx, dNdx);
//	strain_rate.yy = dot(grid_vely, dNdy);
//
//	strain_rate.xy += (g_vel[node.x].x * dNdy.x + g_vel[node.x].y * dNdx.x) / 2.0;
//	strain_rate.xy += (g_vel[node.y].x * dNdy.y + g_vel[node.y].y * dNdx.y) / 2.0;
//	strain_rate.xy += (g_vel[node.z].x * dNdy.z + g_vel[node.z].y * dNdx.z) / 2.0;
//	strain_rate.xy += (g_vel[node.w].x * dNdy.w + g_vel[node.w].y * dNdx.w) / 2.0;
//
//	strain_rate.yx = strain_rate.xy;

	real4 grid_velx = R4(g_vel[node.w].x, g_vel[node.x].x, g_vel[node.y].x, g_vel[node.z].x);
	real4 grid_vely = R4(g_vel[node.w].y, g_vel[node.x].y, g_vel[node.y].y, g_vel[node.z].y);

	realtensor2D F;
	F.xx = dot(grid_velx, dNdx) + 1;
	F.xy = dot(grid_velx, dNdy);
	F.yx = dot(grid_vely, dNdx);
	F.yy = dot(grid_vely, dNdy) + 1;
	strain_rate.xx = F.xx - 1 / 2.0;
	strain_rate.xy = F.xy / 2.0 + F.yx / 2.0;
	strain_rate.yx = F.xy / 2.0 + F.yx / 2.0;
	strain_rate.yy = F.yy - 1 / 2.0;

}

realtensor2D computeGreenStrain(realtensor2D & strain_rate, int4 node, real4 dNdx, real4 dNdy) {

	realtensor2D F = GetDefGrad(node, dNdx, dNdy);

	strain_rate.xx = F.xx * F.xx / 0.2e1 + F.yx * F.yx / 0.2e1 - 0.1e1 / 0.2e1;
	strain_rate.xy = F.xx * F.xy / 0.2e1 + F.yx * F.yy / 0.2e1;
	strain_rate.yx = F.xx * F.xy / 0.2e1 + F.yx * F.yy / 0.2e1;
	strain_rate.yy = F.xy * F.xy / 0.2e1 + F.yy * F.yy / 0.2e1 - 0.1e1 / 0.2e1;

//	realtensor2D R, S;
//	PolarDecomp(F, R, S);
//	strain_rate.xx = S.xx * S.xx / 0.2e1 + S.xy * S.yx / 0.2e1 - 0.1e1 / 0.2e1;
//	strain_rate.xy = S.xx * S.xy / 0.2e1 + S.xy * S.yy / 0.2e1;
//	strain_rate.yx = S.yx * S.xx / 0.2e1 + S.yy * S.yx / 0.2e1;
//	strain_rate.yy = S.xy * S.yx / 0.2e1 + S.yy * S.yy / 0.2e1 - 0.1e1 / 0.2e1;

	return F;

}
void linearStress(int i, int4 node, real4 dNdx, real4 dNdy) {
	realtensor2D strain_rate;
	strain_rate.xx = strain_rate.xy = strain_rate.yx = strain_rate.yy = 0;

	computeInfinitesimalStrain(strain_rate, node, dNdx, dNdy);

	real mu = p_E[i] / (2 * (1 + nu));
	real lambda = p_E[i] * nu / ((1 + nu) * (1 - 2 * nu));

	p_s[i].xx += (2 * mu * strain_rate.xx + lambda * (strain_rate.xx + strain_rate.yy)) * dt;
	p_s[i].xy += (2 * mu * strain_rate.xy) * dt;
	p_s[i].yx += (2 * mu * strain_rate.yx) * dt;
	p_s[i].yy += (2 * mu * strain_rate.yy + lambda * (strain_rate.xx + strain_rate.yy)) * dt;

}

void StVK(int i, int4 node, real4 dNdx, real4 dNdy) {
	realtensor2D E;
	E.xx = E.xy = E.yx = E.yy = 0;
	realtensor2D F = computeGreenStrain(E, node, dNdx, dNdy);
	real mu = p_E[i] / (2 * (1 + nu));
	real lambda = p_E[i] * nu / ((1 + nu) * (1 - 2 * nu));

	p_s[i].xx += (F.xx * (2 * mu * E.xx + lambda * (E.xx + E.yy)) + 2 * F.xy * mu * E.yx);
	p_s[i].xy += (2 * F.xx * mu * E.xy + F.xy * (2 * mu * E.yy + lambda * (E.xx + E.yy)));
	p_s[i].yx += (F.yx * (2 * mu * E.xx + lambda * (E.xx + E.yy)) + 2 * F.yy * mu * E.yx);
	p_s[i].yy += (2 * F.yx * mu * E.xy + F.yy * (2 * mu * E.yy + lambda * (E.xx + E.yy)));

//	p_s[i].xx += (lambda * (E.xx + E.yy) + 2 * mu * E.xx) * dt;
//	p_s[i].xy += (2 * mu * E.xy) * dt;
//	p_s[i].yx += (2 * mu * E.yx) * dt;
//	p_s[i].yy += (lambda * (E.xx + E.yy) + 2 * mu * E.yy) * dt;

}

void NeoHookean(int i, int4 node, real4 dNdx, real4 dNdy) {
	real mu = p_E[i] / (2 * (1 + nu));
	real lambda = p_E[i] * nu / ((1 + nu) * (1 - 2 * nu));
	realtensor2D F = GetDefGrad(node, dNdx, dNdy);

	real J = F.xx * F.yy - F.xy * F.yx;
	if (J <= 0) {
		cout << "J error Neo "<<J << endl;
		exit(1);
	}
	p_s[i].xx += (mu * (F.xx - mu * F.yy / (F.xx * F.yy - F.xy * F.yx)) + lambda * log(J) * F.yy / (F.xx * F.yy - F.xy * F.yx));
	p_s[i].xy += (mu * (F.xy + mu * F.yx / (F.xx * F.yy - F.xy * F.yx)) - lambda * log(J) * F.yx / (F.xx * F.yy - F.xy * F.yx));
	p_s[i].yx += (mu * (F.yx + mu * F.xy / (F.xx * F.yy - F.xy * F.yx)) - lambda * log(J) * F.xy / (F.xx * F.yy - F.xy * F.yx));
	p_s[i].yy += (mu * (F.yy - mu * F.xx / (F.xx * F.yy - F.xy * F.yx)) + lambda * log(J) * F.xx / (F.xx * F.yy - F.xy * F.yx));

}

void CoRotated(int i, int4 node, real4 dNdx, real4 dNdy) {
	real mu = p_E[i] / (2 * (1 + nu));
	real lambda = p_E[i] * nu / ((1 + nu) * (1 - 2 * nu));
	realtensor2D F = GetDefGrad(node, dNdx, dNdy);
	realtensor2D R, S;

	PolarDecomp(F, R, S);

	p_s[i].xx += (2 * mu * (F.xx - R.xx) + lambda * (R.xx * F.xx + R.yx * F.yx - 2 + R.xy * F.xy + R.yy * F.yy) * R.xx);
	p_s[i].xy += (2 * mu * (F.xy - R.xy) + lambda * (R.xx * F.xx + R.yx * F.yx - 2 + R.xy * F.xy + R.yy * F.yy) * R.xy);
	p_s[i].yx += (2 * mu * (F.yx - R.yx) + lambda * (R.xx * F.xx + R.yx * F.yx - 2 + R.xy * F.xy + R.yy * F.yy) * R.yx);
	p_s[i].yy += (2 * mu * (F.yy - R.yy) + lambda * (R.xx * F.xx + R.yx * F.yx - 2 + R.xy * F.xy + R.yy * F.yy) * R.yy);


}

