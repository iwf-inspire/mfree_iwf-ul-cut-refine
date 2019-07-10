/*
 * ================================================
 * 					COPYRIGHT:
 * Institute of Machine Tools & Manufacturing (IWF)
 * Department of Mechanical & Process Engineering
 * 					ETH ZURICH
 * ================================================
 *
 *  This file is part of "mfree_iwf-ul-cut-refine".
 *
 * 	mfree_iwf is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	mfree_iwf-ul-cut-refine is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *  along with mfree_iwf-ul-cut-refine.  If not, see <http://www.gnu.org/licenses/>.
 *
 *	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *  This is the source code used to produce the results
 *  of the metal cutting simulation presented in:
 *
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  "Meshfree Simulation of Metal Cutting:
 *  An Updated Lagrangian Approach with Dynamic Refinement"
 *
 * 	Authored by:
 * 	Mohamadreza Afrasiabi
 * 	Dr. Matthias Roethlin
 * 	Hagen Klippel
 * 	Prof. Dr. Konrad Wegener
 *
 * 	Published by:
 * 	International Journal of Mechanical Sciences
 * 	28 June 2019
 * 	https://doi.org/10.1016/j.ijmecsci.2019.06.045
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * 	For further descriptions, you may refer to the manuscript
 * 	or the previous works of the same research group
 * 	at IWF, ETH Zurich.
 *
 */

#include "correctors.h"

static double stress_angle(double sxx, double sxy, double syy, double eps) {
	double numer = 2.*sxy;
	double denom = sxx - syy + eps;
	return 0.5*atan2(numer,denom);
}

// source: https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int signum(T val) {
    return (T(0) < val) - (val < T(0));
}

// 2D Artificial Monaghan Viscosity
void correctors_mghn_artificial_viscosity(body &b) {
	const double alpha = b.get_sim_data().get_correction_constants().get_art_visc_const().artvisc_alpha();
	const double beta  = b.get_sim_data().get_correction_constants().get_art_visc_const().artvisc_beta();
	const double eta   = b.get_sim_data().get_correction_constants().get_art_visc_const().artvisc_eta();

	const double K = b.get_sim_data().get_physical_constants().K();

	std::vector<particle> &particles = b.get_particles();
	unsigned int n = b.get_num_part();

	for (unsigned int i = 0; i < n; i++) {
		const double xi   = particles[i].x;
		const double yi   = particles[i].y;
		const double vxi  = particles[i].vx;
		const double vyi  = particles[i].vy;

		for (unsigned int j = 0; j < particles[i].num_nbh; j++) {
			unsigned int jdx = particles[i].nbh[j];

			const double xj   = particles[jdx].x;
			const double yj   = particles[jdx].y;
			const double vxj  = particles[jdx].vx;
			const double vyj  = particles[jdx].vy;

			const double vijx = vxi - vxj;
			const double vijy = vyi - vyj;

			const double xij = xi - xj;
			const double yij = yi - yj;

			const double vijposij = vijx*xij + vijy*yij;

			if (vijposij >= 0.) continue;

			const kernel_result w = particles[i].w[j];

			const double hi   = particles[i].h;
			const double rhoi = particles[i].rho;

			if (rhoi < 0.) {
				printf("neg density!\n");
				exit(-1);
			}

			assert(rhoi > 0.);
			const double ci   = sqrt(K/rhoi);

			const double hj   = particles[jdx].h;
			const double mj   = particles[jdx].m;
			const double rhoj = particles[jdx].rho;

			if (rhoj < 0.) {
				printf("neg density!\n");
				exit(-1);
			}

			assert(rhoj > 0.);
			const double cj   = sqrt(K/rhoj);

			const double cij = 0.5*(ci+cj);
			const double hij = 0.5*(hi+hj);
			const double rhoij = 0.5*(rhoi+rhoj);
			const double r2ij = xij*xij + yij*yij;
			const double muij = (hij*vijposij)/(r2ij + eta*eta*hij*hij);
			const double piij = (-alpha*cij*muij + beta*muij*muij)/rhoij;

			particles[i].vx_t += -mj*piij*w.w_x;
			particles[i].vy_t += -mj*piij*w.w_y;
		}
	}
}

// 2D Artificial Monaghan Stress
void correctors_mghn_artificial_stress(body &b) {
	std::vector<particle> &particles  = b.get_particles();
	unsigned int n = b.get_num_part();

	const double eps = b.get_sim_data().get_correction_constants().get_monaghan_const().mghn_eps();

	for (unsigned int i = 0; i < n; i++) {
		double rhoi = particles[i].rho;
		double rhoi21 = 1./(rhoi*rhoi);

		double sxx = particles[i].Sxx - particles[i].p;
		double syy = particles[i].Syy - particles[i].p;

		// off diag
		double sxy = particles[i].Sxy;

		double theta = stress_angle(sxx,sxy,syy,0.);

		double cos_theta = cos(theta);
		double sin_theta = sin(theta);

		double cos_theta2 = cos_theta*cos_theta;
		double sin_theta2 = sin_theta*sin_theta;

		double rot_sxx = cos_theta2*sxx + 2.0*cos_theta*sin_theta*sxy + sin_theta2*syy;
		double rot_syy = sin_theta2*sxx - 2.0*cos_theta*sin_theta*sxy + cos_theta2*syy;

		double rot_rxx = 0.;
		double rot_ryy = 0.;
		if (rot_sxx > 0) rot_rxx = -eps*rot_sxx*rhoi21;
		if (rot_syy > 0) rot_ryy = -eps*rot_syy*rhoi21;

		particles[i].Rxx = cos_theta2*rot_rxx + sin_theta2*rot_ryy;
		particles[i].Rxy = cos_theta*sin_theta*(rot_rxx - rot_ryy);
		particles[i].Ryy = sin_theta2*rot_rxx + cos_theta2*rot_ryy;
	}
}

// 2D XSPH
void correctors_xsph(body &b) {
	const double eps = b.get_sim_data().get_correction_constants().xsph_eps();

	std::vector<particle> &particles = b.get_particles();
	unsigned int n = b.get_num_part();

	for (unsigned int i = 0; i < n; i++) {
		const double vxi = particles[i].vx;
		const double vyi = particles[i].vy;
		const double rhoi = particles[i].rho;

		for (unsigned int j = 0; j < particles[i].num_nbh; j++) {
			const unsigned int jdx = particles[i].nbh[j];
			const kernel_result w = particles[i].w[j];

			const double vxj  = particles[jdx].vx;
			const double vyj  = particles[jdx].vy;
			const double rhoj = particles[jdx].rho;

			const double vijx = vxi - vxj;
			const double vijy = vyi - vyj;

			const double rhoij = 0.5*(rhoi + rhoj);
			const double fac = -eps*w.w*particles[i].m/rhoij;

			particles[i].x_t += fac*vijx;
			particles[i].y_t += fac*vijy;
		}
	}
}
