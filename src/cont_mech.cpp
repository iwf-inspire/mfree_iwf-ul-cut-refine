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

#include "cont_mech.h"

void contmech_continuity(body &b) {
	std::vector<particle> &particles = b.get_particles();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		const double rho  = particles[i].rho;
		const double vx_x = particles[i].vx_x;
		const double vy_y = particles[i].vy_y;

		particles[i].rho_t -= rho*(vx_x + vy_y);
	}
}


void contmech_momentum(body &b) {
	std::vector<particle> &particles = b.get_particles();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		const double Sxx_x = particles[i].Sxx_x;
		const double Sxy_y = particles[i].Sxy_y;
		const double Sxy_x = particles[i].Sxy_x;
		const double Syy_y = particles[i].Syy_y;

		const double rho = particles[i].rho;

		particles[i].vx_t += 1./rho*(Sxx_x + Sxy_y) + particles[i].fcx / particles[i].m + particles[i].ftx / particles[i].m;
		particles[i].vy_t += 1./rho*(Sxy_x + Syy_y) + particles[i].fcy / particles[i].m + particles[i].fty / particles[i].m;
	}
}

void contmech_advection(body &b) {
	std::vector<particle> &particles = b.get_particles();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		particles[i].x_t += particles[i].vx;
		particles[i].y_t += particles[i].vy;
	}
}

void do_boundary_conditions(body &b) {
	std::vector<particle> &particles = b.get_particles();

	// this enforces the fixed boundary conditions
	// demonstrated in Fig. 10 of the manuscript

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		if(particles[i].fixed) {
			particles[i].x    = particles[i].X;
			particles[i].y    = particles[i].Y;
			particles[i].x_t  = 0.;
			particles[i].y_t  = 0.;
			particles[i].vx   = 0.;
			particles[i].vy   = 0.;
			particles[i].vx_t = 0.;
			particles[i].vy_t = 0.;
		}
	}
}
