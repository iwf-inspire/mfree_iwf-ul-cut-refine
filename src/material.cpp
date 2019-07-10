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

#include "material.h"

void material_eos(body &b) {
	std::vector<particle> &particles = b.get_particles();
	double K    = b.get_sim_data().get_physical_constants().K();
	double rho0 = b.get_sim_data().get_physical_constants().rho0();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		const double c0 = sqrt(K/rho0);
		particles[i].p = c0*c0*(particles[i].rho - rho0);
	}
}

void material_stress_rate_jaumann(body &b) {
	std::vector<particle> &particles = b.get_particles();
	double G = b.get_sim_data().get_physical_constants().G();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		const glm::dmat3x3 epsdot = glm::dmat3x3(particles[i].vx_x, 0.5*(particles[i].vx_y + particles[i].vy_x), 0., 0.5*(particles[i].vx_y + particles[i].vy_x), particles[i].vy_y, 0., 0., 0., 0.);
		const glm::dmat3x3 omega  = glm::dmat3x3(0.     , 0.5*(particles[i].vy_x - particles[i].vx_y), 0., 0.5*(particles[i].vx_y - particles[i].vy_x), 0., 0., 0., 0., 0.);
		const glm::dmat3x3 S      = glm::dmat3x3(particles[i].Sxx, particles[i].Sxy, 0., particles[i].Sxy, particles[i].Syy, 0., 0., 0., particles[i].Szz);
		const glm::dmat3x3 I      = glm::dmat3x3(1.);

		const double trace_epsdot = epsdot[0][0] + epsdot[1][1] + epsdot[2][2];

		const glm::dmat3x3 S_t = 2*G*(epsdot - 1./3.*trace_epsdot*I) + omega*S + S*glm::transpose(omega);	//Belytschko (3.7.9)

		particles[i].Sxx_t += S_t[0][0];
		particles[i].Sxy_t += S_t[0][1];
		particles[i].Syy_t += S_t[1][1];
		particles[i].Szz_t += S_t[2][2];
	}
}
