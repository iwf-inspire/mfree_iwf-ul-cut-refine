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

#include "thermal.h"

#include "body.h"

void thermal::heat_conduction_pse(body &b) const {
	std::vector<particle> &particles = b.get_particles();
	unsigned int num_part = b.get_num_part();

	for (unsigned int i = 0; i < num_part; i++) {
		const double Ti = particles[i].T;
		const double xi = particles[i].x;
		const double yi = particles[i].y;
		const double hi = particles[i].h;
		const double hi2 = hi*hi;

		double T_lapl = 0.;

		for (unsigned int j = 0; j < particles[i].num_nbh; j++) {
			unsigned int jdx = particles[i].nbh[j];

			const double Tj = particles[jdx].T;
			const double xj = particles[jdx].x;
			const double yj = particles[jdx].y;
			const double mj = particles[jdx].m;
			const double rhoj = particles[jdx].rho;

			const double xij = xi-xj;
			const double yij = yi-yj;

			const double r = sqrt(xij*xij + yij*yij);
			const double w_pse = 4.0/(hi2*M_PI)*exp(-r*r/(hi2));
			T_lapl += (Tj-Ti)*w_pse*mj/rhoj/(hi2);
		}

		particles[i].T_t += m_alpha*T_lapl;
	}
}

void thermal::heat_conduction_brookshaw(body &b) const {
	std::vector<particle> &particles = b.get_particles();
	unsigned int num_part = b.get_num_part();

	for (unsigned int i = 0; i < num_part; i++) {
		double Ti = particles[i].T;
		double xi = particles[i].x;
		double yi = particles[i].y;

		double T_lapl = 0.;

		for (unsigned int j = 0; j < particles[i].num_nbh; j++) {
			unsigned int jdx = particles[i].nbh[j];
			kernel_result w = particles[i].w[j];

			if (particles[i].idx == particles[jdx].idx) {
				continue;
			}

			double Tj = particles[jdx].T;
			double xj = particles[jdx].x;
			double yj = particles[jdx].y;
			double mj = particles[jdx].m;
			double rhoj = particles[jdx].rho;

			double xij = xi-xj;
			double yij = yi-yj;

			double rij = sqrt(xij*xij + yij*yij);
			double eijx = xij/rij;
			double eijy = yij/rij;

			if (rij < 1e-12) continue;

			double rij1 = 1./rij;

			T_lapl += 2.0*(mj/rhoj)*(Ti-Tj)*rij1*(eijx*w.w_x + eijy*w.w_y);
		}

		particles[i].T_t += m_alpha*T_lapl;
	}
}


void thermal::set_method(thermal_solver solver) {
	m_thermal_solver = solver;
}

void thermal::conduction(body &body) const {

	switch (m_thermal_solver) {
	case thermal_pse:
		heat_conduction_pse(body);
		break;
	case thermal_brookshaw:
		heat_conduction_brookshaw(body);
		break;
	}
}

thermal::thermal(physical_constants pc) {
	assert(pc.tc().k() != 0.);
	m_alpha = pc.tc().k()/(pc.rho0()*pc.tc().cp());
}
