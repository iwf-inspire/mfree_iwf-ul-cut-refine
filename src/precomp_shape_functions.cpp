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

#include "precomp_shape_functions.h"

void precomp_sph(std::vector<particle> &particles, unsigned int n) {

	for (unsigned int i = 0; i < n; i++) {
		double xi = particles[i].x;
		double yi = particles[i].y;

		double hi = particles[i].h;

		for (unsigned int j = 0; j < particles[i].num_nbh; j++) {
			unsigned int jdx = particles[i].nbh[j];

			double xj = particles[jdx].x;
			double yj = particles[jdx].y;

			particles[i].w[j] = cubic_spline(xi, yi, xj, yj, hi);
		}
	}
}

void precomp_cspm(std::vector<particle> &particles, unsigned int n) {

	for (unsigned int i = 0; i < n; i++) {
		double xi = particles[i].x;
		double yi = particles[i].y;

		double hi = particles[i].h;
		glm::dmat2x2 B(0.); // This is in fact the symbol "A_i" in Eq. (29) of the paper

		for (unsigned int j = 0; j < particles[i].num_nbh; j++) {
			unsigned int jdx = particles[i].nbh[j];

			double xj = particles[jdx].x;
			double yj = particles[jdx].y;

			double volj = particles[jdx].m/particles[jdx].rho;

			kernel_result w = cubic_spline(xi, yi, xj, yj, hi);

			B[0][0] += (xj - xi)*w.w_x*volj;
			B[1][0] += (xj - xi)*w.w_y*volj;
			B[0][1] += (yj - yi)*w.w_x*volj;
			B[1][1] += (yj - yi)*w.w_y*volj;
		}

		glm::dmat2x2 invB = glm::inverse(B);

		for (unsigned int j = 0; j < particles[i].num_nbh; j++) {
			unsigned int jdx = particles[i].nbh[j];

			double xj = particles[jdx].x;
			double yj = particles[jdx].y;

			kernel_result w = cubic_spline(xi, yi, xj, yj, hi);

			particles[i].w[j].w = w.w;
			particles[i].w[j].w_x = (w.w_x*invB[0][0] + w.w_y*invB[1][0]);
			particles[i].w[j].w_y = (w.w_x*invB[0][1] + w.w_y*invB[1][1]);
		}
	}
}
