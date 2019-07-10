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

#include "derivatives.h"

void derive_velocity(body &b) {
	std::vector<particle> &particles = b.get_particles();
	unsigned int n = b.get_num_part();

	for (unsigned int i = 0; i < n; i++) {
		double vxi = particles[i].vx;
		double vyi = particles[i].vy;

		double vx_x = 0.;
		double vx_y = 0.;
		double vy_x = 0.;
		double vy_y = 0.;

		for (unsigned int j = 0; j < particles[i].num_nbh; j++) {
			unsigned int jdx = particles[i].nbh[j];
			kernel_result w = particles[i].w[j];

			double vxj = particles[jdx].vx;
			double vyj = particles[jdx].vy;

			double quad_weight = particles[jdx].m/particles[jdx].rho;

			vx_x += (vxj-vxi)*w.w_x*quad_weight;
			vx_y += (vxj-vxi)*w.w_y*quad_weight;
			vy_x += (vyj-vyi)*w.w_x*quad_weight;
			vy_y += (vyj-vyi)*w.w_y*quad_weight;
		}

		particles[i].vx_x = vx_x;
		particles[i].vx_y = vx_y;
		particles[i].vy_x = vy_x;
		particles[i].vy_y = vy_y;
	}
}

void derive_stress_monaghan(body &b) {
	const double wdeltap  = b.get_sim_data().get_correction_constants().get_monaghan_const().mghn_wdeltap();
	const unsigned int corr_exp = b.get_sim_data().get_correction_constants().get_monaghan_const().mghn_corr_exp();

	std::vector<particle> &particles = b.get_particles();
	unsigned int n = b.get_num_part();

	for (unsigned int i = 0; i < n; i++) {
		double Sxxi = particles[i].Sxx - particles[i].p;
		double Sxyi = particles[i].Sxy;
		double Syyi = particles[i].Syy - particles[i].p;

		double Rxxi = particles[i].Rxx;
		double Rxyi = particles[i].Rxy;
		double Ryyi = particles[i].Ryy;

		double rhoi = particles[i].rho;
		double rhoi21 = 1./(rhoi*rhoi);

		double Sxx_x = 0.;
		double Sxy_y = 0.;
		double Sxy_x = 0.;
		double Syy_y = 0.;

		for (unsigned int j = 0; j < particles[i].num_nbh; j++) {
			unsigned int jdx = particles[i].nbh[j];
			kernel_result w = particles[i].w[j];

			double Sxxj = particles[jdx].Sxx - particles[jdx].p;
			double Sxyj = particles[jdx].Sxy;
			double Syyj = particles[jdx].Syy - particles[jdx].p;

			double Rxxj = particles[jdx].Rxx;
			double Rxyj = particles[jdx].Rxy;
			double Ryyj = particles[jdx].Ryy;

			double mj     = particles[jdx].m;
			double rhoj21 = 1./(particles[jdx].rho*particles[jdx].rho);

			double Rxx = 0;
			double Rxy = 0.;
			double Ryy = 0.;

			if (wdeltap > 0 && particles[i].idx != particles[jdx].idx) {
				double fab = w.w/wdeltap;
//				fab = pow(fab,corr_exp);	//dramatically increase performance by for loop!
				double t = 1.;
				for (unsigned int powi = 0; powi < corr_exp; powi++) {
					t = t*fab;
				}
				fab = t;

				Rxx = fab*(Rxxi + Rxxj);
				Rxy = fab*(Rxyi + Rxyj);
				Ryy = fab*(Ryyi + Ryyj);
			}

			Sxx_x += mj*(Sxxi*rhoi21 + Sxxj*rhoj21 + Rxx)*w.w_x;
			Sxy_y += mj*(Sxyi*rhoi21 + Sxyj*rhoj21 + Rxy)*w.w_y;
			Sxy_x += mj*(Sxyi*rhoi21 + Sxyj*rhoj21 + Rxy)*w.w_x;
			Syy_y += mj*(Syyi*rhoi21 + Syyj*rhoj21 + Ryy)*w.w_y;
		}

		particles[i].Sxx_x = Sxx_x*rhoi;
		particles[i].Sxy_y = Sxy_y*rhoi;
		particles[i].Sxy_x = Sxy_x*rhoi;
		particles[i].Syy_y = Syy_y*rhoi;
	}
}
