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

#include "test_density.h"

static body *test_bench_setup_refine_density(unsigned int nbox) {
	// material constants
	double rho0 = 1.0;

	//problem dimensions
	double L = 1.0;
	double dx = L/(nbox-1);
	double hdx = 1.;
	double dt = 0.;

	particle *particles = new particle[nbox*nbox];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = i*dx;
			double py = j*dx;
			particles[part_iter] = particle(part_iter);
			particles[part_iter].x = px;
			particles[part_iter].y = py;
			particles[part_iter].X = px;
			particles[part_iter].Y = py;

			part_iter++;
		}
	}
	unsigned int n = part_iter;

	for (unsigned int i = 0; i < n; i++) {
		particles[i].rho = rho0;
		particles[i].h = hdx*dx;
		particles[i].m = dx*dx*rho0;
		particles[i].refine_step = 0;
		particles[i].split = false;
		particles[i].merge = false;
	}

	//correction constants (monaghan & gray)
	double alpha = 1.;
	double beta  = 1.;
	double eta   = 0.1;

	double art_stress_eps = 0.3;
	kernel_result w = cubic_spline(0, 0, dx, 0, hdx*dx);
	double wdeltap = w.w;
	double stress_exponent = 4.;
	double xsph_eps = 0.5;

	physical_constants physical_constants(0, 1, rho0);

	correction_constants correction_constants(constants_monaghan(wdeltap, stress_exponent, art_stress_eps),
			constants_artificial_viscosity(alpha, beta, eta), xsph_eps);

	simulation_data sim_data(physical_constants, correction_constants);

	body *b = new body(particles, n, sim_data);

	simulation_time *time = &simulation_time::getInstance();
	time->set_t_final(dt);
	time->set_dt(dt);

	global_logger = new logger("square_density");

	return b;
}

void run_test_density_refinement_error(unsigned int nbox, adaptivity::pattern pattern) {
	// reproducing the global refinement test by J. Feldman & J. Bonet [2006]
	// NOTE:
	// comment in the following ONLY for reconstruction of preliminary test
	// Fig. 7 and Fig. 8 from the paper.

	body *b = test_bench_setup_refine_density(nbox);

	// default settings +-+-++-+-+-+-+-+-+-+-+-+-+-
	double alpha_dx = 0.40;
	double beta_h = 0.60;
	double v_cr = 0.40;
	double div_v_cr = 2e+5;
	double SvM_cr = 1e+7;
	double eps_cr = 110;
	double T_cr = 700.;
	glm::dvec2 xy_min = {0.25, 0.25};
	glm::dvec2 xy_max = {0.75, 0.75};
	double frame_width =  0.000350;
	double frame_height = 0.000060;
	unsigned int n_nbh = 10;
	double l_eff = +DBL_MAX;
	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

	adaptivity *adapt = new adaptivity(alpha_dx, beta_h, v_cr, div_v_cr, SvM_cr, eps_cr,
									   T_cr, xy_min, xy_max, frame_width, frame_height,
									   n_nbh, l_eff, false);

	adapt->set_refine_criterion(adaptivity::refine_criteria::position);
	adapt->set_refine_pattern(pattern);

	printf("#particles before refinement = %d\n", (*b).get_num_part());
	global_logger->log(*b, 0);

	adapt->adapt_resolution(*b);

	adapt->dens_before_approx_N2(*b);
	adapt->dens_after_approx_N2(*b);

	printf("#particles after refinement = %d\n", (*b).get_num_part());
	global_logger->log(*b, 1);

	printf("<<< DENSITY TEST COMPLETED! >>> \n");
	exit(0);
}
