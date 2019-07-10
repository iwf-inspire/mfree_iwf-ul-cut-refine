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

#include "test_benches.h"

body *test_bench_setup_rings(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions (monaghan & gray)
	double ri = 0.03;
	double ro = 0.04;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.7;

	double vel_rings =  180.;

	double c0 = physical_constants.c0();
	double dt = 0.1*dx*hdx/(vel_rings + c0);

	printf("using timestep %e\n", dt);

	particle *particles = new particle[nbox*nbox];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = -ro+i*dx; double py = -ro+j*dx;
			double dist = sqrt(px*px + py*py);
			if (dist < ro && dist >= ri) {
				particles[part_iter] = particle(part_iter);
				particles[part_iter].x = px-spacing;
				particles[part_iter].y = py;
				part_iter++;
				particles[part_iter] = particle(part_iter);
				particles[part_iter].x = px+spacing;
				particles[part_iter].y = py;
				part_iter++;
			}
		}
	}

	unsigned int n = part_iter;

	for (unsigned int i = 0; i < n; i++) {
		particles[i].rho = rho0;
		particles[i].h = hdx*dx;
		particles[i].m = dx*dx*rho0;
		particles[i].vx = (particles[i].x < 0.) ? vel_rings : -vel_rings;
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

	correction_constants correction_constants(constants_monaghan(wdeltap, stress_exponent, art_stress_eps),
			constants_artificial_viscosity(alpha, beta, eta), xsph_eps);

	simulation_data sim_data(physical_constants, correction_constants);

	body *b = new body(particles, n, sim_data);

	simulation_time *time = &simulation_time::getInstance();
	time->set_t_final(6e3*dt);
	time->set_dt(dt);

	global_logger = new logger("rings");

	return b;
}

body *test_bench_setup_ring_contact(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions (monaghan & gray)
	double ri = 0.03;
	double ro = 0.04;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.7;

	double vel_rings =  180.;

	double c0 = physical_constants.c0();
	double dt = 0.1*dx*hdx/(vel_rings + c0);

	printf("using timestep %e\n", dt);

	particle *particles = new particle[nbox*nbox];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = -ro+i*dx; double py = -ro+j*dx;
			double dist = sqrt(px*px + py*py);
			if (dist < ro && dist >= ri) {
				particles[part_iter] = particle(part_iter);
				particles[part_iter].x = px-spacing;
				particles[part_iter].y = py;
				part_iter++;
			}
		}
	}

	unsigned int n = part_iter;

	for (unsigned int i = 0; i < n; i++) {
		particles[i].rho = rho0;
		particles[i].h = hdx*dx;
		particles[i].m = dx*dx*rho0;
		particles[i].vx = vel_rings;
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

	correction_constants correction_constants(constants_monaghan(wdeltap, stress_exponent, art_stress_eps),
			constants_artificial_viscosity(alpha, beta, eta), xsph_eps);

	simulation_data sim_data(physical_constants, correction_constants);

	body *b = new body(particles, n, sim_data);

	simulation_time *time = &simulation_time::getInstance();
	time->set_t_final(6e3*dt);
	time->set_dt(dt);

	glm::dvec2 bl(0.,    -2*ro);
	glm::dvec2 br(10*dx, -2*ro);
	glm::dvec2 tr(10*dx, +2*ro);
	glm::dvec2 tl(0.,    +2*ro);

	tool *t = new tool(tl, tr, br, bl, 0.);
	b->set_tool(t);

	global_logger = new logger("rings");

	return b;
}

body *test_bench_setup_disk_impact(unsigned int nbox) {
	physical_constants physical_constants = matlib_tial6v4_Sima_tanh2010_SI();
	double rho0 = physical_constants.rho0();

	//problem dimensions (monaghan & gray)
	double ro = 0.04;
	double spacing = ro + ro/40;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.7;

	double vel_rings = 180.;

	double c0 = physical_constants.c0();
	double dt = 0.1*dx*hdx/(vel_rings + c0);

	printf("using timestep %e\n", dt);

	particle *particles = new particle[nbox*nbox];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = -ro+i*dx; double py = -ro+j*dx;
			double dist = sqrt(px*px + py*py);
			if (dist < ro) {
				particles[part_iter] = particle(part_iter);
				particles[part_iter].x = px-spacing;
				particles[part_iter].y = py;
				part_iter++;
			}
		}
	}

	unsigned int n = part_iter;

	for (unsigned int i = 0; i < n; i++) {
		particles[i].rho = rho0;
		particles[i].h = hdx*dx;
		particles[i].m = dx*dx*rho0;
		particles[i].vx = vel_rings;
		particles[i].T      = physical_constants.jc().Tref();
		particles[i].T_init = physical_constants.jc().Tref();
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

	correction_constants correction_constants(constants_monaghan(wdeltap, stress_exponent, art_stress_eps),
			constants_artificial_viscosity(alpha, beta, eta), xsph_eps);

	simulation_data sim_data(physical_constants, correction_constants);

	body *b = new body(particles, n, sim_data);

	simulation_time *time = &simulation_time::getInstance();
	time->set_t_final(6e3*dt);
	time->set_dt(dt);

	glm::dvec2 bl(0.,    -2*ro);
	glm::dvec2 br(10*dx, -2*ro);
	glm::dvec2 tr(10*dx, +2*ro);
	glm::dvec2 tl(0.,    +2*ro);

	tool *t = new tool(tl, tr, br, bl, 0.);
	b->set_tool(t);

	plasticity *plast = new plasticity(new johnson_cook_Sima_2010(sim_data.get_physical_constants()));
	b->set_plasticity(plast);

	global_logger = new logger("rings");

	return b;
}

body *test_bench_setup_thermal(unsigned int nbox) {
	physical_constants physical_constants = matlib_thermal_synthetic();
	double rho0 = physical_constants.rho0();

	//problem dimensions (monaghan & gray)
	double ri = 0.03;
	double ro = 0.04;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.7;

	double dt = 1e-6;

	printf("using timestep %e\n", dt);

	particle *particles = new particle[nbox*nbox];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = -ro+i*dx; double py = -ro+j*dx;
			double dist = sqrt(px*px + py*py);
			if (dist < ro) {
				particles[part_iter] = particle(part_iter);
				particles[part_iter].x = px;
				particles[part_iter].y = py;
				part_iter++;
			}
		}
	}

	unsigned int n = part_iter;

	for (unsigned int i = 0; i < n; i++) {
		particles[i].rho = rho0;
		particles[i].h = hdx*dx;
		particles[i].m = dx*dx*rho0;
		particles[i].T      = 273;
		particles[i].T_init = 273;

		if (particles[i].x*particles[i].x + particles[i].y*particles[i].y < ri*ri) {
			particles[i].T      = 1e3;
			particles[i].T_init = 1e3;
		}
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

	correction_constants correction_constants(constants_monaghan(wdeltap, stress_exponent, art_stress_eps),
			constants_artificial_viscosity(alpha, beta, eta), xsph_eps);

	simulation_data sim_data(physical_constants, correction_constants);

	body *b = new body(particles, n, sim_data);

	simulation_time *time = &simulation_time::getInstance();
	time->set_t_final(1e3*dt);
	time->set_dt(dt);

	thermal *trml = new thermal(sim_data.get_physical_constants());
	trml->set_method(thermal::thermal_solver::thermal_brookshaw);
//	trml->set_method(thermal::thermal_solver::thermal_pse);
	b->set_thermal(trml);

	global_logger = new logger("rings");

	return b;
}
