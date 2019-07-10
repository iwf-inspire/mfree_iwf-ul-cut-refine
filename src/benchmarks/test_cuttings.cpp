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

#include "test_cuttings.h"

 body *cutting_ref_mr(unsigned int ny) {
	physical_constants physical_constants = matlib_tial6v4_Sima_tanh2010_cm_musec_g();

	double speed = 83.333328*1e-5;
	double mu_fric = 0.35;

	double hi_x = 0.100; double hi_y =  0.060;
	double lo_x = 0.000; double lo_y =  0.030;

	double dy = (hi_y-lo_y)/(ny-1);
	double dx = dy;
	unsigned int nx = (hi_x-lo_x)/dx;
	double hdx = 1.5;

	double c0 = physical_constants.c0();
	double dt = 0.1*dx*hdx/(speed + c0);
	double t_final =  0.1/speed*0.5; // 1mm of cut

	simulation_time *time = &simulation_time::getInstance();
	time->set_t_final(t_final);
	time->set_dt(dt);

	printf("using timestep %e with %d particles\n", dt, nx*ny);

	particle *particles = new particle[nx*ny];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nx; i++) {
		for (unsigned int j = 0; j < ny; j++) {
			double px = i*dx; double py = j*dx;

			particles[part_iter] = particle(part_iter);
			particles[part_iter].x = px + lo_x;
			particles[part_iter].y = py + lo_y;

			particles[part_iter].X = px + lo_x;
			particles[part_iter].Y = py + lo_y;

			part_iter++;
		}
	}

	double rho0 = physical_constants.rho0();
	double T0 = physical_constants.jc().Tref();

	unsigned int n = nx*ny;
	for (unsigned int i = 0; i < n; i++) {
		particles[i].rho = rho0;
		particles[i].h = hdx*dx;
		particles[i].m = dx*dx*rho0;
		particles[i].T = T0;
		particles[i].T_init = T0;

		// fix bottom
		particles[i].fixed = (particles[i].y < lo_y + 0.5*dy) ? true : false;
	}

	// correction constants
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

	// set simulation data
	simulation_data sim_data(physical_constants, correction_constants);

	// plasticity model
	plasticity *plast = new plasticity(new johnson_cook_Sima_2010(physical_constants)); // JC-tanh Sima / Özel 2010
	plast->set_tolerance(1e-6);
	plast->set_dissipation_considered(true);

	// create the body
	body *b = new body(particles, n, sim_data);
	b->set_plasticity(plast);

	double nudge = -0.001;

	glm::dvec2 tl(-0.05      + nudge, 0.0986074);
	glm::dvec2 tr(-0.0086824 + nudge, 0.0986074);
	glm::dvec2 br(-0.0000    + nudge, 0.0486);
	glm::dvec2 bl(-0.05      + nudge, 0.0555074);

	tool *t = new tool(tl, tr, br, bl, 20e-6*100, mu_fric);
	b->set_tool(t);
	t->set_vel(glm::dvec2(speed,0.));

	global_logger = new logger("cutting");
	global_logger->set_tool(t);
	global_logger->set_log_vtk(true);

	return b;
}

 body *cutting_ref_single_resol(unsigned int nbox) {
	/*
	 * ===========================================================
	 * according to (6.3) simulation by Sima & Ozel 2010 -> p. 955
	 * ===========================================================
	 */

	// MODEL 1 & 4 from the paper
	// ----------------------------------------------------------
	// choose your desired material model as the following:
	// ----------------------------------------------------------
	physical_constants pc = matlib_tial6v4_Sima_tanh2010_SI();
	// ----------------------------------------------------------

	bool thermal_conduction = true;
	double hdx = 1.5;
	double rho0 = pc.rho0();
	double T0 = 300.0;
	double thermal_diffusivity = pc.tc().k()/(rho0*pc.tc().cp());

	// workpiece dimensions SI
	double lo_x = 0.00000; double lo_y = 0.00030;
	double hi_x = 0.00200; double hi_y = 0.00060;

	unsigned int ny = nbox;
	double dy = (hi_y-lo_y)/(ny-1);
	double dx = dy;
	unsigned int nx = (hi_x-lo_x)/dx + 1;
	double vc = 500./60.;		// m/min -> m/s
	double nudge = -dx;

	// time settings
	double lc = 1e-3; // 1mm of cut
	double t_final =  lc/vc;
	double dt_empirical = (nbox < 35) ? 1.0e-9 : 5.0e-10;
	double mech_CFL = 0.5*hdx*dx/(pc.c0() + vc);
	double heat_CFL = 0.4*dx*dx/(thermal_diffusivity);
	double dt_mech = fmin(dt_empirical, 0.50*mech_CFL);
	double dt_heat = fmin(dt_empirical, 0.50*heat_CFL);
	double dt = fmin(dt_mech,dt_heat);

	simulation_time *time = &simulation_time::getInstance();
	time->set_t_final(t_final);
	time->set_dt(dt);

	particle *particles  = new particle[nx*ny];

	srand(0);
	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nx; i++) {
		for (unsigned int j = 0; j < ny; j++) {
			double px = i*dx; double py = j*dx;

			particles[part_iter] = particle(part_iter);
			particles[part_iter].x = px + lo_x;
			particles[part_iter].y = py + lo_y;

			particles[part_iter].X = px + lo_x;
			particles[part_iter].Y = py + lo_y;

			part_iter++;
		}
	}

	unsigned int n = nx*ny;
	for (unsigned int i = 0; i < n; i++) {
		particles[i].rho = rho0;
		particles[i].h = hdx*dx;
		particles[i].m = dx*dx*rho0;
		particles[i].T = T0;
		particles[i].T_init = T0;

		// fixtures
		particles[i].fixed = (particles[i].y < lo_y + 0.5*dy) ? true : false;
		particles[i].fixed = particles[i].fixed || (particles[i].x > hi_x - 0.5*dx);
	}

	// correction constants
	double alpha = 1.;
	double beta  = 1.;
	double eta   = 0.1;

	double art_stress_eps = 0.3;
	kernel_result w = cubic_spline(0, 0, dx, 0, hdx*dx);
	double wdeltap = w.w;
	double stress_exponent = 4.;

	double xsph_eps = 0.5;

	correction_constants cs(constants_monaghan(wdeltap, stress_exponent, art_stress_eps),
			constants_artificial_viscosity(alpha, beta, eta), xsph_eps);

	// set simulation data
	simulation_data sim_data(pc, cs);

	// plasticity model
	plasticity *plast = new plasticity(new johnson_cook_Sima_2010(pc)); // JC-tanh Sima / Özel 2010
	plast->set_tolerance(1e-6);
	plast->set_dissipation_considered(true);

	// create the body
	body *b = new body(particles, n, sim_data);

	// save thermal setting to body
	thermal *trml = new thermal(pc);
	trml->set_method(thermal::thermal_solver::thermal_pse);  // optional: thermal_brookshaw

	// tool settings
	float_t rake = 0.00001;
	float_t clear = 11.;
	glm::dvec2 tl(-0.000410  + nudge, 0.000986074);
	float_t length_tool = -0.000086824 - -0.000500000;
	float_t height_tool =  0.000986074 -  0.000555074;

	double mu_friction = 0.35;
	double fillet_radius = 5e-6;
	tool *t = new tool(tl, length_tool, height_tool, rake, clear, fillet_radius, mu_friction);

	double target_feed = 1e-4;	// 0.1 mm
	double current_feed = hi_y - t->low();
	double dist_to_target_feed = fabs(current_feed - target_feed);
	double correction_time = dist_to_target_feed / vc;
	double sign = (current_feed > target_feed) ? 1 : -1.;
	t->set_vel(glm::dvec2(0.,vc));
	t->update_tool(correction_time*sign);
	t->set_vel(glm::dvec2(vc, 0.)); // set actual velocity

	// save settings to body
	b->set_plasticity(plast);
	if (thermal_conduction) b->set_thermal(trml);
	b->set_tool(t);

	global_logger = new logger("cutting");
	global_logger->set_tool(t);
	global_logger->set_log_vtk(true);

	printf("feed: %f, dt %e, num_part %d\n", hi_y - t->low(), dt, nx*ny);
	printf("<<< single-resolution simulation >>>\n");

	return b;
}

 body *cutting_ref_multi_resol_apriori(unsigned int nbox) {
	/*
	 * ===========================================================
	 * according to (6.3) simulation by Sima & Ozel 2010 -> p. 955
	 * ===========================================================
	 */

	// MODEL 3 from the paper
	// ----------------------------------------------------------
	// choose your desired material model as the following:
	// ----------------------------------------------------------
	physical_constants pc = matlib_tial6v4_Sima_tanh2010_SI();
	// ----------------------------------------------------------

	bool thermal_conduction = true;
	double hdx = 1.5;
	double rho0 = pc.rho0();
	double T0 = 300.0;
	double thermal_diffusivity = pc.tc().k()/(rho0*pc.tc().cp());

	// workpiece dimensions SI
	double lo_x = 0.00000; double lo_y = 0.00030;
	double hi_x = 0.00200; double hi_y = 0.00060;
	double lx = hi_x - lo_x;
	double ly = hi_y - lo_y;

	unsigned int ny = nbox;
	double dy = ly/(ny-1);
	double dx = dy;
	unsigned int nx = lx/dx + 1;
	double vc = 500./60.;		// m/min -> m/s
	double nudge = -dx;

	// multi-resolution setup
	double resol_ratio = 2.0;
	double py_split = 0.5*ly + lo_y;
	double dxh = dx;
	double dxl = dxh * resol_ratio;
	unsigned int nxh = nx;
	unsigned int nyh = ny;
	unsigned int nxl = lx/dxl + 1;
	unsigned int nyl = ly/dxl + 1;
	double dVl = dxl *dxl;
	double dVh = dxh *dxh;
	double h0l = hdx * dxl;
	double h0h = hdx * dxh;

	// time settings
	double lc = 1e-3; // 1mm of cut
	double t_final =  lc/vc;
	double dt_empirical = (nbox < 35) ? 1.0e-9 : 5.0e-10;
	double mech_CFL = 0.5*hdx*dx/(pc.c0() + vc);
	double heat_CFL = 0.4*dx*dx/(thermal_diffusivity);
	double dt_mech = fmin(dt_empirical, 0.50*mech_CFL);
	double dt_heat = fmin(dt_empirical, 0.50*heat_CFL);
	double dt = fmin(dt_mech,dt_heat);

	simulation_time *time = &simulation_time::getInstance();
	time->set_t_final(t_final);
	time->set_dt(dt);

	particle *particles  = new particle[nxh*nyh];

	srand(0);
	unsigned int part_iter = 0;

	// 1. create the high resolution region
	for (unsigned int i = 0; i < nxh; i++) {
		for (unsigned int j = 0; j < nyh; j++) {
			double pxh = i*dxh;
			double pyh = j*dxh;

			if ((pyh+lo_y)<(py_split-1.1*dxh)) continue;

			particles[part_iter] = particle(part_iter);
			particles[part_iter].x = pxh + lo_x;
			particles[part_iter].y = pyh + lo_y;
			particles[part_iter].X = pxh + lo_x;
			particles[part_iter].Y = pyh + lo_y;

			particles[part_iter].refine_step = 1;

			/*
			 * high res
			 *
			   +++++++++++++++++
			   +               +
			   +-------CL------+

			 */

			part_iter++;
		}
	}

	// 2. create the low resolution region
	for (unsigned int i = 0; i < nxl; i++) {
		for (unsigned int j = 0; j < nyl; j++) {
			double pxl = i*dxl;
			double pyl = j*dxl;

			if ((pyl+lo_y)>=py_split) continue;

			particles[part_iter] = particle(part_iter);
			particles[part_iter].x = pxl + lo_x;
			particles[part_iter].y = pyl + lo_y;
			particles[part_iter].X = pxl + lo_x;
			particles[part_iter].Y = pyl + lo_y;

			particles[part_iter].refine_step = 0;

			/*
			 * low res
			 *

			   +-------CL------+
			   +               +
			   +++++++++++++++++

			 */

			part_iter++;
		}
	}

	// total #particles
	unsigned int n = part_iter;

	printf("n_single_resolution=%d   n_multi_resolution=%d   \n",nxh*nyh,n);

	for (unsigned int i = 0; i < n; i++) {
		particles[i].rho = rho0;
		particles[i].T = T0;
		particles[i].T_init = T0;
		particles[i].h = (particles[i].refine_step!=0) ? h0h : h0l;
		particles[i].m = (particles[i].refine_step!=0) ? dVh*rho0 : dVl*rho0;
		particles[i].split = false;
		particles[i].merge = false;

		// fixtures
		particles[i].fixed = (particles[i].y < lo_y + 0.5*dxl) ? true : false;
		particles[i].fixed = particles[i].fixed || (particles[i].x > hi_x - 0.5*dxh);
	}

	// correction constants
	double alpha = 1.;
	double beta  = 1.;
	double eta   = 0.1;

	double art_stress_eps = 0.3;
	kernel_result w = cubic_spline(0, 0, dx, 0, hdx*dx);
	double wdeltap = w.w;
	double stress_exponent = 4.;

	double xsph_eps = 0.5;

	correction_constants cs(constants_monaghan(wdeltap, stress_exponent, art_stress_eps),
			constants_artificial_viscosity(alpha, beta, eta), xsph_eps);

	// set simulation data
	simulation_data sim_data(pc, cs);

	// plasticity model
	plasticity *plast = new plasticity(new johnson_cook_Sima_2010(pc)); // JC-tanh Sima / Özel 2010
	plast->set_tolerance(1e-6);
	plast->set_dissipation_considered(true);

	// create the body
	body *b = new body(particles, n, sim_data);

	// save thermal setting to body
	thermal *trml = new thermal(pc);
	trml->set_method(thermal::thermal_solver::thermal_pse);  // optional: thermal_brookshaw

	// tool settings
	float_t rake = 0.00001;
	float_t clear = 11.;
	glm::dvec2 tl(-0.000410  + nudge, 0.000986074);
	float_t length_tool = -0.000086824 - -0.000500000;
	float_t height_tool =  0.000986074 -  0.000555074;

	double mu_friction = 0.35;
	double fillet_radius = 5e-6;
	tool *t = new tool(tl, length_tool, height_tool, rake, clear, fillet_radius, mu_friction);

	double target_feed = 1e-4;	// 0.1 mm
	double current_feed = hi_y - t->low();
	double dist_to_target_feed = fabs(current_feed - target_feed);
	double correction_time = dist_to_target_feed / vc;
	double sign = (current_feed > target_feed) ? 1 : -1.;
	t->set_vel(glm::dvec2(0.,vc));
	t->update_tool(correction_time*sign);
	t->set_vel(glm::dvec2(vc, 0.)); // set actual velocity

	// save settings to body
	b->set_plasticity(plast);
	if (thermal_conduction) b->set_thermal(trml);
	b->set_tool(t);

	global_logger = new logger("cutting");
	global_logger->set_tool(t);
	global_logger->set_log_vtk(true);

	printf("feed: %f, dt %e, num_part %d\n", hi_y - t->low(), dt, n);
	printf("<<< a-priori refinement model >>>\n");

	return b;
}

 body *cutting_ref_multi_resol_dynamic(unsigned int nbox) {
	/*
	 * ===========================================================
	 * according to (6.3) simulation by Sima & Ozel 2010 -> p. 955
	 * ===========================================================
	 */

	// MODEL 2 from the paper
	// ----------------------------------------------------------
	// choose your desired material model as the following:
	// ----------------------------------------------------------
	physical_constants pc = matlib_tial6v4_Sima_tanh2010_SI();
	// ----------------------------------------------------------

	bool thermal_conduction = true;
	double hdx = 1.5;
	double rho0 = pc.rho0();
	double T0 = 300.0;
	double thermal_diffusivity = pc.tc().k()/(rho0*pc.tc().cp());

	// workpiece dimensions SI
	double lo_x = 0.00000; double lo_y = 0.00030;
	double hi_x = 0.00200; double hi_y = 0.00060;
	double lx = hi_x - lo_x;
	double ly = hi_y - lo_y;

	unsigned int ny = nbox;
	double dy = ly/(ny-1);
	double dx = dy;
	unsigned int nx = lx/dx + 1;
	double vc = 500./60.;		// m/min -> m/s
	double nudge = -dx;

	// multi-resolution setup
	double resol_ratio = 2.0;
	double py_split = 0.5*ly + lo_y;
	double dxh = dx;
	double dxl = dxh * resol_ratio;
	unsigned int nxh = nx;
	unsigned int nyh = ny;
	unsigned int nxl = lx/dxl + 1;
	unsigned int nyl = ly/dxl + 1;
	double dVl = dxl *dxl;
	double dVh = dxh *dxh;
	double h0l = hdx * dxl;
	double h0h = hdx * dxh;

	// time settings
	double lc = 1e-3; // 1mm of cut
	double t_final =  lc/vc;
	double dt_empirical = (nbox < 35) ? 1.0e-9 : 5.0e-10;
	double mech_CFL = 0.5*hdx*dx/(pc.c0() + vc);
	double heat_CFL = 0.4*dx*dx/(thermal_diffusivity);
	double dt_mech = fmin(dt_empirical, 0.50*mech_CFL);
	double dt_heat = fmin(dt_empirical, 0.50*heat_CFL);
	double dt = fmin(dt_mech,dt_heat);

	simulation_time *time = &simulation_time::getInstance();
	time->set_t_final(t_final);
	time->set_dt(dt);

	particle *particles  = new particle[nxh*nyh];

	srand(0);
	unsigned int part_iter = 0;

	// 1. create the high resolution region
	for (unsigned int i = 0; i < nxh; i++) {
		for (unsigned int j = 0; j < nyh; j++) {
			double pxh = i*dxh;
			double pyh = j*dxh;

			if ((pyh+lo_y)<(py_split-1.9*dxh) || pxh>0.000117) continue;

			particles[part_iter] = particle(part_iter);
			particles[part_iter].x = pxh + lo_x;
			particles[part_iter].y = pyh + lo_y;
			particles[part_iter].X = pxh + lo_x;
			particles[part_iter].Y = pyh + lo_y;

			particles[part_iter].refine_step = 1;

			/*
			 * high res
			 *
			   +---------------+
			   + HR |          +
			   +-----    LR    +
			   +               +
			   +---------------+

			 */

			part_iter++;
		}
	}

	// 2. create the low resolution region
	for (unsigned int i = 0; i < nxl; i++) {
		for (unsigned int j = 0; j < nyl; j++) {
			double pxl = i*dxl;
			double pyl = j*dxl;

			if ((pyl+lo_y)>=(py_split-1.9*dxh) && pxl<=0.000117) continue;

			particles[part_iter] = particle(part_iter);
			particles[part_iter].x = pxl + lo_x;
			particles[part_iter].y = pyl + lo_y;
			particles[part_iter].X = pxl + lo_x;
			particles[part_iter].Y = pyl + lo_y;

			particles[part_iter].refine_step = 0;

			part_iter++;
		}
	}

	// total #particles
	unsigned int n = part_iter;

	// slight modification for reserved CHILD particles!
	for (unsigned int i = n; i < nxh*nyh; i++) {
		particles[part_iter] = particle(part_iter);
	}

	printf("n_single_resolution=%d   n_current=%d  \n",nxh*nyh,n);

	for (unsigned int i = 0; i < n; i++) {
		particles[i].rho = rho0;
		particles[i].T = T0;
		particles[i].T_init = T0;
		particles[i].h = (particles[i].refine_step!=0) ? h0h : h0l;
		particles[i].m = (particles[i].refine_step!=0) ? dVh*rho0 : dVl*rho0;
		particles[i].split = false;
		particles[i].merge = false;

		// fixtures
		particles[i].fixed = (particles[i].y < lo_y + 0.5*dxl) ? true : false;
		particles[i].fixed = particles[i].fixed || (particles[i].x > hi_x - 0.5*dxh);
	}

	// correction constants
	double alpha = 1.;
	double beta  = 1.;
	double eta   = 0.1;

	double art_stress_eps = 0.3;
	kernel_result w = cubic_spline(0, 0, dx, 0, hdx*dx);
	double wdeltap = w.w;
	double stress_exponent = 4.;

	double xsph_eps = 0.5;

	correction_constants cs(constants_monaghan(wdeltap, stress_exponent, art_stress_eps),
			constants_artificial_viscosity(alpha, beta, eta), xsph_eps);

	// set simulation data
	simulation_data sim_data(pc, cs);

	// plasticity model
	plasticity *plast = new plasticity(new johnson_cook_Sima_2010(pc)); // JC-tanh Sima / Özel 2010
	plast->set_tolerance(1e-6);
	plast->set_dissipation_considered(true);

	// create the body
	body *b = new body(particles, n, sim_data);

	// save thermal setting to body
	thermal *trml = new thermal(pc);
	trml->set_method(thermal::thermal_solver::thermal_pse);  // optional: thermal_brookshaw

	// adaptivity settings

	// default settings +-+-++-+-+-+-+-+-+-+-+-+-+-
	double alpha_dx = 0.50;
	double beta_h = 0.50;
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
	double l_eff = lc + 0.1*lx;
	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

	adaptivity *adapt = new adaptivity(alpha_dx, beta_h, v_cr, div_v_cr, SvM_cr, eps_cr,
									   T_cr, xy_min, xy_max, frame_width, frame_height,
									   n_nbh, l_eff, true);

	adapt->set_refine_criterion(adaptivity::refine_criteria::moving_frame);
	adapt->set_refine_pattern(adaptivity::pattern::cubic_basic);

	// tool settings
	float_t rake = 0.00001;
	float_t clear = 11.;
	glm::dvec2 tl(-0.000410  + nudge, 0.000986074);
	float_t length_tool = -0.000086824 - -0.000500000;
	float_t height_tool =  0.000986074 -  0.000555074;

	double mu_friction = 0.35;
	double fillet_radius = 5e-6;
	tool *t = new tool(tl, length_tool, height_tool, rake, clear, fillet_radius, mu_friction);

	double target_feed = 1e-4;	// 0.1 mm
	double current_feed = hi_y - t->low();
	double dist_to_target_feed = fabs(current_feed - target_feed);
	double correction_time = dist_to_target_feed / vc;
	double sign = (current_feed > target_feed) ? 1 : -1.;
	t->set_vel(glm::dvec2(0.,vc));
	t->update_tool(correction_time*sign);
	t->set_vel(glm::dvec2(vc, 0.)); // set actual velocity
	t->set_edge_coord(glm::dvec2(0., hi_y - target_feed)); // set coord edge

	// save settings to body
	b->set_plasticity(plast);
	if (thermal_conduction) b->set_thermal(trml);
	b->set_tool(t);
	b->set_adaptivity(adapt);

	global_logger = new logger("cutting");
	global_logger->set_tool(t);
	global_logger->set_log_vtk(true);

	printf("feed: %f, dt %e, num_part %d\n", hi_y - t->low(), dt, n);
	printf("<<< dynamic refinement model >>>\n");

	return b;
}
