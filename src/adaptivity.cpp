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

#include "adaptivity.h"
#include "body.h"

static bool inside_bounding_box(glm::dvec2 xlim, glm::dvec2 ylim, glm::dvec2 pos) {
	double xmin = xlim.x; double xmax = xlim.y;
	double ymin = ylim.x; double ymax = ylim.y;
	double x = pos.x;
	double y = pos.y;
	return (x>xmin) && (x<xmax) && (y>ymin) && (y<ymax);
}

void adaptivity::set_refine_criterion(refine_criteria crit) {
	m_refine_criteria = crit;
}

void adaptivity::set_refine_pattern(pattern patt) {
	m_pattern = patt;

	switch (m_pattern) {
	case triangular:
		m_num_child = (plus_one) ?  4 : 3;
		break;
	case cubic_basic:
		m_num_child = 4;
		break;
	case cubic:
		m_num_child = (plus_one) ?  5 : 4;
		break;
	case hexagonal:
		m_num_child = (plus_one) ?  7 : 6;
		break;
	}
}

void adaptivity::flag_reset(body &b) const{
	std::vector<particle> &particles = b.get_particles();
	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		particles[i].split = false;
		particles[i].merge = false;
	}
}

void copy_dad_to_son(const particle &dad, particle &son) {
	son.rho     = dad.rho;
	son.T       = dad.T;
	son.T_init  = dad.T_init;
	son.vx      = dad.vx;
	son.vy      = dad.vy;
}

void adaptivity::dens_before_approx_N2(body &b) const {
	std::vector<particle> &particles = b.get_particles();

	glm::dvec2 xlim(m_xy_min.x, m_xy_max.x);
	glm::dvec2 ylim(m_xy_min.y, m_xy_max.y);

	for (unsigned int i = 0; i < b.get_num_part(); i++) {

		double xi = particles[i].x;
		double yi = particles[i].y;
		double rho0_i = 0.;

		for (unsigned int j = 0; j < b.get_num_part(); j++) {
			if(particles[j].refine_step!=0 && !particles[j].split) continue;

			double xj = particles[j].x;
			double yj = particles[j].y;
			double hj = (particles[j].split) ? (1./m_beta)*particles[j].h : particles[j].h;
			double mj = (particles[j].split) ? m_num_child*particles[j].m : particles[j].m;

			kernel_result w = cubic_spline(xi, yi, xj, yj, hj);

			rho0_i += mj*w.w;
		}

		glm::dvec2 pos(xi,yi);

		particles[i].rho_init = (!inside_bounding_box(xlim,ylim,pos)) ? rho0_i : 1.0;
	}
}

void adaptivity::dens_after_approx_N2(body &b) const {
	std::vector<particle> &particles = b.get_particles();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {

		double xi = particles[i].x;
		double yi = particles[i].y;
		double rho_i = 0.;

		for (unsigned int j = 0; j < b.get_num_part(); j++) {

			double xj = particles[j].x;
			double yj = particles[j].y;
			double hj = particles[j].h;
			double mj = particles[j].m;

			kernel_result w  = cubic_spline(xi, yi, xj, yj, hj);

			rho_i += mj*w.w;
		}
		particles[i].rho = rho_i;
	}
}

void adaptivity::my_extrapol_SPH_values(std::vector<particle> &particles, unsigned int my_idx, unsigned int my_dady) {
	// who are YOU?
	particle pi = particles[my_idx];
	// who's your DAD?
	particle pd = particles[my_dady];

	// SON
	double xi = particles[my_idx].x;
	double yi = particles[my_idx].y;
	double hi = particles[my_idx].h;

	// all the state variables you wish to extrapolate
	double p_i   = 0.;
	double Sxx_i = 0.;
	double Sxy_i = 0.;
	double Syy_i = 0.;
	double Szz_i = 0.;
	double eps_plxx_i = 0.;
	double eps_plxy_i = 0.;
	double eps_plyy_i = 0.;
	double eps_plzz_i = 0.;
	double eps_pleq_i = 0.;

	double denom = 0.;

	// 1ST CONTRIBUTION: neighbors of your DAD
	//*********************
	for (unsigned int j = 0; j < particles[my_dady].num_nbh; j++) {
		unsigned int jdx = particles[my_dady].nbh[j];

		double xj = particles[jdx].x;
		double yj = particles[jdx].y;
		kernel_result w = cubic_spline(xi, yi, xj, yj, hi);

		denom += w.w;

		double p_j   = particles[jdx].p;
		double Sxx_j = particles[jdx].Sxx;
		double Sxy_j = particles[jdx].Sxy;
		double Syy_j = particles[jdx].Syy;
		double Szz_j = particles[jdx].Szz;
		double eps_plxx_j = particles[jdx].eps_plxx;
		double eps_plxy_j = particles[jdx].eps_plxy;
		double eps_plyy_j = particles[jdx].eps_plyy;
		double eps_plzz_j = particles[jdx].eps_plzz;
		double eps_pleq_j = particles[jdx].eps_pl_equiv;

		p_i   += p_j*w.w;
		Sxx_i += Sxx_j*w.w;
		Sxy_i += Sxy_j*w.w;
		Syy_i += Syy_j*w.w;
		Szz_i += Szz_j*w.w;
		eps_plxx_i += eps_plxx_j*w.w;
		eps_plxy_i += eps_plxy_j*w.w;
		eps_plyy_i += eps_plyy_j*w.w;
		eps_plzz_i += eps_plzz_j*w.w;
		eps_pleq_i += eps_pleq_j*w.w;
	}

	// 2ND CONTRIBUTION: your DAD himself
	//******************
	// DAD
	double xk = particles[my_dady].x;
	double yk = particles[my_dady].y;
	kernel_result w = cubic_spline(xi, yi, xk, yk, hi);

	denom += w.w;

	double p_k   = particles[my_dady].p;
	double Sxx_k = particles[my_dady].Sxx;
	double Sxy_k = particles[my_dady].Sxy;
	double Syy_k = particles[my_dady].Syy;
	double Szz_k = particles[my_dady].Szz;
	double eps_plxx_k = particles[my_dady].eps_plxx;
	double eps_plxy_k = particles[my_dady].eps_plxy;
	double eps_plyy_k = particles[my_dady].eps_plyy;
	double eps_plzz_k = particles[my_dady].eps_plzz;
	double eps_pleq_k = particles[my_dady].eps_pl_equiv;

	p_i   += p_k*w.w;
	Sxx_i += Sxx_k*w.w;
	Sxy_i += Sxy_k*w.w;
	Syy_i += Syy_k*w.w;
	Szz_i += Szz_k*w.w;
	eps_plxx_i += eps_plxx_k*w.w;
	eps_plxy_i += eps_plxy_k*w.w;
	eps_plyy_i += eps_plyy_k*w.w;
	eps_plzz_i += eps_plzz_k*w.w;
	eps_pleq_i += eps_pleq_k*w.w;

	assert(denom > 1e-8);

	// save
	particles[my_idx].p   = p_i/denom;
	particles[my_idx].Sxx = Sxx_i/denom;
	particles[my_idx].Sxy = Sxy_i/denom;
	particles[my_idx].Syy = Syy_i/denom;
	particles[my_idx].Szz = Szz_i/denom;
	particles[my_idx].eps_plxx = eps_plxx_i/denom;
	particles[my_idx].eps_plxy = eps_plxy_i/denom;
	particles[my_idx].eps_plyy = eps_plyy_i/denom;
	particles[my_idx].eps_plzz = eps_plzz_i/denom;
	particles[my_idx].eps_pl_equiv = eps_pleq_i/denom;
}

int adaptivity::scan_mark_div_velocity_based(body &b) const {
	unsigned int iter = 0;
	std::vector<particle> &particles = b.get_particles();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {

		double vx_x = particles[i].vx_x;
		double vy_y = particles[i].vy_y;

		// SCAN
		double div_vi = sqrt(vx_x*vx_x + vy_y*vy_y);
		if(div_vi < m_div_v_threshold || particles[i].refine_step >= MAX_REFINE_STEP) continue;

		// MARK
		particles[i].split = true;
		iter++;
	}
	return iter;
}

int adaptivity::scan_mark_velocity_based(body &b) const {
	unsigned int iter = 0;
	std::vector<particle> &particles = b.get_particles();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {

		double vx = particles[i].vx;
		double vy = particles[i].vy;

		// SCAN
		double vi = sqrt(vx*vx + vy*vy);
		if(vi < m_v_threshold || particles[i].refine_step >= MAX_REFINE_STEP) continue;

		// MARK
		particles[i].split = true;
		iter++;
	}
	return iter;
}

int adaptivity::scan_mark_vonmises_based(body &b) const {
	unsigned int iter = 0;
	std::vector<particle> &particles = b.get_particles();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {

		// SCAN ---
		double sxx = particles[i].Sxx - particles[i].p;
		double sxy = particles[i].Sxy;
		double sxz = 0.;
		double syy = particles[i].Syy - particles[i].p;
		double syz = 0.;
		double szz = 0.;

		double svm = sqrt((sxx*sxx + syy*syy + szz*szz) - sxx * syy - sxx * szz - syy * szz + 3.0 * (sxy*sxy + syz*syz + sxz*sxz));

		if(svm < m_SvM_threshold || particles[i].refine_step >= MAX_REFINE_STEP) continue;

		// MARK
		particles[i].split = true;
		iter++;
	}
	return iter;
}

int adaptivity::scan_mark_strain_based(body &b) const {
	unsigned int iter = 0;
	std::vector<particle> &particles = b.get_particles();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {

		// SCAN
		// calculate strain gradient
		double normPL = particles[i].eps_plxx*particles[i].eps_plxx + 2.0*particles[i].eps_plxy*particles[i].eps_plxy + particles[i].eps_plyy*particles[i].eps_plyy;
		double val = sqrt((2./3.)*normPL);

		if(val<m_eps_threshold || particles[i].refine_step>=MAX_REFINE_STEP) continue;

		// MARK
		particles[i].split = true;
		iter++;
	}
	return iter;
}

int adaptivity::scan_mark_temperature_based(body &b) const {
	unsigned int iter = 0;
	std::vector<particle> &particles = b.get_particles();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {

		// SCAN
		double Ti = particles[i].T;
		if (Ti < m_T_threshold || particles[i].refine_step >= MAX_REFINE_STEP) {
			continue;
		}

		// MARK
		particles[i].split = true;
		iter++;
	}
	return iter;
}

int adaptivity::scan_mark_position_based(body &b) const {
	unsigned int iter = 0;
	std::vector<particle> &particles = b.get_particles();

	glm::dvec2 xlim(m_xy_min.x, m_xy_max.x);
	glm::dvec2 ylim(m_xy_min.y, m_xy_max.y);

	for (unsigned int i = 0; i < b.get_num_part(); i++) {

		if(particles[i].refine_step>=MAX_REFINE_STEP) continue;

		// SCAN
		double xi = particles[i].x;
		double yi = particles[i].y;

		glm::dvec2 pos(xi,yi);
		if(!inside_bounding_box(xlim,ylim,pos)) continue;

		// MARK
		particles[i].split = true;
		iter++;
	}
	return iter;
}

int adaptivity::scan_mark_moving_frame(body &b) const {
	std::vector<particle> &particles = b.get_particles();

	// find the maximum "x" of the refined-particles
	double xm = -DBL_MAX;
	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		if(particles[i].refine_step == 0) continue;
		xm = fmax(particles[i].x,xm);
	}

	//======================================================
	simulation_time *time = &simulation_time::getInstance();
	double t = time->get_time();

	// call the initial location of the tool tip
	glm::dvec2 v_tool = b.speed_tool();
	glm::dvec2 tool_tip = b.edge_tool();

	double vc_x = v_tool.x;
	double vc_y = v_tool.y;
	double x0   = tool_tip.x;
	double y0   = tool_tip.y;

	// compute the current location of the tool tip
	double xtool = t*vc_x + x0;
	double ytool = t*vc_y + y0;

	double width_adapt = m_width;
	double nudge = 1e-5;

	glm::dvec2 ylim(ytool - m_height + nudge, +DBL_MAX);
	glm::dvec2 xlim(-DBL_MAX, xtool + width_adapt + nudge);
	//======================================================

	unsigned int iter = 0;
	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		if((particles[i].x > m_l_eff) || particles[i].refine_step>=MAX_REFINE_STEP) continue;

		// SCAN
		glm::dvec2 pos(particles[i].x,particles[i].y);
		if(!inside_bounding_box(xlim,ylim,pos)) continue;

		// MARK
		particles[i].split = true;
		iter++;
	}
	return iter;
}

int adaptivity::scan_mark_neighbor_based(body &b) const {
	unsigned int iter = 0;
	std::vector<particle> &particles = b.get_particles();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {

		// SCAN
		unsigned int num_nbh = particles[i].num_nbh;

		if(num_nbh < m_num_nbh_threshold || particles[i].refine_step >= MAX_REFINE_STEP) continue;

		// MARK
		particles[i].split = true;
		iter++;
	}
	return iter;
}

void adaptivity::perform_split_triangular(body &b) const {

	simulation_time *time = &simulation_time::getInstance();
	unsigned int step = time->get_step();

	std::vector<particle> &particles = b.get_particles();

	// how many "SON" do you have in 2D? ---> 3
	const unsigned int num_SON2D = m_num_child-1;
	double coeff_md = (1./(num_SON2D+1));
	double coeff_m0 = (1./(num_SON2D+1));

	// 1. consider a vector of SON particles
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	std::vector<particle> sons;

	// 2. give birth to SON + assign values to SON
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int son_idx = particles.size();
	for (unsigned int i = 0; i < particles.size(); i++) {

		// only append new "SON" particles if their "DAD" particle:
		// 		A) has been marked as split/merge candidate
		// 		B) has NOT already reached the MAX_REFINE_STEP
		//---------------------------------------------------------
		unsigned int delta_st = step - particles[i].last_refine_at;
		assert(delta_st>=0);
		//if(particles[i].split && particles[i].refine_step<MAX_REFINE_STEP && delta_st>MIN_REFINE_DIFF) {
		if(particles[i].split && particles[i].refine_step<MAX_REFINE_STEP) {

			double x_SON [num_SON2D];
			double y_SON [num_SON2D];
			double h_SON [num_SON2D];
			double m_SON [num_SON2D];

			// 0. call your DAD
			double dx = sqrt(particles[i].m/particles[i].rho);
			double x_DAD = particles[i].x;
			double y_DAD = particles[i].y;
			double h_DAD = particles[i].h;
			double m_DAD = particles[i].m;

			// 1. increase the refinement step of "DAD"
			particles[i].last_refine_at = step;
			particles[i].refine_step ++;


			// 2. give birth to new "SON" - following the given refinement pattern
			/*
			 *               SON[1]
			 *                 |
			 *                 |
			 *                 |
			 * 	  		      DAD
			 * 	            /	  \
			 *             /       \
			 *       SON[0]         SON[2]
			 */

			// SON[0]
			x_SON[0] = x_DAD - cos(M_PI/6.)*m_alpha*dx;
			y_SON[0] = y_DAD - sin(M_PI/6.)*m_alpha*dx;
			h_SON[0] = m_beta*h_DAD;
			m_SON[0] = (plus_one) ? coeff_md*m_DAD : (1./num_SON2D)*m_DAD;

			// SON[1]
			x_SON[1] = x_DAD;
			y_SON[1] = y_DAD + m_alpha*dx;
			h_SON[1] = m_beta*h_DAD;
			m_SON[1] = (plus_one) ? coeff_md*m_DAD : (1./num_SON2D)*m_DAD;

			// SON[2]
			x_SON[2] = x_DAD + cos(M_PI/6.)*m_alpha*dx;
			y_SON[2] = y_DAD - sin(M_PI/6.)*m_alpha*dx;
			h_SON[2] = m_beta*h_DAD;
			m_SON[2] = (plus_one) ? coeff_md*m_DAD : (1./num_SON2D)*m_DAD;

			// 3. append new "SON" to the current particle array
			for (unsigned int ii = 0; ii < num_SON2D; ii++) {
				unsigned int id_DAD = i;

				particle son_particle(son_idx);
				son_idx++;

				son_particle.x = x_SON[ii];
				son_particle.y = y_SON[ii];
				son_particle.X = son_particle.x;
				son_particle.Y = son_particle.y;
				son_particle.m = m_SON[ii];
				son_particle.h = h_SON[ii];
				son_particle.refine_step = 1;
				son_particle.last_refine_at = step;
				son_particle.split = false;
				son_particle.merge = false;

				//copy_dad_to_son(particles[id_DAD], particles[id_SON]);
				copy_dad_to_son(particles[id_DAD], son_particle);
				sons.push_back(son_particle);

				//my_extrapol_SPH_values(b.get_particles(),id_SON,id_DAD);
			}

			// 4. slight modification: DAD himself becomes a SON from now on! "+1" approach
			particles[i].m = (plus_one) ? coeff_m0*m_DAD : 1e-16;
			particles[i].h = m_beta*h_DAD;

			// sanity check - total mass conservation after refinement
			double sum_mass = particles[i].m + m_SON[0] + m_SON[1] + m_SON[2];
			assert(fabs(m_DAD-sum_mass) < 1e-12);
		}
	}

	b.insert_particles(sons);
}

void adaptivity::perform_split_cubic_basic(body &b) const {
	simulation_time *time = &simulation_time::getInstance();
	unsigned int step = time->get_step();

	std::vector<particle> &particles = b.get_particles();

	// how many "SON" do you have in 2D? ---> 3
	const unsigned int num_SON2D = m_num_child-1;

	std::vector<particle> sons;

	// 2. give birth to SON + assign values to SON
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int son_idx = particles.size();
	for (unsigned int i = 0; i < particles.size(); i++) {

		// only append new "SON" particles if their "DAD" particle:
		// 		A) has been marked as split/merge candidate
		// 		B) has NOT already reached the MAX_REFINE_STEP
		//---------------------------------------------------------
		unsigned int delta_st = step - particles[i].last_refine_at;
		assert(delta_st>=0);

		if(particles[i].split && particles[i].refine_step<MAX_REFINE_STEP) {
			double x_SON [num_SON2D];
			double y_SON [num_SON2D];
			double h_SON [num_SON2D];
			double m_SON [num_SON2D];

			// 0. call your DAD
			double dx = sqrt(particles[i].m/particles[i].rho);
			double x_DAD = particles[i].X; // or x if ~regular arrangement
			double y_DAD = particles[i].Y; // or y if ~regular arrangement
			double h_DAD = particles[i].h;
			double m_DAD = particles[i].m;

			// 1. increase the refinement step of "DAD"
			particles[i].last_refine_at = step;
			particles[i].refine_step ++;


			// 2. give birth to new "SON" - following the given refinement pattern
			/*
			 *    DAD -------- SON[0]
			 *     |		     |
			 * 	   |			 |
			 * 	   |			 |
			 *   SON[2] ------ SON[1]
			 */

			if(m_eccentric_cubic) {
				// SON[0]
				x_SON[0] = x_DAD + m_alpha*dx;
				y_SON[0] = y_DAD;
				h_SON[0] = m_beta*h_DAD;
				m_SON[0] = (plus_one) ? (1./(num_SON2D+1))*m_DAD : (1./num_SON2D)*m_DAD;

				// SON[1]
				x_SON[1] = x_DAD + m_alpha*dx;
				y_SON[1] = y_DAD - m_alpha*dx;
				h_SON[1] = m_beta*h_DAD;
				m_SON[1] = (plus_one) ? (1./(num_SON2D+1))*m_DAD : (1./num_SON2D)*m_DAD;

				// SON[2]
				x_SON[2] = x_DAD;
				y_SON[2] = y_DAD - m_alpha*dx;
				h_SON[2] = m_beta*h_DAD;
				m_SON[2] = (plus_one) ? (1./(num_SON2D+1))*m_DAD : (1./num_SON2D)*m_DAD;

			} else {

				// SON[0]
				x_SON[0] = x_DAD - cos(M_PI/4.)*m_alpha*dx;
				y_SON[0] = y_DAD - sin(M_PI/4.)*m_alpha*dx;
				h_SON[0] = m_beta*h_DAD;
				m_SON[0] = (plus_one) ? (1./(num_SON2D+1))*m_DAD : (1./num_SON2D)*m_DAD;

				// SON[1]
				x_SON[1] = x_DAD - cos(M_PI/4.)*m_alpha*dx;
				y_SON[1] = y_DAD + sin(M_PI/4.)*m_alpha*dx;
				h_SON[1] = m_beta*h_DAD;
				m_SON[1] = (plus_one) ? (1./(num_SON2D+1))*m_DAD : (1./num_SON2D)*m_DAD;

				// SON[2]
				x_SON[2] = x_DAD + cos(M_PI/4.)*m_alpha*dx;
				y_SON[2] = y_DAD + sin(M_PI/4.)*m_alpha*dx;
				h_SON[2] = m_beta*h_DAD;
				m_SON[2] = (plus_one) ? (1./(num_SON2D+1))*m_DAD : (1./num_SON2D)*m_DAD;
			}

			// 3. append new "SON" to the current particle array
			for (unsigned int ii = 0; ii < num_SON2D; ii++) {
				unsigned int id_DAD = i;

				particle son_particle(son_idx);
				son_idx++;

				son_particle.x = x_SON[ii];
				son_particle.y = y_SON[ii];
				son_particle.X = son_particle.x;
				son_particle.Y = son_particle.y;
				son_particle.m = m_SON[ii];
				son_particle.h = h_SON[ii];
				son_particle.refine_step = 1;
				son_particle.last_refine_at = step;
				son_particle.split = false;
				son_particle.merge = false;

				//copy_dad_to_son(particles[id_DAD], particles[id_SON]);
				copy_dad_to_son(particles[id_DAD], son_particle);
				sons.push_back(son_particle);

			}

			// 4. slight modification: DAD himself becomes a SON from now on! "+1" approach
			particles[i].x = (m_eccentric_cubic) ? x_DAD : x_DAD + cos(M_PI/4.)*m_alpha*dx;
			particles[i].y = (m_eccentric_cubic) ? y_DAD : y_DAD - sin(M_PI/4.)*m_alpha*dx;
			particles[i].m = (plus_one) ? m_SON[0] : 1e-16;
			particles[i].h = m_beta*h_DAD;
			particles[i].split = false;

			// sanity check - total mass conservation after refinement
			double sum_mass = particles[i].m + m_SON[0] + m_SON[1] + m_SON[2];
			assert(fabs(m_DAD-sum_mass) < 1e-12);
		}
	}

	b.insert_particles(sons);
}

void adaptivity::perform_split_cubic(body &b) const {
	simulation_time *time = &simulation_time::getInstance();
	unsigned int step = time->get_step();

	std::vector<particle> &particles = b.get_particles();

	// how many "SON" do you have in 2D? ---> 4
	const unsigned int num_SON2D = m_num_child-1;

	// 1. consider a vector of SON particles
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	std::vector<particle> sons;

	// 2. give birth to SON + assign values to SON
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int son_idx = particles.size();
	for (unsigned int i = 0; i < particles.size(); i++) {

		// only append new "SON" particles if their "DAD" particle:
		// 		A) has been marked as split/merge candidate
		// 		B) has NOT already reached the MAX_REFINE_STEP
		//---------------------------------------------------------
		unsigned int delta_st = step - particles[i].last_refine_at;
		assert(delta_st>=0);

		if(particles[i].split && particles[i].refine_step < MAX_REFINE_STEP) {
			double x_SON [num_SON2D];
			double y_SON [num_SON2D];
			double h_SON [num_SON2D];
			double m_SON [num_SON2D];

			// 0. call your DAD
			double dx = sqrt(particles[i].m/particles[i].rho);
			double x_DAD = particles[i].x;
			double y_DAD = particles[i].y;
			double h_DAD = particles[i].h;
			double m_DAD = particles[i].m;

			// 1. increase the refinement step of "DAD"
			particles[i].last_refine_at = step;
			particles[i].refine_step ++;


			// 2. give birth to new "SON" - following the given refinement pattern
			/*
			 *    SON[1] --------------- SON[2]
			 *     |						|
			 * 	   |		  DAD			|
			 * 	   |						|
			 *    SON[0] --------------- SON[3]
			 */

			// SON[0]
			x_SON[0] = x_DAD - cos(M_PI/4.)*m_alpha*dx;
			y_SON[0] = y_DAD - sin(M_PI/4.)*m_alpha*dx;
			h_SON[0] = m_beta*h_DAD;
			m_SON[0] = (plus_one) ? (1./(num_SON2D+1))*m_DAD : (1./num_SON2D)*m_DAD;

			// SON[1]
			x_SON[1] = x_DAD - cos(M_PI/4.)*m_alpha*dx;
			y_SON[1] = y_DAD + sin(M_PI/4.)*m_alpha*dx;
			h_SON[1] = m_beta*h_DAD;
			m_SON[1] = (plus_one) ? (1./(num_SON2D+1))*m_DAD : (1./num_SON2D)*m_DAD;

			// SON[2]
			x_SON[2] = x_DAD + cos(M_PI/4.)*m_alpha*dx;
			y_SON[2] = y_DAD + sin(M_PI/4.)*m_alpha*dx;
			h_SON[2] = m_beta*h_DAD;
			m_SON[2] = (plus_one) ? (1./(num_SON2D+1))*m_DAD : (1./num_SON2D)*m_DAD;

			// SON[3]
			x_SON[3] = x_DAD + cos(M_PI/4.)*m_alpha*dx;
			y_SON[3] = y_DAD - sin(M_PI/4.)*m_alpha*dx;
			h_SON[3] = m_beta*h_DAD;
			m_SON[3] = (plus_one) ? (1./(num_SON2D+1))*m_DAD : (1./num_SON2D)*m_DAD;

			// 3. append new "SON" to the current particle array
			for (unsigned int ii = 0; ii < num_SON2D; ii++) {
				unsigned int id_DAD = i;

				particle son_particle(son_idx);

				son_particle.x = x_SON[ii];
				son_particle.y = y_SON[ii];
				son_particle.X = son_particle.x;
				son_particle.Y = son_particle.y;
				son_particle.m = m_SON[ii];
				son_particle.h = h_SON[ii];
				son_particle.refine_step = 1;
				son_particle.last_refine_at = step;
				son_particle.split = false;
				son_particle.merge = false;

				//copy_dad_to_son(particles[id_DAD], particles[id_SON]);
				copy_dad_to_son(particles[id_DAD], son_particle);
				sons.push_back(son_particle);

			}

			// 4. slight modification: DAD himself becomes a SON from now on! "+1" approach
			particles[i].m = (plus_one) ? m_SON[0] : 1e-16;
			particles[i].h = m_beta*h_DAD;

			// sanity check - total mass conservation after refinement
			double sum_mass = particles[i].m + m_SON[0] + m_SON[1] + m_SON[2] + m_SON[3];
			assert(fabs(m_DAD-sum_mass) < 1e-12);

			// 5. update locally for next SONs
		}
	}

	b.insert_particles(sons);
}

void adaptivity::perform_split_hexagonal(body &b)  const {
	simulation_time *time = &simulation_time::getInstance();
	unsigned int step = time->get_step();

	std::vector<particle> &particles = b.get_particles();

	// how many "SON" do you have in 2D? ---> 6
	const unsigned int num_SON2D = m_num_child-1;

	// given by J. Feldman & J. Bonet after minimizing the global density error [2006] with (alpha=60%,beta=60%)
	double coeff_md = (1./(num_SON2D+1)); //0.102181;
	double coeff_m0 = (1./(num_SON2D+1)); //0.386914;

	// 1. consider a vector of SON particles
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	std::vector<particle> sons;

	// 2. give birth to SON + assign values to SON
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int son_idx = particles.size();
	for (unsigned int i = 0; i < particles.size(); i++) {
		// only append new "SON" particles if their "DAD" particle:
		// 		A) has been marked as split/merge candidate
		// 		B) has NOT already reached the MAX_REFINE_STEP
		//---------------------------------------------------------
		unsigned int delta_st = step - particles[i].last_refine_at;
		assert(delta_st>=0);
		//if(particles[i].split && particles[i].refine_step<MAX_REFINE_STEP && delta_st>MIN_REFINE_DIFF) {
		if(particles[i].split && particles[i].refine_step<MAX_REFINE_STEP) {

			double x_SON [num_SON2D];
			double y_SON [num_SON2D];
			double h_SON [num_SON2D];
			double m_SON [num_SON2D];

			// 0. call your DAD
			double dx = sqrt(particles[i].m/particles[i].rho);
			double x_DAD = particles[i].x;
			double y_DAD = particles[i].y;
			double h_DAD = particles[i].h;
			double m_DAD = particles[i].m;

			// 1. increase the refinement step of "DAD"
			particles[i].last_refine_at = step;
			particles[i].refine_step ++;


			// 2. give birth to new "SON" - following the given refinement pattern
			/*
			 *        SON[2]       SON[3]
			 *              \     /
			 *               \   /
			 * 	  SON[1]------DAD------SON[4]
			 * 	             /	 \
			 *              /     \
			 *        SON[0]       SON[5]
			 */

			// SON[0]
			x_SON[0] = x_DAD - cos(M_PI/3.)*m_alpha*dx;
			y_SON[0] = y_DAD - sin(M_PI/3.)*m_alpha*dx;
			h_SON[0] = m_beta*h_DAD;
			m_SON[0] = (plus_one) ? coeff_md*m_DAD : (1./num_SON2D)*m_DAD;

			// SON[1]
			x_SON[1] = x_DAD - m_alpha*dx;
			y_SON[1] = y_DAD;
			h_SON[1] = m_beta*h_DAD;
			m_SON[1] = (plus_one) ? coeff_md*m_DAD : (1./num_SON2D)*m_DAD;

			// SON[2]
			x_SON[2] = x_DAD - cos(M_PI/3.)*m_alpha*dx;
			y_SON[2] = y_DAD + sin(M_PI/3.)*m_alpha*dx;
			h_SON[2] = m_beta*h_DAD;
			m_SON[2] = (plus_one) ? coeff_md*m_DAD : (1./num_SON2D)*m_DAD;

			// SON[3]
			x_SON[3] = x_DAD + cos(M_PI/3.)*m_alpha*dx;
			y_SON[3] = y_DAD + sin(M_PI/3.)*m_alpha*dx;
			h_SON[3] = m_beta*h_DAD;
			m_SON[3] = (plus_one) ? coeff_md*m_DAD : (1./num_SON2D)*m_DAD;

			// SON[4]
			x_SON[4] = x_DAD + m_alpha*dx;
			y_SON[4] = y_DAD;
			h_SON[4] = m_beta*h_DAD;
			m_SON[4] = (plus_one) ? coeff_md*m_DAD : (1./num_SON2D)*m_DAD;

			// SON[5]
			x_SON[5] = x_DAD + cos(M_PI/3.)*m_alpha*dx;
			y_SON[5] = y_DAD - sin(M_PI/3.)*m_alpha*dx;
			h_SON[5] = m_beta*h_DAD;
			m_SON[5] = (plus_one) ? coeff_md*m_DAD : (1./num_SON2D)*m_DAD;


			// 3. append new "SON" to the current particle array
			for (unsigned int ii = 0; ii < num_SON2D; ii++) {
				unsigned int id_DAD = i;

				particle son_particle(son_idx);

				son_particle.x = x_SON[ii];
				son_particle.y = y_SON[ii];
				son_particle.X = son_particle.x;
				son_particle.Y = son_particle.y;
				son_particle.m = m_SON[ii];
				son_particle.h = h_SON[ii];
				son_particle.refine_step = 1;
				son_particle.last_refine_at = step;
				son_particle.split = false;
				son_particle.merge = false;

				//copy_dad_to_son(particles[id_DAD], particles[id_SON]);
				copy_dad_to_son(particles[id_DAD], son_particle);
				sons.push_back(son_particle);

				//my_extrapol_SPH_values(b.get_particles(),id_SON,id_DAD);
			}

			// 4. slight modification: DAD himself becomes a SON from now on! "+1" approach
			particles[i].m = (plus_one) ? coeff_m0*m_DAD : 1e-16;
			particles[i].h = m_beta*h_DAD;

			// sanity check - total mass conservation after refinement
			double sum_mass = particles[i].m + m_SON[0] + m_SON[1] + m_SON[2] + m_SON[3] + m_SON[4] + m_SON[5];
			assert(fabs(m_DAD - sum_mass) < 1e-10);
		}
	}

	b.insert_particles(sons);
}

int adaptivity::evaluate_refinement(body &b) const {
	switch (m_refine_criteria) {
	case velocity:
		return scan_mark_velocity_based(b);
		break;
	case div_velo:
		return scan_mark_div_velocity_based(b);
		break;
	case vonmises:
		return scan_mark_vonmises_based(b);
		break;
	case plastic_strain:
		return scan_mark_strain_based(b);
		break;
	case temperature:
		return scan_mark_temperature_based(b);
		break;
	case position:
		return scan_mark_position_based(b);
		break;
	case moving_frame:
		return scan_mark_moving_frame(b);
		break;
	case num_nbh:
		return scan_mark_neighbor_based(b);
		break;
	}

	//not reachable
	return -1;
}

void adaptivity::do_split(body &b) const {
	switch (m_pattern) {
	case triangular:
		perform_split_triangular(b);
		break;
	case cubic_basic:
		perform_split_cubic_basic(b);
		break;
	case cubic:
		perform_split_cubic(b);
		break;
	case hexagonal:
		perform_split_hexagonal(b);
		break;
	}
}

void adaptivity::adapt_resolution(body &b) const {
	unsigned int n_before = b.get_num_part();

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	flag_reset(b);
	int num_dad = evaluate_refinement(b);

	if (num_dad == 0) return;

	do_split(b);
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	unsigned int n_after = b.get_num_part();

	if(n_after-n_before > (n_before/2))
		printf("WARNING!!! --- too many sons were born, beware the proliferation. \n");
}

adaptivity::adaptivity(double alpha_dx, double beta_h, double v, double div_v, double SvM, double epsPl,
			   double T, glm::dvec2 xy_min, glm::dvec2 xy_max,
			   double frm_width, double frm_height, unsigned int num_nbh, double l_eff, bool eccentric) {

	m_alpha = alpha_dx;
	m_beta = beta_h;
	m_v_threshold = v;
	m_div_v_threshold = div_v;
	m_SvM_threshold = SvM;
	m_eps_threshold = epsPl;
	m_T_threshold = T;
	m_xy_min = xy_min;
	m_xy_max = xy_max;
	m_width  = frm_width;
	m_height = frm_height;
	m_num_nbh_threshold = num_nbh;
	m_l_eff = l_eff;
	m_eccentric_cubic = eccentric;
}
