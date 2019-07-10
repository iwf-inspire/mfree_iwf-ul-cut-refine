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

#include "leap_frog.h"

void leap_frog::init(body &body) {
	std::vector<particle> &particles  = body.get_particles();

	// this can only grow the vector since now coarsening strategy is not applied.
	// resize should not cause re-allocation since the constructor reserves twice the initial particle number,
	// so while this growth strategy is not perfect it should not or rarely have a runtime impact
	if (m_init.size() < particles.size()) {
		m_init.resize(particles.size());
	}

	for (unsigned int i = 0; i < body.get_num_part(); i++) {
		m_init[i] = particles[i];
	}
}

void leap_frog::predict(body &body) const {
	std::vector<particle> &particles  = body.get_particles();
	simulation_time *time = &simulation_time::getInstance();
	double dt = time->get_dt();

	for (unsigned int i = 0; i < body.get_num_part(); i++) {
		particles[i].x   = m_init[i].x   + 0.5*dt*particles[i].x_t;
		particles[i].y   = m_init[i].y   + 0.5*dt*particles[i].y_t;
		particles[i].rho = m_init[i].rho + 0.5*dt*particles[i].rho_t;
		particles[i].h   = m_init[i].h   + 0.5*dt*particles[i].h_t;
		particles[i].vx  = m_init[i].vx  + 0.5*dt*particles[i].vx_t;
		particles[i].vy  = m_init[i].vy  + 0.5*dt*particles[i].vy_t;
		particles[i].Sxx = m_init[i].Sxx + 0.5*dt*particles[i].Sxx_t;
		particles[i].Sxy = m_init[i].Sxy + 0.5*dt*particles[i].Sxy_t;
		particles[i].Syy = m_init[i].Syy + 0.5*dt*particles[i].Syy_t;
		particles[i].Szz = m_init[i].Szz + 0.5*dt*particles[i].Szz_t;
		particles[i].T   = m_init[i].T   + 0.5*dt*particles[i].T_t;
	}
}

void leap_frog::correct(body &body) const {
	std::vector<particle> &particles = body.get_particles();
	simulation_time *time = &simulation_time::getInstance();
	double dt = time->get_dt();

	for (unsigned int i = 0; i < body.get_num_part(); i++) {
		particles[i].x   = m_init[i].x   + dt*particles[i].x_t;
		particles[i].y   = m_init[i].y   + dt*particles[i].y_t;
		particles[i].rho = m_init[i].rho + dt*particles[i].rho_t;
		particles[i].h   = m_init[i].h   + dt*particles[i].h_t;
		particles[i].vx  = m_init[i].vx  + dt*particles[i].vx_t;
		particles[i].vy  = m_init[i].vy  + dt*particles[i].vy_t;
		particles[i].Sxx = m_init[i].Sxx + dt*particles[i].Sxx_t;
		particles[i].Sxy = m_init[i].Sxy + dt*particles[i].Sxy_t;
		particles[i].Syy = m_init[i].Syy + dt*particles[i].Syy_t;
		particles[i].Szz = m_init[i].Szz + dt*particles[i].Szz_t;
		particles[i].T   = m_init[i].T   + dt*particles[i].T_t;
	}
}

void leap_frog::step(body &body) {

	// Update the neighbors by spatial hashing
	body.construct_verlet_lists();

	// Leapfrog predictor step
	init(body);
	predict(body);

	// compute temporal derivatives
	for (unsigned int i = 0; i < body.get_num_part(); i++) {
		body.get_particles()[i].reset();
	}

	// move the tool & do penalty contact
	body.apply_contact();
	body.move_tool();

	// Compute time derivatives of physical properties
	// =================================================
	material_eos(body);
	correctors_mghn_artificial_stress(body);
	derive_stress_monaghan(body);
	derive_velocity(body);
	correctors_mghn_artificial_viscosity(body);
	correctors_xsph(body);
	material_stress_rate_jaumann(body);
	contmech_continuity(body);
	contmech_momentum(body);
	contmech_advection(body);
	// =================================================

	// Solve heat equation
	body.apply_thermal_conduction();

	// Leapfrog corrector step
	correct(body);

	// Perform plasticity by radial return
	body.apply_plasticity();

	// boundary conditions
	do_boundary_conditions(body);

	// restore particles into their original order
	// this step is not necessary for correctness but may be useful for debugging purposes
	// deactivated for performance reasons
    body.restore_order();

    body.apply_adaptivity();
}

leap_frog::leap_frog(unsigned int num_part)  {
	m_init.reserve(2*num_part);
}
