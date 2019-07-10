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

#include "body.h"

void body::apply_plasticity() {
	if (m_plast == 0) return;
	m_plast->plastic_state_by_radial_return(*this);
}

void body::apply_thermal_conduction() {
	if (m_thermal == 0) return;
	m_thermal->conduction(*this);
}

void body::apply_contact() {
	if (m_tool == 0) return;
	contact_apply_tool_to_body_2d(m_tool, *this);
}

void body::apply_adaptivity() {
	if (m_adapt == 0) return;
	m_adapt->adapt_resolution(*this);
}

void body::set_tool(tool *tool) {
	m_tool = tool;
}

void body::move_tool() {
	if (m_tool == 0) return;
	m_tool->update_tool();
}

glm::dvec2 body::speed_tool() {
	return m_tool->get_vel();
}

glm::dvec2 body::edge_tool() {
	return m_tool->get_edge_coord();
}

void body::set_plasticity(plasticity *plasticity) {
	m_plast = plasticity;
}

void body::set_thermal(thermal *thermal) {
	m_thermal = thermal;
}

void body::set_adaptivity(adaptivity *adaptivity) {
	m_adapt = adaptivity;
}

void body::construct_verlet_lists() {
	const unsigned int num_part = m_particles.size();

	m_grid.update_geometry(m_particles, num_part, 2.);
	m_grid.assign_hashes(m_particles, num_part);

	std::sort(m_particles.begin(), m_particles.end(),
			[](const particle &a, const particle &b) {return a.hash < b.hash;});

	m_grid.construct_verlet_lists(m_particles, num_part, 2.);

	m_basis_fun(m_particles, num_part);
}

void body::insert_particles(const std::vector<particle>& additional_particles) {
	m_particles.insert(m_particles.end(), additional_particles.begin(), additional_particles.end());
}

void body::restore_order() {
	std::sort(m_particles.begin(), m_particles.end(),
			[](const particle &a, const particle &b) {return a.idx < b.idx;});
}

void body::set_basis_fun(void (*basis_fun)(std::vector<particle> &particles , unsigned int)) {
	m_basis_fun = basis_fun;
}

simulation_data body::get_sim_data()  const {
	return m_simulation_data;
}

std::vector<particle> &body::get_particles() {
	return m_particles;
}

const std::vector<particle> &body::get_particles() const {
	return m_particles;
}

unsigned int body::get_num_part() const {
	return m_particles.size();
}

body::body(particle* particles, unsigned int n, simulation_data data) :
		m_simulation_data(data) {

	m_particles.resize(n);
	std::copy(particles, particles+n, m_particles.begin());
}

body::body() {}
