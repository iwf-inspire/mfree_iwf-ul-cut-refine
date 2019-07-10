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

#ifndef BODY_H_
#define BODY_H_

#include <type_traits>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "particle.h"
#include "simulation_data.h"
#include "plasticity.h"
#include "precomp_shape_functions.h"
#include "thermal.h"
#include "contact.h"

#include "adaptivity.h"

class body {

private:
	/*
	 "BODY" is comprised of the following encapsulated items:
	*/

	grid m_grid;             			// spatial hashing
	plasticity *m_plast = 0; 			// plasticity algorithm (if any)
	thermal *m_thermal = 0;  			// thermal algorithm (if any)
	adaptivity *m_adapt = 0;  			// adaptivity algorithm (if any)
	tool *m_tool = 0;					// tool (potentially) in contact with this body (if any)
	std::vector<particle> m_particles;  // workpiece particles
	simulation_data m_simulation_data;  // all physical constants
	void (*m_basis_fun)(std::vector<particle> &particles, unsigned int) = &precomp_sph; // basis function chosen SPH

public:
	void set_plasticity(plasticity *plasticity);
	void set_thermal(thermal *thermal);
	void set_tool(tool *tool);
	void set_adaptivity(adaptivity *adaptivity);

	void apply_plasticity();
	void apply_thermal_conduction();
	void apply_contact();
	void move_tool();
	void apply_adaptivity();

	glm::dvec2 speed_tool();
	glm::dvec2 edge_tool();

	void construct_verlet_lists();
	void restore_order();

	void set_basis_fun(void (*m_basis_fun)(std::vector<particle> &particles, unsigned int));

	simulation_data get_sim_data()  const;
	std::vector<particle> &get_particles();
	const std::vector<particle> &get_particles() const;
	unsigned int get_num_part() const;

	void insert_particles(const std::vector<particle>& additional_particles);

	body(particle* particles, unsigned int n, simulation_data data);

	// do not allow copying a body
	body(const body &copy) = delete;
	body& operator= (const body &fraction) = delete;

	body();
};

#endif /* BODY_H_ */
