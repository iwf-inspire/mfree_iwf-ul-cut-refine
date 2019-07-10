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

#ifndef PLASTICITY_H_
#define PLASTICITY_H_

#include <assert.h>
#include <math.h>
#include <vector>

#include "solver.h"
#include "simulation_time.h"
#include "johnson_cook_Sima_2010.h"
#include "particle.h"

/*
 	 The purpose of this implementation:
 	   - specification of plastic state by radial return algorithm.
 	   - "SECANT" scheme is used for root finding procedure of radial return
 	   - This choice is made motivated by the performance of "SECANT" (faster & safer)
*/

class body;

class plasticity {

private:
	double m_tol = 1e-6;
	johnson_cook_Sima_2010 *m_plasticity_model = 0;
	bool m_consider_dissipation = true;
	void print_debug(const std::vector<particle> &particles, unsigned int num_part, unsigned int fail_idx);
	void do_radial_return(std::vector<particle> &particles, unsigned int num_part, simulation_data data);

public:
	void plastic_state_by_radial_return(body &b);
	void set_tolerance(double tol);
	void set_dissipation_considered(bool consider);
	plasticity(johnson_cook_Sima_2010 *plasticity_model);
	plasticity();
};

#endif /* PLASTICITY_H_ */
