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

#ifndef THERMAL_H_
#define THERMAL_H_

#include <math.h>
#include <assert.h>
#include <glm/glm.hpp>

#include "grid.h"
#include "kernel.h"
#include "particle.h"
#include "simulation_data.h"

/*
 This implements heat conduction using either of the desired methods:
 	 1- the particle strength exchange (PSE) method.
 	 2- the Brookshaw-SPH method.

 	 > both schemes discretize the heat equation in a Finite-Difference like approach.
 	 > both schemes are energy conservative. (anti-symmetric form)
 	 > both schemes are numerically efficient and capable of handling adiabatic boundary condition without dummy particles.

 	 For further details, please refer to the following publications:

 	 1- "A general deterministic treatment of derivatives in particle methods."
 	    	By: J. Eldredge et al.
 	    	Journal of Computational Physics 180.2 (2002): 686-709.

 	 2- "A method of calculating radiative heat diffusion in particle simulations”
 	 	 	 By: L. Brookshaw
 	 	     Proceedings of the Astronomical Society of Australia, vol. 6, pp. 207–210, 1985"
*/


class body;

class thermal {
public:
	enum thermal_solver {
		thermal_pse,
		thermal_brookshaw,
	};

	void set_method(thermal_solver solver);
	void conduction(body &body) const;
	thermal(physical_constants pc);

private:
	double m_alpha = 0.;
	thermal_solver m_thermal_solver = thermal_pse;

	void heat_conduction_pse(body &b) const;
	void heat_conduction_brookshaw(body &b) const;
};

#endif /* THERMAL_H_ */
