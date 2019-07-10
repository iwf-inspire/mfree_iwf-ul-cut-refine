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

#ifndef CONTACT_H_
#define CONTACT_H_

#include <glm/glm.hpp>

#include "tool.h"
#include "body.h"

/*
 This is the implementation of the penalty contact algorithm illustrated in Fig. 12 in Section 6.1.

 This function is called as soon as  is detected using the parametrization of the tool.

 - The [TOOL] is rigid and parametrized, as can be found in "tool.h".
 - [TOOL <--> WORKPIECE] contact penetration work as follows.
 	 1- The contact algorithm is established for a query particle "p" as soon as it is found inside the bounding box of the tool.
 	 2- The min distance of "p" from the tool segments is calculated and named as "PENETRATION DEPTH"
 	 3- A penalty contact force is computed proportional to said "PENETRATION DEPTH"
 	 4- Since this location in not admissible, particle "p" is pushed out by exerting this penalty contact force onto it.
 	 5- The resultant Fc (cutting) and Ft (thrust) forces acting on the tool is computed by summing the x and y components of contact and friction forces.

*/

void contact_apply_tool_to_body_2d(const tool *master, body &slave);

#endif /* CONTACT_H_ */
