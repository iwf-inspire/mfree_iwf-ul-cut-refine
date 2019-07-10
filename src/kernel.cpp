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

#include "kernel.h"

#include <stdio.h>

kernel_result cubic_spline(double xi, double yi, double xj, double yj, double h) {
	const double h1 = 1./h;
	const double xij = xi-xj;
	const double yij = yi-yj;
	const double rij = sqrt(xij*xij + yij*yij);
	const double q = rij*h1;

	kernel_result w;

	if (q >= 2.) {
		return w;
	} else if (q >= 1.) {
		const double twominq = (2-q);
		const double fac = 10*(M_1_PI)/7.0*h1*h1;
		w.w =  fac*(0.25*twominq*twominq*twominq);

		const double der = -0.75*twominq*twominq * h1/rij;
		w.w_x = xij*der*fac;
		w.w_y = yij*der*fac;

		return w;
	} else {
		const double fac = 10*(M_1_PI)/7.0*h1*h1;
		w.w = fac*(1 - 1.5*q*q*(1-0.5*q));
		if (rij > 1e-12) {
			const double der = -3.0*q*(1-0.75*q) * h1/rij;
			w.w_x = xij*der*fac;
			w.w_y = yij*der*fac;
		}
		return w;
	}
}
