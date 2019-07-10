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

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "kernel.h"
#include <string.h>
#include <assert.h>
#include <stdio.h>

#define MAX_NBH 600

// "particle" data structure

class particle {

public:
	particle();
	particle(unsigned int idx);
	particle(const particle &other);
	particle &operator=(const particle &other);

	virtual ~particle();

	// reset all time derivatives
	virtual void reset();
	virtual void copy_into(particle &p) const;

	// index + hash for spatial hashing
	unsigned int idx = 0, hash = 0;

	double x = 0., y = 0.;		  /* Current particle location */
	double X = 0.; double Y = 0.; /* Initial particle location */
	double m = 0., rho = 0.;      /* Initial density */
	double rho_init = 0.;         /* Mass and density */
	double quad_weight = 0.;      /* Quadrature weight (equal to m/rho in mfree algos) */
	double h = 0.;				  /* Smoothing length */
	double p = 0.;			      /* Hydrostatic pressure */

	double vx = 0., vy = 0.;							   /* Particle velocity */
	double Sxx = 0., Sxy = 0., Syy = 0., Szz = 0.;		   /* Deviatoric(!) stress tensor components */
	double Sxx_x = 0., Sxy_y = 0., Sxy_x = 0., Syy_y = 0.; /* Stress gradient */
	double vx_x = 0., vx_y = 0., vy_x = 0., vy_y = 0.;     /* Velocity gradient */

	// time derivatives
	double x_t = 0., y_t = 0.;
	double rho_t = 0.;
	double h_t = 0.;
	double vx_t = 0., vy_t = 0.;
	double Sxx_t = 0., Sxy_t = 0., Syy_t = 0., Szz_t = 0.;

	// components of contact (_c) & tangential (_t) forces
	double fcx = 0.; double fcy = 0.;
	double ftx = 0.; double fty = 0.;

	// artificial stress components
	double Rxx = 0., Rxy = 0., Ryy = 0.;

	double T		= 0.;		/* Current particle temperature */
	double T_t		= 0.;		/* Particle temperature derivative wrt time */
	double T_init	= 0.;		/* Initial particle temperature */

	// plastic strain tensor
	double eps_plxx = 0.;
	double eps_plxy = 0.;
	double eps_plyy = 0.;
	double eps_plzz = 0.;

	double eps_pl_equiv = 0.;
	double eps_pl_equiv_dot = 0.;

	// adaptivity
	unsigned int refine_step = 0;
	unsigned int last_refine_at = 0;
	bool split = false;
	bool merge = false;

	// Verlet list stuff
	unsigned int num_nbh = 0;
	unsigned int nbh[MAX_NBH];

	// kernel functions computed
	kernel_result w[MAX_NBH];

	// flag for fixed BC particles
	bool fixed = false;
};

#endif /* PARTICLE_H_ */
