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

#include "particle.h"

particle::particle() {
//	memset(nbh, 0, sizeof(unsigned int)*MAX_NBH);
//	memset(w,   0, sizeof(kernel_result)*MAX_NBH);
};

particle::particle(unsigned int idx) {
	this->idx = idx;

//	memset(nbh, 0, sizeof(unsigned int)*MAX_NBH);
//	memset(w,   0, sizeof(kernel_result)*MAX_NBH);
}

particle::~particle() {};

void particle::reset() {
	x_t  = 0.;
	y_t  = 0.;
	rho_t = 0.;
	h_t = 0.;
	vx_t = 0.;
	vy_t = 0.;
	Sxx_t = 0.;
	Sxy_t = 0.;
	Syy_t = 0.;
	Szz_t = 0.;
	T_t = 0.;
}

particle::particle(const particle &other) {
	idx = other.idx;
	hash = other.hash;

	m = other.m;
	rho = other.rho;
	rho_init = other.rho_init;
	h = other.h;

	x = other.x;
	y = other.y;

	X = other.X;
	Y = other.Y;

	vx = other.vx;
	vy = other.vy;

	Sxx = other.Sxx;
	Sxy = other.Sxy;
	Syy = other.Syy;
	Szz = other.Szz;

	eps_pl_equiv     = other.eps_pl_equiv;
	eps_pl_equiv_dot = other.eps_pl_equiv_dot;

	T = other.T;

	last_refine_at = other.last_refine_at;
	refine_step = other.refine_step;
	split = other.split;
	merge = other.merge;

	fixed = other.fixed;
}

particle &particle::operator=(const particle &other) {
	if (this == &other)
		return *this;

	idx = other.idx;
	hash = other.hash;

	m = other.m;
	rho = other.rho;
	rho_init = other.rho_init;
	h = other.h;

	x = other.x;
	y = other.y;

	X = other.X;
	Y = other.Y;

	vx = other.vx;
	vy = other.vy;

	Sxx = other.Sxx;
	Sxy = other.Sxy;
	Syy = other.Syy;
	Szz = other.Szz;

	eps_pl_equiv     = other.eps_pl_equiv;
	eps_pl_equiv_dot = other.eps_pl_equiv_dot;

	T = other.T;

	last_refine_at = other.last_refine_at;
	refine_step = other.refine_step;
	split = other.split;
	merge = other.merge;

	fixed = other.fixed;

	return *this;

}

void particle::copy_into(particle &p) const {
	p.idx = idx;
	p.hash = hash;

	p.m = m;
	p.rho = rho;
	p.rho_init = rho_init;
	p.h = h;

	p.x = x;
	p.y = y;

	p.X = X;
	p.Y = Y;

	p.vx = vx;
	p.vy = vy;

	p.Sxx = Sxx;
	p.Sxy = Sxy;
	p.Syy = Syy;
	p.Szz = Szz;

	p.eps_pl_equiv = eps_pl_equiv;
	p.eps_pl_equiv_dot = eps_pl_equiv_dot;

	p.eps_plxx = eps_plxx;
	p.eps_plxy = eps_plxy;
	p.eps_plyy = eps_plyy;
	p.eps_plzz = eps_plzz;

	p.T = T;
	p.T_init = T_init;

	p.fcx = fcx;
	p.fcy = fcy;

	p.fixed = fixed;

	p.last_refine_at = last_refine_at;
	p.refine_step = refine_step;
	p.split = split;
	p.merge = merge;

	p.num_nbh = num_nbh;

	memcpy(p.nbh, nbh, sizeof(unsigned int)*num_nbh);
	memcpy(p.w,   w,   sizeof(kernel_result)*num_nbh);
}
