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

#ifndef ADAPTIVITY_H_
#define ADAPTIVITY_H_

#include <math.h>
#include <assert.h>
#include <algorithm>
#include <glm/glm.hpp>
#include <stdio.h>
#include <type_traits>

#include "particle.h"

#define MAX_REFINE_STEP 1
#define MIN_REFINE_DIFF 10

class body;

class adaptivity {
public:
	enum refine_criteria {
		velocity,
		div_velo,
		vonmises,
		plastic_strain,
		temperature,
		position,
		moving_frame,
		num_nbh
	};

	enum pattern {
		triangular,
		cubic_basic,
		cubic,
		hexagonal
	};

	void dens_before_approx_N2(body &b) const;
	void dens_after_approx_N2(body &b) const;

	void set_refine_criterion(refine_criteria crit);
	void set_refine_pattern(pattern patt);
	void adapt_resolution(body &b) const;

	adaptivity(double alpha_dx, double beta_h, double v, double div_v, double SvM, double epsPl,
			   double T, glm::dvec2 xy_min, glm::dvec2 xy_max, double frm_width, double frm_height,
			   unsigned int num_nbh, double l_eff, bool eccentric);

private:
	bool plus_one = true;
	unsigned int m_num_child = 0;

	// settings +-+-++-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-
	double m_alpha = 0.;                 // equal to \kappa in paper: dispersion coefficient
	double m_beta = 0.;                  // equal to \xi in paper: smoothing ratio
	double m_v_threshold = 0.;           // single threshold magnitude for velocity
	double m_div_v_threshold = 0.;       // single threshold magnitude for divergence of velocity
	double m_SvM_threshold = 0.;         // single threshold magnitude for von Mises stress
	double m_eps_threshold = 0.;         // single threshold magnitude for equivalent plastic strain
	double m_T_threshold = 0.;           // single threshold magnitude for temperature
	glm::dvec2 m_xy_min = {0., 0.};      // x_min & y_min of the refinement zone
	glm::dvec2 m_xy_max = {0., 0.};      // x_max & y_max of the refinement zone
	double m_width  = 0.;                // width of the moving refinement frame
	double m_height = 0.;                // width of the moving refinement frame
	unsigned int m_num_nbh_threshold = 0;// min number of neighboring particles
	double m_l_eff = 0.;                 // do not refine after this length
	bool m_eccentric_cubic = true;       // a modifier if cubic_basic pattern is used for cutting test
	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-

	refine_criteria m_refine_criteria = moving_frame;
	pattern m_pattern = cubic_basic;

	void flag_reset(body &b) const;
	void my_extrapol_SPH_values(std::vector<particle> &particles, unsigned int my_idx, unsigned int my_dady);

	int scan_mark_div_velocity_based(body &b) const;
	int scan_mark_velocity_based    (body &b) const;
	int scan_mark_vonmises_based    (body &b) const;
	int scan_mark_strain_based      (body &b) const;
	int scan_mark_temperature_based (body &b) const;
	int scan_mark_position_based    (body &b) const;
	int scan_mark_moving_frame      (body &b) const;
	int scan_mark_neighbor_based    (body &b) const;

	void perform_split_triangular(body &b) const;
	void perform_split_cubic_basic(body &b) const;
	void perform_split_cubic(body &b) const;
	void perform_split_hexagonal(body &b) const;

	int evaluate_refinement(body &b) const;
	void do_split(body &b) const;
};

#endif /* ADAPTIVITY_H_ */
