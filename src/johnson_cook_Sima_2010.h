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

#ifndef JOHNSON_COOK_SIMA_2010_H_
#define JOHNSON_COOK_SIMA_2010_H_

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>

#include "simulation_data.h"
#include "simulation_time.h"

/*
 This file implements the modified Johnson-Cook (1985) plasticity model.
 This material model captures damage through an additional term that softens the material with increasing strain.
 stress state is return to yield surface using a radial return algorithm.

 The expression of this model can be found in Eqs. (18)-(19) of our paper.
 The model includes 5 extra material parameters in addition to the usual 5 parameters
 needed for the standard Johnson-Cook model.

	Please consider the original publication:
	------------------------------------------------------------------------------
	"Modified material constitutive models for serrated chip formation simulations
	 and experimental validation in machining of titanium alloy ti–6al–4v"
						     By: M. Sima, T. Özel
		Int. J. of Machine Tools & Manufacture, 50 (11) (2010), pp. 943-960
*/


class johnson_cook_Sima_2010 {

private:
	double m_A = 0.;
	double m_B = 0.;
	double m_C = 0.;

	double m_m = 0.;
	double m_n = 0.;

	double m_Tmelt = 0.;
	double m_Tref  = 0.;

	double m_eps_dot_ref = 0.;
	double m_tanh_a = 0.0;
	double m_tanh_b = 0.0;
	double m_tanh_c = 0.0;
	double m_tanh_d = 0.0;

	double m_Sima2010_s = 0.0;

	double m_mu = 0.;
	double m_eps_pl_equiv_init = 0.;
	double m_norm_Strial = 0.;

	double m_t = 0.;
public:

	void set_norm_s_trial(double norm_s_trial);
	void set_eps_init(double eps_init);
	void set_temp(double theta);

	double sigma_yield(double eps_pl, double eps_dot_pl);
	double sigma_yield(double eps_pl, double eps_dot_pl, double theta);

	double operator() (double delta_lambda);
	double f(double delta_lambda);

	johnson_cook_Sima_2010(physical_constants pc);
	johnson_cook_Sima_2010();
};

#endif
