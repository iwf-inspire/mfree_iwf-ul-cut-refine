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

#include "johnson_cook_Sima_2010.h"

void johnson_cook_Sima_2010::set_norm_s_trial(double norm_s_trial) {
	m_norm_Strial = norm_s_trial;
}

void johnson_cook_Sima_2010::set_eps_init(double eps_init) {
	m_eps_pl_equiv_init = eps_init;
}

void johnson_cook_Sima_2010::set_temp(double t) {
	m_t = t;
}

double johnson_cook_Sima_2010::operator() (double delta_lambda) {
	simulation_time *time = &simulation_time::getInstance();
	double delta_t = time->get_dt();

	double eps_pl_equiv = m_eps_pl_equiv_init + sqrt(2.0/3.0) * fmax(delta_lambda,0.);
	double eps_pl_equiv_dot = sqrt(2.0/3.0) *  fmax(delta_lambda,0.) / delta_t;

	double sigmaY = sigma_yield(eps_pl_equiv, eps_pl_equiv_dot);

	return m_norm_Strial - 2*m_mu*delta_lambda - sqrt(2.0/3.0)*sigmaY;
}

double johnson_cook_Sima_2010::sigma_yield(double eps_pl, double eps_dot_pl) {
	return sigma_yield(eps_pl, eps_dot_pl, m_t);
}

double johnson_cook_Sima_2010::sigma_yield(double eps_pl, double eps_dot_pl, double t) {
	if (t < m_Tref) {
		printf("Temp (%e) smaller than JC-Reference Temp (%e) - adapted\n", t, m_Tref);
		t = m_Tref;
	}

	double theta = (t - m_Tref)/(m_Tmelt - m_Tref);

	assert (eps_pl >= 0.);

	double corr_B = (m_tanh_a > 0.) ? 1.0/exp(pow(eps_pl, m_tanh_a)) : 1;		// slightly dirty hack to fall back to standard JC
																				// if additional parameters are not defined
	double Term_A = m_A + m_B * pow(eps_pl, m_n) * corr_B;
	double Term_B = 1.0;

	double eps_dot = eps_dot_pl / m_eps_dot_ref;

	if (eps_dot > 1.0) {// Catch potential error at analysis start when plastic strain rate is still 0
		assert(eps_dot >= 0.);
		Term_B = 1.0 + m_C * log(eps_dot);
	} else {
		assert(1.0 + eps_dot >= 0.);
		Term_B = pow((1.0 + eps_dot),m_C);
	}

	assert(theta >= 0.);
	double Term_C = 1.0 - pow(theta, m_m);

	double tanh_D = 1.0 - pow((t/m_Tmelt), m_tanh_d);
	double tanh_S =	pow(t/m_Tmelt, m_tanh_b);

	double Term_D = tanh_D + (1.0 - tanh_D) * pow(tanh(1.0/(pow(eps_pl + tanh_S, m_tanh_c))),m_Sima2010_s);

	return Term_A * Term_B * Term_C * Term_D;
}

double johnson_cook_Sima_2010::f(double delta_lambda) {
	return this->operator ()(delta_lambda);
}

johnson_cook_Sima_2010::johnson_cook_Sima_2010(physical_constants pc) {
	m_A = pc.jc().A();
	m_B = pc.jc().B();
	m_C = pc.jc().C();

	m_m = pc.jc().m();
	m_n = pc.jc().n();

	m_Tref  = pc.jc().Tref();
	m_Tmelt = pc.jc().Tmelt();

	m_tanh_a = pc.jc().a();
	m_tanh_b = pc.jc().b();
	m_tanh_c = pc.jc().c();
	m_tanh_d = pc.jc().d();

	m_Sima2010_s	= pc.jc().s();
	m_eps_dot_ref = pc.jc().eps_ref();

	m_mu = pc.G();
}

johnson_cook_Sima_2010::johnson_cook_Sima_2010() {}
