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

#include "material_library.h"

physical_constants matlib_steel4430() {	// correct name: Steel 4340
	double E    = 200e9;
	double nu   = 0.29;
	double rho0 = 7830.0;

	double JC_A		= 792.0e6;
	double JC_B		= 510.0e6;
	double JC_C		= 0.014;
	double JC_m		= 1.03;
	double JC_n		= 0.26;
	double Tref		= 273.0;
	double Tmelt	= 1793.0;
	double eps_dot_ref = 1;

	double cp =  477.;
	double tq =  .9;
	double k  =  50.;

	johnson_cook_constants jc(JC_A, JC_B, JC_C, JC_m, JC_n, Tmelt, Tref, eps_dot_ref);
	thermal_constants tc(cp, tq, k);
	return physical_constants(nu, E, rho0, jc, tc);
}

physical_constants matlib_ARMCO_iron() {
	double E = 207e9;
	double nu = 0.29;
	double rho0 = 7890.0;

	double JC_A		= 175.0e6;
	double JC_B		= 380.0e6;
	double JC_C		= 0.06;
	double JC_m		= 0.55;
	double JC_n		= 0.32;
	double Tref		= 273.0;
	double Tmelt	= 1811.0;
	double eps_dot_ref = 1;

	double cp = 452.0;
	double tq = .9;

	johnson_cook_constants jc(JC_A, JC_B, JC_C, JC_m, JC_n, Tmelt, Tref, eps_dot_ref);
	thermal_constants tc(cp, tq);
	return physical_constants(nu, E, rho0, jc, tc);
}

physical_constants matlib_OFHC_copper() {
	double E = 124e9;
	double nu = 0.34;
	double rho0 = 8960.0;

	double JC_A		= 90.0e6;
	double JC_B		= 292.0e6;
	double JC_C		= 0.025;
	double JC_m		= 1.09;
	double JC_n		= 0.31;
	double Tref		= 273.0;
	double Tmelt	= 1356.0;
	double eps_dot_ref = 1;

	double cp = 383.0;
	double tq = .9;

	johnson_cook_constants jc(JC_A, JC_B, JC_C, JC_m, JC_n, Tmelt, Tref, eps_dot_ref);
	thermal_constants tc(cp, tq);
	return physical_constants(nu, E, rho0, jc, tc);
}

physical_constants matlib_AISI1045() {
	double E = 205e9;
	double nu = 0.29;
	double rho0 = 7850.0;

	//	Johnson Cook parameters
	double JC_A		= 615.8e6;
	double JC_B		= 667.7e6;
	double JC_C		= 0.0134;
	double JC_m		= 1.078;
	double JC_n		= 0.255;
	double Tref		= 25.0 + 273.0;
	double Tmelt	= 1350.0 + 273.0;
	double eps_dot_ref = 1;

	double cp = 486.0;
	double tq = .9;

	johnson_cook_constants jc(JC_A, JC_B, JC_C, JC_m, JC_n, Tmelt, Tref, eps_dot_ref);
	thermal_constants tc(cp, tq);
	return physical_constants(nu, E, rho0, jc, tc);
}

physical_constants matlib_rubber() {
	double E = 1e7;
	double nu = 0.4;
	double rho0 = 1;

	return physical_constants(nu, E, rho0, johnson_cook_constants());
}

physical_constants matlib_rubber_real() {
	double E = 1e7;
	double nu = 0.4;
	double rho0 = 1.2e3;

	return physical_constants(nu, E, rho0, johnson_cook_constants());
}

physical_constants matlib_thermal_synthetic() {
	double E = 0;
	double nu = 0.;
	double rho0 = 1.;

	double cp = 1.;
	double tq = 0.;
	double k  = 1.;

	return physical_constants(nu, E, rho0, johnson_cook_constants(), thermal_constants(cp, tq, k));
}

physical_constants matlib_tial6v4_lesuer() {
	double E = 110e9;
	double nu = 0.35;
	double rho0 = 4430.0;

	//	Johnson Cook parameters
	double JC_A		= 1098e6;
	double JC_B		= 1092e6;
	double JC_C		= 0.01;
	double JC_m		= 1.1;
	double JC_n		= 0.93;
	double Tref		= 25.0 + 273.0;
	double Tmelt	= 1678.0 + 273.0;
	double eps_dot_ref = 1;

	double cp = 526.0;
	double tq = .9;
	double k  = 6.8;

	johnson_cook_constants jc(JC_A, JC_B, JC_C, JC_m, JC_n, Tmelt, Tref, eps_dot_ref);
	thermal_constants tc(cp, tq, k);
	return physical_constants(nu, E, rho0, jc, tc);
}

physical_constants matlib_tial6v4_johnson_SI() {
	double E = 110e9;
	double nu = 0.35;
	double rho0 = 4430.0;

	//	Johnson Cook parameters
	double JC_A		= 862e6;
	double JC_B		= 331e6;
	double JC_C		= 0.012;
	double JC_m		= 0.8;
	double JC_n		= 0.34;
	double Tref		= 25.0 + 273.0;
	double Tmelt	= 1678.0 + 273.0;
	double eps_dot_ref = 1;

	double cp = 526.0;
	double tq = .9;
	double k  = 6.8;

	thermal_constants tc(cp, tq, k);

	johnson_cook_constants jc(JC_A, JC_B, JC_C, JC_m, JC_n, Tmelt, Tref, eps_dot_ref);
	return physical_constants(nu, E, rho0, jc, tc);
}

physical_constants matlib_tial6v4_johnson_cm_musec_g() {
	//taken from niks simulation
	double E = 1.100000;
	double nu = 0.35;
	double rho0 = 4.4299998;

	//	Johnson Cook parameters
	double JC_A		= 0.0086200;
	double JC_B		= 0.0033100;
	double JC_C		= 0.0100000;
	double JC_m		= 0.8;
	double JC_n		= 0.34;
	double Tref		= 300.;
	double Tmelt	= 1836.0000;
	double eps_dot_ref = 1e-6;

	double cp = 5.260e-06;
	double tq = .9;
	double k  = 6.8e-13;

	johnson_cook_constants jc(JC_A, JC_B, JC_C, JC_m, JC_n, Tmelt, Tref, eps_dot_ref);
	thermal_constants tc(cp, tq, k);
	return physical_constants(nu, E, rho0, jc, tc);
}

physical_constants matlib_tial6v4_Sima_tanh2010_SI() {	// JC-tanh-Model, Modifikation nach Sima/Özel 2010
	// Unit system: L=[m], t=[s], m=[kg]

	double E = 113.80e09;
	double nu = 0.35;
	double rho0 = 4430.0;

	//	Johnson Cook parameters
	double JC_A		= 724.7e6;		/*!< Tabelle 1 -> Lee-Lin[32] - Sima 2010 */
	double JC_B		= 683.1e6;		/*!< Tabelle 1 -> Lee-Lin[32] - Sima 2010 */
	double JC_C		= 0.03500;		/*!< Tabelle 1 -> Lee-Lin[32] - Sima 2010 */
	double JC_m		= 1.0;			/*!< Tabelle 1 -> Lee-Lin[32] - Sima 2010 */
	double JC_n		= 0.47;			/*!< Tabelle 1 -> Lee-Lin[32] - Sima 2010 */

	double Tref		= 298.0;		/*reference temperature*/
	double Tmelt	= 1878.00;		/*melting temperature*/
	double eps_dot_ref = 1.0;

	double JC_tanh_a = 2.0;			   /*!< Tabelle 3, Sima 2010, Modell 3, fett markierte Werte */
	double JC_tanh_b = 5.0;			   /*!< Tabelle 3, Sima 2010, Modell 3, fett markierte Werte */
	double JC_tanh_c = 2.0;			   /*!< Tabelle 3, Sima 2010, Modell 3, fett markierte Werte */
	double JC_tanh_d = 1.0;			   /*!< Tabelle 3, Sima 2010, Modell 3, fett markierte Werte */
	double JC_tanh_Sima2010_s = 0.05;  /*!< Tabelle 3, Sima 2010, Modell 3, fett markierte Werte */

	double cp = 580.0;				/*specific heat capacity*/
	double tq = .9;					/*taylor-quinny coefficient*/
	double k  = 7.3;				/*heat conductivity*/

	johnson_cook_constants jc(JC_A, JC_B, JC_C, JC_m, JC_n, JC_tanh_a, JC_tanh_b, JC_tanh_c, JC_tanh_d, JC_tanh_Sima2010_s, Tmelt, Tref, eps_dot_ref);
	thermal_constants tc(cp, tq, k);
	return physical_constants(nu, E, rho0, jc, tc);
}

physical_constants matlib_tial6v4_Sima_tanh2010_cm_musec_g() {	// JC-tanh-Model, Modifikation nach Sima/Özel 2010
	double E = 1.1380;
	double nu = 0.35;
	double rho0 = 4.430;

	//	Johnson Cook parameters
	double JC_A		= 0.007247;		/*!< Tabelle 1 -> Lee-Lin[32] - Sima 2010 */
	double JC_B		= 0.006831;		/*!< Tabelle 1 -> Lee-Lin[32] - Sima 2010 */
	double JC_C		= 0.03500;		/*!< Tabelle 1 -> Lee-Lin[32] - Sima 2010 */
	double JC_m		= 1.0;			/*!< Tabelle 1 -> Lee-Lin[32] - Sima 2010 */
	double JC_n		= 0.47;			/*!< Tabelle 1 -> Lee-Lin[32] - Sima 2010 */

	double Tref		= 298.0;		/*reference temperature*/
	double Tmelt	= 1878.00;		/*melting temperature*/
	double eps_dot_ref = 1e-6;

	double JC_tanh_a = 2.0;			   /*!< Tabelle 3, Sima 2010, Modell 3, fett markierte Werte */
	double JC_tanh_b = 5.0;			   /*!< Tabelle 3, Sima 2010, Modell 3, fett markierte Werte */
	double JC_tanh_c = 2.0;			   /*!< Tabelle 3, Sima 2010, Modell 3, fett markierte Werte */
	double JC_tanh_d = 1.0;			   /*!< Tabelle 3, Sima 2010, Modell 3, fett markierte Werte */
	double JC_tanh_Sima2010_s = 0.05;  /*!< Tabelle 3, Sima 2010, Modell 3, fett markierte Werte */

	double cp = 580.0e-6;				/*specific heat capacity*/
	double tq = .9;					/*taylor-quinny coefficient*/
	double k  = 7.3e-13;				/*heat conductivity*/

	johnson_cook_constants jc(JC_A, JC_B, JC_C, JC_m, JC_n, JC_tanh_a, JC_tanh_b, JC_tanh_c, JC_tanh_d, JC_tanh_Sima2010_s, Tmelt, Tref, eps_dot_ref);
	thermal_constants tc(cp, tq, k);
	return physical_constants(nu, E, rho0, jc, tc);
}


physical_constants matlib_dummy() {
	double E = 1;
	double nu = 0.25;
	double rho0 = 1;

	return physical_constants(nu, E, rho0, johnson_cook_constants());
}

physical_constants matlib_a2024t351() {
	double E = 73e9;
	double nu = 0.33;
	double rho0 = 2700;

	double JC_A		= 352e6;
	double JC_B		= 440e6;
	double JC_C		= 0.09;
	double JC_m		= 1.03;
	double JC_n		= 0.42;
	double Tref		= 300;
	double Tmelt	= 1793;
	double eps_dot_ref = 1;

	double cp = 875;
	double tq = .9;
	double k = 0.;

	johnson_cook_constants jc(JC_A, JC_B, JC_C, JC_m, JC_n, Tmelt, Tref, eps_dot_ref);
	thermal_constants tc(cp, tq, k);
	return physical_constants(nu, E, rho0, jc, tc);
}
