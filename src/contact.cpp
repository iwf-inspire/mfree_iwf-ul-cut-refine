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

#include "contact.h"

static glm::dvec2 compute_contact_force_nianfei(const tool *master, double pen_depth, glm::dvec2 surf_norm, double alpha, double ms, double dt) {
	// friction force according to
	// "3D adaptive RKPM method for contact problems with elastic–plastic dynamic
	// large deformation" - Nianfei, Guangyao, Shuyao

	const glm::dvec2 n = surf_norm;
	const double gN = pen_depth;

	double dt2 = dt*dt;
	glm::dvec2 fN = -ms*gN*n/dt2*alpha;

	return fN;
}

static glm::dvec2 compute_friction_ldyna(const tool *master, glm::dvec2 fN, glm::dvec2 n, glm::dvec2 vs, glm::dvec2 fricold, double alpha, double ms, double dt, double mu) {
	if (mu == 0.) return glm::dvec2(0.);

	glm::dvec2 vm = master->get_vel();
	glm::dvec2 v = vs-vm;
	glm::dvec2 vr = v - v*n;

	glm::dvec2 kdeltae = alpha*ms*vr/dt;
	double fy = mu*glm::length(fN);
	glm::dvec2 fstar = fricold - kdeltae;

	if (glm::length(fstar) > fy) {
		return fy*fstar/glm::length(fstar);
	} else {
		return fstar;
	}
}

void contact_apply_tool_to_body_2d(const tool *master, body &slave) {
	simulation_time *time = &simulation_time::getInstance();
	double dt = time->get_dt();

	std::vector<particle> &particles = slave.get_particles();

	const double alpha = 0.1;

	for (unsigned int i = 0; i < slave.get_num_part(); i++) {

		double qx = particles[i].x;
		double qy = particles[i].y;

		glm::dvec2 xcntct(0.);
		glm::dvec2 surf_norm(0.);
		glm::dvec2 xslave(qx, qy);

		bool inside = master->contact(xslave, xcntct, surf_norm);

		/*
		 both CONTACT & TANGENTIAL forces are ZERO
		 for particles which are outside the tool
		*/
		if (!inside) {
			particles[i].fcx = 0.;
			particles[i].fcy = 0.;

			particles[i].ftx = 0.;
			particles[i].fty = 0.;

			continue;
		}

		double pen_depth = glm::dot((xslave-xcntct), surf_norm);
		glm::dvec2 fricold(particles[i].ftx, particles[i].fty);

		double ms   = particles[i].m;

		glm::dvec2 vs(particles[i].vx, particles[i].vy);

		glm::dvec2 cntc = compute_contact_force_nianfei(master, pen_depth, surf_norm, alpha, ms, dt);
		glm::dvec2 fric = compute_friction_ldyna(master, cntc, surf_norm, vs, fricold, alpha, ms, dt, master->mu());

		// X and Y components of the contact force
		particles[i].fcx = cntc.x;
		particles[i].fcy = cntc.y;
		// X and Y components of the tangential force
		particles[i].ftx = fric.x;
		particles[i].fty = fric.y;
	}

}
