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

#include "logger.h"

void logger::close() {
	fclose(m_fp_forces);
}

void logger::set_tool(tool *t) {
	m_t = t;
}

void logger::set_log_vtk(bool log_vtk) {
	m_emit_vtk = log_vtk;
}

void logger::set_log_forces(bool log_forces) {
	m_log_forces = log_forces;
}

void logger::add_tracer_particle(unsigned int tracer_idx) {
	m_trace_p.push_back(tracer_idx);
}

void logger::set_folder(const char* folder) {
	strcpy(m_folder, folder);

	fclose(m_fp_forces);

	char buf[256];
	sprintf(buf, "./%s/%s_forces", m_folder, m_case_name);
	m_fp_forces = fopen(buf, "w+");
}

void logger::log(const body &b, unsigned int step) {

	//log forces (if desired)
	if (m_log_forces) {
		double fx = 0.;
		double fy = 0.;

		// sum of X and Y components of both contact & tangential forces
		for (unsigned int i = 0; i < b.get_num_part(); i++) {
			fx += b.get_particles()[i].fcx + b.get_particles()[i].ftx;
			fy += b.get_particles()[i].fcy + b.get_particles()[i].fty;
		}

		simulation_time *time = &simulation_time::getInstance();
		double cur_time = time->get_time();

		fprintf(m_fp_forces, "%e %f %f\n", cur_time, fx, fy);
		fflush(m_fp_forces);
	}

	//trace particles to be traced
	for (const auto it : m_trace_p) {
		fprintf(m_fp_trace, "%f %f ", b.get_particles()[it].x, b.get_particles()[it].y);
	}
	if (m_trace_p.size() != 0) {
		fprintf(m_fp_trace, "\n");
	}

	if (m_emit_vtk) {
		vtk_writer_write(b.get_particles(), step, m_folder);
		if (m_t) {
			vtk_writer_write(m_t, step, m_folder);
		}
	}
}

logger::logger(const char *case_name, const char *foldername) {
	char buf[256];
	sprintf(buf, "./%s/%s_forces", foldername, case_name);
	m_fp_forces = fopen(buf, "w+");
	sprintf(buf, "./%s/trace.txt", foldername);
	m_fp_trace = fopen(buf, "w+");
	strcpy(m_folder, foldername);
	strcpy(m_case_name, case_name);
}
