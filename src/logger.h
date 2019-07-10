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

#ifndef LOGGER_H_
#define LOGGER_H_

#include "tool.h"
#include "body.h"
#include "vtk_writer.h"

#include <vector>

/*
  Logging for visualization purposes
  ------------------------------------------------
  The logger file supports:
    1. simple text representation of tool
	2. forces on tool
	3. textual "vtk" files for particle attributes
  ------------------------------------------------
*/

class logger {

private:
	bool m_log_forces  = true;
	bool m_emit_vtk    = true;

	tool *m_t = 0;
	FILE *m_fp_forces = 0;
	FILE *m_fp_trace = 0;
	std::vector<unsigned int> m_trace_p;
	char m_folder[256] = "results";
	char m_case_name[256] = "case";

public:
	logger(const char *case_name, const char *foldername = "results");
	void close();

	void set_tool(tool *t);
	void set_log_forces(bool log_forces);
	void set_log_vtk(bool log_vtk);
	void add_tracer_particle(unsigned int tracer_idx);
	void set_folder(const char* folder);

	void log(const body &body, unsigned int step);
};

#endif /* LOGGER_H_ */
