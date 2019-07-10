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

#ifndef GRID_H_
#define GRID_H_

#include <glm/glm.hpp>

#include <algorithm>
#include <vector>
#include <float.h>
#include <math.h>
#include <assert.h>

#include "particle.h"
#include "tool.h"

/*
 * This is the implementation of the neighboring search algorithm used.
 * The implementation relies on a simple spatial hashing algorithm & the cell-list data structure.
 * A short explanation in this regard can be found in the paper, right before Section 4.3.

 In short:
 =========
   1- A cell list structure is constructed (std::vector<int> cells from get_cells())
   2- These cell lists are then used to produce verlet lists (list of neighbors saved on each particle)
 */

class grid {

	friend class body;

private:

	std::vector<int> m_cells;

	void assign_hashes(std::vector<particle> &particles, unsigned int n) const;
	bool in_bbox(glm::dvec3 qp) const ;
	void get_bbox(glm::dvec3 &bbmin, glm::dvec3 &bbmax) const;
	void unhash(int idx, unsigned int &i, unsigned int &j) const;
	const std::vector<int> &get_cells(const std::vector<particle> &particles, unsigned int n);
	void update_geometry(const std::vector<particle> &particles, unsigned int n, double kernel_width);
	void debug_print() const;

	// constructs verlet lists (particles[i]->nbh). indices are such that they point
	// into array sorted by hash
	void construct_verlet_lists(std::vector<particle> &particles, unsigned int n, double kernel_width = 2.);

	unsigned int nx() const;
	unsigned int ny() const;

	double bbmin_x() const;
	double bbmin_y() const;

	double dx() const;

	void dbg_print_bbox() const;


	double m_dx;					/*!< Box dimension (edge length of cube) */
	double m_lx, m_ly, m_lz;		/*!< Size in x/y/z direction */
	unsigned int m_nx, m_ny, m_nz;	/*!< Number of boxes in x/y/z direction */
	unsigned int m_num_cell;		/*!< Total number of boxes */

	double m_bbmin_x, m_bbmax_x;	/*!< Min/Max coordinate in x- direction */
	double m_bbmin_y, m_bbmax_y;	/*!< Min/Max coordinate in y- direction */
	double m_bbmin_z, m_bbmax_z;	/*!< Min/Max coordinate in z- direction */
};

#endif /* GRID_H_ */
