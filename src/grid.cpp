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

#include "grid.h"

void grid::assign_hashes(std::vector<particle> &particles , unsigned int n) const {

	for (unsigned int i = 0; i < n; i++) {
		const unsigned int ix = (unsigned int) ((particles[i].x - m_bbmin_x)/m_dx);
		const unsigned int iy = (unsigned int) ((particles[i].y - m_bbmin_y)/m_dx);
		particles[i].hash = ix*m_ny + iy;
	}
}

bool grid::in_bbox(glm::dvec3 qp) const {
	bool in_x = m_bbmin_x <= qp.x && qp.x <= m_bbmax_x;
	bool in_y = m_bbmin_y <= qp.y && qp.y <= m_bbmax_y;
	bool in_z = m_bbmin_z <= qp.z && qp.z <= m_bbmax_z;
	return in_x && in_y && in_z;
}

void grid::get_bbox(glm::dvec3 &bbmin, glm::dvec3 &bbmax) const {
	bbmin.x = m_bbmin_x;
	bbmin.y = m_bbmin_y;
	bbmin.z = m_bbmin_z;

	bbmax.x = m_bbmax_x;
	bbmax.y = m_bbmax_y;
	bbmax.z = m_bbmax_z;
}

void grid::unhash(int idx, unsigned int &i, unsigned int &j) const {
	i = idx/(m_ny);
	j = idx-(i)*(m_ny);
}

const std::vector<int> &grid::get_cells(const std::vector<particle> &particles, unsigned int n)  {

	//needs to be sorted for this to work
	for (unsigned int i = 0; i < n-1; i++) {
		assert(particles[i].hash <= particles[i+1].hash);
	}

	if (m_cells.capacity() < m_num_cell+1) {
		m_cells.resize(std::max((unsigned int) (2*m_cells.size()), 2*(m_num_cell+1)));
	}

	std::fill(m_cells.begin(), m_cells.end(), -1);

	m_cells[particles[0].hash] = 0;

	for (unsigned int i = 0; i < n-1; i++) {
		if (particles[i].hash != particles[i+1].hash) {
			m_cells[particles[i+1].hash] = i+1;
		}
	}
	m_cells[particles[n-1].hash+1] = n;

	//empty boxes are now set to -1
	//in order to iterate through a cell by [cells(cell_index),...,cells(cell_index+1)[
	//those need to be fixed by propagating a "fix" value from the right
	//(such that the above range will just be empty)

	unsigned int fix = n;
	for (auto it = m_cells.rbegin(); it != m_cells.rend(); ++it) {
		if (*it == -1) {
			*it = fix;
		} else {
			fix = *it;
		}
	}

	return m_cells;
}

void grid::update_geometry(const std::vector<particle> &particles, unsigned int n, double kernel_width) {
	double h_max = -DBL_MAX;

	double minx = +DBL_MAX;
	double maxx = -DBL_MAX;
	double miny = +DBL_MAX;
	double maxy = -DBL_MAX;

	for (unsigned int i = 0; i < n; i++) {
		h_max = fmax(particles[i].h, h_max);

		minx = fmin(particles[i].x, minx);
		miny = fmin(particles[i].y, miny);
		maxx = fmax(particles[i].x, maxx);
		maxy = fmax(particles[i].y, maxy);

	}

	//some nudging to prevent round off errors
	m_bbmin_x = minx - 1e-6;
	m_bbmax_x = maxx + 1e-6;
	m_bbmin_y = miny - 1e-6;
	m_bbmax_y = maxy + 1e-6;


	m_dx = h_max*kernel_width;

	m_lx = m_bbmax_x - m_bbmin_x;
	m_ly = m_bbmax_y - m_bbmin_y;

	m_nx = ceil(m_lx/m_dx);
	m_ny = ceil(m_ly/m_dx);
	m_num_cell = m_nx*m_ny;

	assert(m_num_cell < n);
}

void grid::debug_print() const {
	FILE *fp = fopen("grid.txt", "w+");
	for (unsigned int i = 0; i < m_nx; i++) {
		for (unsigned int j = 0; j < m_ny; j++) {
			double x_lo = m_bbmin_x + i*m_dx;
			double x_hi = m_bbmin_x + (i+1)*m_dx;

			double y_lo = m_bbmin_y + j*m_dx;
			double y_hi = m_bbmin_y + (j+1)*m_dx;

			fprintf(fp, "%f %f %f %f\n", x_lo, x_hi, y_lo, y_hi);
		}
	}
	fclose(fp);
}

// constructs verlet lists (particles[i]->nbh). indices are such that they point
// into array sorted by hash
void grid::construct_verlet_lists(std::vector<particle> &particles, unsigned int n, double kernel_width) {

	std::vector<int> cells = get_cells(particles, n);

	for (unsigned int b = 0; b < m_num_cell; b++) {
		unsigned int gi = 0; unsigned int gj = 0;

		unhash((int) b, gi, gj);

		const int low_i  = (int) gi-1 < 0 ? 0 : (int) gi-1;
		const int low_j  = (int) gj-1 < 0 ? 0 : (int) gj-1;
		const int high_i = gi+2 > m_nx ? m_nx : gi+2;
		const int high_j = gj+2 > m_ny ? m_ny : gj+2;

		for (int i = cells[b]; i < cells[b+1]; i++) {

			unsigned int nbh_iter = 0;

			const double hi = particles[i].h;
			const double xi = particles[i].x;
			const double yi = particles[i].y;

			double radius2 = hi*hi*2*2;

			for (int ni = low_i; ni < high_i; ni++) {
				for (int nj = low_j; nj < high_j; nj++) {

					for (int j = cells[ni*m_ny+nj]; j < cells[ni*m_ny+nj+1]; j++) {

						const double xj = particles[j].x;
						const double yj = particles[j].y;

						const double xij = xi-xj;
						const double yij = yi-yj;

						const double r2 = xij*xij + yij*yij;

						if (r2 <= radius2) {
							particles[i].nbh[nbh_iter] = j;
							nbh_iter++;
						}

					}

				}
			}

			assert(nbh_iter < MAX_NBH);

			if (nbh_iter == 0) {
				printf("alarm, particle with no neighbors found!\n");
			}

			particles[i].num_nbh = nbh_iter;
		}
	}
}

unsigned int grid::nx() const {
	return m_nx;
}

unsigned int grid::ny() const {
	return m_ny;
}

double grid::bbmin_x() const {
	return m_bbmin_x;
}

double grid::bbmin_y() const {
	return m_bbmin_y;
}

double grid::dx() const {
	return m_dx;
}

void grid::dbg_print_bbox() const {
	printf("%f %f %f\n", m_bbmin_x, m_bbmin_y, m_bbmin_z);
	printf("%f %f %f\n", m_bbmax_x, m_bbmax_y, m_bbmax_z);
}
