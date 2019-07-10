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

#include "vtk_writer.h"

void vtk_writer_write(const std::vector<particle> &particles, unsigned int step, const char *folder) {
	char buf[256];
	sprintf(buf, "%s/out_%06d.vtk", folder, step);
	FILE *fp = fopen(buf, "w+");

	unsigned int np = particles.size();

	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "mfree iwf\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "\n");

	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");		// Particle positions
	fprintf(fp, "POINTS %d float\n", np);
	for (unsigned int i = 0; i < np; i++) {
		fprintf(fp, "%f %f %f\n", particles[i].x, particles[i].y, 0.);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELLS %d %d\n", np, 2*np);
	for (unsigned int i = 0; i < np; i++) {
		fprintf(fp, "%d %d\n", 1, i);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n", np);
	for (unsigned int i = 0; i < np; i++) {
		fprintf(fp, "%d\n", 1);
	}
	fprintf(fp, "\n");

	fprintf(fp, "POINT_DATA %d\n", np);

	fprintf(fp, "SCALARS density float 1\n");		// Current particle density
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (unsigned int i = 0; i < np; i++) {
		fprintf(fp, "%f\n", particles[i].rho);
	}
	fprintf(fp, "\n");

	fprintf(fp, "SCALARS temperature float 1\n");    // Particle temperature
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (unsigned int i = 0; i < np; i++) {
		fprintf(fp, "%f\n", particles[i].T);
	}
	fprintf(fp, "\n");

	fprintf(fp, "SCALARS Svm float 1\n");        // Particle Von Mises stress
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (unsigned int i = 0; i < np; i++) {
		double sxx = particles[i].Sxx - particles[i].p;
		double sxy = particles[i].Sxy;
		double syy = particles[i].Syy - particles[i].p;
		double szz = particles[i].Szz - particles[i].p;

		double svm = sqrt(fabs((sxx*sxx + syy*syy + szz*szz) - sxx * syy - sxx * szz - syy * szz + 3.0 * (sxy*sxy)));
		fprintf(fp, "%f\n", svm);
	}
	fprintf(fp, "\n");

	fprintf(fp, "SCALARS equiv_plastic_strain float 1\n");		// Current particle's equivalent plastic strain
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (unsigned int i = 0; i < np; i++) {
		fprintf(fp, "%f\n", particles[i].eps_pl_equiv);
	}
	fprintf(fp, "\n");

	fprintf(fp, "VECTORS velocity float\n");		// Particle velocities
	for (unsigned int i = 0; i < np; i++) {
		fprintf(fp, "%f %f %f\n", particles[i].vx, particles[i].vy, 0.);
	}
	fprintf(fp, "\n");

	fprintf(fp, "SCALARS glob_density_err double 1\n");  // global density error acc. to Feldman 2006
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (unsigned int i = 0; i < np; i++) {
		fprintf(fp, "%e\n", (particles[i].rho - particles[i].rho_init)*(particles[i].rho - particles[i].rho_init));
	}
	fprintf(fp, "\n");

	fprintf(fp, "SCALARS mass double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (unsigned int i = 0; i < np; i++) {
		fprintf(fp, "%e\n", particles[i].m);
	}
	fprintf(fp, "\n");

	fclose(fp);
}

void vtk_writer_write(const tool* tool, unsigned int step, const char *folder) {
	auto segments = tool->get_segments();
	if (segments.size() == 0) return;

	assert(segments.size() == 4 || segments.size() == 5);

	struct triangle {
		glm::dvec2 p1, p2, p3;
		triangle(glm::dvec2 p1, glm::dvec2 p2, glm::dvec2 p3) : p1(p1), p2(p2), p3(p3) {}
	};

	std::vector<triangle> triangles;

	//mesh tool "body"
	if (segments.size() == 4) {
        triangles.push_back(triangle(segments[0].left, segments[0].right, segments[1].right));
        triangles.push_back(triangle(segments[2].left, segments[2].right, segments[3].right));
	} else if (segments.size() == 5) {
        triangles.push_back(triangle(segments[0].left, segments[0].right, segments[2].right));
        triangles.push_back(triangle(segments[1].left, segments[1].right, segments[2].right));
        triangles.push_back(triangle(segments[3].left, segments[3].right, segments[4].right));
	}

	//mesh fillet
	if (tool->get_fillet() != 0) {
		const int num_discr = 20;
		auto fillet = tool->get_fillet();
		double t1 = fmin(fillet->t1, fillet->t2);
		double t2 = fmax(fillet->t1, fillet->t2);

		double lo = t1 - 0.1*t1;

		double d_angle = (t2-t1)/(num_discr-1);

		double r = fillet->r;

		for (int i = 0; i < num_discr-1; i++) {
			double angle_1 = lo + (i+0)*d_angle;
			double angle_2 = lo + (i+1)*d_angle;

			glm::dvec2 p1 = glm::dvec2(fillet->p.x, fillet->p.y);
			glm::dvec2 p2 = glm::dvec2(p1.x + r*sin(angle_1), p1.y + r*cos(angle_1));
			glm::dvec2 p3 = glm::dvec2(p1.x + r*sin(angle_2), p1.y + r*cos(angle_2));
			triangles.push_back(triangle(p1, p2, p3));
		}
	}

	int num_tri = triangles.size();

	char buf[256];
	sprintf(buf, "%s/tool_%06d.vtk", folder, step);
	FILE *fp = fopen(buf, "w+");

	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "mfree iwf\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp, "POINTS %d float\n", 3*num_tri);

	for (auto it : triangles) {
		fprintf(fp, "%f %f %f\n", it.p1.x, it.p1.y, 0.);
		fprintf(fp, "%f %f %f\n", it.p2.x, it.p2.y, 0.);
		fprintf(fp, "%f %f %f\n", it.p3.x, it.p3.y, 0.);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELLS %d %d\n", num_tri, 3*num_tri + num_tri);
	for (int i = 0; i < num_tri; i++) {
		fprintf(fp, "3 %d %d %d\n", 3*i+0, 3*i+1, 3*i+2);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n", num_tri);
	for (int i = 0; i < num_tri; i++) {
		fprintf(fp, "5\n");
	}

	fclose(fp);

}
