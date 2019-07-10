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

#ifndef TOOL_H_
#define TOOL_H_

#include <glm/glm.hpp>
#include <vector>
#include <stdio.h>
#include <assert.h>

#include "kernel.h"
#include "simulation_time.h"

/*
 Please take into consideration Fig. (1), Fig. (9), and Fig. (12) of the paper.

 This file shall include/support:
 	 - how tool is parametrized using 4 line segments and possibly a fillet between rake face and clearance.
 	 - construction either by 4 points or top left point + 2 angles + length and height
 	 - queries if a point is inside and contact query (additionally returns penetration depth and direction to closest point on tool)
 	 - the rigid body translation of the tool given the cutting speed
*/

class tool {

private:
	struct line {
		double a = 0.;
		double b = 0.;
		bool vertical = false;

		//return points closest to xq on this line
		glm::dvec2 closest_point(glm::dvec2 xq) const;

		//return intersection point
		glm::dvec2 intersect(line l) const;

		line(double a, double b, bool vertical);
		line();
		// construct line from two points
		line(glm::dvec2 p1, glm::dvec2 p2);
	};

	struct segment {
		glm::dvec2 left;		// left end
		glm::dvec2 right;		// right end
		tool::line l;			// line representation
		glm::dvec2 n;			// normal

		segment(glm::dvec2 left, glm::dvec2 right);
		segment();
		double length() const;
	};

	struct circle_segment {
		double r  = 0.;				//radius
		double t1 = 0.;				//starting angle
		double t2 = 0.;				//end angle
		glm::dvec2 p;			    //center

		// returns distance to circle segment if closest
		// point falls between t1, t2
		// return DBL_MAX otherwise
		double distance(glm::dvec2 qp) const;

		// intersect circle segment with line segment p1,p2
		// returns: 0 if no intersection point falls between t1, t2 or p1,p2 misses segment completely
		//          1 if one intersection point falls between t1, t2. i1 is set
		//          2 if both intersection points fall between t1, t2. i1, i2 is set
		unsigned int intersect(glm::dvec2 p1, glm::dvec2 p2, glm::dvec2 &i1, glm::dvec2 &i2);

		circle_segment(double r, double t1, double t2, glm::dvec2 p);
		circle_segment(glm::dvec2 p1, glm::dvec2 p2, glm::dvec2 p3);
		circle_segment();
	};

private:
	// fit a fillet on line lm to l1 with radius r
	// alternatively: find a point on lm with perpendicular distance r to l1
	glm::dvec2 fit_fillet(double r, line lm, line l1) const;

	// construct segments from list of points
	void construct_segments(std::vector<glm::dvec2> list_p);

	// intersect line tr to br and line br to bl with fillet of radius r
	// (4 points -> 5 points) and set fillet
	std::vector<glm::dvec2> construct_points_and_fillet(glm::dvec2 tl, glm::dvec2 tr, glm::dvec2 br, glm::dvec2 bl, double r);

	// velocity of tool
	glm::dvec2 m_velocity;

	// coordinates of cutting edge of tool
	glm::dvec2 m_edge_coord;

	// friction coefficient
	double m_mu = 0.;

	// components
	circle_segment *m_fillet = 0;
	std::vector<segment> m_segments;

	// experimental
	std::vector<glm::dvec2> m_boundary_particles;
	std::vector<glm::dvec2> m_boundary_normals;
	std::vector<glm::dvec2> m_boundary_weights;

	// if true tool consists of circle segment only
	bool m_chamfer_debug = false;

public:

	struct bbox {
		double bbmin_x = 0.;
		double bbmax_x = 0.;
		double bbmin_y = 0.;
		double bbmax_y = 0.;

		bool in(glm::dvec2 qp);

		// checks if bounding-box (bbox) actually spans a reasonable area
		bool valid() const;

		bbox();
		bbox(glm::dvec2 p1, glm::dvec2 p2);
		bbox(double bbmin_x, double bbmax_x, double bbmin_y, double bbmax_y);
	};

	const std::vector<segment> &get_segments() const;
	const circle_segment *get_fillet() const;

	// project qp onto the tool
	//		qp has to be inside of tool
	glm::dvec2 project(glm::dvec2 qp) const;

	tool::bbox safe_bb(double safety = 0.011) const;

	// return lowest point of tool
	// 		convenient to measure depth of cut (aka feed)
	double low() const;

	// return friction coefficient
	double mu() const;

	// returns distance from qp to tool if qp is inside tool
	// returns -1 otherwise
	// mostly for debugging
	double inside(glm::dvec2 qp) const;

	// visibility test across tool
	bool intersect(glm::dvec2 p1, glm::dvec2 p2) const;
	bool intersect(glm::dvec2 p1, glm::dvec2 p2, double &r) const;

	// establishes contact between tool and query point qp
	// returns false if qp is outside of tool
	// returns true and sets out parameters depth and dir if
	// qp is inside of tool
	bool contact(glm::dvec2 qp, double &depth, glm::dvec2 &dir) const;

	// establishes contact between tool and query point qp
	// returns false if qp is outside of tool
	// returns true and sets out parameters closest point (cp)
	// and normal (n) if qp is inside of tool
	bool contact(glm::dvec2 qp, glm::dvec2 &cp, glm::dvec2 &n) const;

	// move tool with it's velocity
	void update_tool();
	void update_tool(double dt);

	// set velocity of tool
	void set_vel(glm::dvec2 vel);

	// get velocity of tool
	glm::dvec2 get_vel() const;

	// set velocity of tool
	void set_edge_coord(glm::dvec2 coord_tip);

	// get velocity of tool
	glm::dvec2 get_edge_coord() const;

	// returns center (in the sense of center of gravity) of the tool
	//		used for debugging purposes
	glm::dvec2 center() const;

	//chamfer debugging
	void get_chamfer_data(glm::dvec2 &p, double &r) const;
	void set_chamfer(glm::dvec2 cp, double r, double t1, double t2);
	void set_chamfer_debug(bool chamfer_debug);

	// print to file
	void print(FILE *fp);
	void print(unsigned int step, const char *folder = "results");

	// construct tool given by four points and a fillet radius
	tool(glm::dvec2 tl, glm::dvec2 tr, glm::dvec2 br, glm::dvec2 bl, double r, double mu_fric);

	// construct tool given by four points (tool is perfectly sharp)
	tool(glm::dvec2 tl, glm::dvec2 tr, glm::dvec2 br, glm::dvec2 bl, double mu_fric);

	// construct tool given by a reference point, length and height
	// as well as rake and clearance angle (measured from vertically
	// downwards and horizontally leftwards, respectively), and
	// fillet radius r
	// angles are in degrees
	tool(glm::dvec2 tl, double length, double height,
			double rake_angle, double clearance_angle,
			double r, double mu_fric);

	// construct tool given by a reference point, length and height
	// as well as rake and clearance angle (measured from vertically
	// downwards and horizontally leftwards, respectively), the tool
	// is percectly sharp
	// angles are in degrees
	tool(glm::dvec2 tl, double length, double height,
			double rake_angle, double clearance_angle, double mu_fric);

	tool();
};

#endif /* TOOL_H_ */

