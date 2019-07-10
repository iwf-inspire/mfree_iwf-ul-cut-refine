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

#include <iostream>
#include <stdlib.h>
#include <fenv.h>
#include <sys/time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>

#include "particle.h"
#include "contact.h"
#include "vtk_writer.h"

#include "benchmarks/test_density.h"
#include "benchmarks/test_benches.h"
#include "benchmarks/test_cuttings.h"

#include "tool.h"
#include "logger.h"
#include "body.h"

logger *global_logger;

#include <algorithm>
#include <set>
#include <iterator>

#ifdef __FAST_MATH__
#error "Do NOT compile using -ffast-math"
#endif

int main(int argc, char * argv[]) {
	feenableexcept(FE_INVALID | FE_OVERFLOW);

	// clear the "results" directory for a fresh run
	const char *folder = "results";
	struct stat st = {0};
	if (stat(folder, &st) == -1) {
		mkdir(folder, 0777);
	}
	int ret;
	ret = system("rm results/*.txt");
	ret = system("rm results/*.vtk");
	if(ret!=0) printf("results folder was empty!\n");

	// inputs
	opterr = 0;
	int c;

	int model = 1;

	while ((c = getopt (argc, argv, "n:r:m:h:d:")) != -1) {
		switch (c)
		{
		case 'm':
			sscanf(optarg, "%d", &model);
			printf("running model %d\n", model);
			break;
		}
	}

	int nx = 31;
	assert(model >= 1);
	assert(model <= 4);

	/*
	 ==========================
	 *  set up chosen benchmark
	 *  	this runs model 1-4 in the paper
	 *  	other preliminary simulations are available in test_benches.h
	 *  	density reapproximation tests are aviable in test_density.h
	 ==========================
	 */
	body *b = 0;
	switch (model) {
	case 1:
		b = cutting_ref_single_resol(nx);
		break;
	case 2:
		nx = 61;
		b = cutting_ref_multi_resol_apriori(nx);
		break;
	case 3:
		nx = 61;
		b = cutting_ref_multi_resol_dynamic(nx);
		break;
	case 4:
		nx = 61;
		b = cutting_ref_single_resol(nx);
		break;
	}

	/*
	  ==================================
	 settings of the printout
	 at least [num_print] frames are written out
	 ===================================
	 */
	simulation_time *time = &simulation_time::getInstance();
	unsigned int num_step = time->get_t_final()/time->get_dt();
	int num_print = 150;
	unsigned int freq = num_step / num_print;
	unsigned int print_iter = 0;
	struct timeval begin, end, intermediate;
	gettimeofday(&begin, NULL);

	freq = std::max(1, (int) freq);

	/*
	  ========================
	  (2nd-order) LeapFrog scheme is used
	  for the explicit time integration.
	  ========================
	 */
	leap_frog stepper((*b).get_num_part());

	/*
	 * This is the implementation of the main time-loop,
	 * also illustrated by the following flowchart in the paper:
	 * ---------------------------------------------------------
	 * Section 4:
	 * Fig. 5. Flowchart of the model logic for each time-step.
	 *
	 */
	while(!time->finished()) {

		// plot with given frequency
		if (time->get_step() % freq == 0) {

			if (global_logger) {

				/* Write out the results in the desired format (*.txt, *.vtk)
				 * to be read in Matlab or ParaView
				 */
				global_logger->log(*b, print_iter);

				// Report the time left to finish
				gettimeofday(&intermediate, NULL);
				double seconds_so_far = (intermediate.tv_sec - begin.tv_sec) + ((intermediate.tv_usec - begin.tv_usec)/1000000.0);

				double percent_done = 100*time->get_step()/((double) num_step);
				double time_left = seconds_so_far/percent_done*100;

				printf("%06d: #increments %06d, cur time %e, pctg done %f, seconds left: %f\n", print_iter, time->get_step(), time->get_dt()*time->get_step(), percent_done, time_left-seconds_so_far);
				print_iter++;
			}
		}

		/* Carry out the time-stepper:
		 * this is to update the system by evolving the variables
		 * over time using the LeapFrog time stepping
		 */
		stepper.step(*b);

		time->increment_step();
		time->increment_time();
	}

	gettimeofday(&end, NULL);
	double elapsed = (end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0);
	printf("Runtime: %f\n", elapsed);

	return EXIT_SUCCESS;
}
