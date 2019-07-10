# mfree_iwf-ul_cut-refine

This public repository provides the source code to the [publication](https://www.sciencedirect.com/science/article/pii/S0020740319317023) **Meshfree Simulation of Metal Cutting: An Updated Lagrangian Approach with Dynamic Refinement** published by International Journal of Mechanical Sciences ([IJMS](https://www.sciencedirect.com/journal/international-journal-of-mechanical-sciences)). The purpose of this package lies in the effective tailoring of the spatial refinement in the updated Lagrangian frame, tuned for the meshfree simulation of an orthogonal metal cutting problem. Towards this end, a multiplicity of numerical algorithms/stabilizations, as well as a mature physical modeling was employed from the state of the art which are not available in commercial meshfree toolkits like the ones available in LSDYNA or ABAQUS. Some of which include:

* A modified Johnson-Cook constitutive modeling for TiAl6V4, through which the strain softening phenomenon resulting from damage in the machining of titanium alloys can be addressed. This flow rule was first introduced by Sima, M., and Özel, T., "Modified material constitutive models for serrated chip formation simulations and experimental validation in machining of titanium alloy Ti–6Al–4V." International Journal of Machine Tools and Manufacture 50.11 (2010): 943-960.
* Stabilization of the solution using the techniques presented by Gray, J. P., J. J. Monaghan, and R. P. Swift. "SPH elastic dynamics." Computer methods in applied mechanics and engineering 190.49-50 (2001): 6641-6662.
* Thermal coupling by solving the heat conduction in the workpiece with either Particle Strength Exchange or the Brookshaw's scheme.
* Increasing the efficiency of the simulation by tailoring the dynamic refinement algorithm via particle splitting, according to the procedure described by Feldman, J., and J. Bonet. "Dynamic refinement and boundary contact forces in SPH with applications in fluid flow problems." International Journal for Numerical Methods in Engineering 72.3 (2007): 295-324. and Vacondio, R., et al. "Variable resolution for SPH: a dynamic particle coalescing and splitting scheme." Computer Methods in Applied Mechanics and Engineering 256 (2013): 132-148.

As a preliminary benchmark study, 4 refinement patterns (i.e., triangular, cubic, extended cubic, and hexagonal) with uniform mass distribution were first cross-compared in a density approximation error analysis. In a unit square discretized by 11x11 particles, the density error introduced by particle splitting was demonstrated:

![density](https://raw.githubusercontent.com/mroethli/mfree_iwf-ul-cut-refine/master/img/density.png) 

Sketch of the cutting geometry at hand:

![cutting_sketch](https://raw.githubusercontent.com/mroethli/mfree_iwf-ul-cut-refine/master/img/cutting_sketch2.png)

Our initial investigation revealed that the correct chip morphology in this cutting application can be observed by increasing the resolution. Red using ~6200 and blue using ~24000 particles, in the cutting simulation of TiAl6V4 after 1 mm of cut:

![superimposed](https://raw.githubusercontent.com/mroethli/mfree_iwf-ul-cut-refine/master/img/superimposed.png)

Therefore, the cubic pattern together with a moving refinement frame was chosen as the settings for particle refinement in the metal cutting test. By saving up to ~70% of the computational cost using dynamic refinement, this approach allows for remarkable runtime optimization compared to conventional single-resolution simulations. Color depicts the equivalent plastic strain, limited to 100%. Models from left to right: single low-resolution, dynamic refinement, _a priori_ refined configuration, single high-resolution: 

![all_cuts](https://raw.githubusercontent.com/mroethli/mfree_iwf-ul-cut-refine/master/img/all_cuts.png)

An overview of the exemplary results for 1 mm of cut at a cutting speed of 500 m/min is as follows. 

**Benchmarking Overview**

|Model 	                | Resolution at start | Resolution at end | Runtime |
| ---------------------:| -------------------:|------------------:| -------:|
|  single low-resolution|               ~6'200|             ~6'200|       59|
|     dynamic refinement|               ~6'800|            ~11'900|      182|
|  _a priori_ refinement|              ~15'800|            ~15'800|      362|
| single high-resolution|              ~24'400|            ~24'400|      540|



Runtimes are measured and reported in CPU minutes, implemented in C++14, taken on a single core of Intel Core i5-4690 at 3.50 GHz.

Result frames presented above can be viewed using [ParaView](https://www.paraview.org/) using the legacy VTK format.

**mfree_iwf-ul_cut-refine** was tested on various versions of Ubuntu Linux. The only dependency is [GLM](https://glm.g-truc.net/0.9.9/index.html). Make files for both a Release version and a Debug build are provided. **mfree_iwf-ul_cut-refine** was developed at _IWF_ [ETHZ](www.ethz.ch) by the following authors:

* Mohamadreza Afrasiabi, afrasiabi@ethz.ch
* Matthias Röthlin, mroethli@ethz.ch
* Hagen Klippel, hklippel@ethz.ch

**mfree_iwf-ul_cut-refine** is free software and licensed under GPLv3.
