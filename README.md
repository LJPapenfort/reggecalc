reggecalc
=========

Numerical library written in C++, using Boost and Eigen, to simulate the time evolution of spatial hypersurfaces by Regge calculus.

The header-only library relies on C++11 and is multi-threaded, so it should be compiled with the appropriate compiler flags (e.g. on Linux: g++ -m64 -march=native -pthread -std=c++11). It expects the Boost and Eigen libraries in the lib path.

It consists of serveral specialized template classes, implementing structures for the triangulation, the equations and the time evolution. To start the time evolution a hypersurface triangulation of a closed 3-manifold, the initial squard edge lenghts on and between two hyersurfaces and optional additional terms to the Regge equations are needed.

An example of how to use the objects is shown in the main.cpp. A class to construct a QPL hypersurface is given, along with initial edge lengths of the Kasner and Lambda-vacuum universe.

The core classes are
- triangulator_t: 
Once an initial hypersurface triangulation and the edge length are given, the class automatically creates a three-surface Sorkin triangulation, on which the time evolution is carried out. Additionally it defines all involved substructures and simplices, as well as convient methods to gain information about all quantities in the complex.
- regge_functor_t:
This class implements the interface to introduce the Regge equations and customized additional terms. A specialization is given by a cosmological constant term coupled to the Regge equations.
- sorkin_stepper_t:
Takes the triangulation and terms to compute a step-wise time evolution of the hypersurface. Internally the decoupled sets of equations emerging in the Sorkin triangulation are solved for the unknown squared edge lengths, using paralellized threads of Newton-Raphson iterations.

It was written in the course of the master thesis "Cosmology on Simplicial Complexes" at the university of Bonn (GER) and is released under the MPL2.
