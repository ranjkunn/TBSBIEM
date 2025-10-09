.. TBSBIEM documentation master file, created by
   sphinx-quickstart on Thu Aug 28 15:51:46 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TBSBIEM
############

TBSBIEM is an open-source Fortran code for the simulation of spontaneous 3D dynamic rupture propagation on a planar homogeneous or bi-material interface based on 3D Spectral Boundary Integral Equation Method by `Gupta and Ranjith (2024) <https://doi.org/10.1002/nag.3632>`_. This method is based on traction formulation unlike the displacement formulation of `Geubelle and Rice (1995) <http://www.sciencedirect.com/science/article/pii/002250969500043I>`_. The key advantage of this formulation is that it requires half the convolution integrals to solve bi-material interface dynamic ruture problems when compared to that of displacement formultaion. In addition, all convolution kernels are expressed in closed-form for the bi-material interface problem, whereas some kernels are evaluated numerically in the displacement formulation.

..
.. ################
.. Document Title
.. ################

.. Section Heading
..    ***************

.. Subsection Heading
.. ------------------
.. Subsection Heading
.. ==================
.. Sub-subsection Heading
.. ~~~~~~~~~~~~~~~~~~~~~~


.. toctree::
   :maxdepth: 4
   :caption: Table of Contents:

   ./quickstart
   ./Code_Usage
   ./Formulation
   ./examples
   ./publications
   ./citing
   ./authors
