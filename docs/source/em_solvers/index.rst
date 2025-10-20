EM Solvers
==========

**MaxwellLink** ships with two electromagnetic backends. The `Meep <https://meep.readthedocs.io/en/latest/>`_ interface runs a
full finite-difference time-domain (FDTD) grid and streams polarization sources
from **MaxwellLink** molecules, while the single-mode cavity emulator integrates a
single harmonic oscillator entirely in atomic units. Use the pages below for each
EM solver.

.. toctree::
   :maxdepth: 1

   meep
   single_mode_cavity
