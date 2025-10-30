Drivers
=======

**MaxwellLink** ships with several Python molecular drivers that can be accessed either
through the ``mxl_driver`` console script (socket mode) or directly via
:class:`~maxwelllink.molecule.molecule.Molecule` in non-socket mode. These Python drivers inherit from
:class:`~maxwelllink.mxl_drivers.python.models.dummy_model.DummyModel` and provide the unified driver API.
Each driver page documents its parameters, expected outputs, and any runtime prerequisites.

.. toctree::
   :maxdepth: 1

   tls
   qutip
   rttddft
   rtehrenfest
   ase

Additionally, **MaxwellLink** provides a direct connection to `LAMMPS <https://www.lammps.org/>`_ via the custom
``fix mxl``. This C++ driver only supports socket mode.

.. toctree::
   :maxdepth: 1

   lammps
