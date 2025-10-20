Contributing
============

This page highlights the key extension points in **MaxwellLink** and the patterns we
recommend when adding new molecular drivers or electromagnetics (EM) solvers. The
codebase follows a “thin core, pluggable back-ends” design: the socket protocol,
shared unit helpers, and :class:`maxwelllink.Molecule` abstraction hide most of the
coupling logic so new components only need to implement domain-specific details.

Source Layout
-------------

- ``src/maxwelllink/mxl_drivers/python/models``: Python molecular drivers (TLS,
  QuTiP, RTTDDFT, ASE, …) that inherit from
  :class:`maxwelllink.mxl_drivers.python.models.dummy_model.DummyModel`.
- ``src/maxwelllink/mxl_drivers/lammps``: C++ socket client 
  (``fix_maxwelllink``) for LAMMPS, inspired by the `i-PI <https://docs.ipi-code.org/>`_ project.
- ``src/maxwelllink/em_solvers``: EM back-ends such as the Meep wrapper and the
  single-mode cavity solver. Each solver ships its own unit system and molecule
  wrapper.
- ``tests/``: Pytest suites that exercise both socket and embedded modes.
- ``docs/source``: User and developer documentation. Please document any new public
  feature here.

Adding a Python Molecular Driver
--------------------------------

Python drivers encapsulate the quantum or classical model that responds to the EM
field. They are loaded either via the ``mxl_driver`` CLI (socket mode) or
instantiated directly through :class:`maxwelllink.Molecule` (embedded mode).

Driver skeleton
~~~~~~~~~~~~~~~

1. Create ``src/maxwelllink/mxl_drivers/python/models/<name>_model.py``.
2. Subclass :class:`~maxwelllink.mxl_drivers.python.models.dummy_model.DummyModel`.

   - Override :meth:`initialize` if extra setup is required once the hub assigns a time step and molecule id.

   - Override :meth:`propagate` to advance the molecular state under the effective
     electric field expressed in **atomic units**.

   - Override :meth:`calc_amp_vector` to return ``dμ/dt`` in **atomic units**. The base
     class handles the stage/commit protocol used by the socket driver.

   - Override :meth:`append_additional_data` if you wish to stream diagnostics back
     to the EM solver (they appear in ``Molecule.additional_data_history``).

   - Implement ``_dump_to_checkpoint`` / ``_reset_from_checkpoint`` when
     checkpoint/restart is desirable.

3. Expose the driver by updating ``__drivers__`` in
   ``src/maxwelllink/mxl_drivers/python/models/__init__.py``. The key becomes the
   ``--model`` value accepted by ``mxl_driver``.

4. If the driver should be importable from ``maxwelllink`` (e.g., ``mxl.QuTiPModel``),
   add a lazy import branch in ``src/maxwelllink/__init__.py``.

5. Add or update docs under ``docs/source/drivers`` describing user-facing
   parameters.

Example template:

.. code-block:: python

   import numpy as np
   from .dummy_model import DummyModel

   class AwesomeModel(DummyModel):
       def __init__(self, omega, damping, verbose=False, **kwargs):
           super().__init__(verbose=verbose, **kwargs)
           self.omega = float(omega)
           self.damping = float(damping)
           self._current_amp = np.zeros(3)

       def initialize(self, dt_new, molecule_id):
           super().initialize(dt_new, molecule_id)
           # Pre-compute any propagators or allocate arrays here.

       def propagate(self, effective_efield_vec):
           # Update internal state using self.dt (atomic units) and the field.
           self._current_amp = self._solve_heisenberg(effective_efield_vec)

       def calc_amp_vector(self):
           return self._current_amp

Ensure every new parameter is documented and exposed through ``__init__`` so the
driver can be constructed from ``--param`` strings (see
:func:`maxwelllink.mxl_drivers.python.mxl_driver.read_args_kwargs`).

Testing tips
~~~~~~~~~~~~

- Add unit tests in ``tests/`` that run the driver in embedded mode (instantiate
  :class:`maxwelllink.Molecule` with ``driver="<name>"``) and, if possible, through
  the socket handshake using ``SocketHub``. The TLS and ASE tests provide patterns.

- Run ``pytest tests/<area>`` before opening a pull request. Socket tests can be
  slower; scope them narrowly around new functionality.

Implementing a New EM Solver
----------------------------

EM solvers orchestrate the Maxwell-time stepping, query molecules for their source
amplitudes, and convert between native units and atomic units. Existing back-ends
(``meep.py`` and ``single_mode_cavity.py``) demonstrate both a grid-based solver and
an ordinary differential equation (ODE) toy model.

Core building blocks
~~~~~~~~~~~~~~~~~~~~

- :class:`~maxwelllink.em_solvers.dummy_em.DummyEMUnits` stores conversion routines.
  Subclass it to translate native electric fields, source amplitudes, and time steps
  into atomic units.

- :class:`~maxwelllink.em_solvers.dummy_em.MoleculeDummyWrapper` wraps
  :class:`maxwelllink.Molecule` instances so the solver can treat socket and embedded
  drivers uniformly.
  
- :class:`~maxwelllink.em_solvers.dummy_em.DummyEMSimulation` offers a minimal
  container; real solvers typically extend it with solver-specific state and a
  :meth:`run` loop.

- :class:`maxwelllink.sockets.SocketHub` handles multiplexing of socket-mode drivers
  and exposes :meth:`register_molecule_return_id` plus the barrier used during
  time-stepping.

Solver skeleton
~~~~~~~~~~~~~~~

1. Create ``src/maxwelllink/em_solvers/<name>.py``.
2. Implement a units helper:

   .. code-block:: python

      from .dummy_em import DummyEMUnits

      class AwesomeUnits(DummyEMUnits):
          def __init__(self, length_unit_nm, time_unit_fs):
              super().__init__()
              self.length_unit_nm = length_unit_nm
              self.time_unit_fs = time_unit_fs

          def efield_em_to_au(self, field_vec3):
              return np.asarray(field_vec3) * self._field_scale

          def source_amp_au_to_em(self, amp_vec3):
              return np.asarray(amp_vec3) / self._field_scale

          def time_em_to_au(self, time_em):
              return time_em * self.time_unit_fs * FS_TO_AU

3. Author a molecule wrapper that derives from
   :class:`~maxwelllink.em_solvers.dummy_em.MoleculeDummyWrapper`. Stamp solver
   metadata onto the molecule (e.g., spatial kernels, native units) and decide how
   to build EM sources in native data structures. Reuse the molecule’s
   ``init_payload`` so socket drivers receive any solver-specific hints.

4. Implement ``AwesomeSimulation`` that wires everything together. Common steps:

   - Normalize the time step, mesh spacing, and :class:`Molecule` wrappers.

   - Split molecules by mode (socket vs. non-socket) and call
     ``m.initialize_driver`` for embedded drivers.

   - During each time step:

     * Gather fields at molecular sites, convert them to atomic units with
       ``AwesomeUnits.efield_em_to_au``, and call ``propagate`` on each wrapper.
     * After the solver advances its state, request source amplitudes (socket mode
       via :class:`SocketHub`, non-socket via ``calc_amp_vector``) and inject them
       back into the EM update.
     * Push any additional per-molecule diagnostics into
       ``Molecule.additional_data_history``.

5. Export the solver by adding a lazy import branch in ``src/maxwelllink/__init__.py``
   so users can instantiate ``mxl.AwesomeSimulation``. Update
   ``docs/source/em_solvers/index.rst`` with a new page describing runtime knobs.

Prefer mirroring the structure in ``MeepSimulation.run``: insert the MaxwellLink
barrier as the first step function so user-provided callbacks still run each step.

Connecting C++ or External Packages
-----------------------------------

External MD or quantum codes written in C++/Fortran can talk to MaxwellLink through
the socket protocol. The LAMMPS driver (``fix_maxwelllink.cpp``) serves as a
production-ready reference that mirrors the i-PI client workflow. Experienced developers
can modify the LAMMPS driver to connect production-level codes to MaxwellLink.


Testing and Documentation
-------------------------

- Extend ``tests/`` with regression coverage for new solvers/drivers. For socket
  clients, prefer lightweight smoke tests that execute a few time steps.

- Update Sphinx docs so users can discover the feature (driver guide, EM solver
  page, release notes if applicable).

- Run ``make -C docs html`` locally to ensure the documentation builds cleanly.
  Address any warnings related to the new content before submitting a pull request.
