Single-Mode Cavity
==================

The :mod:`maxwelllink.em_solvers.single_mode_cavity` module provides a
lightweight electromagnetic solver that replaces a full FDTD grid with a single
damped harmonic oscillator. It evolves entirely in atomic units while
supporting the same :class:`maxwelllink.Molecule` abstraction used by the Meep
backend, making it well suited for rapid prototyping, regression tests, and
workflows focused on one classical cavity mode.

.. note::

  The cavity replaces the full EM grid with canonical coordinate :math:`q_c` and momentum :math:`p_c` obeying

  .. math::

     \dot{q}_c = p_c, \qquad
     \dot{p}_c = -\omega_c^{2} q_c - \kappa\, p_c + g \sum_{m} \frac{d\mu_{m}}{dt} + D(t),

  where :math:`\omega_c` is ``frequency_au``, :math:`\kappa` is ``damping_au``, :math:`g` is ``coupling_strength``, and :math:`D(t)` is the optional external drive. The sum runs over the selected dipole component of each coupled molecule. The effective electric field returned to the drivers is

  .. math::

     E(t) = - g\, p_c(t),

  ensuring energy balance between the cavity mode and the molecular dipole currents.

Requirements
------------

- No additional dependencies beyond the **MaxwellLink** Python stack are required.

Usage
-----

Socket mode
^^^^^^^^^^^

.. code-block:: python

   import maxwelllink as mxl

   hub = mxl.SocketHub(host="127.0.0.1", port=31415, timeout=10.0, latency=1e-5)
   molecule = mxl.Molecule(hub=hub)

   sim = mxl.SingleModeSimulation(
       dt_au=0.5,
       frequency_au=0.242,
       damping_au=0.0,
       molecules=[molecule],
       coupling_strength=1e-4,
       qc_initial=1e-5,
       hub=hub,
       record_history=True,
   )

   sim.run(steps=4000)

Launch :command:`mxl_driver --model tls --port 31415 ...` (or another driver)
after running this script.

Non-socket mode
^^^^^^^^^^^^^^^

.. code-block:: python

   import maxwelllink as mxl

   tls = mxl.Molecule(
       driver="tls",
       driver_kwargs=dict(
           omega=0.242,
           mu12=187.0,
           orientation=2,
           pe_initial=1e-4,
       ),
   )

   sim = mxl.SingleModeSimulation(
       dt_au=0.5,
       frequency_au=0.242,
       damping_au=0.0,
       molecules=[tls],
       coupling_strength=1e-4,
       qc_initial=1e-5,
       record_history=True,
   )

   sim.run(steps=4000)

Parameters
----------

.. list-table::
   :header-rows: 1

   * - Name
     - Description
   * - ``dt_au``
     - Integration time step in atomic units. Must be positive.
   * - ``frequency_au``
     - Cavity angular frequency :math:`\omega_c` (a.u.).
   * - ``damping_au``
     - Linear damping coefficient :math:`\kappa` applied to the cavity momentum (a.u.).
   * - ``molecules``
     - Iterable of :class:`~maxwelllink.Molecule` instances. Socket and non-socket
       molecules can be mixed in the same simulation.
   * - ``drive``
     - Constant float or callable ``drive(time_au)`` returning the external drive term.
       Defaults to ``0.0``.
   * - ``coupling_strength``
     - Scalar prefactor :math:`g` for the molecular polarization feedback. Default: ``1.0``.
   * - ``coupling_axis``
     - Field component coupled to the molecules: ``0/1/2`` or ``"x"/"y"/"z"``. Default: ``2`` (``z``).
   * - ``hub``
     - Optional :class:`~maxwelllink.SocketHub` shared by all socket-mode molecules.
       The simulation infers the hub from the first socket molecule when omitted.
   * - ``qc_initial``
     - Initial cavity coordinate :math:`q_c(0)` (a.u.). Default: ``0.0``.
   * - ``pc_initial``
     - Initial cavity momentum :math:`\dot{q}_c(0)` (a.u.). Default: ``0.0``.
   * - ``record_history``
     - When ``True`` store histories for time, field, momentum, drive, and net molecular response.
       Default: ``True``.

Returned data
-------------

Calling :class:`~maxwelllink.SingleModeSimulation` with ``record_history=True``
populates:

- :attr:`SingleModeSimulation.time_history` – time stamps in atomic units.
- :attr:`SingleModeSimulation.qc_history` – cavity coordinate :math:`q_c(t)`.
- :attr:`SingleModeSimulation.pc_history` – cavity momentum :math:`\dot{q}_c(t)`.
- :attr:`SingleModeSimulation.drive_history` – external drive values.
- :attr:`SingleModeSimulation.molecule_response_history` – summed molecular response along ``coupling_axis``.

Each :class:`~maxwelllink.Molecule` keeps
:attr:`~maxwelllink.Molecule.additional_data_history`, which records driver
data (e.g., TLS populations, energies, timestamps). Socket drivers may
append extra JSON payloads through the hub; embedded drivers populate the same
history via :meth:`maxwelllink.Molecule.append_additional_data`.

Notes
-----

- Provide either ``steps`` or ``until`` to :meth:`SingleModeSimulation.run`, not both.
- Socket-mode molecules must all bind to the same :class:`~maxwelllink.SocketHub`;
  the simulation waits until every driver acknowledges initialization.
- Setting ``record_history=False`` avoids list allocations for throughput-critical runs.
