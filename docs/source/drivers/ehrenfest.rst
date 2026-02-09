Ehrenfest driver
================

The Ehrenfest driver implements a the Ehrenfest
dynamics model with two electronic states and one classical nuclear coordinate.
It is intended as a lightweight reference implementation for coupling
electron–nuclear nonadiabatic dynamics to explicit electromagnetic field
degrees of freedom in **MaxwellLink**. It can also be extended to include more 
electronic states and nuclear DoFs to deal with more complicated systems.

The driver is provided by
:class:`maxwelllink.mxl_drivers.python.models.Ehrenfest2State1DModel`.

Overview
--------

The model propagates

- a **2×2 electronic density matrix** :math:`\hat{\rho}(t)` in a diabatic
  electronic basis :math:`\{|0\rangle,|1\rangle\}`,
- a **1D classical nuclear coordinate** :math:`R(t)` and momentum :math:`P(t)`,

self-consistently in the mean-field (Ehrenfest) approximation.

Theory
------

Electronic Hamiltonian
^^^^^^^^^^^^^^^^^^^^^^

The molecular electronic Hamiltonian is a 2×2 Hermitian matrix :math:`\hat{H}_0(R)`.
Under the effective local electric field
component :math:`\widetilde{E}(t)`, the total electronic Hamiltonian reads

.. math::

   \hat{H}(R,t) = \hat{H}_0(R) - \widetilde{E}(t)\,\hat{\mu}(R).
:math:`\hat{\mu}(R)` is the dipole operator along the selected polarization axis, 
which also depends on :math:`R`.

Electronic density matrix propagation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The electronic subsystem evolves according to the Liouville–von Neumann equation

.. math::

   \frac{d}{dt}\hat{\rho}(t) = -\frac{i}{\hbar}\Bigl[\hat{H}(R(t),t),\,\hat{\rho}(t)\Bigr].

Density matrix :math:`\hat{\rho}` is propagated using a short-time
unitary propagator :math:`\hat{U}=\exp\!\left(-\frac{i}{\hbar}\hat{H}\Delta t\right)`.

Classical nuclear propagation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The nuclear coordinate is propagated classically,

.. math::

   \dot{R}(t)=\frac{P(t)}{M}, \qquad M\dot{P}(t) = F(R(t),t),

with the force

.. math::

   F(R,t) = -\Big\langle \partial_R \hat{H}(R,t) \Big\rangle
          = -\mathrm{Tr}\!\left[\hat{\rho}(t)\,\partial_R \hat{H}(R,t)\right],

where

.. math::

   \partial_R \hat{H}(R,t) = \partial_R \hat{H}_0(R) - \widetilde{E}(t)\,\partial_R \hat{\mu}(R).

Standard velocity-Verlet algorithm is used to update :math:`(R,P)`.

Dipole moment 
^^^^^^^^^^^^^

The classical dipole (along the selected polarization axis) is calculated as

.. math::

   \mu(t) = \langle \hat{\mu}(R(t))\rangle
          = \mathrm{Tr}\!\left[\hat{\rho}(t)\,\hat{\mu}(R(t))\right].

The gradient of the dipole moment with respect to time passed to Maxwell's equations is evaluated as

.. math::

   \frac{d}{dt}\mu(t)
   =
   \frac{i}{\hbar}\,\mathrm{Tr}\!\left[\hat{\rho}(t)\,[\hat{H}(R(t),t),\hat{\mu}(R(t))]\right]
   + \mathrm{Tr}\!\left[\hat{\rho}(t)\,\partial_R \hat{\mu}(R(t))\right]\dot{R}(t),

with :math:`\dot{R}(t)=P(t)/M`.

This expression includes both the electronic contribution (commutator term) and
the nuclear contribution (explicit :math:`R`-dependence of :math:`\hat{\mu}`).

Requirements
------------

- No additional packages are required beyond **MaxwellLink**'s dependencies.

Usage
-----

Socket mode
^^^^^^^^^^^

.. code-block:: bash

   mxl_driver --model ehrenfest --port 31415 \
     --param "mass=2000, R_initial=2.0, P_initial=0.0, orientation=2, pe_initial=0.0, \
              H0_func=H0_demo, dH0_dR_func=dH0_dR_demo, mu_func=mu_demo, dmu_dR_func=dmu_dR_demo, \
              checkpoint=false, restart=false"

Non-socket mode
^^^^^^^^^^^^^^^

.. code-block:: python

   molecule = mxl.Molecule(
       driver="ehrenfest",
       driver_kwargs={
           "mass": 2000,
           "R_initial": 2.0,
           "P_initial": 0.0,
           "pe_initial": 0.0,
           "orientation": 2,
           "H0_func": H0_demo,
           "dH0_dR_func": dH0_dR_demo,
           "mu_func": mu_demo,
           "dmu_dR_func": dmu_dR_demo,
       },
        # ...
   )

.. note::

   User-defined functions for :math:`\hat{H}_0(R)`, :math:`\partial_R\hat{H}_0(R)`,
   :math:`\hat{\mu}(R)`, and :math:`\partial_R\hat{\mu}(R)` must be made available to the driver
   (e.g., via an importable Python module). See the example in the test directory.

Parameters
----------

.. list-table::
   :header-rows: 1

   * - Name
     - Description
   * - ``mass``
     - Nuclear mass :math:`M` for the classical coordinate :math:`R` (a.u.). Default: ``2000``.
   * - ``R_initial``
     - Initial nuclear coordinate :math:`R(0)` (a.u.). Default: ``0.0``.
   * - ``P_initial``
     - Initial nuclear momentum :math:`P(0)` (a.u.). Default: ``0.0``.
   * - ``orientation``
     - Dipole orientation: ``0`` couples to ``E_x``, ``1`` to ``E_y``, ``2`` to ``E_z``. Default: ``2``.
   * - ``pe_initial``
     - Initial excited-state population. Default: ``0.0``.
   * - ``rho_initial``
     - Optional initial 2×2 complex density matrix overriding ``pe_initial``.
   * - ``H0_func``
     - Callable ``H0_func(R)->(2,2)`` returning the field-free Hamiltonian matrix :math:`\hat{H}_0(R)`.
   * - ``dH0_dR_func``
     - Callable ``dH0_dR_func(R)->(2,2)`` returning :math:`\partial_R\hat{H}_0(R)`.
   * - ``mu_func``
     - Callable ``mu_func(R)->(2,2)`` returning the dipole operator matrix :math:`\hat{\mu}(R)` along the chosen axis.
   * - ``dmu_dR_func``
     - Callable ``dmu_dR_func(R)->(2,2)`` returning :math:`\partial_R\hat{\mu}(R)`.
   * - ``checkpoint``
     - When ``True`` write ``ehrenfest_checkpoint_id_<n>.npz`` after each step. Default: ``False``.
   * - ``restart``
     - When ``True`` resume from the latest checkpoint if present. Default: ``False``.
   * - ``verbose``
     - When ``True`` print diagnostics each step. Default: ``False``.

Returned data
-------------

The driver returns a dictionary (via :meth:`append_additional_data`) with the following:

- ``time_au`` – Simulation time in atomic units.
- ``R_au`` / ``P_au`` – Nuclear coordinate and momentum in atomic units.
- ``force_au`` – Ehrenfest mean-field force in atomic units.
- ``energy_el_au`` – Electronic energy in atomic units.
- ``energy_tot_au`` – Molecular system energy (nuclear kinetic + electronic) in atomic units.
- ``mux_au``, ``muy_au``, ``muz_au`` – Dipole vector components (non-zero along the selected orientation) in atomic units.
- ``Pe`` / ``Pg`` – Excited and ground-state populations.
- ``Pge_real`` / ``Pge_imag`` – Real and imaginary parts of the coherence.

Notes
-----

- The reference implementation is formulated in a diabatic representation.
- For stable long-time propagation, the implementation enforces Hermiticity and unit trace of
  the density matrix after each short-time update.
- When coupling to Maxwell/FDTD, the :math:`d\mu/dt` term is evaluated consistently with
  the instantaneous Hamiltonian and the explicit :math:`R`-dependence of :math:`\hat{\mu}(R)`.