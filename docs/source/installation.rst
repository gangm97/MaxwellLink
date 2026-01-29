Installation
============

The recommended way to install **MaxwellLink** is through a conda environment that
provides third-party EM solvers and molecular drivers alongside Python 3.9+.

Recommended Installation Option
---------------------------------

Prerequisites
~~~~~~~~~~~~~~~

- A recent Python (``>=3.9``). **MaxwellLink** is tested with CPython.
- A working MPI stack (e.g. MPICH or OpenMPI) whenever you plan to run Meep or
  **MaxwellLink** under ``mpirun``.
- C/Fortran toolchains supplied by your platform if you plan to build drivers
  such as LAMMPS from source.

Create a conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Create a fresh conda environment for MaxwellLink
   CONDA_ENV="mxl"
   conda create -n "$CONDA_ENV" python=3.13
   conda activate "$CONDA_ENV"

Install **MaxwellLink**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install maxwelllink

This installs the Python package together with the ``mxl_driver`` console entry
point used to launch molecular drivers.

Optional EM solvers
~~~~~~~~~~~~~~~~~~~~~

Currently, **MaxwellLink** ships with `Meep <https://meep.readthedocs.io/en/latest/>`_ and a simple single-mode cavity solver. While the built-in
single-mode cavity is mainly for prototyping and debugging, `Meep <https://meep.readthedocs.io/en/latest/>`_ is a full-featured FDTD package
suitable for production simulations. One can install `Meep <https://meep.readthedocs.io/en/latest/>`_ with MPI support via conda:

.. code-block:: bash

   conda install -c conda-forge pymeep="*=mpi_mpich_*"

.. note::
   Users can also install `Meep <https://meep.readthedocs.io/en/latest/>`_  from source following the
   `official instructions <https://meep.readthedocs.io/en/latest/Installation/>`_.

Optional driver dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A two-level system (TLS) model ships with **MaxwellLink** and does not require extra packages.
Beyond this lightweight driver, **MaxwellLink** supports several molecular drivers that depend on
third-party packages. Install any molecular drivers below that you want to use from the command line. Each molecular driver can
be pulled into the same conda environment.

.. code-block:: bash

   # Model Hamiltonians via QuTiP
   conda install -c conda-forge qutip

   # Psi4 for RT-TDDFT and Ehrenfest drivers
   conda install -c conda-forge psi4

   # Atomic Simulation Environment (ASE) for classical/MM drivers
   conda install -c conda-forge ase

   # LAMMPS driver with fix mxl (installs a custom binary lmp_mxl)
   mxl_install_lammps

The `LAMMPS <https://www.lammps.org/>`_ helper downloads, patches, and builds a `LAMMPS <https://www.lammps.org/>`_ executable that contains
``fix mxl``; alternatively copy the provided ``fix_maxwelllink.cpp`` and ``fix_maxwelllink.h`` files in source code (src/maxwelllink/mxl_drivers/lammps/) into
your existing `LAMMPS <https://www.lammps.org/>`_ build and recompile.


Install **MaxwellLink** from source
--------------------------------------

In the same conda environment, if you are interested in chatting with or contributing to **MaxwellLink**, clone the repository and install the package:

.. code-block:: bash

   git clone https://github.com/TaoELi/MaxwellLink.git
   cd MaxwellLink
   pip install .


Enter the developer mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install developer dependencies to set up a development environment:

.. code-block:: bash

   pip install ".[dev]"


Testing and Validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After contributing to the code, add a test file in the ``tests/`` folder and run the test suite to ensure everything works as expected.

Use ``pytest`` to run the tests. From the root directory of the repository, run:

.. code-block:: bash

   # Run the fast tests (TLS under socket/non-socket modes)
   pytest -m core -v

   # Run the optional Psi4 and QuTiP tests when the dependencies are available
   pytest -m "optional or slow" -v

If the optional tests are skipped, confirm that the corresponding packages were
installed or adjust your ``PYTHONPATH`` so that the drivers can import them.


Build the Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The documentation can be built locally. Use the optional `docs` extra to install the Sphinx toolchain:

.. code-block:: bash

   pip install ".[docs]"


Then generate the API and HTML pages using Sphinx and open them in your default browser:

.. code-block:: bash

   make doc html


Vibe simulations
--------------------------------------

Users can run ``vibe simulations``, i.e., using natural language prompts to set up and run
**MaxwellLink** simulations. To enable this feature, install either `VS Code IDE + Codex extension <https://developers.openai.com/codex/ide>`_ or `Codex CLI <https://developers.openai.com/codex/cli/>`_.

Then, open Codex at the root directory of the **MaxwellLink** repository (requiring installing **MaxwellLink** from source) and start to chat:

.. code-block:: text

   In my local machine, run an initially weakly excited two-level system coupled to 2d vacuum using meep fdtd and plot the excited-state population dynamics

It supports ``vibe simulations`` on both local machines and HPC clusters. See :doc:`agent_skills` for more details.

.. note::

   For users running ``vibe simulations`` on local machines, usually the agent runs on a sandboxed environment. For ``vibe simulations`` on local machines (**not on HPC**), this might conflict with the MPI environment used by `Meep <https://meep.readthedocs.io/en/latest/>`_, leading to failed simulations.
   To resolve this issue, consider installing the serial version of `Meep <https://meep.readthedocs.io/en/latest/>`_ (i.e., without MPI support) instead:

   .. code-block:: bash

      conda install -c conda-forge pymeep="*=nompi_*"

Testing and Validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, ``pytest`` does not run unit tests for ``vibe simulations`` that require Codex. To run these tests, first install ``codex`` (either the VS Code extension or the CLI tool), 
then make sure you have logged in with your OpenAI API key or account credentials. Finally, run the unit tests with the ``agent`` marker:

.. code-block:: bash

   # run agent tests, default run on local machines
   pytest -m agent -v
   # run agent tests on local machines
   pytest -m agent --codex-prompts local -v
   # run agent tests on hpc machines
   pytest -m agent --codex-prompts hpc -v
   # run agent tests by skipping the passed tests in the last run
   pytest -m agent -v --codex-resume

.. note::

   A **passed** simulation here means that the agent was able to generate a working simulation script and get some results. It does **not guarantee** that the simulation results are **physically accurate**.


Running a batch of prompts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users can also use this unit test feature to automatically run a set of Codex prompts (each line corresponds to one independent prompt) stored in a text file:

.. code-block:: bash

   # run agent tests, default run on local machines
   pytest -m agent -v --codex-prompts-file path/to/your/prompts.txt

The output is stored at ``tests/test_agents/your_prompt_file_name/`` by default.