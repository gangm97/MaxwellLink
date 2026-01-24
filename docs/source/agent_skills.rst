Agent Skills
============

Agent Skills for **MaxwellLink** provide a simple way to get started with this code. Using natural language prompts, users can easily create the input files and run jobs in both local machines and HPC systems. Users can then inspect the input files and modify them as needed for more customized simulations.

Prerequisites
~~~~~~~~~~~~~~

This feature requires to clone the MaxwellLink Github repository to access the skills folder. First, either in your local machines or HPC clusters, clone the repository and install **MaxwellLink**:

.. code-block:: bash

   git clone https://github.com/TaoELi/MaxwellLink.git
   cd MaxwellLink
   pip install .

Then, follow :doc:`installation` to install the third-party EM solver (MEEP FDTD) and molecular drivers.

After installation, you can use MaxwellLink's agent skills to scaffold and run light-matter simulations automatically.

Automatic light-matter simulations on local machines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please watch the following `walkthrough video <https://www.youtube.com/watch?v=ttAvvNByMLg>`_ for an introduction to using MaxwellLink's agent skills with VS Code and the Codex extension:

.. youtube:: ttAvvNByMLg

To use MaxwellLink's agent skills on your local machine via VS Code and the Codex extension, follow these steps:

1. Open VS Code -> ``File`` -> ``Open Folder...`` -> select ``path/to/MaxwellLink``.
2. Install/enable the **Codex** extension (from Marketplace). Make sure the extension has access to the workspace.
3. Open the Codex chat panel (usually the side activity bar) and provide your prompt. The agent will load the relevant skills at ``skills/`` and try to accomplish your request.
4. When prompted, let the agent run the suggested terminal commands in VS Code's integrated terminal; it will create ``projects/YYYY-MM-DD-<name>/`` and create input files for **MaxwellLink**.

The above video tutorial uses the following input prompt:

.. code-block:: text

   In my local machine, run an initially weakly excited two-level system coupled to 2d vacuum using meep fdtd and plot the excited-state population dynamics


