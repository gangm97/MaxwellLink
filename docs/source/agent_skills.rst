Agent Skills
============

Automatic light-matter simulations on local machines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please watch the following `walkthrough video <https://www.youtube.com/watch?v=ttAvvNByMLg>`_ for an introduction to using MaxwellLink's agent skills with VS Code and the Codex extension:

.. youtube:: ttAvvNByMLg

To use MaxwellLink's agent skills on your local machine via VS Code and the Codex extension, follow these steps:

1. Open VS Code -> ``File`` -> ``Open Folder...`` -> select ``path/to/MaxwellLink``.
2. Install/enable the **Codex** extension (from Marketplace). Make sure the extension has access to the workspace.
3. Open the Codex chat panel (usually the side activity bar) and say what you want, e.g. "Use the MaxwellLink skills to scaffold a Meep + TLS socket job." The agent will load the relevant skills at ``skills/``.
4. When prompted, let the agent run the suggested terminal commands in VS Code's integrated terminal; it will create ``projects/YYYY-MM-DD-<name>/`` and fill template files.

The above video tutorial uses the following input prompt:

.. code-block:: text

   In my local machine, run an initially weakly excited two-level system coupled to 2d vacuum using meep fdtd and plot the excited-state population dynamics


