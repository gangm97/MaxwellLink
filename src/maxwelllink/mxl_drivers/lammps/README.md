# LAMMPS driver connecting with FDTD engines via MaxwellLink

To use this functionality, please copy **fix_maxwelllink.h** and **fix_maxwelllink.cpp** to the LAMMPS source code (lammps/src/MISC/) and then recompile LAMMPS.

Then the recompiled LAMMPS code can connect with MaxwellLink SocketHub via the following fix (similar as **fix ipi**):
```bash
fix 1 all mxl host port
```
Here, **host** is the IP address of the machine where the FDTD engine is running (such as **localhost**), and **port** should match the port number in **SocketHub**. Please check [../README.md](../README.md) for details regarding setting the **SocketHub**.

