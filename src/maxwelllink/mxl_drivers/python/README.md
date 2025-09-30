# Python drivers connecting with FDTD engines via MaxwellLink

A few python drivers for the coupled Maxwell-molecular dynamics simulations using the external FDTD engine are provided here.

**All parameters in python drivers are in atomic units**. MaxwellLink converts the FDTD units to atomic units when sending data to Python drivers. Similarly, Python drivers send data to MaxwellLink in atomic units, and MaxwellLink converts the data to FDTD units for propagating Maxwell's equations.

## **electronic two-level system (tls)** 

A sample bash input for setting up the **tls** model is as follows: 

```bash
mxl_driver.py --model tls --port 31415 --param "omega=0.242, mu12=187, orientation=2, pe_initial=0.01" --verbose
```

1. The **--port** number should match that in **SocketHub**; see [../README.md](../README.md) for the initialization of **SocketHub** in MaxwellLink. 

2. The **--param** input provides the initialization of the Python class of the two-level system: [models/tls_model.py](./models/tls_model.py). 

3. **orientation=2** means setting the dipole vector to orient along the **z** (0-x, 1-y, 2-z) direction. 

4. **pe_initial=0.01** sets the initial electronic excited-state population.

5. If **checkpoint=true**, after each step, the necessary checkpoint data will be written to disk.

6. If **restart=true** and **checkpoint=true**, the driver code can restart from the checkpoint file in disk and resume the simulation. This setting is necessary when many  drivers are connected to the FDTD engine at the same time. If one driver is terminated in one machine (by different reasons), all the other drivers will pause the simulation and wait for the restart of this driver. 

7. If the driver and the FDTD code are running in different HPC nodes, please also set **--address <FDTD_CODE_IP_OR_DNS>** so that this driver can connect with FDTD.


## **real-time time-dependent density functional theory (rttddft)** 

A sample bash input for setting up the **rttddft** model is as follows: 

```bash
mxl_driver.py --model rttddft --port 31415 --param "molecule_xyz=PATH_TO_MaxwellLink/tests/data/hcn.xyz, checkpoint=false, restart=false" --verbose
```

Similar as the **tls**, here the **--param** input provides the initialization of the Python class of the **rttddft** system: [models/rttddft_model.py](./models/rttddft_model.py). Currently, the Psi4 electronic structure engine is used to provide necessary quantities for propagating RT-TDDFT. 

We can freely change the molecular geometry, level of functional, basis, time step in RT-TDDFT, required memory, and number of CPU threads:
```bash
--param "molecule_xyz=PATH_TO_MaxwellLink/tests/data/hcn.xyz, functional="SCF",basis="sto-3g", dt_rttddft_au=0.04, memory="2GB",num_threads=1, checkpoint=false, restart=false"
```

A few tips when running RT-TDDFT with MaxwellLink:

1. The second line in the xyz file should be "charge multiplicity" (such as **0 1**). 

2. The time step in RT-TDDFT **dt_rttddft_au** (in atomic units) should be the same or smaller than the FDTD time step. If the RT-TDDFT time step is smaller than the FDTD time step, **dt_rttddft_au** will be adjusted automatically so an integer number of RT-TDDFT steps are propagated per FDTD step.

3. If **checkpoint** is **true**, after each FDTD step (or a few RT-TDDFT steps), the necessary checkpoint data in RT-TDDFT will be written to disk.

4. If **restart=true** and **checkpoint=true**, the driver code can restart from the checkpoint file in disk and resume the simulation. This setting is necessary when many  drivers are connected to MaxwellLink at the same time. If one driver is terminated in one machine (by different reasons), all the other drivers will pause the simulation and wait for the restart of this driver. 

5. If the driver and the FDTD code are running in different HPC nodes, please also set **--address <FDTD_CODE_IP_OR_DNS>** so that this driver can connect with FDTD.



## **Write your own driver** 

In [models/](./models/), each driver class inherits the **DummyModel** class. To implement a new Python driver, please add a new class file in [models/](./models/) and provide the implemenation of the class methods in **DummyModel** class.