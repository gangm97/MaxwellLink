# MaxwellLink
A flexible framework for self-consistent EM-molecular simulations via a socket protocol

## Install from source
Assuming Anaconda is installed, MaxwellLink can be installed from source to access the latest features:
```bash
# create a new conda environment
CONDA_ENV="mxl_build"
# install MEEP for the FDTD engine (the only supported FDTD engine for now)
conda create -n $CONDA_ENV -c conda-forge pymeep

# install MaxwellLink
git clone git@github.com:TaoELi/MaxwellLink.git
cd PATH_TO_MAXWELLLINK/
pip install .

# [optional] install Psi4 quantum chemistry code for RT-TDDFT driver
conda install conda-forge::psi4

# [optional] install modified LAMMPS code (with fix mxl support) for classical MD driver
# the command below works after installing MaxwellLink with pip install .
mxl_install_lammps
# If the above command fails, please try the bash code below instead for installing lammps [working for Linux and MacOS]
# bash ./src/maxwelllink/mxl_drivers/lammps/mxl_install_lammps.sh
```

## Uninstall 
```bash
pip uninstall maxwelllink
```

## Test
```bash
pytest -v
```
