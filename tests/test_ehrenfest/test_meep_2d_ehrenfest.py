#--------------------------------------------------------------------------------------#
# Copyright (c) 2026 MaxwellLink                                                       #
# This file is part of MaxwellLink. Repository: https://github.com/TaoELi/MaxwellLink  #
# If you use this code, always credit and cite arXiv:2512.06173.                       #
# See AGENTS.md and README.md for details.                                             #
#--------------------------------------------------------------------------------------#


import numpy as np
import maxwelllink as mxl
from maxwelllink import sockets as mxs
import meep as mp
from eh_hamiltonian import *    # import 1D model Hamiltonian

def test_meep_2d_ehrenfest(plotting=False):

    molecule = mxl.Molecule(
        center=mp.Vector3(0, 0, 0),
        size=mp.Vector3(1, 1, 1),
        sigma=0.1,
        dimensions=2,
        rescaling_factor=1,
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
            "verbose": False,
    },
)


    sim = mxl.MeepSimulation(
        molecules=[molecule],
        cell_size=mp.Vector3(8, 8, 0),
        boundary_layers=[mp.PML(3.0)],
        resolution=10,
        # fix a units system
        time_units_fs=0.1,
    )


    total_steps = 10000
    sim.run(
        steps=total_steps,
        )
    
    if plotting:
        import matplotlib.pyplot as plt
        time_au = np.array([entry["time_au"] for entry in molecule.additional_data_history])
        R_au = np.array([entry["R_au"] for entry in molecule.additional_data_history])
        P_au = np.array([entry["P_au"] for entry in molecule.additional_data_history])

        fig, ax1 = plt.subplots(figsize=(6.00, 2.36))
        ax1.plot(time_au, R_au, color='blue', label="position") 
        ax1.set_xlabel(r"$t$ [a.u.]")
        ax1.set_ylabel(r"$R$ [a.u.]")
        ax2 = ax1.twinx()
        ax2.plot(time_au, P_au, color='green', label="momentum") 
        ax2.set_ylabel(r"$P$ [a.u.]")
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2)
        fig.tight_layout()
        plt.savefig("trajectory.png", dpi=400)
        plt.show()
if __name__ == "__main__":
    test_meep_2d_ehrenfest(plotting=True)
