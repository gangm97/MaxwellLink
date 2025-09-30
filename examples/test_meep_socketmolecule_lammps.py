import meep as mp
import maxwelllink as mxl

cell = mp.Vector3(8, 8, 0)
geometry = []
sources_non_molecule = []
pml_layers = [mp.PML(3.0)]
resolution = 10

# define a TLS molecule
dipole_moment = 1e-1
frequency = 1.0

hub = mxl.SocketHub(
    host="localhost", port=31889, unixsocket="lmp_test", timeout=10.0, latency=1e-5
)

molecule1 = mxl.SocketMolecule(
    hub=hub,
    molecule_id=0,
    resolution=resolution,
    center=mp.Vector3(0, 0, 0),
    size=mp.Vector3(1, 1, 1),
    sigma=0.1,
    dimensions=2,
    time_units_fs=10,
)

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources_non_molecule,
    boundary_layers=pml_layers,
    resolution=resolution,
)

sim.run(
    mxl.update_molecules(hub=hub, sources_non_molecule=[], molecules=[molecule1]),
    until=90,
)

hub.stop()

"""
# plot the TLS population relaxation dynamics
dmux_dt = np.array([additional_data["dmudt_au"][0] for additional_data in molecule1.additional_data_history])
time_fs = np.array([np.real(additional_data["t_fs"]) for additional_data in molecule1.additional_data_history])

# if we are running this script directly, plot the results
if __name__ == "__main__" and mp.am_master():
    plt.figure()
    plt.plot(time_fs, dmux_dt, label="Numerical")
    plt.xlabel("Time [fs]")
    plt.ylabel("dmux/dt")
    plt.title("TLS Population Relaxation Dynamics")
    plt.legend()
    plt.grid()
    plt.show()
"""
