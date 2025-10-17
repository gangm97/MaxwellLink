# molecule_abstract.py
"""
EM-agnostic molecular classes for coupling quantum molecular dynamics with EM solvers.

This module intentionally avoids importing any specific EM solver (e.g., MEEP).
Backends (e.g., em_solver/meep.py) provide the EM-specific pieces:
  * building/updating sources for injection
  * computing regularized E-field integrals
  * step functions for the solver main loop
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Dict, Optional

from .sockets import SocketHub, am_master, mpi_bcast_from_master
from .units import FS_TO_AU
from .mxl_drivers.python.models import __drivers__


# ---------------------------------------------------------------------
# Minimal Vector3 so we don't depend on mp.Vector3 in the abstract layer
# ---------------------------------------------------------------------
@dataclass
class Vector3:
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0


# ---------------------------------------------------------------------
# A novel molecule class which can support socket and non-socket modes
# This is a future extension and is not yet integrated with the rest of the codebase.
# ---------------------------------------------------------------------


class Molecule:
    """
    A molecule class which can support both socket and non-socket modes.
    EM-specific source creation and field-integral computation are provided by backends.
    """

    def __init__(
        self,
        hub: Optional[SocketHub] = None,
        driver: Optional[str] = None,
        center: Optional[Vector3] = None,
        size: Optional[Vector3] = None,
        dimensions: Optional[int] = None,
        sigma: Optional[float] = None,
        resolution: Optional[int] = None,
        init_payload: Optional[Dict] = None,
        driver_kwargs: Optional[Dict] = None,
        rescaling_factor: float = 1.0,
    ):
        self.hub = hub
        self.driver = driver
        self.center = center
        self.size = size
        self.dimensions = dimensions
        self.sigma = sigma
        # optional setting resolution for building real-space sources in FDTD code
        self.resolution = resolution
        # optional setting init_payload for socket communication with molecular drivers
        self.init_payload = init_payload if init_payload is not None else {}
        self.rescaling_factor = rescaling_factor

        self.molecule_id = -1  # Default ID for non-socket mode
        self.time_units_fs = 0.0  # Optional time unit when connecting to FDTD code
        self.init_payload["dt_au"] = 0.0  # Placeholder for dt_au in init_payload

        # if resolution is provided, we also compute dx and dt
        self.dx = 1.0 / resolution if resolution is not None else 0.0
        self.dt = 0.5 / resolution if resolution is not None else 0.0

        # reserve for sources and additional data history
        self.sources = []
        self.additional_data_history = []

        if self.dimensions not in (1, 2, 3) and self.dimensions is not None:
            raise ValueError("Molecule only supports 1D, 2D and 3D simulations.")

        # identify the spatial polarization kernel function
        self.polarization_fingerprint = {
            "dimensions": self.dimensions,
            "sigma": self.sigma,
            "size": [self.size.x, self.size.y, self.size.z],
            "rescaling_factor": self.rescaling_factor,
            "center": [self.center.x, self.center.y, self.center.z],
        }
        self.polarization_fingerprint_hash = hash(
            json.dumps(self.polarization_fingerprint)
        )

        self.mode = "none"
        if hub is not None and driver is None:
            self.mode = "socket"
        elif hub is None and driver is not None:
            self.mode = "non-socket"
        else:
            raise ValueError("Either hub or driver must be provided, but not both.")

        if self.mode == "socket":
            # register with the hub on rank 0
            if am_master():
                self.molecule_id = self.hub.register_molecule_return_id()
                print(
                    f"[Init Molecule] Under socket mode, registered molecule with ID {self.molecule_id}"
                )
            # if using mpi, we also need to broadcast the molecule_id to other ranks
            self.molecule_id = mpi_bcast_from_master(self.molecule_id)
        elif self.mode == "non-socket":
            self.driver = str(driver).lower()
            if self.driver not in list(__drivers__.keys()):
                raise ValueError(
                    f"[Init Molecule] Unsupported driver: {self.driver}, only supports {list(__drivers__.keys())}"
                )
            print(
                f"[Init Molecule] Operating in non-socket mode, using driver: {self.driver}"
            )
            # now we initialize the driver
            try:
                self.d_f = __drivers__[self.driver](**driver_kwargs)
            except ImportError:
                # specific errors have already been triggered
                raise
            except Exception as err:
                print(
                    f"Error setting up molecular dynamics model {self.driver} with args {driver_kwargs}"
                )
                print(__drivers__[self.driver].__doc__)
                print("Error trace: ")
                raise err

    def _refresh_time_units(self, time_units_fs):
        self.time_units_fs = time_units_fs
        self.dt_au = self.dt * time_units_fs * FS_TO_AU
        # also refresh in init_payload for drivers
        self.init_payload["dt_au"] = self.dt_au

    def _refresh_time_step(self, dt_em):
        self.dt = dt_em
        self.dt_au = self.dt * self.time_units_fs * FS_TO_AU
        # also refresh in init_payload for drivers
        self.init_payload["dt_au"] = self.dt_au

    def initialize_driver(self, dt_au, molecule_id):
        if self.mode == "socket":
            raise NotImplementedError(
                "Socket mode driver initialization is not implemented yet."
            )
        elif self.mode == "non-socket":
            # use the driver to initialize
            self.dt_au = float(dt_au)
            self.molecule_id = int(molecule_id)
            self.d_f.initialize(self.dt_au, self.molecule_id)
        else:
            raise RuntimeError(
                "Molecule is not properly initialized in socket or non-socket mode."
            )

    def propagate(self, efield_vec3):
        if self.mode == "socket":
            raise NotImplementedError("Socket mode propagation is not implemented yet.")
        elif self.mode == "non-socket":
            # use the driver to propagate
            self.d_f.propagate(efield_vec3)
        else:
            raise RuntimeError(
                "Molecule is not properly initialized in socket or non-socket mode."
            )

    def calc_amp_vector(self):
        if self.mode == "socket":
            raise NotImplementedError(
                "Socket mode amplitude calculation is not implemented yet."
            )
        elif self.mode == "non-socket":
            # use the driver to calculate amplitude vector
            return self.d_f.calc_amp_vector()
        else:
            raise RuntimeError(
                "Molecule is not properly initialized in socket or non-socket mode."
            )
