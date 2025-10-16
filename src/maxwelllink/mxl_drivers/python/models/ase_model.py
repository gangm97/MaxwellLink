import os
import numpy as np
from typing import Optional, Sequence, Union, Dict
import ast

try:
    from .dummy_model import DummyModel
except:
    from dummy_model import DummyModel

try:
    from ase import units
except ImportError as e:
    raise ImportError(
        "ASE package is required for ASEModel. Please install it via 'conda install conda_forge::ase'."
    ) from e

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

# from ase.md.langevin import Langevin
from ase.io import read as ase_read

# ----- Unit constants -----
# 1 a.u. of time = 0.02418884326505 fs  => fs_to_au = 1 fs in a.u.
FS_TO_AU = 41.341373335  # you already use this in your code
# 1 Å = BOHR_PER_ANG * a0
BOHR_PER_ANG = 1.889726124565062
# E(a.u.) = 5.142206747e11 V/m; F(eV/Å) for q=1 is E(V/m) * 1e-10
# So F(eV/Å) = q * E(a.u.) * 51.422067476  (≈ 5.1422e11 * 1e-10)
FORCE_PER_EFIELD_AU_EV_PER_ANG = 51.422067476


def _parse_kwargs_string(s: str) -> Dict:
    """
    Parse a compact 'k1=v1,k2=v2' into dict with numbers/bools auto-cast.
    (Used for preset params and for passing kwargs to user spec.)

    + **`s`** (str): Input string.
    """
    if not s:
        return {}

    def _strip_quotes(s):
        if not isinstance(s, str):
            return s
        s = s.strip()
        if (len(s) >= 2) and ((s[0] == s[-1]) and s[0] in ("'", '"')):
            return s[1:-1]
        return s

    s = _strip_quotes(s)
    out = {}
    for token in s.split(","):
        if not token.strip():
            continue
        if "=" not in token:
            # allow bare flags as True
            out[token.strip()] = True
            continue
        k, v = token.split("=", 1)
        k = k.strip()
        v = v.strip()
        # try literal (int/float/bool/list/tuple) -> else string
        try:
            out[k] = ast.literal_eval(v)
        except Exception:
            # also allow "true"/"false"
            lv = v.lower()
            if lv in ("true", "false"):
                out[k] = lv == "true"
            else:
                out[k] = v
    return out


def _build_calculator(name: str, **kwargs) -> Calculator:
    """
    Minimal factory for common ASE calculators.

    Examples:
    - name='psi4' -> Psi4 via ase.calculators.psi4 (Psi4 binary required)
    - name='dftb'  -> pip install dftbplus (and set DFTB+ binary)
    - name='orca'  -> ORCA via ase.calculators.orca (ORCA binary required)

    + **`name`** (str): Name of the calculator.
    + **`kwargs`**: Additional kwargs passed to the calculator constructor.
    """
    n = (name or "").strip().lower()

    if n == "psi4":
        from ase.calculators.psi4 import Psi4

        return Psi4(**kwargs)

    if n == "orca":
        from ase.calculators.orca import ORCA

        return ORCA(**kwargs)

    if n == "dftb":
        from ase.calculators.dftb import Dftb

        return Dftb(**kwargs)

    # Fallback: try to import path "ase.calculators.<name>"
    try:
        mod = __import__(f"ase.calculators.{n}", fromlist=["*"])
        # look for a Calculator subclass with same-cased name, else first subclass
        for attr in dir(mod):
            cls = getattr(mod, attr)
            try:
                if issubclass(cls, Calculator) and cls is not Calculator:
                    return cls(**kwargs)
            except Exception:
                pass
    except Exception as e:
        raise ImportError(
            f"Unknown calculator '{name}'. Install or extend _build_calculator()."
        ) from e

    raise ImportError(
        f"Failed to construct calculator for '{name}'. Extend _build_calculator()."
    )


class ForceAugmenter(Calculator):
    """
    ASE Calculator wrapper that adds an external uniform E-field force F_i^ext = q_i E
    to each atom i, where q_i are per-atom charges (either fixed or recomputed each step).
    """

    implemented_properties = ("energy", "forces")

    def __init__(self, base, charges=None, recompute_charges=False, verbose=False):
        """
        + **`base`** (Calculator): An ASE Calculator instance to wrap.
        + **`charges`** (array-like or None): Per-atom charges in |e| units. If None, set recompute_charges=True.
        + **`recompute_charges`** (bool): If True, query charges each step (e.g. Mulliken). Default is False.
        + **`verbose`** (bool): Whether to print verbose output. Default is False.
        """
        super().__init__()
        self.base = base
        self.charges = None if charges is None else np.asarray(charges, float).copy()
        self.recompute_charges = bool(recompute_charges)
        self.verbose = bool(verbose)
        self._E_au = np.zeros(3, float)

        # small caches for saving the computed results given a molecular geometry
        self._cache_key = None
        self._cache_energy = None
        self._cache_forces = None

    def set_field_au(self, Evec3_au):
        """
        Set the external uniform E-field vector in atomic units (a.u.).

        + **`Evec3_au`** (array-like): 3-element array-like representing the E-field vector in a.u.
        """
        self._E_au = np.asarray(Evec3_au, float).reshape(3)

    def _geom_key(self, atoms: Atoms):
        """
        Build a key for the current molecular geometry.

        + **`atoms`** (Atoms): An ASE Atoms object.

        Returns:
        - A hashable key representing the current geometry.
        """
        # Positions + cell + numbers
        pos = atoms.get_positions()
        cell = atoms.cell.array
        Z = atoms.get_atomic_numbers()
        return (
            pos.shape,
            pos.tobytes(),
            cell.shape,
            cell.tobytes(),
            Z.shape,
            Z.tobytes(),
        )

    def calculation_required(self, atoms, properties):
        """
        Determine whether a recalculation is required based on changes in the atomic configuration.

        + **`atoms`** (Atoms): An ASE Atoms object.
        + **`properties`** (tuple): A tuple of properties to be calculated (e.g., 'energy', 'forces').

        Returns:
        - A boolean indicating whether a recalculation is required.
        """
        # If we have a cached forces/energy for the current molecular geometry, no recalc needed.
        key_now = self._geom_key(atoms)
        if self._cache_key is not None and key_now == self._cache_key:
            # We can satisfy any subset of ('energy','forces') from cache
            have_forces = self._cache_forces is not None
            have_energy = self._cache_energy is not None
            need_forces = "forces" in properties
            need_energy = "energy" in properties
            if (not need_forces or have_forces) and (not need_energy or have_energy):
                return False
        # Otherwise, defer to the base (positions/cell/numbers changes) -> recalc
        if hasattr(self.base, "calculation_required"):
            return self.base.calculation_required(atoms, properties)
        return super().calculation_required(atoms, properties)

    def calculate_external_force(self, atoms):
        """
        Calculate the external force on each atom due to the uniform E-field.

        + **`atoms`** (Atoms): An ASE Atoms object.

        Returns:
        - A numpy array of shape (N, 3) representing the external forces on each atom in eV/Angstrom (internal units of ASE forces).
        """
        if True:
            # Resolve charges
            if self.recompute_charges:
                q = None
                getq = getattr(self.base, "get_charges", None)
                if callable(getq):
                    try:
                        # print("[ForceAugmenter] calculate() charges called")
                        q = np.asarray(getq(atoms), float)
                        self.charges = q
                    except Exception:
                        q = None
                if q is None:
                    if self.charges is None:
                        raise RuntimeError(
                            "recompute_charges=True but base has no get_charges(); pass fixed 'charges=' instead."
                        )
                    q = self.charges
            else:
                if self.charges is None:
                    raise RuntimeError(
                        "Need 'charges=' or recompute_charges=True with a supporting calculator."
                    )
                q = self.charges

            # Add uniform-field force in eV/Angstrom
            Fext = q.reshape(-1, 1) * (self._E_au * FORCE_PER_EFIELD_AU_EV_PER_ANG)
            return Fext

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        """
        Calculate the requested properties for the given atomic configuration.

        + **`atoms`** (Atoms): An ASE Atoms object.
        + **`properties`** (tuple): A tuple of properties to be calculated (e.g., 'energy', 'forces').
        + **`system_changes`** (list): A list of changes in the atomic configuration.
        """
        key_now = self._geom_key(atoms)

        # If cache matches and covers requested properties, serve from cache and return
        if self._cache_key is not None and key_now == self._cache_key:
            if "energy" in properties and self._cache_energy is not None:
                self.results["energy"] = self._cache_energy
            if "forces" in properties and self._cache_forces is not None:
                self.results["forces"] = (
                    self._cache_forces + self.calculate_external_force(atoms)
                )
            if (("energy" not in properties) or (self._cache_energy is not None)) and (
                ("forces" not in properties) or (self._cache_forces is not None)
            ):
                return

        # Ask base ONLY for what was requested
        props_for_base = tuple(
            p
            for p in properties
            if p in getattr(self.base, "implemented_properties", ())
        )
        if not props_for_base:
            # Fallback: many base calculators accept energy-only
            props_for_base = ("energy",) if "energy" in properties else properties

        # Clear results for this call
        self.results.clear()

        # Compute on base
        self.base.calculate(atoms, props_for_base, system_changes)

        # Copy energy if requested
        if "energy" in properties and "energy" in self.base.results:
            self.results["energy"] = float(self.base.results["energy"])
        elif "energy" in properties:
            # Not all calculators compute energy on a forces-only request
            self.results["energy"] = None

        # Forces: only if requested now
        if "forces" in properties:
            f_base = self.base.results.get("forces", None)
            if f_base is None:
                raise RuntimeError(
                    "Base calculator did not provide forces when requested."
                )
            f = np.array(f_base, dtype=float, copy=True)

            # Update caches
            self._cache_key = key_now
            self._cache_energy = self.results.get("energy", None)
            self._cache_forces = f.copy()

            Fext = self.calculate_external_force(atoms)
            f += Fext
            self.results["forces"] = f


class ASEModel(DummyModel):
    r"""
    General BOMD (Born-Oppenheimer MD) driver using ASE.

    This model provides the two key functionalities like other Models supported by MaxwellLink:

    1. **MD coupled to E-field**: Injects F_i^ext = q_i E (uniform field) to the molecular forces
    in MD simulations, where q_i are per-atom charges (constant user-supplied or calculator-reported each step).

    2. **Return source amplitude**: $\dot{\mu} = \sum_i q_i v_i$ (converted to atomic units)
    """

    def __init__(
        self,
        atoms: Union[str, Atoms],
        calculator: str = "xtb",
        calc_kwargs: Optional[dict] = None,
        charges: Optional[Sequence[float]] = None,
        recompute_charges: bool = False,
        n_substeps: int = 1,
        temperature_K: float = 0.0,
        verbose: bool = False,
        checkpoint: bool = False,
        restart: bool = False,
        **extra,
    ):
        """
        + **`atoms`**: either an ASE Atoms object or a path to structure file (e.g., xyz) readable by ASE.
        + **`calculator`**: name of ASE calculator ('psi4', 'dftb', 'orca', ...)
        + **`calc_kwargs`**: dict of kwargs passed to calculator constructor
        + **`charges`**: a string "[-1.0 1.0]" representing an array of per-atom charges (in |e|), separated by **space** (not **comma**).
        If None, need to set recompute_charges=True
        + **`recompute_charges`**: if True, query charges each step (e.g., Mulliken). Default is False.
        + **`n_substeps`**: number of MD sub-steps per MEEP step. Default is 1.
        + **`temperature_K`**: initial temperature in Kelvin for Maxwell-Boltzmann distribution. Default is 0.0 K.
        + **`verbose`**: whether to print verbose output. Default is False.
        + **`checkpoint`**: whether to enable checkpointing. Default is False.
        + **`restart`**: whether to restart from a checkpoint if available. Default is False.
        """
        super().__init__(verbose=verbose, checkpoint=checkpoint, restart=restart)

        # atoms
        if isinstance(atoms, Atoms):
            self.atoms = atoms.copy()
        else:
            # treat as file path
            self.atoms = ase_read(str(atoms))

        self.calc_name = calculator
        # the input for calc_kwargs is as follows: --params 'xx=xx, calc_kwargs=k1=v1,k2=v2, yy=yy'
        # recover the first hit (k1=v1)
        self.calc_kwargs = _parse_kwargs_string(calc_kwargs)
        # all the other extra kwargs go into calc_kwargs (such as k2=v2)
        _own_keys = {
            "atoms",
            "calculator",
            "calc_kwargs",
            "charges",
            "recompute_charges",
            "n_substeps",
            "temperature_K",
            "verbose",
            "checkpoint",
            "restart",
        }
        for k in list(extra.keys()):
            if k not in _own_keys:
                self.calc_kwargs[k] = extra.pop(k)

        print("[ASEModel] calculator name =", self.calc_name)
        print("[ASEModel] calculator kwargs =", self.calc_kwargs)

        if charges is None:
            self.user_charges = None
        else:
            try:
                self.user_charges = np.fromstring(charges.strip("[]"), sep=" ")
            except Exception as e:
                raise ValueError(
                    "Failed to parse 'charges' string into array; use format like '[0.1 -0.2 0.0 ...]'"
                ) from e
        self.recompute_charges = bool(recompute_charges)
        self.n_substeps = max(int(n_substeps), 1)
        self.temperature_K = temperature_K

        print("[ASEModel] user_charges =", self.user_charges)
        print("[ASEModel] recompute_charges =", self.recompute_charges)
        print("[ASEModel] n_substeps =", self.n_substeps)
        print("[ASEModel] temperature_K =", self.temperature_K)
        if self.user_charges is None and not self.recompute_charges:
            raise RuntimeError(
                "ASEModel needs charges: pass 'charges=' or set 'recompute_charges=True' with calculator support."
            )

        self.integrator = None
        self.forcewrap = None
        self._last_amp = np.zeros(3)

        # cached for dmu/dt
        self._charges = None
        self._vel_angs_per_fs = None

    # -------------- heavy-load initialization (at INIT) --------------

    def initialize(self, dt_new, molecule_id):
        """
        Initialize the model with the new time step and molecule ID.

        + **`dt_new`** (float): The new time step in atomic units (a.u.).
        + **`molecule_id`** (int): The ID of the molecule.
        """
        self.dt = float(dt_new)
        self.molecule_id = int(molecule_id)

        self.checkpoint_filename = f"ase_checkpoint_id_{self.molecule_id}.npz"

        # MD step in fs
        dt_fs = (self.dt / FS_TO_AU) / self.n_substeps
        if dt_fs <= 0.0:
            raise ValueError("Non-positive dt_fs computed.")

        # Base calculator and force wrapper
        base_calc = _build_calculator(self.calc_name, **self.calc_kwargs)
        self.forcewrap = ForceAugmenter(
            base=base_calc,
            charges=self.user_charges,
            recompute_charges=self.recompute_charges,
            verbose=self.verbose,
        )
        self.atoms.calc = self.forcewrap

        # Initialize velocities
        if self.temperature_K >= 0.0:
            MaxwellBoltzmannDistribution(
                self.atoms, temperature_K=float(self.temperature_K)
            )

        if self.checkpoint and self.restart:
            self._reset_from_checkpoint()
            self.restarted = True

        # Choose the VelocityVerlet integrator
        self.integrator = VelocityVerlet(
            self.atoms, timestep=dt_fs * units.fs, logfile=None, loginterval=0
        )

        if self.verbose:
            print(
                f"[ASEModel {self.molecule_id}] dt={self.dt:.6e} a.u. "
                f"-> dt_fs={dt_fs:.6e} fs; substeps={self.n_substeps}; "
                f"calculator={self.calc_name}({self.calc_kwargs})"
            )

    # -------------- one FDTD step under E-field --------------

    def propagate(self, effective_efield_vec):
        """
        Propagate the BO molecular dynamics given the effective electric field vector.

        + **`effective_efield_vec`**: Effective electric field vector in the form [Ex, Ey, Ez].
        """
        # 1. set the field for the wrapper (a.u.)
        self.E_vec = np.asarray(effective_efield_vec, float).reshape(3)

        if self.verbose:
            print(
                f"[ASEModel {self.molecule_id}] t={self.t:.6f} a.u., E={self.E_vec[2]:.10e} a.u."
            )

        self.forcewrap.set_field_au(self.E_vec)

        # 2. do n_substeps
        self.integrator.run(self.n_substeps)

        # 3. cache per-atom charges & velocities at the end of the step
        # charges (either recomputed or fixed)
        if self.recompute_charges:
            self._charges = self.forcewrap.charges
        else:
            self._charges = self.user_charges

        # velocities (Angstrom/fs)
        vel = self.atoms.get_velocities()
        if vel is None:
            vel = np.zeros((len(self.atoms), 3), float)
        self._vel_angs_per_fs = np.asarray(vel, float)

        if self.verbose:
            print(
                f"[ASEModel {self.molecule_id}] t={self.t:.6e} au, E_au={self.forcewrap._E_au[2]:.10e} a.u.,"
                f"q={self._charges}, v_angs_per_fs={self._vel_angs_per_fs}"
            )

        # advance model time in a.u.
        self.t += self.dt

    def calc_amp_vector(self):
        """
        Return the amplitude vector (dmu/dt) for the current time step in atomic units.

        In classical MD: dmu/dt = sum_i q_i v_i
        """
        if self._vel_angs_per_fs is None or self._charges is None:
            return np.zeros(3, float)

        v_au = self._vel_angs_per_fs * (BOHR_PER_ANG / FS_TO_AU)
        amp = (self._charges.reshape(-1, 1) * v_au).sum(axis=0)

        return amp

    # ------------ optional operation / checkpoint --------------

    def append_additional_data(self):
        """
        Append additional data to be sent back to MaxwellLink, which can be retrieved by the user
        via the Python interface: maxwelllink.SocketMolecule.additional_data_history, where
        additional_data_history is a list of dictionaries.

        Returns:
        - A dictionary containing additional data.
        """
        d = {
            "time_au": float(self.t),
            "temperature_K": float(self.atoms.get_temperature()),
        }
        return d

    def _dump_to_checkpoint(self):
        """
        Dump the internal state of the model to a checkpoint.
        """
        np.savez(
            self.checkpoint_filename,
            t=self.t,
            positions=self.atoms.get_positions(),
            velocities=self.atoms.get_velocities(),
        )

    def _reset_from_checkpoint(self):
        """
        Reset the internal state of the model from a checkpoint.
        """
        if not os.path.exists(self.checkpoint_filename):
            # No checkpoint file found means this driver has not been paused or terminated abnormally
            # so we just start fresh.
            if self.verbose:
                print(
                    f"[ASEModel] No checkpoint file found for id={self.molecule_id}, starting fresh."
                )
        else:
            data = np.load(self.checkpoint_filename)
            self.t = float(data["t"])
            pos = np.asarray(data["positions"], float)
            vel = np.asarray(data["velocities"], float)
            self.atoms.set_positions(pos)
            self.atoms.set_velocities(vel)
            if self.verbose:
                print(f"[ASEModel] Restarted from checkpoint for id={self.molecule_id}")

    def _snapshot(self):
        """
        Return a snapshot of the internal state for propagation. Deep copy the arrays to avoid mutation issues.
        """
        return {
            "time": self.t,
            "positions": self.atoms.get_positions().copy(),
            "velocities": self.atoms.get_velocities().copy(),
        }

    def _restore(self, snapshot):
        """
        Restore the internal state from a snapshot.
        """
        self.t = snapshot["time"]
        self.atoms.set_positions(snapshot["positions"])
        self.atoms.set_velocities(snapshot["velocities"])
