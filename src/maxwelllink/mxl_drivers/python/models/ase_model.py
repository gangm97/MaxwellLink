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


# ---------- small helper: parse k=v,k2=v2 string into dict ----------
def _strip_quotes(s):
    if not isinstance(s, str):
        return s
    s = s.strip()
    if (len(s) >= 2) and ((s[0] == s[-1]) and s[0] in ("'", '"')):
        return s[1:-1]
    return s


def _parse_kwargs_string(s: str) -> Dict:
    """
    Parse a compact 'k1=v1,k2=v2' into dict with numbers/bools auto-cast.
    (Used for preset params and for passing kwargs to user spec.)

    + **`s`** (str): Input string.
    """
    if not s:
        return {}
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


# ---------- small helper: build a calculator by name ----------
def _build_calculator(name: str, **kwargs) -> Calculator:
    """
    Minimal factory for common ASE calculators.
    Examples:
    - name='dftb'  -> pip install dftbplus (and set DFTB+ binary)
    - name='orca'  -> ORCA via ase.calculators.orca (ORCA binary required)
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


# ---------- wrapper to add qE force each step ----------
class ForceAugmenter(Calculator):
    implemented_properties = ("energy", "forces")

    def __init__(self, base, charges=None, recompute_charges=False, verbose=False):
        super().__init__()
        self.base = base
        self.charges = None if charges is None else np.asarray(charges, float).copy()
        self.recompute_charges = bool(recompute_charges)
        self.verbose = bool(verbose)
        self._E_au = np.zeros(3, float)

        # --- small cache for a single geometry ---
        self._cache_key = None
        self._cache_energy = None
        self._cache_forces = None

    def set_field_au(self, Evec3_au):
        self._E_au = np.asarray(Evec3_au, float).reshape(3)

    def _geom_key(self, atoms: Atoms):
        """Build a cheap, robust key for current 'state' we care about."""
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
        # If we have a cached forces/energy for the current geometry+field, no recalc needed.
        key_now = self._geom_key(atoms)
        if self._cache_key is not None and key_now == self._cache_key:
            # We can satisfy any subset of ('energy','forces') from cache
            have_forces = self._cache_forces is not None
            have_energy = self._cache_energy is not None
            need_forces = "forces" in properties
            need_energy = "energy" in properties
            if (not need_forces or have_forces) and (not need_energy or have_energy):
                return False  # cache is good
        # Otherwise, defer to the base (positions/cell/numbers changes) -> recalc
        if hasattr(self.base, "calculation_required"):
            return self.base.calculation_required(atoms, properties)
        return super().calculation_required(atoms, properties)

    def calculate_external_force(self, atoms):
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

            # Add uniform-field force in eV/Å
            Fext = q.reshape(-1, 1) * (self._E_au * FORCE_PER_EFIELD_AU_EV_PER_ANG)
            return Fext

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        # print("[ForceAugmenter] calculate() called")

        key_now = self._geom_key(atoms)

        # If cache matches and covers requested properties, serve from cache and return
        if self._cache_key is not None and key_now == self._cache_key:
            if "energy" in properties and self._cache_energy is not None:
                self.results["energy"] = self._cache_energy
            if "forces" in properties and self._cache_forces is not None:
                self.results["forces"] = (
                    self._cache_forces + self.calculate_external_force(atoms)
                )
                # print("[ForceAugmenter] served forces from cache")
            if (("energy" not in properties) or (self._cache_energy is not None)) and (
                ("forces" not in properties) or (self._cache_forces is not None)
            ):
                return

        # Ask base ONLY for what was requested (don’t trigger extra work)
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
            # print("[ForceAugmenter] calculate() energy called")
        elif "energy" in properties:
            # Not all calculators compute energy on a forces-only request
            self.results["energy"] = None

        # Forces: only if requested now
        if "forces" in properties:
            # print("[ForceAugmenter] calculate() forces called")
            f_base = self.base.results.get("forces", None)
            if f_base is None:
                raise RuntimeError(
                    "Base calculator did not provide forces when requested."
                )
            f = np.array(f_base, dtype=float, copy=True)

            # Update cache (store what we actually have now)
            self._cache_key = key_now
            self._cache_energy = self.results.get("energy", None)
            self._cache_forces = f.copy()  # cache base forces

            # print("[ForceAugmenter] base forces =", f)

            Fext = self.calculate_external_force(atoms)
            f += Fext
            self.results["forces"] = f


class ASEModel(DummyModel):
    """
    General BOMD (Born–Oppenheimer MD) driver using ASE.

    Coupling to E-field:
        - injects F_i^ext = q_i E (uniform field), where q_i are per-atom charges
          (constant user-supplied or calculator-reported each step).

    Returned source amplitude:
        - \dot{mu} = sum_i q_i v_i  (converted to atomic units)

    Key options:
        atoms: either an ASE Atoms object or a path to structure file readable by ASE.
        calculator: name of ASE calculator ('xtb', 'dftb', 'gpaw', 'orca', ...)
        calc_kwargs: dict of kwargs passed to calculator constructor
        charges: array of per-atom charges (in e). If None, set recompute_charges=True
                 and use a calculator that exposes get_charges(atoms).
        recompute_charges: if True, query charges each step (e.g. Mulliken).
        n_substeps: number of MD substeps per MEEP tick.
        dt_scale: optional scale on dt (e.g. do smaller MD steps than dt_au from MEEP).
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
        self.calc_kwargs = _parse_kwargs_string(_strip_quotes(calc_kwargs))
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

        self.integrator = None  # ASE MD object
        self.forcewrap = None  # ForceAugmenter
        self._last_amp = np.zeros(3)  # last \dot{mu}

        # cached for dmu/dt
        self._charges = None  # current charges (e)
        self._vel_angs_per_fs = None  # (N,3) Å/fs

    # -------------- initialization after we know dt_au, molid --------------
    def initialize(self, dt_new, molecule_id):
        self.dt = float(dt_new)
        self.molecule_id = int(molecule_id)

        # MD step in fs
        dt_fs = (self.dt / FS_TO_AU) / self.n_substeps
        if dt_fs <= 0.0:
            raise ValueError("Non-positive dt_fs computed; check dt_scale.")

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

    # -------------- one MEEP step under E-field --------------
    def propagate(self, effective_efield_vec):
        """
        One MEEP tick: do n_substeps of MD under constant E (uniform).
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

        # velocities (Å/fs); ensure present
        vel = self.atoms.get_velocities()  # shape (N,3)
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
        \dot{mu} = sum_i q_i v_i
        Convert v from Å/fs to a0/a.u.: v_a0_per_au = v(Å/fs) * (BOHR_PER_ANG / FS_TO_AU)
        Return vector in atomic units.
        """
        if self._vel_angs_per_fs is None or self._charges is None:
            return np.zeros(3, float)

        v_au = self._vel_angs_per_fs * (BOHR_PER_ANG / FS_TO_AU)  # (N,3)
        amp = (self._charges.reshape(-1, 1) * v_au).sum(axis=0)  # (3,)

        """
        if self.verbose:
            ke = self.atoms.get_kinetic_energy()  # eV
            try:
                pe = self.atoms.get_potential_energy()  # eV
            except Exception:
                pe = np.nan
            print(
                f"[ASEModel {self.molecule_id}] t={self.t:.6e} au, E={self.forcewrap._E_au} au, "
                f"amp={amp}, E_kin={ke:.3f} eV, E_pot={pe:.3f} eV"
            )
        """
        return amp

    # -------------- optional extras --------------
    def append_additional_data(self):
        d = {
            "time_au": float(self.t),
            "temperature_K": float(self.atoms.get_temperature()),
            # "Etot_eV": float(self.atoms.get_total_energy()),
        }
        return d
