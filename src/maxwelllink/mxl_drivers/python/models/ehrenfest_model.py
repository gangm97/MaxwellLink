import numpy as np
from scipy.linalg import expm
import os

try:
    from .dummy_model import DummyModel
except:
    from dummy_model import DummyModel


class Ehrenfest2State1DModel(DummyModel):
    """
      - Ehrenfest model linked to MaxwellLink:
      - 2 electronic states (density matrix rho, 2x2)
      - 1 classical nuclear coordinate R with momentum P
      - dipole coupling to effective local E-field (one polarization axis)
    """

    def __init__(
        self,
        # nuclear parameters
        mass: float = 2000.0,
        R_initial: float = 0.0,
        P_initial: float = 0.0,
        # dipole / field coupling
        orientation: int = 2,
        # electronic initial condition
        rho_initial=None,
        pe_initial: float = 0.0,
        # constants
        hbar: float = 1.0,
        # checkpointing / verbose
        checkpoint: bool = False,
        restart: bool = False,
        verbose: bool = False,
        # user-supplied model functions
        H0_func=None,
        dH0_dR_func=None,
        mu_func=None,
        dmu_dR_func=None,
    ):
        super().__init__(verbose, checkpoint, restart)

        # orientation (same convention as TLSModel)
        self.orientation_idx = int(orientation)
        if self.orientation_idx < 0 or self.orientation_idx > 2:
            raise ValueError("Orientation must be 0 (x), 1 (y), or 2 (z).")

        # physical parameters
        self.M = float(mass)
        self.hbar = float(hbar)

        # nuclear state
        self.R = float(R_initial)
        self.P = float(P_initial)

        # model functions
        self._H0_func = H0_func
        self._dH0_dR_func = dH0_dR_func
        self._mu_func = mu_func
        self._dmu_dR_func = dmu_dR_func

        # electronic state
        if rho_initial is not None:
            rho = np.asarray(rho_initial, dtype=np.complex128)
            if rho.shape != (2, 2):
                raise ValueError("rho_initial must be (2,2).")
            self.rho = rho
            # normalize trace just in case
            self.rho /= np.trace(self.rho).real
        else:
            self.rho = np.zeros((2, 2), dtype=np.complex128)
            self._reset_population(pe_initial)

        # bookkeeping
        self.dipole_vec = np.zeros(3, dtype=float)
        self.energy_el = 0.0
        self.energy_tot = 0.0
        self.force = 0.0

        # cache last received E along dipole axis (needed by calc_amp_vector)
        self._last_E_loc = 0.0

        # checkpoint
        self.checkpoint_filename = None
        self.restarted = False

    # -------------- heavy-load initialization (at INIT) --------------
    def initialize(self, dt_new, molecule_id):
        self.dt = float(dt_new)
        self.molecule_id = int(molecule_id)
        self.checkpoint_filename = f"ehrenfest_checkpoint_id_{self.molecule_id}.npz"
        print(
            "init Ehrenfest2State1DModel with dt = %.6f a.u., molecule ID = %d"
            % (self.dt, self.molecule_id)
        )

        if self.restart and self.checkpoint:
            self._reset_from_checkpoint()
            self.restarted = True

    # ---------------- model hooks (user must provide) ----------------

    def H0(self, R: float) -> np.ndarray:
        if self._H0_func is None:
            raise NotImplementedError("Provide H0_func(R) or override H0().")
        H = np.asarray(self._H0_func(R), dtype=np.complex128)
        if H.shape != (2, 2):
            raise ValueError("H0(R) must return (2,2).")
        return H

    def dH0_dR(self, R: float) -> np.ndarray:
        if self._dH0_dR_func is None:
            raise NotImplementedError("Provide dH0_dR_func(R) or override dH0_dR().")
        dH = np.asarray(self._dH0_dR_func(R), dtype=np.complex128)
        if dH.shape != (2, 2):
            raise ValueError("dH0_dR(R) must return (2,2).")
        return dH

    def mu_op(self, R: float) -> np.ndarray:
        if self._mu_func is None:
            raise NotImplementedError("Provide mu_func(R) or override mu_op().")
        mu = np.asarray(self._mu_func(R), dtype=np.complex128)
        if mu.shape != (2, 2):
            raise ValueError("mu_op(R) must return (2,2).")
        return mu

    def dmu_dR(self, R: float) -> np.ndarray:
        if self._dmu_dR_func is None:
            raise NotImplementedError("Provide dmu_dR_func(R) or override dmu_dR().")
        dmu = np.asarray(self._dmu_dR_func(R), dtype=np.complex128)
        if dmu.shape != (2, 2):
            raise ValueError("dmu_dR(R) must return (2,2).")
        return dmu

    # ------------ internal functions -------------
    def _reset_population(self, excited_population: float):
        """Create a pure state with specified excited population and set rho=|psi><psi|."""
        if excited_population < 0 or excited_population > 1:
            raise ValueError("Excited population must be between 0 and 1.")
        c0 = np.sqrt(1.0 - excited_population)
        c1 = np.sqrt(excited_population)
        C = np.array([[c0], [c1]], dtype=np.complex128)
        self.rho = C @ C.conj().T

    def _full_H(self, R: float, E_loc: float) -> np.ndarray:
        """H = H0(R) - mu(R)*E_loc"""
        return self.H0(R) - self.mu_op(R) * E_loc

    def _dH_dR(self, R: float, E_loc: float) -> np.ndarray:
        """dH/dR = dH0/dR - dmu/dR * E_loc"""
        return self.dH0_dR(R) - self.dmu_dR(R) * E_loc

    def _unitary(self, H: np.ndarray, dt: float) -> np.ndarray:
        """U = exp(-i H dt / hbar)"""
        return expm(-1j * H * dt / self.hbar)

    def _ehrenfest_force(self, R: float, rho: np.ndarray, E_loc: float) -> float:
        """F = -Tr(rho dH/dR)"""
        val = np.trace(rho @ self._dH_dR(R, E_loc))
        return float(-np.real(val))

    def _dipole(self, R: float, rho: np.ndarray) -> float:
        """mu = Tr(rho mu(R))"""
        return float(np.real(np.trace(rho @ self.mu_op(R))))

    def _mu_dot(self, R: float, P: float, rho: np.ndarray, E_loc: float) -> float:
        """
        d/dt <mu> = (i/hbar) Tr[rho [H,mu]] + Tr[rho dmu/dR] * (P/M)
        """
        H = self._full_H(R, E_loc)
        mu = self.mu_op(R)
        comm = H @ mu - mu @ H
        electronic = (1j / self.hbar) * np.trace(rho @ comm)
        nuclear = np.trace(rho @ self.dmu_dR(R)) * (P / self.M)
        return float(np.real(electronic + nuclear))

    # ---------------- one step propagation ----------------

    def propagate(self, effective_efield_vec):
        """
        Advance (rho, R, P) by one dt using a symmetric scheme:

          1) electronic half-step at (R_n, E_n)
          2) nuclear full-step (velocity-Verlet) using rho_half
          3) electronic half-step at (R_{n+1}, E_n)
        """
        E_loc = float(effective_efield_vec[self.orientation_idx])
        self._last_E_loc = E_loc  # cache for calc_amp_vector

        if self.verbose:
            print(
                f"[molecule ID {self.molecule_id}] t={self.t:.6f} a.u. "
                f"R={self.R:.6f} P={self.P:.6f} E_loc={E_loc:.6e}"
            )

        # --- (1) electronic half step ---
        H_n = self._full_H(self.R, E_loc)
        U_half = self._unitary(H_n, 0.5 * self.dt)
        rho_half = U_half @ self.rho @ U_half.conj().T
        rho_half = 0.5 * (rho_half + rho_half.conj().T)
        rho_half /= np.trace(rho_half).real

        # --- (2) nuclear full step: velocity-Verlet with mean-field force ---
        F_n = self._ehrenfest_force(self.R, rho_half, E_loc)
        P_half = self.P + 0.5 * self.dt * F_n
        R_new = self.R + self.dt * (P_half / self.M)

        F_new = self._ehrenfest_force(R_new, rho_half, E_loc)
        P_new = P_half + 0.5 * self.dt * F_new

        # --- (3) electronic half step at new R ---
        H_new = self._full_H(R_new, E_loc)
        U_half2 = self._unitary(H_new, 0.5 * self.dt)
        rho_new = U_half2 @ rho_half @ U_half2.conj().T
        rho_new = 0.5 * (rho_new + rho_new.conj().T)
        rho_new /= np.trace(rho_new).real

        # commit step-local update
        self.R = float(R_new)
        self.P = float(P_new)
        self.rho = rho_new
        self.t += self.dt

        # update dipole vector and energies for logging
        mu_scalar = self._dipole(self.R, self.rho)
        self.dipole_vec = np.zeros(3, dtype=float)
        self.dipole_vec[self.orientation_idx] = mu_scalar

        self.energy_el = float(np.real(np.trace(self.H0(self.R) @ self.rho)))
        T_nuc = 0.5 * (self.P**2) / self.M
        self.energy_tot = float(T_nuc + self.energy_el)
        self.force = float(F_new)

    def calc_amp_vector(self):
        """
        Return [dmu_x/dt, dmu_y/dt, dmu_z/dt].

        This is called right after propagate() in DummyModel.stage_step.
        We use cached self._last_E_loc from the same step.
        """
        mu_dot = self._mu_dot(self.R, self.P, self.rho, float(self._last_E_loc))
        amp_vec = np.zeros(3, dtype=float)
        amp_vec[self.orientation_idx] = mu_dot
        return amp_vec

    def append_additional_data(self):
        """
        Provide diagnostics similar to TLSModel.
        """
        data = {}
        data["time_au"] = self.t
        data["R_au"] = self.R
        data["P_au"] = self.P
        data["force_au"] = self.force
        data["energy_el_au"] = self.energy_el
        data["energy_tot_au"] = self.energy_tot

        data["mux_au"] = float(self.dipole_vec[0])
        data["muy_au"] = float(self.dipole_vec[1])
        data["muz_au"] = float(self.dipole_vec[2])

        data["Pg"] = float(np.real(self.rho[0, 0]))
        data["Pe"] = float(np.real(self.rho[1, 1]))
        data["Pge_real"] = float(np.real(self.rho[0, 1]))
        data["Pge_imag"] = float(np.imag(self.rho[0, 1]))
        return data

    # ---------------- optional: checkpointing ----------------

    def _dump_to_checkpoint(self):
        np.savez(
            self.checkpoint_filename,
            density_matrix=self.rho,
            time=self.t,
            R=self.R,
            P=self.P,
        )

    def _reset_from_checkpoint(self):
        if not os.path.exists(self.checkpoint_filename):
            print(
                "[checkpoint] No checkpoint file found for molecule ID %d, starting fresh."
                % self.molecule_id
            )
        else:
            data = np.load(self.checkpoint_filename)
            self.rho = np.asarray(data["density_matrix"], dtype=np.complex128)
            self.t = float(data["time"])
            self.R = float(data["R"])
            self.P = float(data["P"])

    # ---------------- stage/commit support ----------------

    def _snapshot(self):
        """
        Deep copy of state so DummyModel.stage_step can safely preview and rollback.
        """
        return {
            "time": float(self.t),
            "density_matrix": self.rho.copy(),
            "R": float(self.R),
            "P": float(self.P),
            "dipole_vec": self.dipole_vec.copy(),
            "energy_el": float(self.energy_el),
            "energy_tot": float(self.energy_tot),
            "force": float(self.force),
            "_last_E_loc": float(self._last_E_loc),
        }

    def _restore(self, snapshot):
        self.t = float(snapshot["time"])
        self.rho = snapshot["density_matrix"].copy()
        self.R = float(snapshot["R"])
        self.P = float(snapshot["P"])
        self.dipole_vec = snapshot["dipole_vec"].copy()
        self.energy_el = float(snapshot["energy_el"])
        self.energy_tot = float(snapshot["energy_tot"])
        self.force = float(snapshot["force"])
        self._last_E_loc = float(snapshot["_last_E_loc"])
