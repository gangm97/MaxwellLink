import numpy as np


class DummyModel:
    """
    A dummy quantum dynamics model for demonstration purposes.
    This class serves as a template for implementing specific quantum dynamics models.
    It provides the necessary interface for integration with the SocketMEEP framework.
    """

    def __init__(self, verbose=False, checkpoint=False, restart=False):
        """
        Initialize the necessary parameters for the dummy quantum dynamics model.

        TIPS: The computational load of this step should be minimal.

        + **`verbose`** (bool): Whether to print verbose output. Default is False.
        """
        self.dt = 0.0  # time step in a.u.
        self.molecule_id = -1  # molecule ID
        self.verbose = verbose
        self.t = 0.0  # current time in a.u.

        self._preview = None  # deep-copied model after the proposed step
        self._pending_amp = None  # amplitude from the previewed step
        self._have = False
        self.checkpoint = checkpoint
        self.restart = restart

    def initialize(self, dt_new, molecule_id):
        """
        Set the time step and molecule ID for this quantum dynamics model and provide necessary initialization.
        This function will be called in the driver code after the molecule ID is assigned
        (the INIT stage of socket communication).

        TIPS: The major computational load of initialization should be done here.

        + **`dt_new`** (float): The new time step in atomic units (a.u.).
        + **`molecule_id`** (int): The ID of the molecule.
        """
        self.dt = float(dt_new)
        self.molecule_id = int(molecule_id)
        # perform any additional initialization here as needed

    def propagate(self, effective_efield_vec):
        """
        Propagate the quantum molecular dynamics for one MEEP step given the effective electric field vector.
        This method should be overridden by subclasses to implement specific propagation logic.

        TIPS: One can implement sub-steps (running many steps for the model per MEEP call) or
        macrosteps (running one step for the model per few MEEP calls) within this function as needed.

        + **`effective_efield_vec`**: Effective electric field vector in the form [Ex, Ey, Ez].
        """
        raise NotImplementedError("This method should be overridden by subclasses.")

    def calc_amp_vector(self):
        """
        Update the source amplitude vector after propagating this molecule for one time step.
        This method should be overridden by subclasses to implement specific source update logic.

        Amplitude vector should be calculated by dP/dt, where P is the classical dipole vector of
        the molecule.

        Returns:
        - A numpy array representing the amplitude vector in the form [dPx/dt, dPy/dt, dPz/dt].
        """
        # update the amplitude vector here as needed
        return np.array([0.0, 0.0, 0.0])

    def append_additional_data(self):
        """
        Append additional data to be sent back to MEEP, which can be retrieved by the user
        via the Python interface of MEEP: mp.SocketMolecule.additional_data_history.

        This method can be overridden by subclasses to include specific data.
        Returns:
        - A dictionary containing additional data.
        """
        data = {}
        return data

    def _dump_to_checkpoint(self):
        """
        Dump the internal state of the model to a checkpoint. Please implement this method if checkpointing is needed.
        This method can be overridden by subclasses to implement specific checkpoint logic.
        """
        pass

    def _reset_from_checkpoint(self):
        """
        Reset the internal state of the model from a checkpoint. Please implement this method if one needs to restart
        from a checkpoint (done in the self.initialize() stage).

        This method can be overridden by subclasses to implement specific reset logic.
        """
        pass

    def _snapshot(self):
        """
        Return a snapshot of the internal state for propagation. Deep copy is required.
        """
        snapshot = {
            "time": self.t,
        }
        return snapshot

    def _restore(self, snapshot):
        """
        Restore the internal state from a snapshot.
        """
        self.t = snapshot["time"]

    def stage_step(self, E_vec):
        # 1. work on a deep copy so committed state is untouched
        previous_state = self._snapshot()
        self.propagate(effective_efield_vec=E_vec)
        amp_vec = np.asarray(self.calc_amp_vector(), dtype=float)

        preview = self._snapshot()
        self._restore(previous_state)  # restore previous state

        # 2. stash for commit
        self._preview = preview
        self._pending_amp = amp_vec
        self._have = True

    def have_result(self):
        return self._have

    def commit_step(self):
        """Commit the previewed step and return the staged amplitude."""
        if not self._have or self._preview is None or self._pending_amp is None:
            # Defensive: no result staged – return zeros
            return np.zeros(3, float)

        # Atomically replace the committed model state
        self._restore(self._preview)
        amp = self._pending_amp

        if self.checkpoint:
            self._dump_to_checkpoint(self.molecule_id)

        # Clear staging
        self._preview = None
        self._pending_amp = None
        self._have = False
        return amp
