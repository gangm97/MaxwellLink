import numpy as np


class DummyModel:
    """
    A dummy quantum dynamics model for demonstration purposes.
    This class serves as a template for implementing specific quantum dynamics models.
    It provides the necessary interface for integration with the MaxwellLink framework.
    """

    def __init__(self, verbose=False, checkpoint=False, restart=False):
        """
        Initialize the necessary parameters for the dummy quantum dynamics model.

        TIPS: The computational load of this step should be minimal.

        This method *should be* overridden by subclasses if more member variables need to be initialized.

        + **`verbose`** (bool): Whether to print verbose output. Default is False.
        + **`checkpoint`** (bool): Whether to enable checkpointing. Default is False.
        + **`restart`** (bool): Whether to restart from a checkpoint if available. Default is False.
        """
        self.dt = 0.0  # time step in a.u.
        self.molecule_id = -1  # molecule ID
        self.verbose = verbose
        self.t = 0.0  # current time in a.u.

        self._preview = None  # deep-copied molecular state after the proposed step
        self._pending_amp = None  # amplitude (dmu/dt) from the previewed step
        self._have = False  # whether a step has been finished and amplitude is ready
        self.checkpoint = checkpoint
        self.restart = restart

    def initialize(self, dt_new, molecule_id):
        """
        Set the time step and molecule ID for this quantum dynamics model and provide necessary initialization.
        This function will be called in the driver code after the molecule ID is assigned
        (the INIT stage of socket communication).

        TIPS: The major computational load of initialization should be done here.

        This method *should be* overridden by subclasses if more member variables need to be initialized.

        + **`dt_new`** (float): The new time step in atomic units (a.u.).
        + **`molecule_id`** (int): The ID of the molecule assigned by SocketHub.
        """
        self.dt = float(dt_new)  # reset the time step
        self.molecule_id = int(
            molecule_id
        )  # reset the molecule ID assigned by SocketHub
        # perform any additional initialization here as needed

    def propagate(self, effective_efield_vec):
        """
        Propagate the quantum molecular dynamics for one MEEP step given the effective electric field vector.
        This method should be overridden by subclasses to implement specific propagation logic.

        TIPS: One can implement sub-steps (running many steps for the model per MEEP call) or
        macrosteps (running one step for the model per few MEEP calls) within this function as needed.

        This method *must be* overridden by subclasses.

        + **`effective_efield_vec`**: Effective electric field vector in the form [Ex, Ey, Ez].
        """
        raise NotImplementedError("This method should be overridden by subclasses.")

    def calc_amp_vector(self):
        """
        Update the source amplitude vector after propagating this molecule for one time step.
        This method should be overridden by subclasses to implement specific source update logic.

        Amplitude vector should be calculated by dmu/dt, where mu is the classical dipole vector of
        the molecule.

        This method *must be* overridden by subclasses.

        Returns:
        - A numpy array representing the amplitude vector in the form [dmu_x/dt, dmu_y/dt, dmu_z/dt].
        """
        # update the amplitude vector here as needed
        return np.array([0.0, 0.0, 0.0])

    def append_additional_data(self):
        """
        Append additional data to be sent back to MaxwellLink, which can be retrieved by the user
        via: MaxwellLink.SocketMolecule.additional_data_history, where additional_data_history is a list of
        dictionaries.

        This method can be *optionally* overridden by subclasses to send additional data to MaxwellLink.

        Returns:
        - A dictionary containing additional data.
        """
        data = {}
        return data

    def _dump_to_checkpoint(self, molid):
        """
        Dump the internal state of the model to a checkpoint. Please implement this method if checkpointing is needed.

        This method can be *optionally* overridden by subclasses to implement specific checkpoint logic.

        + **`molid`**: Molecule ID for checkpointing.
        """
        pass

    def _reset_from_checkpoint(self, molid):
        """
        Reset the internal state of the model from a checkpoint. Please implement this method if one needs to restart
        from a checkpoint (done in the self.initialize() stage).

        This method can be *optionally* overridden by subclasses to implement specific reset logic.

        + **`molid`**: Molecule ID for checkpointing.
        """
        pass

    def _snapshot(self):
        """
        Return a snapshot of the internal state for propagation. Deep copy is required.

        This method can be *optionally* overridden by subclasses if the stage-commit protocol is used.

        If not overridden, only the time variable will be snapshotted.

        Returns:
        - A dictionary containing the snapshot of the internal state.
        """
        snapshot = {
            "time": self.t,
        }
        return snapshot

    def _restore(self, snapshot):
        """
        Restore the internal state from a snapshot.

        This method can be *optionally* overridden by subclasses if the stage-commit protocol is used.

        If not overridden, only the time variable will be restored, and the stage-commit protocol defined below
        (self.stage_step and self.commit_step) is not used.
        """
        self.t = snapshot["time"]

    def stage_step(self, E_vec):
        """
        Stage a propagation step with the given effective electric field vector.
        This method performs the propagation and calculates the amplitude vector, but does not commit the changes to
        the internal state. The result can be committed later using the self.commit_step method.

        This method should *not* be overridden by subclasses.

        + **`E_vec`**: Effective electric field vector in the form [Ex, Ey, Ez].
        """
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
        """
        Check if a staged step is ready to be committed.

        This method should *not* be overridden by subclasses.

        Returns:
        - A boolean indicating whether a staged step is ready.
        """
        return self._have

    def commit_step(self):
        """
        Commit the previewed step and return the staged amplitude.
        This method applies the changes from the staged step to the internal state and returns the calculated
        amplitude vector.

        This method should *not* be overridden by subclasses.

        Returns:
        - A numpy array representing the amplitude vector in the form [dmu_x/dt, dmu_y/dt, dmu_z/dt].
        """
        if not self._have or self._preview is None or self._pending_amp is None:
            # No result staged: return zeros
            return np.zeros(3, float)

        # Commit the new molecular state
        self._restore(self._preview)
        amp = self._pending_amp

        # Optionally checkpointing
        if self.checkpoint:
            self._dump_to_checkpoint(self.molecule_id)

        # Clear staging
        self._preview = None
        self._pending_amp = None
        self._have = False
        return amp
