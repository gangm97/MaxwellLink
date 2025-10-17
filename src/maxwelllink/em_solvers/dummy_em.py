import numpy as np
from ..molecule import Molecule


class DummyEMUnits:
    def __init__(self):
        pass

    def efield_em_to_au(self, Emu_vec3):
        return np.asarray(Emu_vec3, dtype=float)

    def source_amp_au_to_em(self, amp_au_vec3):
        return np.asarray(amp_au_vec3, dtype=float)

    def time_em_to_au(self, time_em: float):
        return float(time_em)

    def units_helper(self) -> str:
        return "Dummy EM units (1:1 conversion)"


class MoleculeDummyWrapper:
    def __init__(self, molecule: Molecule):
        self.molecule = molecule
        self.em_units = DummyEMUnits()


class DummyEMSimulation:
    def __init__(self):
        pass

    def run(self, until: float):
        pass
