import numpy as np

def H0_demo(R: float) -> np.ndarray:
    k = 0.02
    R0 = 0.0
    R1 = -2.0
    Delta = 0.04
    Vc = 0.02
    a = 0.5

    V00 = 0.5 * k * (R + R0) ** 2
    V11 = 0.5 * k * (R + R1) ** 2 + Delta
    V01 = Vc * np.exp(-a * R * R)
    return np.array([[V00, V01], [V01, V11]], dtype=np.float64)

def dH0_dR_demo(R: float) -> np.ndarray:
    k = 0.02
    R0 = 0.0
    R1 = -2.0
    Vc = 0.02
    a = 0.5

    dV00 = k * (R + R0)
    dV11 = k * (R + R1)
    dV01 = Vc * np.exp(-a * R * R) * (-2.0 * a * R)
    return np.array([[dV00, dV01], [dV01, dV11]], dtype=np.float64)

def mu_demo(R: float) -> np.ndarray:
    mu01_0 = 5.0
    b = 0.2
    mu01 = mu01_0 * np.exp(-b * R * R)
    mu00 = 0.0
    mu11 = 0.0
    return np.array([[mu00, mu01], [mu01, mu11]], dtype=np.float64)

def dmu_dR_demo(R: float) -> np.ndarray:
    mu01_0 = 5.0
    b = 0.2
    mu01 = mu01_0 * np.exp(-b * R * R)
    dmu01 = mu01 * (-2.0 * b * R)
    return np.array([[0.0, dmu01], [dmu01, 0.0]], dtype=np.float64)
