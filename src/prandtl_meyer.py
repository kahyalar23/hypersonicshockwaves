"""
Prandtl-Meyer Expansion Fan
=============================
Computes the Prandtl-Meyer function ν(M) and isentropic relations for
supersonic expansion turns.

References
----------
Anderson, J. D. (2003). *Modern Compressible Flow*, 3rd ed. McGraw-Hill.
"""

import math
import numpy as np
from scipy.optimize import brentq


def nu(M: float, gamma: float = 1.4) -> float:
    """Prandtl-Meyer function ν(M) in degrees.

    Parameters
    ----------
    M : float
        Mach number (M ≥ 1).
    gamma : float
        Ratio of specific heats (default 1.4).

    Returns
    -------
    float
        Prandtl-Meyer angle in degrees.
    """
    if M < 1.0:
        raise ValueError(f"Mach number must be ≥ 1 (got {M}).")
    g = gamma
    A = math.sqrt((g + 1) / (g - 1))
    B = math.sqrt((g - 1) / (g + 1) * (M ** 2 - 1))
    return math.degrees(A * math.atan(B) - math.atan(math.sqrt(M ** 2 - 1)))


def mach_from_nu(nu_deg: float, gamma: float = 1.4) -> float:
    """Invert the Prandtl-Meyer function: find M given ν (degrees).

    Uses Brent's root-finding method.

    Parameters
    ----------
    nu_deg : float
        Prandtl-Meyer angle in degrees.
    gamma : float
        Ratio of specific heats.

    Returns
    -------
    float
        Corresponding Mach number.
    """
    # Maximum ν for gamma=1.4 is ~130.45°
    nu_max = nu(1e6, gamma)
    if nu_deg < 0 or nu_deg > nu_max:
        raise ValueError(
            f"ν = {nu_deg:.2f}° is outside [0, {nu_max:.2f}°]."
        )
    return brentq(lambda m: nu(m, gamma) - nu_deg, 1.0, 1e5)


class PrandtlMeyer:
    """Prandtl-Meyer expansion from Mach M1 through a turning angle Δθ.

    Parameters
    ----------
    M1 : float
        Upstream Mach number (≥ 1).
    delta_theta : float
        Expansion turning angle in degrees (must be positive).
    gamma : float
        Ratio of specific heats (default 1.4).

    Examples
    --------
    >>> pm = PrandtlMeyer(M1=5, delta_theta=10)
    >>> round(pm.M2, 3)
    6.297
    """

    def __init__(self, M1: float, delta_theta: float, gamma: float = 1.4):
        if M1 < 1.0:
            raise ValueError(f"M1 must be ≥ 1 (got {M1}).")
        if delta_theta <= 0:
            raise ValueError("Expansion turning angle must be positive.")

        self.M1 = float(M1)
        self.delta_theta = float(delta_theta)
        self.gamma = float(gamma)

        nu1 = nu(M1, gamma)
        nu2 = nu1 + delta_theta
        self._nu1 = nu1
        self._nu2 = nu2
        self._M2 = mach_from_nu(nu2, gamma)

        # Isentropic relations (p/p0, T/T0 are functions of M)
        g = gamma
        self._T_ratio = self._isen_T_ratio(M1, g) / self._isen_T_ratio(self._M2, g)
        self._p_ratio = self._isen_p_ratio(M1, g) / self._isen_p_ratio(self._M2, g)
        self._rho_ratio = self._isen_rho_ratio(M1, g) / self._isen_rho_ratio(self._M2, g)

    @staticmethod
    def _isen_T_ratio(M: float, g: float) -> float:
        """T₀/T = 1 + (γ-1)/2 · M²"""
        return 1.0 + (g - 1.0) / 2.0 * M ** 2

    @staticmethod
    def _isen_p_ratio(M: float, g: float) -> float:
        """p₀/p = (1 + (γ-1)/2 · M²)^(γ/(γ-1))"""
        return (1.0 + (g - 1.0) / 2.0 * M ** 2) ** (g / (g - 1.0))

    @staticmethod
    def _isen_rho_ratio(M: float, g: float) -> float:
        """ρ₀/ρ = (1 + (γ-1)/2 · M²)^(1/(γ-1))"""
        return (1.0 + (g - 1.0) / 2.0 * M ** 2) ** (1.0 / (g - 1.0))

    # ------------------------------------------------------------------
    # Mach angle of the leading and trailing Mach waves
    # ------------------------------------------------------------------
    @property
    def mu1(self) -> float:
        """Leading Mach wave angle μ₁ = arcsin(1/M1) in degrees."""
        return math.degrees(math.asin(1.0 / self.M1))

    @property
    def mu2(self) -> float:
        """Trailing Mach wave angle μ₂ = arcsin(1/M2) in degrees."""
        return math.degrees(math.asin(1.0 / self._M2))

    @property
    def M2(self) -> float:
        """Downstream Mach number."""
        return self._M2

    @property
    def T_ratio(self) -> float:
        """Static temperature ratio T₂/T₁."""
        return self._T_ratio

    @property
    def p_ratio(self) -> float:
        """Static pressure ratio p₂/p₁ (< 1 for expansion)."""
        return self._p_ratio

    @property
    def rho_ratio(self) -> float:
        """Density ratio ρ₂/ρ₁."""
        return self._rho_ratio

    def summary(self) -> dict:
        return {
            "M1": self.M1,
            "delta_theta_deg": self.delta_theta,
            "nu1_deg": round(self._nu1, 4),
            "nu2_deg": round(self._nu2, 4),
            "M2": round(self.M2, 6),
            "p2/p1": round(self.p_ratio, 6),
            "T2/T1": round(self.T_ratio, 6),
            "rho2/rho1": round(self.rho_ratio, 6),
            "mu1_deg": round(self.mu1, 4),
            "mu2_deg": round(self.mu2, 4),
        }

    def __repr__(self) -> str:
        return (
            f"PrandtlMeyer(M1={self.M1}, Δθ={self.delta_theta}°)\n"
            f"  ν₁        = {self._nu1:.4f}°\n"
            f"  ν₂        = {self._nu2:.4f}°\n"
            f"  M2        = {self.M2:.4f}\n"
            f"  p2/p1     = {self.p_ratio:.4f}\n"
            f"  T2/T1     = {self.T_ratio:.4f}\n"
            f"  rho2/rho1 = {self.rho_ratio:.4f}"
        )
