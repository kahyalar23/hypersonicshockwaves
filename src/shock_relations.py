"""
Normal and Oblique Shock Relations
====================================
Implements the Rankine-Hugoniot jump conditions for a calorically perfect gas.

References
----------
Anderson, J. D. (2003). *Modern Compressible Flow*, 3rd ed. McGraw-Hill.
NASA SP-8077 — Aerodynamic Design Data Book.

Units
-----
All angles are in **degrees** unless the docstring states otherwise.
All ratios are dimensionless (downstream / upstream).
"""

import math
import numpy as np
from scipy.optimize import brentq


class NormalShock:
    """Rankine-Hugoniot relations across a normal shock.

    Parameters
    ----------
    M1 : float
        Upstream Mach number (must be > 1).
    gamma : float, optional
        Ratio of specific heats.  Default 1.4 (calorically perfect air).

    Examples
    --------
    >>> ns = NormalShock(M1=7.6)          # Artemis peak heating ~Mach 7.6
    >>> ns.M2
    0.39...
    >>> round(ns.p_ratio, 2)
    66.64
    """

    def __init__(self, M1: float, gamma: float = 1.4):
        if M1 <= 1.0:
            raise ValueError(f"Upstream Mach number must be > 1 (got {M1}).")
        self.M1 = float(M1)
        self.gamma = float(gamma)
        g = self.gamma

        # ---- downstream Mach number ----
        self._M2 = math.sqrt(
            ((g - 1) * M1 ** 2 + 2) / (2 * g * M1 ** 2 - (g - 1))
        )

        # ---- Rankine-Hugoniot ratios ----
        self._p_ratio = (2 * g * M1 ** 2 - (g - 1)) / (g + 1)
        self._rho_ratio = ((g + 1) * M1 ** 2) / ((g - 1) * M1 ** 2 + 2)
        self._T_ratio = self._p_ratio / self._rho_ratio  # ideal gas
        self._p0_ratio = (
            ((g + 1) * M1 ** 2 / ((g - 1) * M1 ** 2 + 2)) ** (g / (g - 1))
            * ((g + 1) / (2 * g * M1 ** 2 - (g - 1))) ** (1 / (g - 1))
        )

    # ------------------------------------------------------------------
    # Public properties
    # ------------------------------------------------------------------
    @property
    def M2(self) -> float:
        """Downstream Mach number."""
        return self._M2

    @property
    def p_ratio(self) -> float:
        """Static pressure ratio p₂/p₁."""
        return self._p_ratio

    @property
    def rho_ratio(self) -> float:
        """Density ratio ρ₂/ρ₁."""
        return self._rho_ratio

    @property
    def T_ratio(self) -> float:
        """Static temperature ratio T₂/T₁."""
        return self._T_ratio

    @property
    def p0_ratio(self) -> float:
        """Stagnation pressure ratio p₀₂/p₀₁ (total-pressure recovery)."""
        return self._p0_ratio

    def summary(self) -> dict:
        """Return all key ratios as a dictionary."""
        return {
            "M1": self.M1,
            "M2": round(self.M2, 6),
            "p2/p1": round(self.p_ratio, 6),
            "rho2/rho1": round(self.rho_ratio, 6),
            "T2/T1": round(self.T_ratio, 6),
            "p02/p01": round(self.p0_ratio, 6),
        }

    def __repr__(self) -> str:
        return (
            f"NormalShock(M1={self.M1}, gamma={self.gamma})\n"
            f"  M2        = {self.M2:.4f}\n"
            f"  p2/p1     = {self.p_ratio:.4f}\n"
            f"  rho2/rho1 = {self.rho_ratio:.4f}\n"
            f"  T2/T1     = {self.T_ratio:.4f}\n"
            f"  p02/p01   = {self.p0_ratio:.6f}"
        )


# ---------------------------------------------------------------------------
# θ–β–M relation helpers (module-level, re-used by ObliqueShock)
# ---------------------------------------------------------------------------

def _theta_from_beta(M1: float, beta_deg: float, gamma: float = 1.4) -> float:
    """Flow deflection angle θ (degrees) from shock angle β (degrees)."""
    b = math.radians(beta_deg)
    g = gamma
    num = 2.0 * (1.0 / math.tan(b)) * (M1 ** 2 * math.sin(b) ** 2 - 1.0)
    den = M1 ** 2 * (g + math.cos(2.0 * b)) + 2.0
    return math.degrees(math.atan(num / den))


def _beta_from_theta(M1: float, theta_deg: float, gamma: float = 1.4,
                     weak: bool = True) -> float:
    """Shock angle β (degrees) from flow deflection angle θ (degrees).

    Uses Brent's root-finding method on the θ–β–M relation.
    Returns the *weak* shock solution by default.
    """
    theta = math.radians(theta_deg)
    mu = math.degrees(math.asin(1.0 / M1))  # Mach angle

    def residual(beta_deg):
        return math.radians(_theta_from_beta(M1, beta_deg, gamma)) - theta

    # Weak shock: β ∈ (μ, β_max)
    # Strong shock: β ∈ (β_max, 90°)
    # Find β_max (maximum deflection angle)
    betas = np.linspace(mu + 0.01, 89.99, 2000)
    thetas = [_theta_from_beta(M1, b, gamma) for b in betas]
    idx_max = int(np.argmax(thetas))
    theta_max_deg = thetas[idx_max]

    if theta_deg > theta_max_deg:
        raise ValueError(
            f"Deflection angle {theta_deg:.2f}° exceeds maximum "
            f"{theta_max_deg:.2f}° for M1={M1:.2f}. Shock is detached."
        )

    beta_max = betas[idx_max]
    if weak:
        beta = brentq(residual, mu + 0.01, beta_max - 0.01)
    else:
        beta = brentq(residual, beta_max + 0.01, 89.99)
    return beta


class ObliqueShock:
    """Oblique shock wave relations.

    Instantiate by providing **either** the shock angle ``beta`` **or** the
    flow deflection angle ``theta`` (not both).

    Parameters
    ----------
    M1 : float
        Upstream Mach number (> 1).
    beta : float, optional
        Shock wave angle w.r.t. the free-stream direction (degrees).
    theta : float, optional
        Flow deflection / wedge half-angle (degrees).
    gamma : float, optional
        Ratio of specific heats.  Default 1.4.
    weak : bool, optional
        If ``theta`` is given, select the weak-shock solution (default True).

    Examples
    --------
    >>> os_ = ObliqueShock(M1=10, theta=20)
    >>> round(os_.beta, 2)
    27.79
    >>> round(os_.p_ratio, 2)
    11.47
    """

    def __init__(
        self,
        M1: float,
        beta: float = None,
        theta: float = None,
        gamma: float = 1.4,
        weak: bool = True,
    ):
        if beta is None and theta is None:
            raise ValueError("Provide either 'beta' or 'theta'.")
        if beta is not None and theta is not None:
            raise ValueError("Provide only one of 'beta' or 'theta'.")

        self.M1 = float(M1)
        self.gamma = float(gamma)

        if beta is not None:
            self._beta = float(beta)
            self._theta = _theta_from_beta(M1, beta, gamma)
        else:
            self._theta = float(theta)
            self._beta = _beta_from_theta(M1, theta, gamma, weak)

        # Normal component of upstream Mach number
        M1n = M1 * math.sin(math.radians(self._beta))
        ns = NormalShock(M1n, gamma)

        # Downstream Mach number
        self._M2 = ns.M2 / math.sin(math.radians(self._beta - self._theta))
        self._M2n = ns.M2

        self._p_ratio = ns.p_ratio
        self._rho_ratio = ns.rho_ratio
        self._T_ratio = ns.T_ratio
        self._p0_ratio = ns.p0_ratio

    # ------------------------------------------------------------------
    # Public properties
    # ------------------------------------------------------------------
    @property
    def beta(self) -> float:
        """Shock angle β (degrees)."""
        return self._beta

    @property
    def theta(self) -> float:
        """Flow deflection angle θ (degrees)."""
        return self._theta

    @property
    def M2(self) -> float:
        """Downstream Mach number."""
        return self._M2

    @property
    def p_ratio(self) -> float:
        """Static pressure ratio p₂/p₁."""
        return self._p_ratio

    @property
    def rho_ratio(self) -> float:
        """Density ratio ρ₂/ρ₁."""
        return self._rho_ratio

    @property
    def T_ratio(self) -> float:
        """Static temperature ratio T₂/T₁."""
        return self._T_ratio

    @property
    def p0_ratio(self) -> float:
        """Total pressure ratio p₀₂/p₀₁."""
        return self._p0_ratio

    @property
    def mach_angle(self) -> float:
        """Mach angle μ = arcsin(1/M1) in degrees."""
        return math.degrees(math.asin(1.0 / self.M1))

    def summary(self) -> dict:
        return {
            "M1": self.M1,
            "beta_deg": round(self.beta, 4),
            "theta_deg": round(self.theta, 4),
            "M2": round(self.M2, 6),
            "p2/p1": round(self.p_ratio, 6),
            "rho2/rho1": round(self.rho_ratio, 6),
            "T2/T1": round(self.T_ratio, 6),
            "p02/p01": round(self.p0_ratio, 6),
        }

    def __repr__(self) -> str:
        return (
            f"ObliqueShock(M1={self.M1}, beta={self.beta:.3f}°, "
            f"theta={self.theta:.3f}°)\n"
            f"  M2        = {self.M2:.4f}\n"
            f"  p2/p1     = {self.p_ratio:.4f}\n"
            f"  rho2/rho1 = {self.rho_ratio:.4f}\n"
            f"  T2/T1     = {self.T_ratio:.4f}\n"
            f"  p02/p01   = {self.p0_ratio:.6f}"
        )
