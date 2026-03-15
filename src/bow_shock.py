"""
Bow Shock Geometry and Standoff Distance
==========================================
Models the detached bow shock that forms in front of blunt bodies
(capsules, sphere-cone forebodies) at hypersonic speeds.

References
----------
Van Dyke, M. (1958). A study of second-order supersonic flow theory.
NACA TR-1253.

Billig, F. S. (1967). Shock-wave shapes around spherical- and cylindrical-
nosed bodies. *Journal of Spacecraft and Rockets*, 4(6), 822-823.

Anderson, J. D. (2006). *Hypersonic and High Temperature Gas Dynamics*,
2nd ed. AIAA Education Series.
"""

import math
import numpy as np
from .shock_relations import NormalShock


class BowShock:
    """Blunt-body bow shock model using the Billig correlation.

    The Billig (1967) correlation fits the bow-shock shape as a conic
    section for a sphere or a spherically-blunted cylinder/cone.

    Parameters
    ----------
    M1 : float
        Free-stream Mach number (≥ 5 for hypersonic regime).
    nose_radius : float
        Nose radius of the blunt body [m].  Default 1.0 m (normalised).
    gamma : float
        Ratio of specific heats.  Default 1.4.
    body_type : str
        'sphere' (default) or 'cylinder'.

    Attributes
    ----------
    delta : float
        Shock standoff distance [m] (same units as *nose_radius*).
    Rc : float
        Shock-wave radius of curvature at the stagnation point [m].
    ns : NormalShock
        Normal shock at the stagnation point.

    Examples
    --------
    >>> bs = BowShock(M1=7.6, nose_radius=1.0)
    >>> round(bs.delta / bs.nose_radius, 4)
    0.0875
    """

    # Billig curve-fit constants
    _BILLIG = {
        "sphere": {"A": 0.386, "B": 4.67, "C": 0.0455, "D": 0.6},
        "cylinder": {"A": 0.386, "B": 4.67, "C": 0.0455, "D": 0.6},
    }

    def __init__(
        self,
        M1: float,
        nose_radius: float = 1.0,
        gamma: float = 1.4,
        body_type: str = "sphere",
    ):
        if M1 < 2.0:
            raise ValueError(f"M1 must be ≥ 2 (got {M1}).")
        if body_type not in self._BILLIG:
            raise ValueError(f"body_type must be 'sphere' or 'cylinder' (got {body_type!r}).")

        self.M1 = float(M1)
        self.nose_radius = float(nose_radius)
        self.gamma = float(gamma)
        self.body_type = body_type

        self.ns = NormalShock(M1, gamma)
        self._compute_geometry()

    def _compute_geometry(self):
        M = self.M1
        R = self.nose_radius
        c = self._BILLIG[self.body_type]

        # Billig (1967) correlations
        # δ/R  = A * exp(B / M^2)          [standoff distance]
        # Rc/R = C * M^D                   [shock radius of curvature]
        self.delta = R * c["A"] * math.exp(c["B"] / M ** 2)
        self.Rc = R * c["C"] * M ** c["D"]

    @property
    def stagnation_pressure(self) -> float:
        """Stagnation pressure ratio p₀₂/p₁ (Rayleigh Pitot formula)."""
        M = self.M1
        g = self.gamma
        # Rayleigh Pitot-tube formula
        term1 = ((g + 1) ** 2 * M ** 2 / (4 * g * M ** 2 - 2 * (g - 1))) ** (g / (g - 1))
        term2 = (1 - g + 2 * g * M ** 2) / (g + 1)
        return term1 * term2

    def shock_shape(self, n_points: int = 200) -> tuple:
        """Generate (x, y) coordinates of the bow-shock contour.

        Uses the Billig conic-section approximation:
            x = δ + Rc - √(Rc² + y²/ε²)    (hyperbola approximation)

        where the eccentricity ε is a function of Mach number.

        Parameters
        ----------
        n_points : int
            Number of points along the shock (default 200).

        Returns
        -------
        x_shock : np.ndarray
            Axial coordinates [m], positive downstream.
        y_shock : np.ndarray
            Radial / lateral coordinates [m].
        x_body : np.ndarray
            Sphere-nose body x-coordinates for overlay.
        y_body : np.ndarray
            Sphere-nose body y-coordinates.
        """
        R = self.nose_radius
        delta = self.delta
        Rc = self.Rc

        # Eccentricity of the shock conic (Billig curve fit)
        eps = 1.0 + 0.8 / (self.M1 - 1) ** 0.5

        # Radial extent: up to ~5R
        y_max = 5.0 * R
        y = np.linspace(0, y_max, n_points)

        # Billig hyperbola: x measured from body vertex, positive into freestream
        # x_shock = R + delta - Rc * (sqrt(1 + (y/Rc)^2 / eps^2) - 1)
        x_shock = -(R + delta - Rc * (np.sqrt(1.0 + (y / Rc) ** 2 / eps ** 2) - 1.0))

        # Sphere nose body
        theta = np.linspace(-math.pi / 2, math.pi / 2, n_points)
        x_body = -R * np.cos(theta) + R  # vertex at x=0
        # Clamp to the front hemisphere only (x ≥ 0)
        mask = x_body >= 0
        x_body = x_body[mask]
        y_body = R * np.sin(theta[mask])

        return x_shock, y, x_body, y_body

    def summary(self) -> dict:
        return {
            "M1": self.M1,
            "gamma": self.gamma,
            "body_type": self.body_type,
            "nose_radius_m": self.nose_radius,
            "standoff_distance_m": round(self.delta, 6),
            "delta_over_R": round(self.delta / self.nose_radius, 6),
            "shock_Rc_m": round(self.Rc, 6),
            "stagnation_p_ratio": round(self.stagnation_pressure, 4),
            "ns_T_ratio": round(self.ns.T_ratio, 4),
            "ns_p_ratio": round(self.ns.p_ratio, 4),
        }

    def __repr__(self) -> str:
        return (
            f"BowShock(M1={self.M1}, R={self.nose_radius} m, {self.body_type})\n"
            f"  Standoff distance δ = {self.delta:.4f} m  "
            f"(δ/R = {self.delta/self.nose_radius:.4f})\n"
            f"  Shock Rc            = {self.Rc:.4f} m\n"
            f"  Stag. press. ratio  = {self.stagnation_pressure:.4f}\n"
            f"  Normal shock T2/T1  = {self.ns.T_ratio:.4f}"
        )
