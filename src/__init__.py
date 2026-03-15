"""
Hypersonic Shock Waves — Hipersonik Şok Dalgaları
==================================================
Physics toolkit for Mach 5+ aerodynamics.

Modules
-------
shock_relations  : Normal and oblique shock Rankine-Hugoniot equations
prandtl_meyer    : Prandtl-Meyer expansion fan
bow_shock        : Bow-shock geometry and standoff distance
visualization    : Plotting helpers
"""

from .shock_relations import NormalShock, ObliqueShock
from .prandtl_meyer import PrandtlMeyer
from .bow_shock import BowShock

__version__ = "1.0.0"
__all__ = ["NormalShock", "ObliqueShock", "PrandtlMeyer", "BowShock"]
