"""
Unit tests — Prandtl-Meyer Expansion Fan
"""

import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from src.prandtl_meyer import PrandtlMeyer, nu, mach_from_nu


# ---------------------------------------------------------------------------
# nu(M) function
# ---------------------------------------------------------------------------
class TestNuFunction:
    def test_nu_M1_is_zero(self):
        """ν(1) must be 0."""
        assert abs(nu(1.0)) < 1e-8

    def test_nu_M2(self):
        """ν(2.0) ≈ 26.38° (Anderson table A.5)."""
        assert abs(nu(2.0) - 26.38) < 0.05

    def test_nu_M5(self):
        """ν(5.0) ≈ 76.92° — computed from standard P-M formula."""
        assert abs(nu(5.0) - 76.92) < 0.1

    def test_nu_increases_with_M(self):
        assert nu(3.0) > nu(2.0) > nu(1.5)

    def test_nu_subsonic_raises(self):
        with pytest.raises(ValueError):
            nu(0.8)


# ---------------------------------------------------------------------------
# mach_from_nu (inverse)
# ---------------------------------------------------------------------------
class TestMachFromNu:
    def test_roundtrip_M3(self):
        M_in = 3.0
        nu_val = nu(M_in)
        M_out = mach_from_nu(nu_val)
        assert abs(M_out - M_in) < 1e-4

    def test_roundtrip_M10(self):
        M_in = 10.0
        nu_val = nu(M_in)
        M_out = mach_from_nu(nu_val)
        assert abs(M_out - M_in) < 1e-3

    def test_nu_zero_gives_M1(self):
        M = mach_from_nu(0.0)
        assert abs(M - 1.0) < 1e-4

    def test_negative_nu_raises(self):
        with pytest.raises(ValueError):
            mach_from_nu(-1.0)


# ---------------------------------------------------------------------------
# PrandtlMeyer class
# ---------------------------------------------------------------------------
class TestPrandtlMeyer:
    def test_M2_greater_than_M1(self):
        pm = PrandtlMeyer(M1=5, delta_theta=10)
        assert pm.M2 > pm.M1

    def test_pressure_decreases(self):
        """Expansion → pressure drops."""
        pm = PrandtlMeyer(M1=5, delta_theta=10)
        assert pm.p_ratio < 1.0

    def test_temperature_decreases(self):
        pm = PrandtlMeyer(M1=5, delta_theta=10)
        assert pm.T_ratio < 1.0

    def test_density_decreases(self):
        pm = PrandtlMeyer(M1=5, delta_theta=10)
        assert pm.rho_ratio < 1.0

    def test_M2_known(self):
        """M1=5, Δθ=10° → M2 ≈ 6.297 (from P-M function inversion)."""
        pm = PrandtlMeyer(M1=5, delta_theta=10)
        assert abs(pm.M2 - 6.297) < 0.05

    def test_mu1_correct(self):
        pm = PrandtlMeyer(M1=5, delta_theta=5)
        expected = math.degrees(math.asin(1.0 / 5.0))
        assert abs(pm.mu1 - expected) < 1e-6

    def test_mu2_less_than_mu1(self):
        """Higher M2 → smaller Mach angle."""
        pm = PrandtlMeyer(M1=5, delta_theta=10)
        assert pm.mu2 < pm.mu1

    def test_zero_turn_raises(self):
        with pytest.raises(ValueError):
            PrandtlMeyer(M1=5, delta_theta=0)

    def test_negative_turn_raises(self):
        with pytest.raises(ValueError):
            PrandtlMeyer(M1=5, delta_theta=-5)

    def test_subsonic_M1_raises(self):
        with pytest.raises(ValueError):
            PrandtlMeyer(M1=0.8, delta_theta=5)

    def test_summary_keys(self):
        pm = PrandtlMeyer(M1=5, delta_theta=10)
        s = pm.summary()
        expected_keys = {"M1", "delta_theta_deg", "nu1_deg", "nu2_deg",
                         "M2", "p2/p1", "T2/T1", "rho2/rho1", "mu1_deg", "mu2_deg"}
        assert expected_keys == set(s.keys())

    def test_large_expansion(self):
        """Large turning angle should still converge."""
        pm = PrandtlMeyer(M1=2, delta_theta=30)
        assert pm.M2 > 2.0

    def test_repr_contains_delta_theta(self):
        pm = PrandtlMeyer(M1=5, delta_theta=10)
        assert "10" in repr(pm)
