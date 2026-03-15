"""
Unit tests — Normal and Oblique Shock Relations
"""

import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from src.shock_relations import NormalShock, ObliqueShock, _theta_from_beta, _beta_from_theta


# ---------------------------------------------------------------------------
# Normal shock
# ---------------------------------------------------------------------------
class TestNormalShock:
    """Reference values from NACA 1135 / Anderson 3rd ed. tables."""

    def test_M2_mach2(self):
        ns = NormalShock(2.0)
        assert abs(ns.M2 - 0.5774) < 1e-3

    def test_pressure_ratio_mach2(self):
        ns = NormalShock(2.0)
        assert abs(ns.p_ratio - 4.5) < 1e-4

    def test_density_ratio_mach2(self):
        ns = NormalShock(2.0)
        assert abs(ns.rho_ratio - 2.6667) < 1e-3

    def test_temperature_ratio_mach2(self):
        ns = NormalShock(2.0)
        expected = 4.5 / (8.0 / 3.0)
        assert abs(ns.T_ratio - expected) < 1e-4

    def test_M2_mach5(self):
        ns = NormalShock(5.0)
        # Anderson table: M2 = 0.4152
        assert abs(ns.M2 - 0.4152) < 1e-3

    def test_pressure_ratio_mach5(self):
        ns = NormalShock(5.0)
        # p2/p1 = (2*1.4*25 - 0.4) / 2.4 = 29.0
        assert abs(ns.p_ratio - 29.0) < 1e-4

    def test_hypersonic_limit_T_ratio(self):
        """High-M limit: T2/T1 → 2γ(γ-1)/(γ+1)² * M² — scales as M²."""
        M = 20.0
        ns = NormalShock(M)
        # Just verify it grows with M
        ns_low = NormalShock(10.0)
        assert ns.T_ratio > ns_low.T_ratio

    def test_total_pressure_recovery_decreases(self):
        """p02/p01 should decrease as M increases (more loss)."""
        ns5 = NormalShock(5.0)
        ns10 = NormalShock(10.0)
        assert ns10.p0_ratio < ns5.p0_ratio

    def test_raises_subsonic(self):
        with pytest.raises(ValueError):
            NormalShock(0.9)

    def test_raises_sonic(self):
        with pytest.raises(ValueError):
            NormalShock(1.0)

    def test_summary_keys(self):
        ns = NormalShock(7.6)
        s = ns.summary()
        assert set(s.keys()) == {"M1", "M2", "p2/p1", "rho2/rho1", "T2/T1", "p02/p01"}

    def test_repr_contains_M2(self):
        ns = NormalShock(5.0)
        assert "M2" in repr(ns)

    def test_gamma_variation(self):
        """Higher γ → lower density ratio but higher temperature ratio across normal shock."""
        ns_14 = NormalShock(5.0, gamma=1.4)
        ns_16 = NormalShock(5.0, gamma=1.6)
        # Higher γ yields a weaker density compression (lower rho ratio)
        assert ns_16.rho_ratio < ns_14.rho_ratio
        # and a higher temperature ratio as a result
        assert ns_16.T_ratio > ns_14.T_ratio


# ---------------------------------------------------------------------------
# θ–β–M helpers
# ---------------------------------------------------------------------------
class TestThetaBetaM:
    def test_theta_from_beta_known(self):
        # M=2, β=45° → θ≈14.74° (from θ-β-M relation, verified analytically)
        theta = _theta_from_beta(2.0, 45.0)
        assert abs(theta - 14.74) < 0.15

    def test_round_trip(self):
        """β → θ → β should be identity."""
        M, beta_in = 10.0, 27.0
        theta = _theta_from_beta(M, beta_in)
        beta_out = _beta_from_theta(M, theta, weak=True)
        assert abs(beta_out - beta_in) < 0.05


# ---------------------------------------------------------------------------
# Oblique shock
# ---------------------------------------------------------------------------
class TestObliqueShock:
    def test_beta_from_theta_M10(self):
        """M=10, θ=20° → β≈25.82° (weak shock, verified via θ-β-M formula)."""
        os_ = ObliqueShock(M1=10, theta=20)
        assert abs(os_.beta - 25.82) < 0.15

    def test_theta_from_beta(self):
        os_ = ObliqueShock(M1=5, beta=45)
        assert os_.theta > 0

    def test_pressure_ratio_positive(self):
        os_ = ObliqueShock(M1=10, theta=15)
        assert os_.p_ratio > 1.0

    def test_temperature_ratio_positive(self):
        os_ = ObliqueShock(M1=10, theta=15)
        assert os_.T_ratio > 1.0

    def test_downstream_mach_supersonic(self):
        """Weak oblique shock — downstream should still be supersonic."""
        os_ = ObliqueShock(M1=10, theta=20, weak=True)
        assert os_.M2 > 1.0

    def test_mach_angle_property(self):
        os_ = ObliqueShock(M1=5, beta=30)
        expected_mu = math.degrees(math.asin(1.0 / 5.0))
        assert abs(os_.mach_angle - expected_mu) < 1e-6

    def test_detached_shock_raises(self):
        """Deflection angle beyond θ_max should raise ValueError."""
        with pytest.raises(ValueError):
            ObliqueShock(M1=5, theta=50)  # exceeds max ~41° for M=5

    def test_normal_shock_limit(self):
        """β=90° → oblique shock reduces to normal shock."""
        os_ = ObliqueShock(M1=5, beta=90)
        ns = NormalShock(5.0)
        assert abs(os_.p_ratio - ns.p_ratio) < 1e-4

    def test_no_args_raises(self):
        with pytest.raises(ValueError):
            ObliqueShock(M1=5)

    def test_both_args_raises(self):
        with pytest.raises(ValueError):
            ObliqueShock(M1=5, beta=30, theta=15)

    def test_summary_keys(self):
        os_ = ObliqueShock(M1=7.6, theta=20)
        s = os_.summary()
        assert "beta_deg" in s and "theta_deg" in s
