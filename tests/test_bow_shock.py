"""
Unit tests — Bow Shock geometry and standoff distance
"""

import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
import numpy as np
from src.bow_shock import BowShock


class TestBowShock:
    def test_standoff_distance_positive(self):
        bs = BowShock(M1=5.0, nose_radius=1.0)
        assert bs.delta > 0

    def test_standoff_over_R_decreases_with_mach(self):
        """δ/R decreases as M increases (Billig correlation)."""
        bs5 = BowShock(M1=5.0, nose_radius=1.0)
        bs10 = BowShock(M1=10.0, nose_radius=1.0)
        assert bs10.delta / bs10.nose_radius < bs5.delta / bs5.nose_radius

    def test_nose_radius_scales_linearly(self):
        bs1 = BowShock(M1=7.6, nose_radius=1.0)
        bs2 = BowShock(M1=7.6, nose_radius=2.0)
        assert abs(bs2.delta / bs1.delta - 2.0) < 1e-6

    def test_Rc_positive(self):
        bs = BowShock(M1=7.6, nose_radius=1.0)
        assert bs.Rc > 0

    def test_stagnation_pressure_gt1(self):
        bs = BowShock(M1=5.0, nose_radius=1.0)
        assert bs.stagnation_pressure > 1.0

    def test_shock_shape_returns_correct_shapes(self):
        bs = BowShock(M1=7.6, nose_radius=1.0)
        x_s, y_s, x_b, y_b = bs.shock_shape(100)
        assert len(x_s) == 100
        assert len(y_s) == 100

    def test_shock_shape_y_nonnegative(self):
        bs = BowShock(M1=7.6, nose_radius=1.0)
        _, y_s, _, _ = bs.shock_shape(100)
        assert np.all(y_s >= 0)

    def test_invalid_body_type_raises(self):
        with pytest.raises(ValueError):
            BowShock(M1=7.6, nose_radius=1.0, body_type="cube")

    def test_subsonic_raises(self):
        with pytest.raises(ValueError):
            BowShock(M1=0.8)

    def test_summary_keys(self):
        bs = BowShock(M1=7.6, nose_radius=1.0)
        s = bs.summary()
        assert "standoff_distance_m" in s
        assert "delta_over_R" in s
        assert "stagnation_p_ratio" in s

    def test_normal_shock_property(self):
        bs = BowShock(M1=10.0, nose_radius=1.0)
        assert bs.ns.M2 < 1.0  # subsonic behind stagnation normal shock

    def test_repr_contains_standoff(self):
        bs = BowShock(M1=7.6, nose_radius=1.0)
        assert "Standoff" in repr(bs)
