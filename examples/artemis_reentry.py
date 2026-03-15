"""
NASA Artemis Capsule Re-entry — Hypersonic Shock Analysis
===========================================================
The Orion capsule (Artemis programme) re-enters the atmosphere at
approximately Mach 32 at ~120 km altitude and decelerates through the
peak-heating corridor around Mach 7–10 at ~60–70 km altitude.

This script analyses three representative trajectory points:

  * Point A — Peak dynamic pressure,  Mach 7.6,  h = 65 km
  * Point B — Peak heating,            Mach 10,   h = 75 km
  * Point C — Deep hypersonic entry,   Mach 25,   h = 100 km

For each point the script computes:
  - Normal shock relations at the stagnation point
  - Bow-shock standoff distance (Billig correlation)
  - Stagnation temperature and pressure
  - Saves individual plots to the output/ directory

Usage
-----
    python examples/artemis_reentry.py
"""

import os
import sys
import math

# allow running from repo root without install
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from src.shock_relations import NormalShock, ObliqueShock
from src.bow_shock import BowShock
from src.visualization import (
    plot_bow_shock,
    plot_normal_shock_properties,
    plot_stagnation_heating,
)

# ---------------------------------------------------------------------------
# Orion capsule geometry (approximate)
# ---------------------------------------------------------------------------
ORION_NOSE_RADIUS = 6.0        # metres (spherical heat-shield nose radius)
ORION_CONE_HALF_ANGLE = 32.5   # degrees (forebody half-cone angle)

# ---------------------------------------------------------------------------
# US Standard Atmosphere — simplified values for trajectory points
# ---------------------------------------------------------------------------
TRAJECTORY = [
    {
        "label": "A — Peak Dynamic Pressure",
        "mach": 7.6,
        "altitude_km": 65,
        "T_inf": 233.0,   # K
        "p_inf": 10.9,    # Pa
        "rho_inf": 1.63e-4,  # kg/m³
    },
    {
        "label": "B — Peak Heating",
        "mach": 10.0,
        "altitude_km": 75,
        "T_inf": 208.0,
        "p_inf": 2.4,
        "rho_inf": 4.0e-5,
    },
    {
        "label": "C — Deep Hypersonic Entry",
        "mach": 25.0,
        "altitude_km": 100,
        "T_inf": 195.0,
        "p_inf": 0.032,
        "rho_inf": 5.6e-7,
    },
]

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def analyse_trajectory_point(pt: dict):
    """Print and return key aerothermodynamic properties for one trajectory point."""
    M = pt["mach"]
    T_inf = pt["T_inf"]
    p_inf = pt["p_inf"]
    gamma = 1.4

    print(f"\n{'='*60}")
    print(f"  Trajectory Point  :  {pt['label']}")
    print(f"  Altitude          :  {pt['altitude_km']} km")
    print(f"  Mach number       :  {M}")
    print(f"  T∞                :  {T_inf:.1f} K")
    print(f"  p∞                :  {p_inf:.4f} Pa")
    print(f"{'='*60}")

    # 1. Normal shock at stagnation point
    ns = NormalShock(M, gamma)
    T2 = T_inf * ns.T_ratio
    p2 = p_inf * ns.p_ratio

    print(f"\n--- Normal Shock (stagnation) ---")
    print(f"  M₂                : {ns.M2:.4f}")
    print(f"  p₂/p₁             : {ns.p_ratio:.2f}")
    print(f"  ρ₂/ρ₁             : {ns.rho_ratio:.4f}")
    print(f"  T₂/T₁             : {ns.T_ratio:.2f}")
    print(f"  T₂ (absolute)     : {T2:.1f} K  ({T2 - 273.15:.1f} °C)")
    print(f"  p₂ (absolute)     : {p2:.4f} Pa")
    print(f"  Total press recov : {ns.p0_ratio:.6f}")

    # 2. Stagnation temperature
    T_stag = T_inf * (1 + (gamma - 1) / 2 * M ** 2)
    p_stag = p_inf * (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (gamma - 1))
    print(f"\n--- Stagnation Conditions ---")
    print(f"  T₀                : {T_stag:.1f} K  ({T_stag - 273.15:.1f} °C)")
    print(f"  p₀ (isentropic)   : {p_stag:.4f} Pa")

    # 3. Bow shock
    bs = BowShock(M, ORION_NOSE_RADIUS, gamma, "sphere")
    print(f"\n--- Bow Shock (Billig correlation, R_nose = {ORION_NOSE_RADIUS} m) ---")
    print(f"  Standoff distance : {bs.delta:.4f} m  (δ/R = {bs.delta/ORION_NOSE_RADIUS:.4f})")
    print(f"  Shock Rc          : {bs.Rc:.4f} m")

    # 4. Forebody oblique shock (cone side)
    try:
        os_ = ObliqueShock(M, theta=ORION_CONE_HALF_ANGLE, gamma=gamma)
        print(f"\n--- Oblique Shock (forebody cone, θ = {ORION_CONE_HALF_ANGLE}°) ---")
        print(f"  Shock angle β     : {os_.beta:.2f}°")
        print(f"  Downstream M      : {os_.M2:.4f}")
        print(f"  p₂/p₁             : {os_.p_ratio:.4f}")
        print(f"  T₂/T₁             : {os_.T_ratio:.4f}")
    except ValueError as exc:
        print(f"\n  [Detached oblique shock — {exc}]")

    return ns, bs


def plot_artemis_bow_shocks():
    """Side-by-side bow shock plots for all three trajectory points."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 6), tight_layout=True)
    fig.suptitle(
        "Orion Capsule Bow Shock — Three Trajectory Points",
        fontsize=14, fontweight="bold",
    )

    for ax, pt in zip(axes, TRAJECTORY):
        M = pt["mach"]
        bs = BowShock(M, ORION_NOSE_RADIUS, body_type="sphere")
        x_shock, y_shock, x_body, y_body = bs.shock_shape(300)

        ax.plot(x_shock, y_shock, "r-", lw=2, label="Bow shock")
        ax.plot(x_shock, -y_shock, "r-", lw=2)
        ax.fill_betweenx(y_body, x_body, x_body.max() + ORION_NOSE_RADIUS,
                         alpha=0.3, color="steelblue")
        ax.fill_betweenx(-y_body, x_body, x_body.max() + ORION_NOSE_RADIUS,
                         alpha=0.3, color="steelblue")
        ax.plot(x_body, y_body, "b-", lw=2)
        ax.plot(x_body, -y_body, "b-", lw=2)

        info = (
            f"δ/R = {bs.delta/ORION_NOSE_RADIUS:.4f}\n"
            f"T₂/T₁ = {bs.ns.T_ratio:.1f}\n"
            f"p₂/p₁ = {bs.ns.p_ratio:.1f}"
        )
        ax.text(0.03, 0.97, info, transform=ax.transAxes,
                fontsize=8, va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.85))
        ax.set_title(f"{pt['label']}\nM∞ = {M}", fontsize=9)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_aspect("equal")

    out = os.path.join(OUTPUT_DIR, "artemis_bow_shocks.png")
    fig.savefig(out, bbox_inches="tight")
    print(f"\n[Saved] {out}")
    return fig


def plot_artemis_summary_table():
    """Bar chart summary of stagnation temperatures across trajectory points."""
    labels, T_stag_vals, p_stag_vals = [], [], []
    gamma = 1.4
    for pt in TRAJECTORY:
        M, T_inf, p_inf = pt["mach"], pt["T_inf"], pt["p_inf"]
        T_stag = T_inf * (1 + (gamma - 1) / 2 * M ** 2)
        p_stag = p_inf * (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (gamma - 1))
        labels.append(f"M = {M}")
        T_stag_vals.append(T_stag)
        p_stag_vals.append(p_stag)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), tight_layout=True)
    fig.suptitle("Orion Capsule — Stagnation Properties", fontsize=13, fontweight="bold")

    colors = ["#d62728", "#ff7f0e", "#1f77b4"]
    ax1.bar(labels, T_stag_vals, color=colors, edgecolor="k", linewidth=0.8)
    ax1.axhline(1811, color="orange", ls="--", lw=1.5, label="Steel melt ~1811 K")
    ax1.axhline(3695, color="purple", ls="--", lw=1.5, label="W melt ~3695 K")
    ax1.set_ylabel("Stagnation Temperature T₀ (K)")
    ax1.set_title("Stagnation Temperature")
    ax1.legend(fontsize=8)
    for bar, val in zip(ax1.patches, T_stag_vals):
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 50,
                 f"{val:.0f} K", ha="center", va="bottom", fontsize=9)

    ax2.bar(labels, p_stag_vals, color=colors, edgecolor="k", linewidth=0.8)
    ax2.set_ylabel("Stagnation Pressure p₀ (Pa)")
    ax2.set_title("Stagnation Pressure")
    for bar, val in zip(ax2.patches, p_stag_vals):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.02,
                 f"{val:.3g} Pa", ha="center", va="bottom", fontsize=9)

    out = os.path.join(OUTPUT_DIR, "artemis_summary.png")
    fig.savefig(out, bbox_inches="tight")
    print(f"[Saved] {out}")
    return fig


if __name__ == "__main__":
    print("=" * 60)
    print("  NASA Artemis / Orion — Hypersonic Re-entry Analysis")
    print("=" * 60)

    for pt in TRAJECTORY:
        analyse_trajectory_point(pt)

    plot_artemis_bow_shocks()
    plot_artemis_summary_table()

    # also save the general property sweep
    plot_normal_shock_properties(
        save_path=os.path.join(OUTPUT_DIR, "normal_shock_properties.png")
    )
    plot_stagnation_heating(
        save_path=os.path.join(OUTPUT_DIR, "stagnation_heating.png")
    )

    print("\nAll output figures saved to  output/")
