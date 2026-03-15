"""
Parametric Shock Wave Study — Mach 5 to 25
============================================
Sweeps the Mach number range from 5 to 25 and computes key shock-wave
properties for both normal and oblique shocks, then plots the results.

Usage
-----
    python examples/parametric_study.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.shock_relations import NormalShock, ObliqueShock
from src.bow_shock import BowShock
from src.visualization import (
    plot_theta_beta_mach,
    plot_oblique_shock_heatmap,
)

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

GAMMA = 1.4
M_MIN, M_MAX = 5, 25
N = 200


def normal_shock_sweep():
    """Compute and tabulate normal shock properties for M = 5 to 25."""
    M_arr = np.linspace(M_MIN, M_MAX, N)
    results = {
        "M": M_arr,
        "M2": [], "p_ratio": [], "rho_ratio": [], "T_ratio": [], "p0_ratio": [],
    }
    for M in M_arr:
        ns = NormalShock(M, GAMMA)
        results["M2"].append(ns.M2)
        results["p_ratio"].append(ns.p_ratio)
        results["rho_ratio"].append(ns.rho_ratio)
        results["T_ratio"].append(ns.T_ratio)
        results["p0_ratio"].append(ns.p0_ratio)

    for key in ("M2", "p_ratio", "rho_ratio", "T_ratio", "p0_ratio"):
        results[key] = np.array(results[key])
    return results


def oblique_shock_sweep(theta_values=(15, 25, 35)):
    """Compute oblique shock properties for selected wedge angles."""
    M_arr = np.linspace(M_MIN, M_MAX, N)
    sweeps = {}
    for theta in theta_values:
        p_r, T_r, M2_r = [], [], []
        for M in M_arr:
            try:
                os_ = ObliqueShock(M, theta=theta, gamma=GAMMA)
                p_r.append(os_.p_ratio)
                T_r.append(os_.T_ratio)
                M2_r.append(os_.M2)
            except ValueError:
                p_r.append(np.nan)
                T_r.append(np.nan)
                M2_r.append(np.nan)
        sweeps[theta] = {
            "M": M_arr,
            "p_ratio": np.array(p_r),
            "T_ratio": np.array(T_r),
            "M2": np.array(M2_r),
        }
    return sweeps


def bow_shock_standoff_sweep():
    """Compute standoff distance ratio δ/R for M = 5 to 25."""
    M_arr = np.linspace(M_MIN, M_MAX, 100)
    delta_R = []
    for M in M_arr:
        bs = BowShock(M, nose_radius=1.0, gamma=GAMMA)
        delta_R.append(bs.delta)
    return M_arr, np.array(delta_R)


def plot_parametric_normal(data: dict, save_path: str):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), tight_layout=True)
    fig.suptitle(
        f"Normal Shock Properties — Parametric Study (M = {M_MIN}–{M_MAX})",
        fontsize=13, fontweight="bold",
    )
    M = data["M"]
    sets = [
        (axes[0, 0], data["p_ratio"], "p₂/p₁", "Pressure Ratio", "#d62728"),
        (axes[0, 1], data["rho_ratio"], "ρ₂/ρ₁", "Density Ratio", "#1f77b4"),
        (axes[1, 0], data["T_ratio"], "T₂/T₁", "Temperature Ratio", "#ff7f0e"),
        (axes[1, 1], data["p0_ratio"], "p₀₂/p₀₁", "Total-Pressure Recovery", "#2ca02c"),
    ]
    for ax, vals, ylabel, title, col in sets:
        ax.plot(M, vals, color=col, lw=2)
        ax.set_xlabel("Mach Number M₁")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
    fig.savefig(save_path, bbox_inches="tight")
    print(f"[Saved] {save_path}")
    return fig


def plot_parametric_oblique(sweeps: dict, save_path: str):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), tight_layout=True)
    fig.suptitle(
        "Oblique Shock Properties — Parametric Study",
        fontsize=13, fontweight="bold",
    )
    colors = {15: "#1f77b4", 25: "#ff7f0e", 35: "#d62728"}
    for theta, data in sweeps.items():
        col = colors[theta]
        M = data["M"]
        axes[0].plot(M, data["p_ratio"], color=col, lw=2, label=f"θ = {theta}°")
        axes[1].plot(M, data["T_ratio"], color=col, lw=2, label=f"θ = {theta}°")
        axes[2].plot(M, data["M2"], color=col, lw=2, label=f"θ = {theta}°")

    axes[0].set_ylabel("p₂/p₁")
    axes[0].set_title("Pressure Ratio")
    axes[1].set_ylabel("T₂/T₁")
    axes[1].set_title("Temperature Ratio")
    axes[2].set_ylabel("M₂")
    axes[2].set_title("Downstream Mach")

    for ax in axes:
        ax.set_xlabel("Mach Number M₁")
        ax.legend(fontsize=9)

    fig.savefig(save_path, bbox_inches="tight")
    print(f"[Saved] {save_path}")
    return fig


def plot_parametric_standoff(M_arr, delta_R, save_path: str):
    fig, ax = plt.subplots(figsize=(9, 5), tight_layout=True)
    ax.plot(M_arr, delta_R, "b-", lw=2)
    ax.fill_between(M_arr, 0, delta_R, alpha=0.15, color="blue")
    ax.set_xlabel("Mach Number M∞")
    ax.set_ylabel("Standoff Distance δ / R_nose")
    ax.set_title(
        "Bow Shock Standoff Distance — Billig Correlation",
        fontweight="bold",
    )
    # Mark Artemis point
    M_art = 7.6
    from src.bow_shock import BowShock as _BS
    d_art = _BS(M_art, 1.0).delta
    ax.plot(M_art, d_art, "ro", ms=8, zorder=5, label=f"Artemis M≈{M_art}")
    ax.annotate(
        f"Artemis\nδ/R = {d_art:.4f}",
        xy=(M_art, d_art),
        xytext=(M_art + 2, d_art + 0.005),
        arrowprops=dict(arrowstyle="->", color="red"),
        fontsize=9,
        color="red",
    )
    ax.legend()
    fig.savefig(save_path, bbox_inches="tight")
    print(f"[Saved] {save_path}")
    return fig


if __name__ == "__main__":
    print("Running parametric study (M = 5 – 25)…\n")

    ns_data = normal_shock_sweep()
    plot_parametric_normal(
        ns_data,
        os.path.join(OUTPUT_DIR, "parametric_normal_shock.png"),
    )

    oblique_data = oblique_shock_sweep(theta_values=(15, 25, 35))
    plot_parametric_oblique(
        oblique_data,
        os.path.join(OUTPUT_DIR, "parametric_oblique_shock.png"),
    )

    M_arr, delta_R = bow_shock_standoff_sweep()
    plot_parametric_standoff(
        M_arr, delta_R,
        os.path.join(OUTPUT_DIR, "bow_shock_standoff.png"),
    )

    plot_theta_beta_mach(
        mach_numbers=(5, 8, 10, 15, 20),
        save_path=os.path.join(OUTPUT_DIR, "theta_beta_mach.png"),
    )

    plot_oblique_shock_heatmap(
        save_path=os.path.join(OUTPUT_DIR, "oblique_shock_heatmap.png"),
    )

    print("\nDone.  All plots saved to  output/")
