"""
Visualization Helpers
======================
All plotting functions used by the examples.  Each function accepts
a ``save_path`` argument; when provided the figure is saved to disk
instead of (or in addition to) being displayed.
"""

import math
import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend — safe for servers/CI
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib import cm

from .shock_relations import NormalShock, ObliqueShock, _theta_from_beta
from .prandtl_meyer import PrandtlMeyer
from .bow_shock import BowShock


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
_STYLE = {
    "figure.dpi": 150,
    "axes.grid": True,
    "grid.alpha": 0.35,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 11,
}
plt.rcParams.update(_STYLE)


# ---------------------------------------------------------------------------
# 1.  Normal shock property sweep
# ---------------------------------------------------------------------------
def plot_normal_shock_properties(
    M_range=(1.5, 25),
    gamma: float = 1.4,
    save_path: str = None,
):
    """Plot p₂/p₁, ρ₂/ρ₁, T₂/T₁, and p₀₂/p₀₁ vs upstream Mach number."""
    M_arr = np.linspace(M_range[0], M_range[1], 500)
    p_r, rho_r, T_r, p0_r = [], [], [], []
    for M in M_arr:
        ns = NormalShock(M, gamma)
        p_r.append(ns.p_ratio)
        rho_r.append(ns.rho_ratio)
        T_r.append(ns.T_ratio)
        p0_r.append(ns.p0_ratio)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), tight_layout=True)
    fig.suptitle(
        f"Normal Shock Properties  (γ = {gamma})", fontsize=14, fontweight="bold"
    )

    datasets = [
        (axes[0, 0], p_r, "p₂/p₁", "Static Pressure Ratio", "#e05c5c"),
        (axes[0, 1], rho_r, "ρ₂/ρ₁", "Density Ratio", "#5c8be0"),
        (axes[1, 0], T_r, "T₂/T₁", "Temperature Ratio", "#e09c5c"),
        (axes[1, 1], p0_r, "p₀₂/p₀₁", "Total Pressure Recovery", "#5cb85c"),
    ]
    for ax, data, ylabel, title, color in datasets:
        ax.plot(M_arr, data, color=color, lw=2)
        ax.axvspan(5, M_range[1], alpha=0.07, color="red", label="Hypersonic (M≥5)")
        ax.set_xlabel("Mach Number M₁")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=9)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# 2.  θ–β–M diagram
# ---------------------------------------------------------------------------
def plot_theta_beta_mach(
    mach_numbers=(2, 3, 5, 8, 10, 15, 20),
    gamma: float = 1.4,
    save_path: str = None,
):
    """Classic θ–β–M diagram for oblique shock waves."""
    fig, ax = plt.subplots(figsize=(10, 7), tight_layout=True)
    ax.set_title("θ–β–M Diagram — Oblique Shock Waves", fontsize=14, fontweight="bold")
    ax.set_xlabel("Shock Wave Angle β (degrees)")
    ax.set_ylabel("Flow Deflection Angle θ (degrees)")

    colors = cm.viridis(np.linspace(0.15, 0.9, len(mach_numbers)))
    for M, col in zip(mach_numbers, colors):
        mu = math.degrees(math.asin(1.0 / M))
        betas = np.linspace(mu + 0.1, 89.9, 800)
        thetas = [_theta_from_beta(M, b, gamma) for b in betas]
        ax.plot(betas, thetas, color=col, lw=1.8, label=f"M = {M}")
        # mark maximum deflection
        idx_max = int(np.argmax(thetas))
        ax.plot(betas[idx_max], thetas[idx_max], "o", color=col, ms=5)

    # Sonic line (M2 = 1 downstream)
    ax.axvline(90, color="k", lw=0.8, ls="--", alpha=0.4)
    ax.legend(title="Free-stream Mach", loc="upper right", fontsize=9)
    ax.set_xlim(0, 90)
    ax.set_ylim(0, None)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# 3.  Bow shock shape around a blunt body
# ---------------------------------------------------------------------------
def plot_bow_shock(
    M1: float,
    nose_radius: float = 1.0,
    gamma: float = 1.4,
    body_type: str = "sphere",
    save_path: str = None,
):
    """Visualise the detached bow shock around a blunt nose."""
    bs = BowShock(M1, nose_radius, gamma, body_type)
    x_shock, y_shock, x_body, y_body = bs.shock_shape(300)

    fig, ax = plt.subplots(figsize=(10, 7), tight_layout=True)
    ax.set_title(
        f"Bow Shock — M∞ = {M1},  δ/R = {bs.delta/nose_radius:.4f}",
        fontsize=13,
        fontweight="bold",
    )

    # Shock (top + bottom mirror)
    ax.plot(x_shock, y_shock, "r-", lw=2.2, label="Bow shock")
    ax.plot(x_shock, -y_shock, "r-", lw=2.2)

    # Body
    ax.fill_betweenx(y_body, x_body, x_body.max() + nose_radius,
                     alpha=0.3, color="steelblue")
    ax.fill_betweenx(-y_body, x_body, x_body.max() + nose_radius,
                     alpha=0.3, color="steelblue")
    ax.plot(x_body, y_body, "b-", lw=2, label="Body surface")
    ax.plot(x_body, -y_body, "b-", lw=2)

    # Free-stream arrows
    x_arr = np.full(5, x_shock.min() - nose_radius)
    y_arr = np.linspace(-3 * nose_radius, 3 * nose_radius, 5)
    for xi, yi in zip(x_arr, y_arr):
        ax.annotate(
            "", xy=(xi + 0.6 * nose_radius, yi),
            xytext=(xi, yi),
            arrowprops=dict(arrowstyle="->", color="gray", lw=1.2),
        )

    ax.text(
        x_arr[2] + 0.3 * nose_radius, y_arr[-1] + 0.2 * nose_radius,
        f"M∞ = {M1}", color="gray", fontsize=10,
    )

    ax.set_xlabel("x / R  (axial)")
    ax.set_ylabel("y / R  (radial)")
    ax.set_aspect("equal")
    ax.legend(loc="upper right")

    # Info box
    info = (
        f"δ = {bs.delta:.4f} m\n"
        f"Rc = {bs.Rc:.4f} m\n"
        f"T₂/T₁ = {bs.ns.T_ratio:.2f}\n"
        f"p₂/p₁ = {bs.ns.p_ratio:.2f}"
    )
    ax.text(
        0.02, 0.97, info,
        transform=ax.transAxes,
        fontsize=9, va="top",
        bbox=dict(boxstyle="round,pad=0.4", fc="lightyellow", alpha=0.8),
    )

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# 4.  Mach cone visualisation
# ---------------------------------------------------------------------------
def plot_mach_cone(
    mach_numbers=(5, 10, 15, 20),
    save_path: str = None,
):
    """Illustrate Mach cones for several hypersonic Mach numbers."""
    fig, axes = plt.subplots(1, len(mach_numbers), figsize=(14, 4), tight_layout=True)
    fig.suptitle("Mach Cones at Hypersonic Speeds", fontsize=13, fontweight="bold")

    L = 5.0  # body length (arbitrary)
    colors = cm.plasma(np.linspace(0.2, 0.85, len(mach_numbers)))

    for ax, M, col in zip(axes, mach_numbers, colors):
        mu = math.degrees(math.asin(1.0 / M))
        half_angle = mu  # Mach cone half-angle

        x = np.array([0, L])
        y_top = np.tan(math.radians(half_angle)) * x
        y_bot = -y_top

        ax.fill_between(x, y_bot, y_top, alpha=0.25, color=col)
        ax.plot(x, y_top, color=col, lw=2)
        ax.plot(x, y_bot, color=col, lw=2)
        ax.plot([0], [0], "ko", ms=5)

        # free-stream arrow
        ax.annotate("", xy=(0, 0), xytext=(-1.5, 0),
                    arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))

        ax.set_title(f"M = {M}\nμ = {mu:.1f}°", fontsize=10)
        ax.set_xlim(-1.8, L + 0.5)
        ax.set_ylim(-2, 2)
        ax.set_aspect("equal")
        ax.axis("off")

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# 5.  Stagnation heating vs Mach (simplified)
# ---------------------------------------------------------------------------
def plot_stagnation_heating(
    M_range=(5, 30),
    T_inf: float = 220.0,
    save_path: str = None,
):
    """Adiabatic wall temperature vs Mach number (simplified flat-plate estimate).

    Uses:  T_aw ≈ T∞ · (1 + r·(γ-1)/2 · M²)
    where r = recovery factor ≈ √Pr ≈ 0.845 for turbulent flow.
    """
    gamma = 1.4
    r = 0.845  # turbulent recovery factor
    M_arr = np.linspace(M_range[0], M_range[1], 300)
    T_aw = T_inf * (1.0 + r * (gamma - 1) / 2.0 * M_arr ** 2)
    T_stag = T_inf * (1.0 + (gamma - 1) / 2.0 * M_arr ** 2)

    fig, ax = plt.subplots(figsize=(10, 5), tight_layout=True)
    ax.plot(M_arr, T_stag, "r-", lw=2, label="Stagnation temp T₀ (K)")
    ax.plot(M_arr, T_aw, "b--", lw=2, label="Adiabatic wall temp Taw (K)")
    ax.axhline(1811, color="orange", ls=":", lw=1.5, label="Steel melting ~1811 K")
    ax.axhline(3695, color="purple", ls=":", lw=1.5, label="Tungsten melting ~3695 K")
    ax.axvspan(5, M_range[1], alpha=0.06, color="red")

    # Artemis marker
    M_artemis = 7.6
    T_art = T_inf * (1.0 + (gamma - 1) / 2.0 * M_artemis ** 2)
    ax.annotate(
        f"Artemis M≈{M_artemis}\nT₀≈{T_art:.0f} K",
        xy=(M_artemis, T_art),
        xytext=(M_artemis + 2, T_art + 500),
        arrowprops=dict(arrowstyle="->", color="black"),
        fontsize=9,
    )

    ax.set_xlabel("Mach Number")
    ax.set_ylabel("Temperature (K)")
    ax.set_title(
        f"Stagnation & Adiabatic Wall Temperatures (T∞ = {T_inf} K)",
        fontweight="bold",
    )
    ax.legend(fontsize=9)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# 6.  Oblique shock grid (pressure ratio heat-map)
# ---------------------------------------------------------------------------
def plot_oblique_shock_heatmap(
    M_range=(5, 25),
    theta_range=(5, 45),
    gamma: float = 1.4,
    save_path: str = None,
):
    """Heat-map of p₂/p₁ over (M1, θ) space for oblique shocks."""
    M_arr = np.linspace(M_range[0], M_range[1], 80)
    theta_arr = np.linspace(theta_range[0], theta_range[1], 80)
    Z = np.full((len(theta_arr), len(M_arr)), np.nan)

    for j, M in enumerate(M_arr):
        for i, th in enumerate(theta_arr):
            try:
                os_ = ObliqueShock(M, theta=th, gamma=gamma)
                Z[i, j] = os_.p_ratio
            except ValueError:
                Z[i, j] = np.nan  # detached shock

    fig, ax = plt.subplots(figsize=(10, 6), tight_layout=True)
    pcm = ax.pcolormesh(M_arr, theta_arr, Z, cmap="hot_r", shading="auto")
    fig.colorbar(pcm, ax=ax, label="p₂/p₁")
    ax.set_xlabel("Free-stream Mach Number M₁")
    ax.set_ylabel("Deflection Angle θ (degrees)")
    ax.set_title(
        "Oblique Shock Pressure Ratio p₂/p₁ (γ = 1.4)",
        fontweight="bold",
    )
    ax.contour(M_arr, theta_arr, Z, levels=10, colors="white", linewidths=0.5, alpha=0.5)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig
