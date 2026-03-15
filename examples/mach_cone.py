"""
Mach Cone & Shock Angle Visualisation
=======================================
Illustrates Mach cones and bow-shock shapes for a range of
hypersonic Mach numbers.

Usage
-----
    python examples/mach_cone.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import matplotlib
matplotlib.use("Agg")

from src.visualization import plot_mach_cone, plot_bow_shock

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

if __name__ == "__main__":
    # Mach cones
    fig1 = plot_mach_cone(
        mach_numbers=(5, 10, 15, 20),
        save_path=os.path.join(OUTPUT_DIR, "mach_cone.png"),
    )
    print(f"[Saved] {os.path.join(OUTPUT_DIR, 'mach_cone.png')}")

    # Bow shock at Mach 7.6 (Artemis peak heating)
    fig2 = plot_bow_shock(
        M1=7.6,
        nose_radius=1.0,
        save_path=os.path.join(OUTPUT_DIR, "bow_shock_M7_6.png"),
    )
    print(f"[Saved] {os.path.join(OUTPUT_DIR, 'bow_shock_M7_6.png')}")

    # Bow shock at Mach 25 (deep hypersonic)
    fig3 = plot_bow_shock(
        M1=25.0,
        nose_radius=1.0,
        save_path=os.path.join(OUTPUT_DIR, "bow_shock_M25.png"),
    )
    print(f"[Saved] {os.path.join(OUTPUT_DIR, 'bow_shock_M25.png')}")
