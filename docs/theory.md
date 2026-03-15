# Hypersonic Shock Waves — Theory Reference

## 1. Governing Equations

### 1.1 Normal Shock — Rankine-Hugoniot Relations

For a calorically perfect gas with ratio of specific heats γ, the Rankine-Hugoniot
conditions across a normal shock with upstream Mach number M₁ are:

| Quantity | Formula |
|---|---|
| M₂ | √[ ((γ−1)M₁² + 2) / (2γM₁² − (γ−1)) ] |
| p₂/p₁ | (2γM₁² − (γ−1)) / (γ+1) |
| ρ₂/ρ₁ | (γ+1)M₁² / ((γ−1)M₁² + 2) |
| T₂/T₁ | (p₂/p₁) / (ρ₂/ρ₁) |
| p₀₂/p₀₁ | [ρ₂/ρ₁]^(γ/(γ−1)) · [p₁/p₂]^(1/(γ−1)) |

**Hypersonic limit (M₁ → ∞, γ = 1.4):**

- M₂ → √((γ−1)/(2γ)) ≈ 0.378
- ρ₂/ρ₁ → (γ+1)/(γ−1) = 6
- T₂/T₁ → 2γ(γ−1)/(γ+1)² · M₁² (grows as M²)

---

### 1.2 Oblique Shock — θ–β–M Relation

The relationship between the flow deflection angle θ, shock wave angle β, and
upstream Mach number M₁:

```
tan θ = 2 cot β · [M₁² sin²β − 1] / [M₁²(γ + cos 2β) + 2]
```

The oblique shock is treated as a normal shock in the direction normal to the wave.
The normal component of the upstream Mach number is:

```
M₁ₙ = M₁ sin β
```

All Rankine-Hugoniot ratios are computed using M₁ₙ.  
The downstream Mach number is:

```
M₂ = M₂ₙ / sin(β − θ)
```

**Key features:**
- For M₁ > 1 there is a maximum deflection angle θ_max beyond which the shock detaches.
- Two solutions exist for θ < θ_max: *weak shock* (smaller β) and *strong shock* (larger β).
- In most practical aerodynamic applications the weak-shock solution is observed.

---

### 1.3 Prandtl-Meyer Expansion Fan

When a supersonic flow turns away from itself around a convex corner, an isentropic
expansion fan forms.  The Prandtl-Meyer function is:

```
ν(M) = √((γ+1)/(γ−1)) · arctan(√((γ−1)/(γ+1)(M²−1))) − arctan(√(M²−1))
```

For an expansion through turning angle Δθ:

```
ν₂ = ν₁ + Δθ
```

M₂ is found by inverting ν(M₂) = ν₂.  
All downstream properties follow from isentropic relations.

---

### 1.4 Bow Shock Standoff Distance (Billig Correlation)

For a blunt body at hypersonic speeds, the bow shock detaches and stands off the
nose at a distance δ.  The Billig (1967) empirical correlation gives:

```
δ/R = A · exp(B / M∞²)
Rc/R = C · M∞^D
```

where R is the nose radius of curvature and the constants (A, B, C, D) = (0.386, 4.67,
0.0455, 0.6) for a sphere.  The shock shape is approximated as a hyperbola:

```
x = δ + Rc − √(Rc² + y²/ε²)
```

with eccentricity ε ≈ 1 + 0.8 / √(M∞ − 1).

---

## 2. Stagnation (Total) Conditions

For isentropic flow ahead of the shock:

```
T₀/T = 1 + (γ−1)/2 · M²
p₀/p = [1 + (γ−1)/2 · M²]^(γ/(γ−1))
```

The **Rayleigh Pitot-tube formula** gives the stagnation pressure behind a normal
shock relative to the free-stream static pressure:

```
p₀₂/p₁ = [(γ+1)²M₁²/(4γM₁²−2(γ−1))]^(γ/(γ−1)) · (1−γ+2γM₁²)/(γ+1)
```

---

## 3. Aerodynamic Heating

A simplified estimate of the adiabatic wall temperature for a flat plate in turbulent
flow is:

```
T_aw = T∞ · (1 + r · (γ−1)/2 · M²)
```

where r ≈ √Pr ≈ 0.845 is the turbulent recovery factor for air.

The stagnation-point heat flux is often estimated via the **Fay-Riddell** formula
(not implemented here — requires full real-gas thermochemistry).

---

## 4. Applicability to the NASA Artemis Programme

The Orion capsule (used on Artemis I–IV) re-enters the Earth's atmosphere at
approximately Mach 32 (~11 km/s) after a translunar return trajectory.  The peak
aerothermal environment occurs during the *skip* and *constant-g* phases:

| Phase | Mach | Altitude (km) | T₀ approx |
|---|---|---|---|
| Entry interface | ~32 | 120 | >10 000 K |
| Peak heating | ~10 | 75 | ~3 200 K |
| Peak dynamic pressure | ~7.6 | 65 | ~1 800 K |
| Subsonic parachute deploy | <0.8 | ~8 | — |

At these conditions real-gas effects (dissociation, ionisation, vibration) become
important and the calorically-perfect-gas assumption overestimates the shock
temperature by a factor of ~2–3.  Nonetheless the perfect-gas model gives correct
trends and is the standard first-principles teaching tool.

---

## 5. References

1. Anderson, J. D. (2003). *Modern Compressible Flow*, 3rd ed. McGraw-Hill.
2. Anderson, J. D. (2006). *Hypersonic and High Temperature Gas Dynamics*, 2nd ed. AIAA.
3. Billig, F. S. (1967). Shock-wave shapes around spherical- and cylindrical-nosed bodies.
   *Journal of Spacecraft and Rockets*, 4(6), 822–823.
4. NASA SP-8077 — Aerodynamic Design Data Book.
5. NACA TN 1428 — Equations, tables and charts for compressible flow.
