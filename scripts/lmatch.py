#!/usr/bin/env python
"""
L-network matching calculator.
Matches a real source impedance Z_s (typically 50 Ω) to complex load Z_L at a given frequency.

Topology A  (shunt at source, series at load):  requires R_s >= R_L   — step down
Topology B  (series at source, shunt at load):  requires R_s*R_L <= |Z_L|^2 — step up

Each topology has two solutions (sign choice), giving up to 4 total.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np


def _comp_from_X(X, omega, label):
    f = omega / (2 * np.pi)
    if abs(X) < 1e-15:
        return f"{label}: short  (X=0)"
    if X > 0:
        return f"{label}: L = {X/omega*1e9:.3f} nH  (X=+{X:.3f} Ω, inductive)"
    else:
        return f"{label}: C = {-1/(omega*X)*1e12:.3f} pF  (X={X:.3f} Ω, capacitive)"


def _comp_from_B(B, omega, label):
    if abs(B) < 1e-15:
        return f"{label}: open   (B=0)"
    if B > 0:
        return f"{label}: C = {B/omega*1e12:.3f} pF  (B=+{B:.4f} S, capacitive)"
    else:
        return f"{label}: L = {-1/(omega*B)*1e9:.3f} nH  (B={B:.4f} S, inductive)"


def lmatch(Z_s, Z_L, freq_hz):
    """
    Design L-networks matching Z_s to Z_L at freq_hz.
    Z_s must be real (e.g. 50.0). Z_L may be complex.
    Returns (solutions, omega) where each solution is a dict with keys:
      topology, Q (or None), X_se (series reactance), B_sh (shunt susceptance).
    """
    omega = 2 * np.pi * freq_hz
    R_s   = float(np.real(Z_s))
    R_L   = float(np.real(Z_L))
    X_L   = float(np.imag(Z_L))

    solutions = []

    # --- Topology A: shunt at source, series at load ---
    # Derived from:  R_L / [R_L^2 + (X_se + X_L)^2] = 1/R_s
    # → (X_se + X_L)^2 = R_L*(R_s - R_L);  condition: R_s >= R_L
    disc_A = R_L * (R_s - R_L)
    if disc_A >= 0:
        sqrt_A = np.sqrt(disc_A)
        Q      = sqrt_A / R_L          # = sqrt(R_s/R_L - 1)
        for sign in (+1, -1):
            X_se = -X_L + sign * sqrt_A
            B_sh = sign * sqrt_A / (R_s * R_L)
            solutions.append({'topology': 'A (shunt-source / series-load)', 'Q': Q,
                               'X_se': X_se, 'B_sh': B_sh})

    # --- Topology B: series at source, shunt at load ---
    # Derived from:  G_L / [G_L^2 + (B_sh + B_L)^2] = R_s  (admittance domain)
    # → (B_sh + B_L)^2 = G_L*(1 - R_s*G_L)/R_s;  condition: R_s*G_L <= 1
    Y_L   = 1.0 / complex(Z_L)
    G_L   = float(Y_L.real)
    B_L   = float(Y_L.imag)
    disc_B = G_L * (1.0 - R_s * G_L) / R_s if G_L > 0 else -1.0
    if disc_B >= 0:
        sqrt_B = np.sqrt(disc_B)
        for sign in (+1, -1):
            D    = sign * sqrt_B
            B_sh = -B_L + D
            X_se = R_s * D / G_L
            solutions.append({'topology': 'B (series-source / shunt-load)', 'Q': None,
                               'X_se': X_se, 'B_sh': B_sh})

    return solutions, omega


def print_lmatch(Z_s, Z_L, freq_hz):
    solutions, omega = lmatch(Z_s, Z_L, freq_hz)
    print(f"\nL-network:  Z_s = {float(Z_s):.1f} Ω  →  Z_L = {complex(Z_L):.3f} Ω"
          f"  at {freq_hz/1e9:.3f} GHz")
    if not solutions:
        print("  No real solution (check R_s vs R_L and |Z_L|^2/(R_s*R_L) >= 1).")
        return
    for i, s in enumerate(solutions):
        q_str = f"  Q={s['Q']:.3f}" if s['Q'] is not None else ""
        print(f"\n  Solution {i+1} — {s['topology']}{q_str}")
        print(f"    {_comp_from_B(s['B_sh'], omega, 'Shunt')}")
        print(f"    {_comp_from_X(s['X_se'], omega, 'Series')}")


if __name__ == '__main__':
    freq_hz = 2.4e9
    Z_s     = 50.0       # source (real, Ω)
    Z_L     = 20 - 30j   # load — replace with sp.z_in(gl) or sp.z_out(gs)

    print_lmatch(Z_s, Z_L, freq_hz)
