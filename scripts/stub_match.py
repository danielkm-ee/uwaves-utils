#!/usr/bin/env python
"""
Single-stub tuner: find stub position d and stub length l to match Z_L to Z0.
Outputs two solutions in electrical degrees and fraction of wavelength.
Both short-circuit and open-circuit stub options are given for each solution.

Usage: python stub_match.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np


def _pos_deg(beta):
    """Map arctan result to (0°, 180°)."""
    if beta <= 0:
        beta += np.pi
    return np.degrees(beta)


def stub_match(Z_L, Z0=50.0):
    """
    Single-stub tuner: match Z_L to Z0.

    Returns list of solution dicts with keys:
      d_deg, d_wl          — line section from load to stub junction
      b_d                  — residual normalized susceptance at junction
      l_short_deg/_wl      — short-circuit stub length
      l_open_deg/_wl       — open-circuit stub length
    """
    y_L = Z0 / complex(Z_L)   # normalized load admittance
    g_L = y_L.real
    b_L = y_L.imag

    if g_L < 1e-10:
        print("Warning: load has zero conductance. No shunt-stub solution.")
        return []

    y_L_sq  = g_L**2 + b_L**2
    a_coef  = g_L - y_L_sq
    b_coef  = 2.0 * b_L
    c_coef  = g_L - 1.0

    if abs(a_coef) < 1e-10:
        if abs(b_coef) < 1e-10:
            print("Load is already matched (y_L = 1 + j0). No stub needed.")
            return []
        t_vals = [-c_coef / b_coef]
    else:
        disc = b_coef**2 - 4.0 * a_coef * c_coef
        if disc < 0:
            print(f"No real solution (discriminant = {disc:.4g}).")
            return []
        sq = np.sqrt(disc)
        t_vals = [(-b_coef + sq) / (2.0 * a_coef),
                  (-b_coef - sq) / (2.0 * a_coef)]

    solutions = []
    for t in t_vals:
        beta_d = np.arctan(t)
        if beta_d < 0:
            beta_d += np.pi
        d_deg = np.degrees(beta_d)
        d_wl  = beta_d / (2.0 * np.pi)

        # Residual susceptance at junction: Im{y(d)} given Re{y(d)}=1
        # b_d = [b_L + t*(1 - |y_L|^2) - b_L*t^2] / [g_L*(1+t^2)]
        b_d = (b_L + t * (1.0 - y_L_sq) - b_L * t**2) / (g_L * (1.0 + t**2))

        # Short-circuit stub: b_sc = -cot(β*l) = b_stub = -b_d
        # → cot(β*l) = b_d  →  tan(β*l) = 1/b_d
        if abs(b_d) < 1e-10:
            l_sc_deg = 90.0   # quarter-wave short gives zero susceptance
        else:
            l_sc_deg = _pos_deg(np.arctan(1.0 / b_d))

        # Open-circuit stub: b_oc = tan(β*l) = b_stub = -b_d
        l_oc_deg = _pos_deg(np.arctan(-b_d))

        solutions.append({
            'd_deg':       d_deg,
            'd_wl':        d_wl,
            'b_d':         b_d,
            'l_short_deg': l_sc_deg,
            'l_short_wl':  l_sc_deg / 360.0,
            'l_open_deg':  l_oc_deg,
            'l_open_wl':   l_oc_deg / 360.0,
        })
    return solutions


def print_stub_match(Z_L, Z0=50.0, freq_hz=None):
    print(f"\nSingle-stub tuner:  Z_L = {complex(Z_L):.3f} Ω,  Z0 = {Z0:.1f} Ω")
    if freq_hz:
        print(f"  Frequency: {freq_hz/1e9:.3f} GHz")
    solns = stub_match(Z_L, Z0)
    if not solns:
        return
    for i, s in enumerate(solns):
        print(f"\n  Solution {i+1}:  (b_junction = {s['b_d']:.4f})")
        print(f"    d           = {s['d_deg']:.2f}°  ({s['d_wl']:.4f} λ)  — line from load to stub")
        print(f"    Short stub  l = {s['l_short_deg']:.2f}°  ({s['l_short_wl']:.4f} λ)")
        print(f"    Open  stub  l = {s['l_open_deg']:.2f}°  ({s['l_open_wl']:.4f} λ)")
    if freq_hz:
        lam_mm = 3e8 / freq_hz * 1e3
        print(f"\n  λ = {lam_mm:.2f} mm at {freq_hz/1e9:.3f} GHz (free space)")


if __name__ == '__main__':
    freq_hz = 2.4e9
    Z0      = 50.0
    Z_L     = 25 - 50j   # Ω — replace with sp.z_in(gl) or sp.z_out(gs)

    print_stub_match(Z_L, Z0, freq_hz)
