#!/usr/bin/env python
"""
Gain calculator: G_T, G_P, G_A for a two-port given S-params and reflection coefficients.
Usage: edit variables below or import gain_summary() into another script.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
from utils import s_params
from helpers import polar, db


def gain_summary(sp: s_params, gamma_s, gamma_l):
    g_in  = sp.gamma_in(gamma_l)
    g_out = sp.gamma_out(gamma_s)
    gt = sp.gain_t(gamma_s, gamma_l, gamma_out=g_out, gamma_in=g_in)
    gp = sp.gain_p(gamma_l, gamma_in=g_in)
    ga = sp.gain_a(gamma_s, gamma_out=g_out)

    print(f"  Gamma_S   = {gamma_s:.4f}  |Γ_S|={np.abs(gamma_s):.4f}  ∠{np.angle(gamma_s, deg=True):.1f}°")
    print(f"  Gamma_L   = {gamma_l:.4f}  |Γ_L|={np.abs(gamma_l):.4f}  ∠{np.angle(gamma_l, deg=True):.1f}°")
    print(f"  Gamma_in  = {g_in:.4f}   |Γ_in| ={np.abs(g_in):.4f}")
    print(f"  Gamma_out = {g_out:.4f}  |Γ_out|={np.abs(g_out):.4f}")
    print(f"  G_T = {db(gt):.2f} dB  (transducer)")
    print(f"  G_P = {db(gp):.2f} dB  (operating power)")
    print(f"  G_A = {db(ga):.2f} dB  (available)")
    return gt, gp, ga


if __name__ == '__main__':
    s11 = polar(0.85, -65.0)
    s12 = polar(0.04,  45.0)
    s21 = polar(4.20, 120.0)
    s22 = polar(0.45, -35.0)
    sp = s_params(s11, s12, s21, s22)

    gs, gl = sp.conjugate_match()
    print("--- Gain at simultaneous conjugate match ---")
    gain_summary(sp, gs, gl)
