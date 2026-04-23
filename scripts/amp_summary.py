#!/usr/bin/env python
'''
Stability and gain summary calculator for a two-port transistor amplifier.
Usage: edit the S-parameters below and run directly.
'''
import numpy as np
from utils import s_params
from helpers import polar, db

def summarize(sp: s_params, label: str = ""):
    if label:
        print(f"\n{'='*50}")
        print(f"  {label}")
        print(f"{'='*50}")

    det = sp._det()
    k, _ = sp.rollet()
    mu1, mu2 = sp.mu()
    C_L, r_L, C_S, r_S = sp.stability_circles()

    print(f"\n--- Stability ---")
    print(f"  |Delta|  = {np.abs(det):.4f}")
    print(f"  K        = {k:.4f}")
    print(f"  mu1      = {mu1:.4f}")
    print(f"  mu2      = {mu2:.4f}")
    uc = (mu1 > 1) or (mu2 > 1)
    print(f"  Status   : {'UNCONDITIONALLY STABLE' if uc else 'POTENTIALLY UNSTABLE'}")

    print(f"\n--- Stability Circles ---")
    print(f"  Output (Gamma_L): center = {C_L:.4f}  radius = {r_L:.4f}")
    print(f"  Input  (Gamma_S): center = {C_S:.4f}  radius = {r_S:.4f}")
    print(f"  Stable origin?  : {'Yes (|S11|<1)' if np.abs(sp.s11) < 1 else 'No  (|S11|>=1)'} [output circle]")
    print(f"                    {'Yes (|S22|<1)' if np.abs(sp.s22) < 1 else 'No  (|S22|>=1)'} [input circle]")

    print(f"\n--- Gain ---")
    msg = sp.gmsg()
    mag_val = sp.mag()
    print(f"  MSG      = {db(msg):.2f} dB  ({msg:.4f} linear)")
    if k >= 1:
        print(f"  MAG      = {db(mag_val):.2f} dB  ({mag_val:.4f} linear)")
    else:
        print(f"  MAG      : N/A (potentially unstable, K < 1)")
    print(f"  |S21|^2  = {db(np.abs(sp.s21)**2):.2f} dB")

    if uc:
        gs_opt, gl_opt = sp.conjugate_match()
        print(f"\n--- Simultaneous Conjugate Match ---")
        print(f"  Gamma_S_opt = {gs_opt:.4f}  |Gamma_S| = {np.abs(gs_opt):.4f}  ang = {np.angle(gs_opt, deg=True):.1f} deg")
        print(f"  Gamma_L_opt = {gl_opt:.4f}  |Gamma_L| = {np.abs(gl_opt):.4f}  ang = {np.angle(gl_opt, deg=True):.1f} deg")
        gamma_in  = sp.gamma_in(gl_opt)
        gamma_out = sp.gamma_out(gs_opt)
        gt = sp.gain_t(gs_opt, gl_opt)
        print(f"  G_T (max)   = {db(gt):.2f} dB")
    else:
        print(f"\n--- Conjugate Match ---")
        print(f"  Skipped: device is potentially unstable. Choose Gamma_L/S from stable region.")


if __name__ == "__main__":
    s11 = polar(0.85, -65.0)
    s12 = polar(0.04,  45.0)
    s21 = polar(4.20,  120.0)
    s22 = polar(0.45,  -35.0)

    sp = s_params(s11, s12, s21, s22)
    summarize(sp, label="SAV-541+ @ [freq] GHz, Vds=[V] V, Ids=[mA] mA")
