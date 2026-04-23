#!/usr/bin/env python
"""
Query S-parameters from a .s2p file at a target frequency.
Supports DB, MA, and RI formats (auto-detected from # option line).
Usage: python s2p_query.py <file.s2p> <freq_GHz>
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
from utils import s_params
from helpers import db


def _parse_s2p(path):
    freqs, data, fmt = [], [], 'DB'
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('!'):
                continue
            if line.startswith('#'):
                tokens = line.upper().split()
                for f_tok in ('DB', 'MA', 'RI'):
                    if f_tok in tokens:
                        fmt = f_tok
                continue
            vals = list(map(float, line.split()))
            freqs.append(vals[0])
            data.append(vals[1:])
    return np.array(freqs), np.array(data), fmt


def _to_complex(a, b, fmt):
    if fmt == 'DB':
        return 10**(a / 20.0) * np.exp(1j * np.radians(b))
    if fmt == 'MA':
        return a * np.exp(1j * np.radians(b))
    return a + 1j * b  # RI


def load_s2p(path, freq_hz):
    """Load nearest S-params to freq_hz from a Touchstone .s2p file."""
    freqs, data, fmt = _parse_s2p(path)
    idx = int(np.argmin(np.abs(freqs - freq_hz)))
    r = data[idx]
    # Touchstone column order: S11, S21, S12, S22
    s11 = _to_complex(r[0], r[1], fmt)
    s21 = _to_complex(r[2], r[3], fmt)
    s12 = _to_complex(r[4], r[5], fmt)
    s22 = _to_complex(r[6], r[7], fmt)
    return s_params(s11, s12, s21, s22), freqs[idx]


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: s2p_query.py <file.s2p> <freq_GHz>")
        sys.exit(1)
    path      = sys.argv[1]
    freq_ghz  = float(sys.argv[2])
    sp, actual = load_s2p(path, freq_ghz * 1e9)

    print(f"Nearest: {actual/1e9:.4f} GHz  (requested {freq_ghz:.3f} GHz)")
    print(f"  S11 = {np.abs(sp.s11):.4f} ∠ {np.angle(sp.s11, deg=True):.2f}°")
    print(f"  S12 = {np.abs(sp.s12):.4f} ∠ {np.angle(sp.s12, deg=True):.2f}°")
    print(f"  S21 = {np.abs(sp.s21):.4f} ∠ {np.angle(sp.s21, deg=True):.2f}°")
    print(f"  S22 = {np.abs(sp.s22):.4f} ∠ {np.angle(sp.s22, deg=True):.2f}°")

    k, det = sp.rollet()
    mu1, mu2 = sp.mu()
    uc = (mu1 > 1) or (mu2 > 1)
    print(f"\n  K={k:.3f}  |Δ|={np.abs(det):.4f}  μ1={mu1:.3f}  μ2={mu2:.3f}")
    print(f"  {'UNCONDITIONALLY STABLE' if uc else 'POTENTIALLY UNSTABLE'}")
    print(f"  {'MAG' if k >= 1 else 'MSG'} = {db(sp.mag()):.2f} dB")

    if uc:
        gs, gl   = sp.conjugate_match()
        g_in     = sp.gamma_in(gl)
        g_out    = sp.gamma_out(gs)
        zin      = sp.z_in(gl)
        zout     = sp.z_out(gs)
        gt       = sp.gain_t(gs, gl, gamma_out=g_out, gamma_in=g_in)
        print(f"\n  Γ_S_opt = {np.abs(gs):.4f} ∠ {np.angle(gs, deg=True):.1f}°")
        print(f"  Γ_L_opt = {np.abs(gl):.4f} ∠ {np.angle(gl, deg=True):.1f}°")
        print(f"  Z_in    = {zin.real:.3f} + j{zin.imag:.3f} Ω")
        print(f"  Z_out   = {zout.real:.3f} + j{zout.imag:.3f} Ω")
        print(f"  G_T     = {db(gt):.2f} dB")
    else:
        C_L, r_L, C_S, r_S = sp.stability_circles()
        print(f"\n  Output stability circle: center={C_L:.4f}  r={r_L:.4f}")
        print(f"  Input  stability circle: center={C_S:.4f}  r={r_S:.4f}")
        print(f"  Stable origin? output={'Yes' if np.abs(sp.s11)<1 else 'No'}  input={'Yes' if np.abs(sp.s22)<1 else 'No'}")
