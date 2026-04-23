# Microwave Design Utilities
Library and utility scripts for maching networks, stability analysis, and microwave amplifier design.

---

## Library (`uwaves-utils/`)

### `utils.py` — `s_params` class

Core two-port S-parameter object. Instantiate with four complex values:

```python
from utils import s_params
sp = s_params(s11, s12, s21, s22)   # all complex
```

| Method | Returns | Notes |
|--------|---------|-------|
| `sp._det()` | Δ = S11·S22 − S12·S21 | |
| `sp.rollet()` | `(K, Δ)` | K > 1 and \|Δ\| < 1 → UC stable |
| `sp.mu()` | `(μ1, μ2)` | Edwards–Sinsky; UC stable if either > 1 |
| `sp.stability_circles()` | `(C_L, r_L, C_S, r_S)` | Centers/radii in Γ-plane |
| `sp.conjugate_match()` | `(Γ_S_opt, Γ_L_opt)` | Only valid when UC stable |
| `sp.mag()` | MAG (linear) | Falls back to MSG when K < 1 |
| `sp.gmsg()` | MSG = \|S21/S12\| (linear) | |
| `sp.gamma_in(gamma_l)` | Γ_in | |
| `sp.gamma_out(gamma_s)` | Γ_out | |
| `sp.z_in(gamma_l, z0=50)` | Z_in (Ω) | |
| `sp.z_out(gamma_s, z0=50)` | Z_out (Ω) | |
| `sp.gain_t(gs, gl)` | G_T (linear) | Call `gamma_in`/`gamma_out` first, or pass results as kwargs |
| `sp.gain_p(gl)` | G_P (linear) | Needs `gamma_in` cached |
| `sp.gain_a(gs)` | G_A (linear) | Needs `gamma_out` cached |

### `helpers.py` — Utility functions

```python
from helpers import polar, db, dbv, gamma2z, z2gamma

polar(mag, deg)        # mag∠deg° → complex
db(x)                  # 10·log10|x|  (power dB)
dbv(x)                 # 20·log10|x|  (voltage dB)
gamma2z(gamma, z0=50)  # Γ → Z
z2gamma(z, z0=50)      # Z → Γ
```

---

## Scripts (`./scripts/`)

### `scripts/amp_summary.py` — Full stability + gain summary

Edit S-params at the bottom, run directly. Prints |Δ|, K, μ1/μ2, stability circles, MSG/MAG, and (if UC stable) the simultaneous conjugate match and G_T_max.

```bash
.venv/bin/python calculators/amp_summary.py
```

---

### `scripts/s2p_query.py` — Load S-params from a .s2p file

Auto-detects DB/MA/RI format. Finds the nearest frequency row and prints a full summary (stability, conjugate match, Z_in, Z_out, G_T).

```bash
.venv/bin/python ./scripts/s2p_query.py <file.s2p> <freq_GHz>

# Example
.venv/bin/python ./scripts/s2p_query.py data/measured/SAV_541_2V_30MA.S2P 2.4
```

Also importable:

```python
from scripts.s2p_query import load_s2p
sp, actual_freq = load_s2p("path/to/file.s2p", 2.4e9)
```

---

### `scripts/gain_calc.py` — G_T / G_P / G_A at arbitrary Γ_S, Γ_L

Edit S-params and Γ_S/Γ_L at the bottom. Importable `gain_summary()` returns `(gt, gp, ga)` in linear.

```bash
.venv/bin/python calculators/scripts/gain_calc.py
```

```python
from scripts.gain_calc import gain_summary
gt, gp, ga = gain_summary(sp, gamma_s, gamma_l)
```

---

### `scripts/lmatch.py` — L-network matching

Matches real Z_s (e.g. 50 Ω) to complex Z_L. Returns up to 4 solutions (2 topologies × 2 sign choices) with L/C values in nH/pF.

- **Topology A** (shunt at source, series at load): requires R_s ≥ R_L
- **Topology B** (series at source, shunt at load): requires R_s·G_L ≤ 1

```bash
.venv/bin/python ./scripts/lmatch.py   # edit Z_s, Z_L, freq_hz
```

```python
from scripts.lmatch import print_lmatch
print_lmatch(50.0, 20-30j, 2.4e9)

# For matching to transistor input/output:
gl_opt = sp.conjugate_match()[1]
print_lmatch(50.0, sp.z_in(gl_opt), 2.4e9)
```

---

### `scripts/stub_match.py` — Single-stub tuner

Matches Z_L to Z0 via a shunt stub. Always produces two solutions. Outputs stub position `d` and stub length `l` for both short-circuit and open-circuit stubs, in degrees and wavelength fractions.

```bash
.venv/bin/python ./scripts/stub_match.py   # edit Z_L, Z0, freq_hz
```

```python
from scripts.stub_match import print_stub_match
print_stub_match(25-50j, Z0=50.0, freq_hz=2.4e9)
```

---

## Usage Example

```python
from utils import s_params
from helpers import polar, db
from scripts.gain_calc import gain_summary
from scripts.lmatch import print_lmatch
from scripts.stub_match import print_stub_match

# 1. Build S-param object
sp = s_params(polar(0.85,-65), polar(0.04,45), polar(4.2,120), polar(0.45,-35))

# 2. Check stability
k, det = sp.rollet()
mu1, mu2 = sp.mu()

# 3. Simultaneous conjugate match (UC stable only)
gs, gl = sp.conjugate_match()
zin  = sp.z_in(gl)
zout = sp.z_out(gs)

# 4. Design matching networks
print_lmatch(50.0, zin,  freq_hz=2.4e9)   # input  match
print_lmatch(50.0, zout, freq_hz=2.4e9)   # output match
print_stub_match(zin,  freq_hz=2.4e9)
print_stub_match(zout, freq_hz=2.4e9)
```
