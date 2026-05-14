"""
Compute state-to-state rate coefficients for H2O + H2 collisions using
MQCT cross sections, for a specified symmetry system and initial H2 level.

One output file per (system, j_H2) combination.

Parameters
----------
  --system    Symmetry system (required):
                para_H2O-para_H2 | ortho_H2O-para_H2
                para_H2O-ortho_H2 | ortho_H2O-ortho_H2
  --j_h2     Initial H2 rotational level: para-H2: 8,6,4,2,0 / ortho-H2: 7,5,3,1
  --data_dir  Root directory of MQCT cross section data
                (default: ./MQCT_Cross_Sections)
                Note: unzip the per-system .zip files before running.
                Expected layout: MQCT_Cross_Sections/<system>/<system>/U*.txt
  --out_dir   Output directory for rate coefficient files
                (default: ./MQCT_Rate_Coefficients)
  --A_low     Sub-threshold power-law exponent (default: 1.0)

Usage
-----
  python compute_rates_by_jh2.py --system para_H2O-para_H2 --j_h2 4
  python compute_rates_by_jh2.py --system ortho_H2O-ortho_H2 --j_h2 3 \\
      --data_dir /path/to/MQCT_Cross_Sections --out_dir /path/to/MQCT_Rate_Coefficients

Output columns
--------------
  J1 KA1 KC1 J2(H2)  J1p KA1p KC1p J2p(H2)  k(T1) k(T2) ...
"""

import glob
import os
import argparse
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

# Physical constants
K_TO_WN = 0.6950345740          # cm^-1 / K
SPEED   = 4.84997439836149e-13  # prefactor (cm^3 s^-1 K^-0.5 wn^2)

TEMPERATURES = [10., 15., 20., 25., 30., 40., 50., 60., 70., 80.,
                100., 150., 200., 300., 400., 500., 800., 1000., 1200.,
                1500., 2000.]

SYSTEMS = [
    "para_H2O-para_H2",
    "ortho_H2O-para_H2",
    "para_H2O-ortho_H2",
    "ortho_H2O-ortho_H2",
]

# Column indices (0-based after line.split())
_J_I, _KA_I, _KC_I = 0, 1, 2
_J_F, _KA_F, _KC_F = 3, 4, 5
_H2_I, _H2_F       = 6, 7
_U                  = 12
_SIGU_QU            = 13
_SIGU_EX            = 14
_DE                 = 19


def load_data(data_dir, j_h2_filter):
    """
    Read all U*.txt files and collect transitions where H2 initial level == j_h2_filter.
    """
    u_files = sorted(glob.glob(os.path.join(data_dir, "U*.txt")))
    if not u_files:
        raise FileNotFoundError(f"No U*.txt files in {data_dir}")

    data, dE_map, trans_order = {}, {}, None

    for fpath in u_files:
        file_keys = []
        with open(fpath) as f:
            f.readline()  # skip header
            for line in f:
                if not line.strip():
                    continue
                c = line.split()
                if int(c[_H2_I]) != j_h2_filter:
                    continue
                key = (int(c[_J_I]), int(c[_KA_I]), int(c[_KC_I]), int(c[_H2_I]),
                       int(c[_J_F]), int(c[_KA_F]), int(c[_KC_F]), int(c[_H2_F]))
                row = [float(c[_U]), float(c[_SIGU_QU]), float(c[_SIGU_EX])]
                if key not in data:
                    data[key]   = []
                    dE_map[key] = float(c[_DE])
                data[key].append(row)
                file_keys.append(key)
        if trans_order is None:
            trans_order = file_keys

    for key in data:
        arr = np.array(data[key])
        data[key] = arr[arr[:, 0].argsort()]

    return trans_order, data, dE_map


def build_sigma_func(U_raw, s_raw, U0, A_low=1.0):
    mask = s_raw > 1e-25
    Uf, Sf = U_raw[mask], s_raw[mask]
    if len(Uf) < 2:
        return lambda u: 0.0
    try:
        b = -(Uf[-2] - Uf[-1]) / np.log(Sf[-2] / Sf[-1])
        a = Sf[-1] * np.exp(Uf[-1] / b)
    except Exception:
        a, b = Sf[-1], 500.0
    kind   = 'cubic' if len(Uf) >= 4 else 'linear'
    interp = interp1d(np.log10(Uf), np.log10(Sf), kind=kind, fill_value='extrapolate')
    r_u, r_s = Uf[0], Sf[0]

    def fn(u):
        if u <= U0:
            return 1e-50
        elif u < r_u:
            return r_s * ((u - U0) / (r_u - U0)) ** A_low
        elif u <= Uf[-2]:
            return 10.0 ** interp(np.log10(u))
        else:
            return a * np.exp(-u / b)
    return fn


def rate_at_T(arr, dE, g_i, g_f, T, A_low=1.0):
    """
    Compute rate coefficient k(T) for a transition
    """
    U0 = abs(dE) / 4.0
    kT = K_TO_WN * T

    if dE <= 0.0:
        s_fwd, s_rev = arr[:, 1], arr[:, 2]
    else:
        s_fwd, s_rev = arr[:, 2], arr[:, 1]

    f_fwd = build_sigma_func(arr[:, 0], s_fwd, U0, A_low)
    f_rev = build_sigma_func(arr[:, 0], s_rev, U0, A_low)
    sigma_avg = lambda u: (g_i * f_fwd(u) + g_f * f_rev(u)) / 2.0

    pref = SPEED * np.sqrt(T) / kT**2 * np.exp(-dE / (2.0 * kT))

    def integrand(u):
        eps = abs(dE) / (4.0 * u)
        return sigma_avg(u) * u * (1.0 - eps**2) * np.exp(-u / kT * (1.0 + eps**2))

    val, _ = quad(integrand, U0, 15000.0, limit=200)
    k = pref * val
    if not np.isfinite(k) or k < 0.0 or k > 1.0e-7:
        return 0.0
    return k


def main():
    parser = argparse.ArgumentParser(
        description="Compute H2O + H2 state-to-state rate coefficients from MQCT cross sections."
    )
    parser.add_argument("--system", required=True, choices=SYSTEMS,
                        help="Symmetry system (e.g. para_H2O-ortho_H2)")
    parser.add_argument("--j_h2", type=int, required=True,
                        help="Initial H2 rotational level (para-H2: even 8,6,4,2,0; ortho-H2: odd 7,5,3,1)")
    parser.add_argument("--data_dir", default=None,
                        help="Root directory containing MQCT cross section data (default: ./MQCT_Cross_Sections)")
    parser.add_argument("--out_dir", default=None,
                        help="Output directory for rate files (default: ./MQCT_Rate_Coefficients)")
    parser.add_argument("--A_low", type=float, default=1.0,
                        help="Sub-threshold power-law exponent (default 1.0)")
    args = parser.parse_args()

    # Derive H2O and H2 symmetry labels from system name
    h2o_sym, h2_sym = args.system.split("-")        # e.g. "para_H2O", "ortho_H2"
    h2o_tag = "pH2O" if h2o_sym.startswith("para") else "oH2O"
    h2_tag  = "pH2"  if h2_sym.startswith("para")  else "oH2"

    # Warn if j_H2 parity doesn't match the symmetry species
    is_para_h2 = h2_sym.startswith("para")
    if is_para_h2 and args.j_h2 % 2 != 0:
        print(f"Warning: para-H2 levels are even (0,2,4,...), got j_H2={args.j_h2}", flush=True)
    elif not is_para_h2 and args.j_h2 % 2 == 0:
        print(f"Warning: ortho-H2 levels are odd (1,3,5,...), got j_H2={args.j_h2}", flush=True)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_base  = args.data_dir or os.path.join(script_dir, "MQCT_Cross_Sections")
    out_base   = args.out_dir  or os.path.join(script_dir, "MQCT_Rate_Coefficients")

    data_dir = os.path.join(data_base, args.system, args.system)
    os.makedirs(out_base, exist_ok=True) 

    out_path = os.path.join(out_base, f"rates_{h2o_tag}_{h2_tag}_j{args.j_h2}.txt")

    trans_order, data, dE_map = load_data(data_dir, args.j_h2)

    seen, ordered_keys = set(), []
    for k in trans_order:
        if k not in seen:
            seen.add(k)
            ordered_keys.append(k)

    n_trans = len(ordered_keys)
    n_temps = len(TEMPERATURES)
    print(f"{args.system}  j_H2={args.j_h2} : {n_trans} transitions × {n_temps} temperatures → {out_path}", flush=True)

    results = []
    for key in ordered_keys:
        j_i, ka_i, kc_i, h2_i, j_f, ka_f, kc_f, h2_f = key
        g_i = 2 * j_i + 1
        g_f = 2 * j_f + 1
        ks = []
        for T in TEMPERATURES:
            try:
                k = rate_at_T(data[key], dE_map[key], g_i, g_f, T, args.A_low)
                ks.append(max(k, 0.0))
            except Exception:
                ks.append(0.0)
        results.append((key, ks))

    QW, RW = 6, 15
    qn_hdr   = (f"{'J1':>{QW}}{'KA1':>{QW}}{'KC1':>{QW}}{'J2':>{QW}}"
                f"{'J1p':>{QW}}{'KA1p':>{QW}}{'KC1p':>{QW}}{'J2p':>{QW}}")
    temp_hdr = "".join(f"{'T_' + str(T):>{RW}}" for T in TEMPERATURES)

    with open(out_path, "w") as out:
        out.write(f"number of temperatures :\n{n_temps}\n")
        out.write(f"number of transitions :\n{n_trans}\n")
        out.write(f"system : {args.system}  j_H2_init={args.j_h2}\n")
        out.write("order : J1 KA1 KC1 J2(H2)  J1p KA1p KC1p J2p(H2)\n")
        out.write(qn_hdr + temp_hdr + "\n")
        for (j_i, ka_i, kc_i, h2_i, j_f, ka_f, kc_f, h2_f), ks in results:
            qn  = (f"{j_i:{QW}d}{ka_i:{QW}d}{kc_i:{QW}d}{h2_i:{QW}d}"
                   f"{j_f:{QW}d}{ka_f:{QW}d}{kc_f:{QW}d}{h2_f:{QW}d}")
            rts = "".join(f"{k:{RW}.6E}" for k in ks)
            out.write(qn + rts + "\n")

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
