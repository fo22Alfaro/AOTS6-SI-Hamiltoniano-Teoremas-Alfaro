#!/usr/bin/env python3
"""PROTOCOLO P1: n mínimo y p-valor condicionado a datos reales.

No inventa p-valores. Si no hay CSV, solo imprime la tabla de potencia.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
from math import erfc


def z_ppf(p: float) -> float:
    """Aprox. Acklam de Φ^{-1} para p en (0,1)."""
    if not 0.0 < p < 1.0:
        raise ValueError("p in (0,1)")
    a = [
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577509590705e02,
        -3.066479806614736e01,
        2.506628277459239e00,
    ]
    b = [
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068155288572e01,
    ]
    c = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e00,
        -2.549732539343734e00,
        4.374664141464968e00,
        2.938163982698783e00,
    ]
    d = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e00,
        3.754408661907416e00,
    ]
    plow = 0.02425
    phigh = 1 - plow
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (
            (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
        )
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(
            (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
        )
    q = p - 0.5
    r = q * q
    return (
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
        / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)
    )


def n_per_group(d: float, alpha: float = 0.05, power: float = 0.8) -> int:
    z_a = z_ppf(1 - alpha / 2)
    z_b = z_ppf(power)
    n = 2.0 * (z_a + z_b) ** 2 / (d**2)
    return int(math.ceil(n))


def two_sample_t_p(x: np.ndarray, y: np.ndarray) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        raise ValueError("n>=2 por grupo")
    mx, my = x.mean(), y.mean()
    vx, vy = x.var(ddof=1), y.var(ddof=1)
    sp2 = ((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2)
    se = math.sqrt(sp2 * (1 / nx + 1 / ny))
    t = (mx - my) / se if se > 0 else 0.0
    # p bilateral vía normal (aprox. n grande); para n chico usar scipy.t
    p = erfc(abs(t) / math.sqrt(2.0))
    d = (mx - my) / math.sqrt(sp2) if sp2 > 0 else 0.0
    return {"n_x": nx, "n_y": ny, "t": t, "p_approx_normal": p, "cohen_d": d}


def bonferroni(alpha: float, k: int) -> float:
    return alpha / k


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv-a", type=Path, default=None, help="columna única, grupo A")
    ap.add_argument("--csv-b", type=Path, default=None, help="columna única, grupo B")
    ap.add_argument("--k-tests", type=int, default=99)
    args = ap.parse_args()

    print("PROTOCOLO P1 AOTS6 — no hay p-valor sin datos")
    print("alpha=0.05  power=0.80  dos colas")
    for d in (0.2, 0.5, 0.8):
        print(f"  Cohen d={d:.1f}  n_por_grupo>={n_per_group(d)}")
    print(f"Bonferroni 99 hallazgos: alpha'={bonferroni(0.05, args.k_tests):.6e}")

    if args.csv_a and args.csv_b:
        a = np.loadtxt(args.csv_a, delimiter=",")
        b = np.loadtxt(args.csv_b, delimiter=",")
        print(two_sample_t_p(a, b))
    else:
        print("p-valor: NO CALCULADO (faltan --csv-a y --csv-b)")


if __name__ == "__main__":
    main()
