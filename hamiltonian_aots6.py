#!/usr/bin/env python3
"""H_AOTS6 discreto en una retícula de T^6 (6-toro).

No es un solver de Yang-Mills ni de P vs NP.
Construye el Laplaciano periódico 1D por eje (producto de Kronecker)
más potencial V_A. Unidades SI si el usuario pasa m_i, Lambda_i, V0 en SI.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

HBAR = 1.054571817e-34  # J s


@dataclass
class Scales:
    masses: np.ndarray  # shape (6,) kg
    lambdas: np.ndarray  # shape (6,) SI mixto; ver DEF A1
    V0: float = 0.0  # J
    V1: float = 0.0  # J
    kappa_A: float = 26.3  # hiperparámetro adimensional


def circular_distance(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    d = np.abs(u - v)
    return np.minimum(d, 1.0 - d)


def d_A(a: np.ndarray, b: np.ndarray) -> float:
    delta = circular_distance(np.asarray(a, float), np.asarray(b, float))
    return float(np.sqrt(np.sum(delta**2)))


def diameter() -> float:
    return float(np.sqrt(1.5))


def potential_VA(theta: np.ndarray, sc: Scales) -> float:
    theta = np.asarray(theta, float) % 1.0
    onsite = sc.V0 * np.sum(1.0 - np.cos(2.0 * np.pi * theta))
    couple = 0.0
    for i, j in itertools.combinations(range(6), 2):
        couple += np.cos(2.0 * np.pi * (theta[i] - theta[j]))
    return float(onsite + sc.kappa_A * sc.V1 * couple)


def bound_VA(sc: Scales) -> float:
    """TEOREMA ALFARO I."""
    return 6.0 * abs(sc.V0) + 15.0 * abs(sc.kappa_A * sc.V1)


def laplacian_circle(N: int, mass: float, lam: float) -> np.ndarray:
    """-ħ²/(2 m Λ²) * ∂²/∂θ² en N puntos, θ=k/N."""
    if N < 3:
        raise ValueError("N>=3")
    dtheta = 1.0 / N
    coef = -(HBAR**2) / (2.0 * mass * lam**2 * dtheta**2)
    # matriz de ∂² periódica: (u_{k+1}-2u_k+u_{k-1})
    eye = np.eye(N)
    plus = np.roll(eye, 1, axis=1)
    minus = np.roll(eye, -1, axis=1)
    D2 = plus + minus - 2.0 * eye
    return coef * D2


def hamiltonian_1axis_demo(N: int = 16, axis: int = 0, sc: Scales | None = None) -> np.ndarray:
    """Demo 1D: un solo eje + potencial onsite. 6D pleno es N^6 y no cabe en RAM."""
    if sc is None:
        sc = Scales(
            masses=np.ones(6),
            lambdas=np.ones(6),
            V0=1e-21,
            V1=0.0,
            kappa_A=26.3,
        )
    Hkin = laplacian_circle(N, float(sc.masses[axis]), float(sc.lambdas[axis]))
    thetas = np.arange(N) / N
    V = np.diag([sc.V0 * (1.0 - np.cos(2.0 * np.pi * th)) for th in thetas])
    return Hkin + V


def spectrum_demo(N: int = 16) -> np.ndarray:
    H = hamiltonian_1axis_demo(N=N)
    evals = np.linalg.eigvalsh(H)
    return np.sort(evals.real)


if __name__ == "__main__":
    print("diametro T^6 (Teorema Alfaro IV) =", diameter())
    sc = Scales(masses=np.ones(6), lambdas=np.ones(6), V0=1.0, V1=0.1, kappa_A=26.3)
    print("cota |V_A| (Teorema Alfaro I) =", bound_VA(sc))
    e = spectrum_demo(12)
    print("E0..E3 demo 1D (J, masas=1, Lambda=1, no físico):", e[:4])
