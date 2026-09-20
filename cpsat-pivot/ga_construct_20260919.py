#!/usr/bin/env python3
"""
Gaussian-approximation (GA) density evolution for polar construction on the
binary-input AWGN channel, as a replacement proxy for the BEC Bhattacharyya
recursion in gl_shear_search_20260919.bha_seq.

Why: bha_seq uses  z -> (2z - z^2, z^2), which is EXACT only on the BEC.  On
AWGN the "-" branch is an upper bound (Z^- <= 2Z - Z^2), so the error compounds
over m levels and the resulting Z-sum is not a probability (ours reached 18.3
at m=8).  GA instead tracks the mean LLR mu of a symmetric Gaussian:

    mu^-  = phi^{-1}( 1 - (1 - phi(mu))^2 )        (check-node / worse)
    mu^+  = 2 mu                                    (variable-node / better)
    mu_0  = 2 / sigma^2,   sigma^2 = 1 / (2 R Eb/N0)

with Chung's approximation of phi (Chung, Richardson, Urbanke 2001):

    phi(x) = exp(-0.4527 x^0.86 + 0.0218)                     0 < x < 10
    phi(x) = sqrt(pi/x) exp(-x/4) (1 - 10/(7x))               x >= 10

Per-channel SC error probability and the union bound on SC block error:

    P_e(i) = Q( sqrt(mu_i / 2) ),      P_block <~ sum_{i in J} P_e(i)

The branch order matches bha_seq: (worse, better), so index conventions agree.

CAVEAT: GA is still an SC metric.  It is a better SC proxy than the BEC
Z-sum on AWGN, but it does not model SCL path survival.  For SCL the
honest route is simulation of the top candidates.
"""
from __future__ import annotations
import math
from functools import lru_cache
from typing import List

_A, _B, _C = 0.4527, 0.86, 0.0218


def phi(x: float) -> float:
    if x <= 0.0:
        return 1.0
    if x < 10.0:
        return math.exp(-_A * (x ** _B) + _C)
    return math.sqrt(math.pi / x) * math.exp(-x / 4.0) * (1.0 - 10.0 / (7.0 * x))


def phi_inv(y: float) -> float:
    """Numerical inverse of phi (phi is strictly decreasing on x>0)."""
    if y >= 1.0:
        return 0.0
    if y <= 0.0:
        return float("inf")
    if y > phi(10.0):                       # closed form on the x<10 branch
        return ((_C - math.log(y)) / _A) ** (1.0 / _B)
    lo, hi = 10.0, 20.0
    while phi(hi) > y:
        hi *= 2.0
        if hi > 1e12:
            return hi
    for _ in range(200):                    # bisection
        mid = 0.5 * (lo + hi)
        if phi(mid) > y:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def qfunc(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2.0))


@lru_cache(maxsize=None)
def ga_mu(m: int, ebn0_db: float, rate: float) -> tuple:
    """Mean LLRs of the 2^m synthetic AWGN channels.  Branch order (worse, better)."""
    sigma2 = 1.0 / (2.0 * rate * 10 ** (ebn0_db / 10.0))
    mu = [2.0 / sigma2]
    for _ in range(m):
        nxt = []
        for x in mu:
            p = phi(x)
            nxt.append(phi_inv(1.0 - (1.0 - p) ** 2))   # worse
            nxt.append(2.0 * x)                          # better
        mu = nxt
    return tuple(mu)


@lru_cache(maxsize=None)
def ga_pe(m: int, ebn0_db: float, rate: float) -> tuple:
    """Per-channel SC error probability P_e(i) = Q(sqrt(mu_i/2))."""
    return tuple(qfunc(math.sqrt(x / 2.0)) for x in ga_mu(m, ebn0_db, rate))


def ga_seq(m: int, ebn0_db: float, rate: float) -> List[float]:
    """Drop-in replacement for bha_seq: per-channel score, smaller = better."""
    return list(ga_pe(m, ebn0_db, rate))


def block_error(pe, J) -> float:
    """1 - prod(1 - P_e(i)) over the information set: exact under independence."""
    acc = 1.0
    for i in J:
        acc *= (1.0 - pe[i])
    return 1.0 - acc
