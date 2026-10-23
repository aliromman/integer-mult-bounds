"""Exact sharp-header counts without enumerating compressed matrices.

Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
The field order q is a supplied prime power.
"""
from fractions import Fraction
from math import prod


def need(ok, message):
    if not ok:
        raise ValueError(message)


def injection(q, n, k):
    if k < 0 or n < k:
        return 0
    return prod(q**n - q**i for i in range(k))


def gl(q, n):
    return injection(q, n, n)


def gaussian(q, n, k):
    if k < 0 or k > n:
        return 0
    value, rem = divmod(injection(q, n, k), gl(q, k))
    need(rem == 0, 'integer Gaussian coefficient')
    return value


def joint_nullities(q, d):
    need(type(q) is int and q >= 2 and type(d) is int and d >= 0, 'field order and dimension')
    g = [[gaussian(q, n, k) for k in range(n + 1)] for n in range(d + 1)]
    count = {}
    for a in range(d, -1, -1):
        for b in range(d - a, -1, -1):
            incidence = g[d][a] * g[d - a][b] * q**(a*b + d*(d-a-b))
            larger = sum(g[i][a] * g[j][b] * count[i, j]
                         for i in range(a, d + 1)
                         for j in range(b, d - i + 1) if (i, j) != (a, b))
            count[a, b] = incidence - larger
            need(count[a, b] >= 0, 'nonnegative exact stratum')
    need(sum(count.values()) == q**(d*d), 'all matrices counted')
    return count


def unrestricted(q, r, d, rho):
    return sum(gaussian(q, d-rho, k-rho) * injection(q, r, k) * q**(d*(r-k))
               for k in range(rho, min(d, r) + 1))


def full_rank(q, r, d, rho):
    if r < d:
        return 0
    return injection(q, r, d) * q**(rho*(r-d)) * injection(q, r-d, d-rho)


def headers(q, d, profile):
    need(type(d) is int and d >= 1 and len(profile) == 4 and
         all(type(r) is int and r >= 0 for r in profile), 'header dimensions')
    n = sum(profile)
    need(n >= d, 'pivot fits ambient space')
    strata = joint_nullities(q, d)
    arbitrary_by_rank = {r: [unrestricted(q, r, d, rho) for rho in range(d+1)] for r in profile}
    full_by_rank = {r: [full_rank(q, r, d, rho) for rho in range(d+1)] for r in profile[1:]}
    sharp = common = 0
    for (a, b), multiplicity in strata.items():
        ranks = (d-b, d-b, d-b, d-a)
        arbitrary = [arbitrary_by_rank[r][rank] for r, rank in zip(profile, ranks)]
        good = [full_by_rank[r][rank] for r, rank in zip(profile[1:], ranks[1:])]
        common += multiplicity * prod(arbitrary)
        sharp += multiplicity * arbitrary[0] * prod(good)
    completion = gl(q, n-d)
    total = gl(q, n)
    return dict(sharp=sharp*completion, common=common*completion, total=total,
                sharp_probability=Fraction(sharp*completion, total),
                conditional_probability=Fraction(sharp, common) if common else None)
