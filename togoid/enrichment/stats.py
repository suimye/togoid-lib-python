"""
Statistics for over-representation (enrichment) analysis.

This module deliberately depends only on the Python standard library. The
hypergeometric survival function is evaluated in log space with ``math.lgamma``
so that the enrichment core of ``togoid`` stays installable without SciPy or
statsmodels, which are heavy dependencies for an ID-conversion library.
"""
from math import lgamma, exp, isfinite
from typing import List, Sequence

__all__ = [
    "log_binomial_coefficient",
    "hypergeometric_pmf",
    "hypergeometric_sf",
    "hypergeometric_pvalue",
    "benjamini_hochberg",
    "fdr",
    "fold_enrichment",
]


def log_binomial_coefficient(n: int, k: int) -> float:
    """
    Return ``log(C(n, k))``.

    Args:
        n: Size of the set to choose from.
        k: Number of elements chosen.

    Returns:
        Natural logarithm of the binomial coefficient, or ``-inf`` when the
        coefficient is zero (i.e. ``k`` outside ``[0, n]``).
    """
    if k < 0 or k > n or n < 0:
        return float("-inf")
    # lgamma(x + 1) == log(x!), and is stable for the large factorials that
    # genome-scale background sets (N ~ 20000) would otherwise overflow.
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def hypergeometric_pmf(k: int, N: int, M: int, n: int) -> float:
    """
    Probability of drawing exactly ``k`` successes.

    The parameter names follow the convention used throughout this package:

    Args:
        k: Observed overlap (successes drawn).
        N: Background size (total population).
        M: Term size (successes in the population).
        n: Query size (number drawn).

    Returns:
        ``P(X = k)`` for ``X ~ Hypergeometric(N, M, n)``.
    """
    log_p = (
        log_binomial_coefficient(M, k)
        + log_binomial_coefficient(N - M, n - k)
        - log_binomial_coefficient(N, n)
    )
    if not isfinite(log_p):
        return 0.0
    return exp(log_p)


def hypergeometric_sf(k: int, N: int, M: int, n: int) -> float:
    """
    Upper-tail probability ``P(X >= k)`` of the hypergeometric distribution.

    This is the one-sided over-representation p-value: the probability of seeing
    an overlap at least as large as the observed one by chance.

    Args:
        k: Observed overlap.
        N: Background size.
        M: Term size (number of background genes annotated with the term).
        n: Query size (number of genes in the cluster/query list).

    Returns:
        P-value in ``[0, 1]``. Returns 1.0 for ``k <= 0`` (an overlap of zero or
        less is always at least as likely as observed).

    Examples:
        5 of a 20-gene query fall in a 50-gene term, out of 2000 genes::

            >>> round(hypergeometric_sf(5, 2000, 50, 20), 8)
            9.316e-05

        An overlap of zero is never surprising::

            >>> hypergeometric_sf(0, 100, 10, 10)
            1.0
    """
    if k <= 0:
        return 1.0
    if N <= 0 or M <= 0 or n <= 0:
        return 1.0

    # The support of the distribution is [max(0, n - (N - M)), min(n, M)].
    upper = min(n, M)
    if k > upper:
        return 0.0

    # Summing from the far tail inwards keeps the smallest terms first, which
    # avoids losing them to rounding when the p-value is extremely small.
    total = 0.0
    for i in range(upper, k - 1, -1):
        total += hypergeometric_pmf(i, N, M, n)

    # Guard against accumulated floating point drift.
    if total > 1.0:
        return 1.0
    if total < 0.0:
        return 0.0
    return total


def benjamini_hochberg(pvalues: Sequence[float]) -> List[float]:
    """
    Benjamini-Hochberg FDR correction.

    Args:
        pvalues: Sequence of raw p-values, in any order.

    Returns:
        List of adjusted p-values in the same order as the input.

    Examples:
        >>> [round(p, 4) for p in benjamini_hochberg([0.01, 0.04, 0.03, 0.005])]
        [0.02, 0.04, 0.04, 0.02]

        An empty input gives an empty result::

            >>> benjamini_hochberg([])
            []
    """
    m = len(pvalues)
    if m == 0:
        return []

    # Sort descending so the cumulative minimum enforces monotonicity in one pass.
    order = sorted(range(m), key=lambda i: pvalues[i], reverse=True)

    adjusted = [0.0] * m
    running_min = 1.0
    for rank_from_top, idx in enumerate(order):
        rank = m - rank_from_top  # 1-based rank in ascending order
        value = pvalues[idx] * m / rank
        running_min = min(running_min, value)
        adjusted[idx] = running_min

    return adjusted


def fold_enrichment(k: int, N: int, M: int, n: int) -> float:
    """
    Observed overlap divided by the overlap expected under independence.

    Args:
        k: Observed overlap.
        N: Background size.
        M: Term size.
        n: Query size.

    Returns:
        Fold enrichment, or 0.0 when the expected overlap is zero.

    Examples:
        The expected overlap is ``n * M / N`` = 10, so 20 observed is 2-fold::

            >>> fold_enrichment(20, 1000, 100, 100)
            2.0
            >>> fold_enrichment(5, 100, 0, 10)
            0.0
    """
    if N <= 0:
        return 0.0
    expected = (n * M) / N
    if expected <= 0:
        return 0.0
    return k / expected


# --------------------------------------------------------------------- #
# Aliases matching the R implementation
# --------------------------------------------------------------------- #
#
# The R package names these after what they produce rather than after the
# distribution function or the method: togoid_hypergeometric_pvalue() and
# togoid_fdr(). Both names are correct, and someone moving between the two
# implementations should not have to learn two vocabularies, so each function is
# reachable under either name. The originals stay the canonical ones.

#: Alias of :func:`hypergeometric_sf`, matching R's ``togoid_hypergeometric_pvalue``.
hypergeometric_pvalue = hypergeometric_sf

#: Alias of :func:`benjamini_hochberg`, matching R's ``togoid_fdr``.
fdr = benjamini_hochberg
