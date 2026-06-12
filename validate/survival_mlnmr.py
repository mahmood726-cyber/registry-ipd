#!/usr/bin/env python
"""EXPERIMENTAL — no external R/flexsurv oracle; validated only against its own Python output (an independent scipy.optimize fit of the same Poisson log-likelihood in harvest/test_phase3c_step5.py); not for primary published estimates without independent confirmation.

A native survival-likelihood ML-NMR engine: piecewise-exponential (Poisson) network meta-regression.

This is the one synthesis engine the portfolio did NOT already have (SYNTHESIS-VISION.md Sec 7): a ML-NMR
with a *native survival likelihood* rather than the RMST-as-Gaussian route of Phase 3c step 3. It fits a
piecewise-exponential (PWE) proportional-hazards model with study x interval baselines, treatment effects,
a prognostic covariate, and a treatment-by-covariate EFFECT-MODIFIER interaction, by maximum Poisson
likelihood (the standard PWE <-> Poisson equivalence: events_si ~ Poisson(exp(eta) * person_time_si),
offset = log person_time). IPD trials enter as fine patient-interval rows; AD (reconstructed-curve) trials
enter as per-arm per-interval aggregate rows -- both are just (study, interval, treatment, x, events,
person_time) tuples at different granularity, so a curve-only survival trial joins natively.

Linear predictor for a row (study j, interval s, treatment k, covariate x):
    eta = mu_{j,s} + d_k * 1[k != ref] + gamma_k * x * 1[k != ref]
where {mu_{j,s}} are nuisance study-interval baselines, {d_k} treatment effects vs reference, and {gamma_k}
the effect-modifier interactions (the estimands). gamma_k lives in the within-study treatment contrast
d_k + gamma_k*x, so a prognostic covariate effect on the baseline hazard is absorbed by mu_{j,s} (shared by
both arms) and need not be modelled separately. Fit by IRLS; covariance = (X' W X)^-1.

Self-contained (numpy only). Cross-validated against an independent scipy.optimize fit of the same Poisson
log-likelihood in harvest/test_phase3c_step5.py.
"""
import math

import numpy as np


def build_design(rows, reference, include_interaction=True):
    """rows: list of dicts {study, interval, treatment, x, events, person_time}.
    Returns (X, y, offset, colnames, index) where index maps parameter roles to column indices.
    Columns: one baseline per (study, interval); d_k for each non-reference treatment; and (if
    include_interaction) gamma_k for each. With include_interaction=False the model is a PH pool that
    ignores the effect modifier -- the comparator that proves the modifier matters.
    """
    studies_intervals = sorted({(r["study"], r["interval"]) for r in rows})
    treatments = sorted({r["treatment"] for r in rows})
    nonref = [t for t in treatments if t != reference]
    si_index = {si: i for i, si in enumerate(studies_intervals)}
    n_base = len(studies_intervals)
    k = len(nonref)
    trt_col = {t: n_base + i for i, t in enumerate(nonref)}            # d_k
    int_col = {t: n_base + k + i for i, t in enumerate(nonref)} if include_interaction else {}
    ncol = n_base + (2 * k if include_interaction else k)

    X = np.zeros((len(rows), ncol), dtype=float)
    y = np.zeros(len(rows), dtype=float)
    offset = np.zeros(len(rows), dtype=float)
    for i, r in enumerate(rows):
        X[i, si_index[(r["study"], r["interval"])]] = 1.0
        if r["treatment"] != reference:
            X[i, trt_col[r["treatment"]]] = 1.0
            if include_interaction:
                X[i, int_col[r["treatment"]]] = r["x"]
        y[i] = r["events"]
        offset[i] = math.log(max(r["person_time"], 1e-12))
    colnames = ([f"mu[{a}|{b}]" for (a, b) in studies_intervals] + [f"d[{t}]" for t in nonref]
                + ([f"gamma[{t}]" for t in nonref] if include_interaction else []))
    index = {"d": trt_col, "gamma": int_col, "nonref": nonref, "n_base": n_base,
             "has_interaction": include_interaction}
    return X, y, offset, colnames, index


def fit_poisson_irls(X, y, offset, max_iter=100, tol=1e-10, ridge=1e-8):
    """Poisson GLM by IRLS. Returns (beta, cov, n_iter, converged). Tiny ridge stabilises the many
    near-collinear nuisance baselines without materially shifting the effect parameters."""
    n, p = X.shape
    b = np.zeros(p, dtype=float)
    R = ridge * np.eye(p)
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        eta = X @ b + offset
        eta = np.clip(eta, -30, 30)
        mu = np.exp(eta)
        W = mu                                   # Poisson IRLS weights = mu
        # working response z = eta - offset + (y - mu)/mu
        z = (eta - offset) + (y - mu) / np.maximum(mu, 1e-12)
        XtW = X.T * W
        A = XtW @ X + R
        rhs = XtW @ z
        b_new = np.linalg.solve(A, rhs)
        if np.max(np.abs(b_new - b)) < tol:
            b = b_new
            converged = True
            break
        b = b_new
    eta = np.clip(X @ b + offset, -30, 30)
    mu = np.exp(eta)
    info = (X.T * mu) @ X + R
    cov = np.linalg.inv(info)
    return b, cov, it, converged


def neg_loglik(b, X, y, offset):
    """Poisson negative log-likelihood (up to constant), for independent cross-validation."""
    eta = np.clip(X @ b + offset, -30, 30)
    mu = np.exp(eta)
    return float(np.sum(mu - y * eta))


class SurvivalMLNMR:
    """Native survival-likelihood ML-NMR (piecewise-exponential Poisson) with an effect modifier."""

    def __init__(self, reference):
        self.reference = reference

    def fit(self, rows, include_interaction=True):
        X, y, offset, colnames, index = build_design(rows, self.reference, include_interaction)
        b, cov, n_iter, converged = fit_poisson_irls(X, y, offset)
        return SurvivalMLNMRFit(b, cov, colnames, index, X, y, offset, n_iter, converged)


class SurvivalMLNMRFit:
    def __init__(self, b, cov, colnames, index, X, y, offset, n_iter, converged):
        self.b, self.cov, self.colnames, self.index = b, cov, colnames, index
        self.X, self.y, self.offset, self.n_iter, self.converged = X, y, offset, n_iter, converged

    def coef(self, role, treatment=None):
        if role == "beta":
            return float(self.b[self.index["beta"]])
        return float(self.b[self.index[role][treatment]])

    def se(self, role, treatment=None):
        j = self.index["beta"] if role == "beta" else self.index[role][treatment]
        return float(math.sqrt(max(self.cov[j, j], 0.0)))

    def gamma(self, treatment):
        """Effect-modifier interaction for a treatment (vs reference)."""
        return self.coef("gamma", treatment), self.se("gamma", treatment)

    def treatment_logHR_at(self, treatment, x):
        """Population-adjusted log-HR for treatment vs reference at covariate value x: d_k + gamma_k * x.
        If the model has no interaction term (modifier-ignorant PH pool), returns the flat d_k."""
        j_d = self.index["d"][treatment]
        c = np.zeros_like(self.b)
        c[j_d] = 1.0
        eff = self.b[j_d]
        if self.index["has_interaction"]:
            j_g = self.index["gamma"][treatment]
            eff = eff + self.b[j_g] * x
            c[j_g] = x
        var = float(c @ self.cov @ c)
        return float(eff), math.sqrt(max(var, 0.0))
