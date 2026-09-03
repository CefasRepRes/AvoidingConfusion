#!/usr/bin/env python3
"""NumPyro transcription of the Spence et al. Section 4 model.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from ml_prediction_results_utilities import load_json_payload

try:
    import jax

    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    import numpyro
    import numpyro.distributions as dist
    from numpyro.infer import MCMC, NUTS
except Exception:
    jax = jnp = numpyro = dist = MCMC = NUTS = None

LOG = logging.getLogger(__name__)

# The published Section 4 analysis used 20 chains, 5,000 warm-up iterations,
# and 5,000 retained iterations per chain.  mc_samples controls retained draws
# per chain here so callers may deliberately request a smaller development run.
NUMPYRO_MCMC_WARMUP = 5000
NCHAINS = 20


def extract_validation_error_model(payload: Any) -> Optional[Dict[str, Any]]:
    """Return the classifier-error block exported by the validation workflow."""
    if not isinstance(payload, dict):
        return None
    for key in ("validation_error_model", "classifier_error_model"):
        value = payload.get(key)
        if isinstance(value, dict):
            return value
    return None


def extract_confusion_counts(
    error_model: Dict[str, Any], class_order: List[str]
) -> np.ndarray:
    """Return the raw true-by-predicted confusion-count matrix."""
    if not isinstance(error_model, dict):
        raise ValueError("Validation error model must be a dictionary")

    matrix = error_model.get("confusion_counts")
    if not isinstance(matrix, list) or not matrix:
        raise ValueError(
            "Section 4 alignment requires raw 'confusion_counts'; "
            "Dirichlet posterior parameters alone are insufficient."
        )

    counts = np.asarray(matrix)
    n = len(class_order)
    if counts.shape != (n, n):
        raise ValueError("Confusion matrix shape does not match class_order")
    if not np.all(np.isfinite(counts)) or np.any(counts < 0):
        raise ValueError("Confusion counts must be finite and non-negative")

    rounded = np.rint(counts).astype(np.int32)
    if not np.allclose(counts, rounded):
        raise ValueError("Confusion counts must be integers")
    if np.any(rounded.sum(axis=1) == 0):
        raise ValueError(
            "Every true class needs at least one validation observation for "
            "the Section 4 confusion-matrix likelihood."
        )
    return rounded


def derive_dirichlet_posterior_parameters(
    error_model: Dict[str, Any], class_order: List[str], prior: float = 1.0
) -> List[List[float]]:
    """Compatibility helper for plotting P based on validation data alone."""
    counts = extract_confusion_counts(error_model, class_order)
    return (counts.astype(float) + prior).tolist()


def _simple(df: pd.DataFrame, classes: List[str]) -> pd.DataFrame:
    """Fallback used only when the NumPyro/JAX stack is unavailable."""
    rows = []
    for _, source in df.iterrows():
        row = {key: source[key] for key in ("timestamp", "run_name", "blob_name")}
        for class_name in classes:
            value = float(source.get(class_name, 0))
            row.update(
                {
                    f"{class_name}_corrected_mean": value,
                    f"{class_name}_corrected_median": value,
                    f"{class_name}_corrected_lower": max(0.0, value - 0.5),
                    f"{class_name}_corrected_upper": value + 0.5,
                }
            )
        rows.append(row)
    return pd.DataFrame(rows, index=df.index)

def allocation_probabilities(q: np.ndarray, p: np.ndarray) -> np.ndarray:
    """
    Return Pr(true=i | predicted=j, q, P).

    Shapes
    ------
    q:
        draw x location x true_class

    p:
        draw x true_class x predicted_class

    returns:
        draw x location x true_class x predicted_class
    """
    q = np.asarray(q, dtype=float)
    p = np.asarray(p, dtype=float)

    # Pr(true=i, predicted=j) = q_i * P_ij
    numerator = q[:, :, :, None] * p[:, None, :, :]

    # For each predicted class j, normalise across possible true classes i.
    denominator = numerator.sum(axis=2, keepdims=True)

    if np.any(~np.isfinite(denominator)) or np.any(denominator <= 0):
        raise ValueError("Invalid allocation-probability denominator")

    allocation = numerator / denominator

    # Each predicted-class column must sum to one over true classes.
    if not np.allclose(allocation.sum(axis=2), 1.0):
        raise AssertionError(
            "Allocation probabilities do not sum to one across true classes"
        )

    return allocation


def reconstruct_latent_true_counts(
    prevalence_samples: np.ndarray,
    confusion_samples: np.ndarray,
    observed_counts: np.ndarray,
    *,
    seed: int,
    return_flow: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Allocate each observed predicted count back to latent true classes."""
    q = np.asarray(prevalence_samples, dtype=float)
    p = np.asarray(confusion_samples, dtype=float)
    w = np.asarray(observed_counts, dtype=np.int64)

    if q.ndim != 3 or p.ndim != 3 or w.ndim != 2:
        raise ValueError(
            "Expected draw x location x true, draw x true x predicted, "
            "and location x predicted arrays"
        )
    draws, locations, true_classes = q.shape
    if p.shape[0] != draws or p.shape[1] != true_classes:
        raise ValueError("Posterior q and P dimensions are incompatible")
    if w.shape != (locations, p.shape[2]):
        raise ValueError("Observed counts and posterior dimensions are incompatible")
    if np.any(w < 0) or not np.all(np.isfinite(q)) or not np.all(np.isfinite(p)):
        raise ValueError("Inputs must be finite and observed counts non-negative")

    allocation = allocation_probabilities(q, p)
    flow = np.zeros((draws, locations, true_classes, p.shape[2]), dtype=np.int64)
    rng = np.random.default_rng(seed)
    for s in range(draws):
        for l in range(locations):
            for j in range(w.shape[1]):
                if w[l, j] > 0:
                    pvals = np.asarray(allocation[s, l, :, j], dtype=np.float64)
                    pvals = np.clip(pvals, 0.0, None)
                    total = pvals.sum(dtype=np.float64)
                    if not np.isfinite(total) or total <= 0.0:
                        raise ValueError(
                            "Invalid allocation probabilities at draw=%s, location=%s, predicted=%s"
                            % (s, l, j)
                        )
                    pvals /= total
                    prefix = pvals[:-1].sum(dtype=np.float64)
                    if prefix >= 1.0:
                        pvals[:-1] *= np.nextafter(1.0, 0.0) / prefix
                        prefix = pvals[:-1].sum(dtype=np.float64)
                    pvals[-1] = 1.0 - prefix
                    flow[s, l, :, j] = rng.multinomial(int(w[l, j]), pvals)

    latent = flow.sum(axis=3)
    expected = (allocation * w[None, :, None, :]).sum(axis=3)
    assert np.array_equal(flow.sum(axis=2), np.broadcast_to(w, (draws, locations, w.shape[1])))
    assert np.array_equal(latent.sum(axis=2), np.broadcast_to(w.sum(axis=1), (draws, locations)))

    return (latent, expected, flow) if return_flow else (latent, expected)

def joint_model(
    v: np.ndarray,
    y: np.ndarray,
    p_prior: float = 1.0,
    q_prior: float = 1.0,
):
    """
    Section 4 joint model.

    P is learned jointly from the raw validation confusion counts and field
    observations. Priors and likelihoods follow ``plankton_model.stan``.
    """
    if not isinstance(v, np.ndarray):
        v = np.asarray(v)

    if not isinstance(y, np.ndarray):
        y = np.asarray(y)

    k = v.shape[0]

    # ------------------------------------------------------------------
    # Confusion matrix
    # ------------------------------------------------------------------

    p = numpyro.sample(
        "confusion_matrix",
        dist.Dirichlet(
            jnp.full(k, p_prior, dtype=jnp.float64)
        ).expand((k,))
    )

    validation_totals = jnp.sum(v, axis=1)

    numpyro.sample(
        "validation_confusion_counts",
        dist.Multinomial(
            total_count=validation_totals,
            probs=p,
        ),
        obs=v,
    )

    # ------------------------------------------------------------------
    # Ecological prevalence model
    # ------------------------------------------------------------------

    ecological_mean = numpyro.sample(
        "ecological_mean",
        dist.Dirichlet(
            jnp.full(k, q_prior, dtype=jnp.float64)
        ),
    )

    # Stan: real<lower=0> b_par; b_par ~ normal(0, 1000).
    # HalfNormal(1000) is the equivalent normal distribution truncated at zero.
    ecological_concentration = numpyro.sample(
        "ecological_concentration",
        dist.HalfNormal(1000.0),
    )

    ecological_alpha = numpyro.deterministic(
        "ecological_alpha",
        ecological_concentration * ecological_mean,
    )

    q = numpyro.sample(
        "prevalence",
        dist.Dirichlet(ecological_alpha).expand((y.shape[0],)),
    )

    # ------------------------------------------------------------------
    # Total-count process
    # ------------------------------------------------------------------

    log_mu = numpyro.sample(
        "log_mu",
        dist.Normal(0.0, 1000.0),
    )

    mu = numpyro.deterministic(
        "mu",
        jnp.exp(log_mu),
    )

    # Stan: real<lower=0> k_nb; k_nb ~ normal(0, 1000).
    kappa = numpyro.sample(
        "k",
        dist.HalfNormal(1000.0),
    )

    numpyro.sample(
        "observed_totals",
        dist.NegativeBinomial2(
            mean=mu,
            concentration=kappa,
        ),
        obs=jnp.sum(y, axis=1),
    )

    # ------------------------------------------------------------------
    # Observation model
    # ------------------------------------------------------------------

    reported_prob = q @ p

    numpyro.sample(
        "observed_counts",
        dist.Multinomial(
            total_count=jnp.sum(y, axis=1),
            probs=reported_prob,
        ),
        obs=y,
    )

def fit_joint(
    v: np.ndarray,
    y: np.ndarray,
    samples: int,
    seed: int,
    warmup: int = NUMPYRO_MCMC_WARMUP,
    chains: int = NCHAINS,
):
    """Fit the Section 4 joint model using NumPyro."""
    if MCMC is None:
        raise RuntimeError("JAX/NumPyro unavailable")

    if numpyro is not None:
        numpyro.set_host_device_count(chains)

    mcmc = MCMC(
        NUTS(joint_model),
        num_warmup=warmup,
        num_samples=max(1, int(samples)),
        num_chains=chains,
        chain_method="parallel",
        progress_bar=True,
    )
    mcmc.run(
        jax.random.PRNGKey(seed),
        v=jnp.asarray(v),
        y=jnp.asarray(y),
    )
    return {k: np.asarray(x) for k, x in mcmc.get_samples().items()}


def validation_only_draws(v: np.ndarray, n: int, seed: int, prior: float = 1.0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    alpha = v + prior
    return np.stack([[rng.dirichlet(row) for row in alpha] for _ in range(n)])


def export_diagnostics(out: str, flow: np.ndarray, ef: np.ndarray, q: np.ndarray, p: np.ndarray, v: np.ndarray, times: List[Any], classes: List[str], target: str, seed: int):
    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)
    fi = classes.index(target)
    if ef.ndim != 4:
        raise ValueError(
            "Expected allocation flow must have shape (draws, locations, true_classes, predicted_classes)"
        )
    rows = []
    for tt, stamp in enumerate(times):
        total = flow[:, tt, :, fi].sum(axis=1)
        exp_total = ef[:, tt, :, fi].sum(axis=1)
        for s, name in enumerate(classes):
            x = flow[:, tt, s, fi]
            ex = ef[:, tt, s, fi]
            rows.append(
                dict(
                    timestamp=stamp,
                    predicted_source=name,
                    target_true_class=target,
                    sampled_mean=x.mean(),
                    sampled_median=np.median(x),
                    sampled_lower=np.quantile(x, 0.025),
                    sampled_upper=np.quantile(x, 0.975),
                    expected_mean=ex.mean(),
                    expected_median=np.median(ex),
                    expected_lower=np.quantile(ex, 0.025),
                    expected_upper=np.quantile(ex, 0.975),
                    observed_total=float(total.mean()),
                    expected_total=float(exp_total.mean()),
                )
            )
    pd.DataFrame(rows).to_csv(out_dir / "target_source_attribution.csv", index=False)

    validation_only = validation_only_draws(v, len(p), seed + 2)
    rows = []
    for j, name in enumerate(classes):
        joint = p[:, fi, j]
        validation = validation_only[:, fi, j]
        rows.append(
            dict(
                true_class=target,
                predicted_class=name,
                joint_mean=joint.mean(),
                joint_median=np.median(joint),
                joint_lower=np.quantile(joint, 0.025),
                joint_upper=np.quantile(joint, 0.975),
                validation_only_mean=validation.mean(),
                validation_only_median=np.median(validation),
                validation_only_lower=np.quantile(validation, 0.025),
                validation_only_upper=np.quantile(validation, 0.975),
            )
        )
    pd.DataFrame(rows).to_csv(out_dir / "target_confusion_joint_vs_validation_only.csv", index=False)
    np.savez_compressed(
        out_dir / "allocation_diagnostic_arrays.npz",
        class_order=np.asarray(classes),
        q=q,
        P=p,
        allocation_flow=flow,
        expected_allocation_flow=ef,
    )


def build_uncertainty_dataframe(
    plot_df: pd.DataFrame,
    validation_json: Any,
    *,
    mc_samples: int = 5000,
    mc_seed: int = 42,
    target_class: str = "fish_larvae",
    diagnostics_dir: Optional[str] = None,
) -> pd.DataFrame:
    """Fit the Section 4 joint model and summarise corrected latent counts."""
    error_model = extract_validation_error_model(validation_json)
    if error_model is None:
        raise ValueError("Validation JSON has no classifier error model")

    classes = [str(x) for x in error_model.get("class_order", []) if str(x)]
    if not classes:
        raise ValueError("Validation error model has no class_order")
    confusion_counts = extract_confusion_counts(error_model, classes)
    n = len(classes)

    data_columns = {
        column
        for column in plot_df.columns
        if column not in {"timestamp", "run_name", "blob_name"}
    }
    unknown = sorted(data_columns.difference(classes))
    if unknown:
        raise ValueError(
            "Observed classes are absent from the validation model: "
            + ", ".join(unknown)
        )

    if any(x is None for x in (jax, jnp, numpyro, dist, MCMC, NUTS)):
        LOG.warning(
            "JAX/NumPyro is unavailable; returning uncorrected fallback intervals"
        )
        return _simple(plot_df, classes)

    observed = plot_df.reindex(columns=classes, fill_value=0).to_numpy(float)
    observed_counts = np.rint(observed).astype(np.int32)
    if not np.allclose(observed, observed_counts) or np.any(observed_counts < 0):
        raise ValueError("Observed counts must be non-negative integers")
    totals = observed_counts.sum(axis=1).astype(np.int32)
    if np.any(totals == 0):
        raise ValueError(
            "Section 4's multinomial likelihood requires a positive total at "
            "every retained location/timestamp"
        )

    posterior = fit_joint(
        confusion_counts,
        observed_counts,
        max(1, int(mc_samples)),
        mc_seed,
        warmup=NUMPYRO_MCMC_WARMUP,
        chains=NCHAINS,
    )

    q = np.asarray(posterior["prevalence"])
    p = np.asarray(posterior["confusion_matrix"])
    allocation = allocation_probabilities(q, p)
    latent, expected, flow = reconstruct_latent_true_counts(
        q, p, observed_counts, seed=mc_seed + 1, return_flow=True
    )
    expected_flow = allocation * observed_counts[None, :, None, :]
    count_samples = latent.astype(float)
    composition = count_samples / totals[None, :, None]

    corrected_mean = expected.mean(axis=0)
    corrected_median = np.median(expected, axis=0)
    corrected_lower = np.percentile(expected, 5, axis=0)
    corrected_upper = np.percentile(expected, 95, axis=0)

    condition_number = float(np.mean([np.linalg.cond(draw.T) for draw in p]))
    ecological_alpha = np.asarray(posterior["ecological_alpha"])
    ecological_mean = np.asarray(posterior["ecological_mean"])
    ecological_concentration = np.asarray(posterior["ecological_concentration"])
    mu = np.asarray(posterior["mu"])
    kappa = np.asarray(posterior["k"])

    rows = []
    for location, (_, source) in enumerate(plot_df.iterrows()):
        row = {
            "timestamp": source["timestamp"],
            "run_name": source["run_name"],
            "blob_name": source["blob_name"],
            "observed_total": float(totals[location]),
            "mean_sampled_condition_number": condition_number,
            "total_process_mu_mean": float(mu.mean()),
            "total_process_mu_lower": float(np.percentile(mu, 5)),
            "total_process_mu_upper": float(np.percentile(mu, 95)),
            "total_process_k_mean": float(kappa.mean()),
            "total_process_k_lower": float(np.percentile(kappa, 5)),
            "total_process_k_upper": float(np.percentile(kappa, 95)),
            "ecological_concentration_mean": float(ecological_concentration.mean()),
            "ecological_concentration_lower": float(
                np.percentile(ecological_concentration, 5)
            ),
            "ecological_concentration_upper": float(
                np.percentile(ecological_concentration, 95)
            ),
        }
        for i, class_name in enumerate(classes):
            row.update(
                {
                    f"{class_name}_shared_alpha_mean": float(
                        ecological_alpha[:, i].mean()
                    ),
                    f"{class_name}_shared_alpha_lower": float(
                        np.percentile(ecological_alpha[:, i], 5)
                    ),
                    f"{class_name}_shared_alpha_upper": float(
                        np.percentile(ecological_alpha[:, i], 95)
                    ),
                    f"{class_name}_ecological_mean": float(
                        ecological_mean[:, i].mean()
                    ),
                    f"{class_name}_corrected_mean": float(
                        corrected_mean[location, i]
                    ),
                    f"{class_name}_corrected_median": float(
                        corrected_median[location, i]
                    ),
                    f"{class_name}_corrected_lower": float(
                        corrected_lower[location, i]
                    ),
                    f"{class_name}_corrected_upper": float(
                        corrected_upper[location, i]
                    ),
                    f"{class_name}_boundary_frequency": float(
                        np.mean(count_samples[:, location, i] == 0)
                    ),
                    f"{class_name}_expected_composition": float(
                        composition[:, location, i].mean()
                    ),
                    f"{class_name}_composition_median": float(
                        np.median(composition[:, location, i])
                    ),
                    f"{class_name}_composition_lower": float(
                        np.percentile(composition[:, location, i], 5)
                    ),
                    f"{class_name}_composition_upper": float(
                        np.percentile(composition[:, location, i], 95)
                    ),
                    f"{class_name}_expected_true_count": float(
                        corrected_mean[location, i]
                    ),
                    f"{class_name}_latent_true_count_median": float(
                        corrected_median[location, i]
                    ),
                    f"{class_name}_latent_true_count_lower": float(
                        corrected_lower[location, i]
                    ),
                    f"{class_name}_latent_true_count_upper": float(
                        corrected_upper[location, i]
                    ),
                }
            )
        rows.append(row)

    if diagnostics_dir:
        export_diagnostics(
            diagnostics_dir,
            flow,
            expected_flow,
            q,
            p,
            confusion_counts,
            plot_df["timestamp"].tolist(),
            classes,
            target_class,
            mc_seed,
        )

    return pd.DataFrame(rows, index=plot_df.index)


def load_validation_uncertainty_dataframe(
    plot_df,
    validation_source,
    *,
    mc_samples,
    mc_seed,
    blob_service_client=None,
    target_class="fish_larvae",
    diagnostics_dir=None,
):
    """Load validation JSON, fit the Section 4 model, and return summaries."""
    if not validation_source:
        return None
    try:
        payload = load_json_payload(
            validation_source, blob_service_client=blob_service_client
        )
        return build_uncertainty_dataframe(
            plot_df,
            payload,
            mc_samples=mc_samples,
            mc_seed=mc_seed,
            target_class=target_class,
            diagnostics_dir=diagnostics_dir,
        )
    except Exception as exc:
        LOG.exception("Skipping validation uncertainty: %s", exc)
        return None