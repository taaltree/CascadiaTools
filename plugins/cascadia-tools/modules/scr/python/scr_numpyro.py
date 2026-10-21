"""ecotools scr: marginal SCR likelihood fitted with NumPyro NUTS.

Ported from ../../../python/fit_numpyro.py (SCR GPU benchmark, approach
"marginal"), extended with per-trap effort, trap-level detection covariates,
density covariates on an irregular habitat mask, and a Poisson (count)
observation model. Activity centres and data-augmentation indicators are
integrated over the habitat mask; only the top-level parameters are sampled.

  psi ~ U(0,1), sigma ~ U(0, sigma_max), bD ~ N(0,2), bP ~ N(0,2)
  binomial: det0 = p0 ~ U(0,1), p0_j = expit(logit(p0) + XP_j bP)
            log p_gj = log p0_j - d2_gj / (2 sigma^2); y_ij ~ Bin(u_j, p_gj)
  poisson:  det0 = lambda0 ~ U(0, lam_max), log lambda0_j = log lambda0 + XP_j bP
            y_ij ~ Poisson(u_j lambda0_j exp(-d2_gj / (2 sigma^2)))   (u_j = total effort)
  pi_g = softmax(XD bD)

Usage: python scr_numpyro.py in.json out.json
in.json keys: d2 (G x J), use (J), XD (G x PD), XP (J x PP), nzi, nzj, nzc (0-based),
n, M, sigma_max, poisson, lam_max, chains, warmup, samples, seed, x64, init (dict).
out.json: {"draws": {name: [[chain0...], [chain1...]]}, "timing": {...}}
"""
import json
import sys
import time

inp, outp = sys.argv[1], sys.argv[2]
with open(inp) as fh:
    d = json.load(fh)

import jax  # noqa: E402

jax.config.update("jax_enable_x64", bool(d.get("x64", True)))
import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402
import numpyro  # noqa: E402
import numpyro.distributions as dist  # noqa: E402
from jax.scipy.special import logsumexp  # noqa: E402
from numpyro.infer import MCMC, NUTS, init_to_value  # noqa: E402

FT = jnp.float64 if d.get("x64", True) else jnp.float32
d2 = jnp.asarray(np.asarray(d["d2"], dtype=float), FT)
G, J = d2.shape
use = jnp.asarray(np.asarray(d["use"], dtype=float), FT)
XD = np.asarray(d["XD"], dtype=float).reshape(G, -1) if d["PD"] > 0 else np.zeros((G, 0))
XP = np.asarray(d["XP"], dtype=float).reshape(J, -1) if d["PP"] > 0 else np.zeros((J, 0))
PD, PP = XD.shape[1], XP.shape[1]
XD, XP = jnp.asarray(XD, FT), jnp.asarray(XP, FT)
n, M = int(d["n"]), float(d["M"])
nzi = jnp.asarray(np.asarray(d["nzi"], dtype=np.int64))
nzj = jnp.asarray(np.asarray(d["nzj"], dtype=np.int64))
nzc = jnp.asarray(np.asarray(d["nzc"], dtype=float), FT)
smax = float(d["sigma_max"])
POIS = bool(d.get("poisson", False))
lam_max = float(d.get("lam_max", 10.0))


def log1mexp(lp):
    return jnp.where(lp > -0.693, jnp.log(-jnp.expm1(lp)), jnp.log1p(-jnp.exp(lp)))


def model():
    psi = numpyro.sample("psi", dist.Uniform(0.0, 1.0))
    if POIS:
        det0 = numpyro.sample("det0", dist.Uniform(0.0, lam_max))
    else:
        det0 = numpyro.sample("det0", dist.Uniform(0.0, 1.0))
    sigma = numpyro.sample("sigma", dist.Uniform(0.0, smax))
    if PD > 0:
        bD = numpyro.sample("bD", dist.Normal(0.0, 2.0).expand([PD]).to_event(1))
        lw = XD @ bD
        lpi = lw - logsumexp(lw)
    else:
        lpi = jnp.full(G, -jnp.log(G), dtype=FT)
    eta = jnp.log(det0) if POIS else jnp.log(det0) - jnp.log1p(-det0)
    if PP > 0:
        bP = numpyro.sample("bP", dist.Normal(0.0, 2.0).expand([PP]).to_event(1))
        eta = eta + XP @ bP
    else:
        eta = jnp.full(J, eta, dtype=FT)
    lq = eta if POIS else jax.nn.log_sigmoid(eta)
    logp = lq[None, :] - d2 / (2 * sigma ** 2)             # G x J
    if POIS:
        base = -(jnp.exp(logp) * use[None, :]).sum(1) + lpi
        contrib = nzc[:, None] * logp[:, nzj].T            # nnz x G (log u terms are constants)
    else:
        log1mp = log1mexp(logp)
        base = (use[None, :] * log1mp).sum(1) + lpi
        contrib = nzc[:, None] * (logp - log1mp)[:, nzj].T
    ll = jnp.zeros((n, G), dtype=FT).at[nzi].add(contrib) + base
    lf_det = logsumexp(ll, axis=1)
    lf0 = logsumexp(base)
    loglik = n * jnp.log(psi) + lf_det.sum() + (M - n) * jnp.logaddexp(jnp.log(psi) + lf0, jnp.log1p(-psi))
    numpyro.factor("loglik", loglik)


init = {k: (jnp.asarray(v, FT) if isinstance(v, list) else v) for k, v in d["init"].items()}
kernel = NUTS(model, init_strategy=init_to_value(values=init), target_accept_prob=0.8)
C = int(d["chains"])
mcmc = MCMC(kernel, num_warmup=int(d["warmup"]), num_samples=int(d["samples"]), num_chains=C,
            chain_method="vectorized", progress_bar=False)
t0 = time.perf_counter()
mcmc.run(jax.random.PRNGKey(int(d["seed"])), extra_fields=("diverging", "num_steps"))
smp = mcmc.get_samples(group_by_chain=True)
_ = float(smp["psi"].sum())
t = time.perf_counter() - t0
out = {}
for k in ("psi", "det0", "sigma"):
    out[k] = np.asarray(smp[k]).tolist()
for nm, P in (("bD", PD), ("bP", PP)):
    if P > 0:
        a = np.asarray(smp[nm])
        for q in range(P):
            out[f"{nm}{q + 1}"] = a[:, :, q].tolist()
ex = mcmc.get_extra_fields(group_by_chain=True)
timing = dict(t_total=t, divergences=int(np.asarray(ex["diverging"]).sum()),
              mean_num_steps=float(np.asarray(ex["num_steps"]).mean()), device=str(jax.devices()[0]))
with open(outp, "w") as fh:
    json.dump({"draws": out, "timing": timing}, fh)
print(json.dumps(timing))
