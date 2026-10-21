"""scr-open NumPyro engine: open-population (Jolly-Seber) SCR, NUTS.

Port of ../../../python/fit_numpyro_open.py (SCR GPU benchmark, marginal
approach) with trap usage per period (K_jt) and any number of density
covariates. Alive states (1 not entered, 2 alive, 3 dead) are summed out with a
forward algorithm in log space and activity centres are summed over habitat-mask
pixels, so NUTS only sees p0, sigma, phi, the entry parameters and the density
betas. Entry uses the superpopulation parameterisation (psi ~ Uniform(0,1),
b ~ Dirichlet(1,..,1); gamma_t derived), matching the Nimble engines.

Usage: python scro_numpyro.py in.json out.json
Derived quantities (N_t, recruits, superpopulation) are computed in R from the
parameter draws (modules/scr-open/R/post.R).
"""
import json
import sys
import time

inp = json.load(open(sys.argv[1]))
import jax  # noqa: E402

jax.config.update("jax_enable_x64", bool(inp.get("x64", True)))
import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402
import numpyro  # noqa: E402
import numpyro.distributions as dist  # noqa: E402
from jax.scipy.special import logsumexp  # noqa: E402
from numpyro.infer import MCMC, NUTS, init_to_value  # noqa: E402

FT = jnp.float64 if inp.get("x64", True) else jnp.float32
NEG = -1e30
d2 = jnp.asarray(np.array(inp["d2"], dtype=float), FT)            # G x J
G, J = d2.shape
Kjt = jnp.asarray(np.array(inp["Kjt"], dtype=float).reshape(J, -1), FT)   # J x T
T = Kjt.shape[1]
PD = int(inp["PD"])
XD = jnp.asarray(np.array(inp["XD"], dtype=float).reshape(G, PD) if PD else np.zeros((G, 0)), FT)
n, M = int(inp["n"]), float(inp["M"])
nzi = jnp.asarray(np.array(inp["nzi"], dtype=np.int64))
nzj = jnp.asarray(np.array(inp["nzj"], dtype=np.int64))
nzt = jnp.asarray(np.array(inp["nzt"], dtype=np.int64))
nzc = jnp.asarray(np.array(inp["nzc"], dtype=float), FT)
seen = jnp.asarray(np.array(inp["seen"], dtype=bool).reshape(n, T))
smax = float(inp["sigma_max"])


def log1mexp(lp):
    return jnp.where(lp > -0.693, jnp.log(-jnp.expm1(lp)), jnp.log1p(-jnp.exp(lp)))


def forward(log_alive, seenb, lgam, lgam1, lphi, lphi1):
    zero = jnp.broadcast_to(jnp.where(seenb, NEG, 0.0), log_alive.shape)
    a = jnp.stack([lgam1[0] + zero[..., 0], lgam[0] + log_alive[..., 0],
                   jnp.full(zero[..., 0].shape, NEG)], axis=-1)
    for t in range(1, T):
        ne = a[..., 0] + lgam1[t] + zero[..., t]
        en = jnp.logaddexp(a[..., 0] + lgam[t], a[..., 1] + lphi) + log_alive[..., t]
        de = jnp.logaddexp(a[..., 1] + lphi1, a[..., 2]) + zero[..., t]
        a = jnp.stack([ne, en, de], axis=-1)
    return logsumexp(a, axis=-1)


def model():
    p0 = numpyro.sample("p0", dist.Uniform(0.0, 1.0))
    sigma = numpyro.sample("sigma", dist.Uniform(0.0, smax))
    phi = numpyro.sample("phi", dist.Uniform(0.0, 1.0))
    psi = numpyro.sample("psi", dist.Uniform(0.0, 1.0))
    b = numpyro.sample("b", dist.Dirichlet(jnp.ones(T, dtype=FT)))
    cb_prev = jnp.concatenate([jnp.zeros(1, dtype=FT), jnp.cumsum(b)[:-1]])
    gam = numpyro.deterministic("gam", jnp.clip(psi * b / (1.0 - psi * cb_prev), 1e-12, 1.0 - 1e-12))
    if PD:
        bD = numpyro.sample("bD", dist.Normal(0.0, 2.0).expand([PD]).to_event(1))
        lw = XD @ bD
        lpi = lw - logsumexp(lw)
    else:
        lpi = jnp.full(G, -jnp.log(G), dtype=FT)
    lgam, lgam1 = jnp.log(gam), jnp.log1p(-gam)
    lphi, lphi1 = jnp.log(phi), jnp.log1p(-phi)
    logp = jnp.log(p0) - d2 / (2 * sigma ** 2)                  # G x J
    l1 = log1mexp(logp)
    base = l1 @ Kjt                                              # G x T
    lg = logp - l1
    col = lg[:, nzj].T * nzc[:, None]                            # nnz x G
    flat = jnp.zeros((n * T, G), dtype=FT).at[nzi * T + nzt].add(col)
    ll = flat.reshape(n, T, G).transpose(0, 2, 1) + base[None, :, :]   # n x G x T
    lf = forward(ll, seen[:, None, :], lgam, lgam1, lphi, lphi1)      # n x G
    lf_det = logsumexp(lpi[None, :] + lf, axis=1)
    lf0_g = forward(base, jnp.zeros((G, T), dtype=bool), lgam, lgam1, lphi, lphi1)
    lf0 = logsumexp(lpi + lf0_g)
    numpyro.factor("loglik", lf_det.sum() + (M - n) * lf0)


ini = inp["init"]
init = dict(p0=ini["p0"], sigma=ini["sigma"], phi=ini["phi"], psi=ini["psi"], b=np.asarray(ini["b"], dtype=float))
if PD:
    init["bD"] = np.asarray(ini["bD"], dtype=float)
C = int(inp["chains"])
kernel = NUTS(model, init_strategy=init_to_value(values=init), target_accept_prob=0.85)
mcmc = MCMC(kernel, num_warmup=int(inp["warmup"]), num_samples=int(inp["samples"]), num_chains=C,
            chain_method="vectorized", progress_bar=False)
t0 = time.perf_counter()
mcmc.run(jax.random.PRNGKey(int(inp["seed"])), extra_fields=("num_steps", "diverging"))
s = mcmc.get_samples(group_by_chain=True)
t_total = time.perf_counter() - t0
ex = mcmc.get_extra_fields(group_by_chain=True)
draws = {k: np.asarray(s[k]).reshape(C, -1).tolist() for k in ["p0", "sigma", "phi"]}
g = np.asarray(s["gam"])
for t in range(T):
    draws[f"gam{t + 1}"] = g[:, :, t].tolist()
if PD:
    b = np.asarray(s["bD"])
    for k in range(PD):
        draws[f"bD{k + 1}"] = b[:, :, k].tolist()
out = dict(draws=draws, timing=dict(t_total=t_total, device=str(jax.devices()[0]),
                                    divergences=int(np.asarray(ex["diverging"]).sum()),
                                    mean_num_steps=float(np.asarray(ex["num_steps"]).mean())))
json.dump(out, open(sys.argv[2], "w"))
print(json.dumps(out["timing"]))
