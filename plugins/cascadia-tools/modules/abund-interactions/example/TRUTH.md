# Truth used to simulate the example (for checking only; real data have no truth)

Model: N ~ Poisson(1.0 * exp(0.5 z(woodland_pct))); r = 0.08 per marten per day;
logit psi = logit(0.6) + 0.5 z(woodland_pct) + 0.3 z(broadleaf_pct) + N (-0.6 + 0.6 z(broadleaf_pct)); p = 0.15 per day.
z() = standardised over the 260 sites (as the module does).
Realised: total marten-site abundance 274; sites with marten 0.592; sites with red squirrel 0.504.
