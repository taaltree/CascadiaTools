# Worked example: genetic diversity and structure of tailed frogs in eight streams

Module: `popgen-structure`. Appendix S1 worked example and website "try an example".
All numbers below come from the runs listed in section (g). The data are
**simulated**, so the true answers are known (see `README_data.md`).

---

## (a) The question, in the researcher's words

> "I swabbed coastal tailed frogs in eight Coast Range streams and had them
> genotyped at 12 microsatellites. I want to know how much genetic diversity
> each stream has, whether any population is inbred or isolated, how
> different the streams are from each other, and whether streams farther
> apart are more different. Is any stream a candidate for genetic rescue?"

## (b) The files the researcher has

`tailed_frog_microsats.csv`: one row per frog, two columns per locus, 0 = missing
(first 5 frogs, first 4 loci shown; the file has 12 loci and 200 frogs):

```
individual,population,Atru01,Atru01.2,Atru02,Atru02.2,Atru03,Atru03.2,Atru04,Atru04.2,...
ATNFSIU-01,NFSIU,158,158,240,244,114,118,250,250,...
ATNFSIU-02,NFSIU,158,158,216,240,100,118,238,250,...
ATNFSIU-03,NFSIU,152,164,240,240,114,118,230,250,...
ATNFSIU-04,NFSIU,0,0,216,240,108,114,238,250,...
ATNFSIU-05,NFSIU,144,168,216,216,100,114,230,250,...
```

`tailed_frog_sites.csv`: one row per stream (first 5 of 8):

```
population,stream,latitude,longitude,elevation_m
NFSIU,North Fork Siuslaw River,44.112,-123.853,210
LOBST,Lobster Creek,44.205,-123.642,340
FIVEM,Five Mile Creek,44.081,-123.917,150
DRIFT,Drift Creek,44.418,-123.985,120
```

## (c) What the researcher types to Claude

> I have microsatellite genotypes for coastal tailed frogs from 8 streams
> (tailed_frog_microsats.csv, two columns per locus, 0 means missing) and the
> stream coordinates (tailed_frog_sites.csv). Please tell me how genetically
> diverse each stream is, whether any population is inbred or isolated, how
> different the streams are and whether that increases with distance, and
> how many genetic groups there are. I don't code, so please give me a report
> with figures and explain what it means for management.

## (d) What the agent does, step by step (eco-analyze)

**1. Project folder.** Creates `projects/<name>/` with the raw files copied
into `data_raw/` (never edited), an analysis plan and a log of every command.
(For this appendix the files live in `modules/popgen-structure/example/` and
the outputs in `results/popgen-structure/si/`.)

**2. Data intake (eco-data-intake).** The agent reads both files and
describes them back:

* 200 frogs in 8 streams: NFSIU 30, LOBST 24, FIVEM 28, DRIFT 18, MARYS 32,
  ALSEA 26, YACHT 22, SMITH 20.
* 12 microsatellite loci, none monomorphic. 2.6% of genotypes are missing.
* All 8 streams have coordinates. They are 6.2 to 80.8 km apart.
* The streams lie in several catchments (a 2-D landscape), so distance
  effects are tested on log distance.
* No duplicated IDs, no half-scored genotypes, no unmatched stream names.

**3. Choose the module.** `popgen-structure` fits: diploid codominant
markers, several sampled populations, questions about diversity, FIS, FST,
isolation by distance and clusters. (Effective size or relatedness would be
`popgen-ne`.)

**4. Plan (eco-analysis-plan).** Written before fitting:

* Estimands: Ho, He and allelic richness per stream; FIS per stream and
  overall; global and pairwise FST (Weir & Cockerham 1984); isolation-by-
  distance slope of FST/(1-FST) on log distance, with a Mantel test; number
  of genetic clusters (k-means + BIC, DAPC).
* Intervals: 95%, jackknife over loci plus the covariance from shared
  individuals, t(11 df).
* Checks: Hardy-Weinberg tests per stream and locus (FDR-corrected); loci
  that deviate in more than half the streams would be flagged as possible
  null alleles.
* Adequate result: planted-truth gate passes at this design; no systematic
  Hardy-Weinberg failure.

**5. Planted truth at this design (eco-planted-truth).** The agent simulates
200 datasets that copy the study design (8 populations, the same 18-32
frogs per population, 12 microsatellites with up to 14 alleles, 3% missing,
one small population of Ne = 30, others Ne 100-150) and refits each one.
Gate: |bias| at most 10% of truth (or 3 Monte Carlo SE) and 95%-interval
coverage at least 0.904.

| quantity | truth (mean) | mean estimate | relative bias | coverage | power | gate |
|---|---|---|---|---|---|---|
| global FST | 0.0757 | 0.0812 | +7.2% | 0.930 | 1.00 | PASS |
| global FIS | -0.0001 | -0.00002 | bias 0.0001 | 0.965 | 0.03 (false-positive rate) | PASS |
| mean He | 0.743 | 0.743 | +0.04% | 0.925 | | PASS |
| mean Ho | 0.743 | 0.743 | +0.01% | 0.940 | | PASS |
| FST, neighbouring pair | 0.0469 | 0.0465 | -0.8% | 0.910 | | PASS |
| FST, most distant pair | 0.0991 | 0.0975 | -1.6% | 0.965 | | PASS |
| IBD slope | 0.0691 | 0.0682 | -1.2% | 0.965 | 0.995 | PASS |

All 7 PASS; 200 of 200 fits ran. In one sentence for the researcher:
*"At your design the method recovers the true values with less than 8% bias,
and its 95% intervals contain the truth 91-97% of the time."*

The module's simulator puts populations on a grid. The example's own
history used the real stream map and an almost isolated Drift Creek. So the
agent also re-ran planted truth with **that exact process** (200 histories;
`dev/popgen-structure/si/si_design_truth.R`) and scored all 57 quantities.
Every one had relative bias under 6%. 52 of 57 passed the coverage bound.
The 5 misses were close: 4 of 28 pairwise FST (coverage 0.890-0.900) and
He in SMITH (0.885). Global FST: bias +0.8%, coverage 0.970. IBD slope:
+1.9%, coverage 0.970.

One thing this check revealed: global FST must be compared with a truth that
weights streams the way the sample does (see section e).

**6. Fit.** The agent runs the module on the real files (`config.yaml`).
It takes under a minute.

**7. Checks (eco-convergence, eco-model-checking).**

* Convergence: not applicable. The estimators are closed-form, not MCMC.
  The fit reports "converged"; 12 loci, 0 monomorphic.
* Interval health: the shared-individual part of the variance is 8% (median).
* Goodness of fit: **0 of 94** stream x locus Hardy-Weinberg tests deviate
  after false-discovery-rate correction. No locus looks like a null allele.
  Global FIS = 0.003 (95% CI -0.028 to 0.034): random mating within streams.

**8. Interpretation.** Excerpts from the report (`popgen-structure_report.html`),
each checked against the estimates table:

> "Overall differentiation among the 8 sampled populations is FST = 0.094
> (95% CI 0.075 to 0.114) ... about 9.4% of the gene diversity lies between
> populations."

> "Pairwise FST ranges from 0.007 (95% CI -0.004 to 0.018) between FIVEM and
> NFSIU to 0.366 (95% CI 0.276 to 0.456) between DRIFT and SMITH; 27 of 28
> pairs have intervals excluding 0."

> "Genetic diversity: mean expected heterozygosity is 0.678 (95% CI 0.627 to
> 0.728). It is lowest in DRIFT (0.256 (95% CI 0.126 to 0.386)) ... small,
> isolated populations lose diversity by drift and are the first candidates
> to monitor for inbreeding depression and to consider for genetic rescue."

> "Isolation by distance: FST/(1-FST) increases by 0.0566 (95% CI 0.0217 to
> 0.0916) per unit of log distance ...; Mantel r = 0.17, permutation P = 0.240.
> The evidence is mixed. ... DRIFT is strongly differentiated from all others
> ... Do not conclude that differentiation increases with distance; re-run
> without DRIFT (or with more populations) to check."

> "Clustering ... found K = 2 genetic clusters for 8 sampled populations."
> (Drift Creek versus the other seven.)

The island-model translation (about 2.4 migrants per generation, 1.9-3.1) is
reported only as an order of magnitude, as the report says.

**Figures** (`results/popgen-structure/si/fit/figures/`, PDF + PNG):

| id | caption (short) |
|---|---|
| `diversity` | Ho, He and allelic richness (32 gene copies) per stream with 95% intervals. Drift Creek stands far below the rest. |
| `fst_matrix` | Pairwise FST heatmap. The bright DRIFT column (0.24-0.37) against 0.01-0.08 elsewhere. |
| `ibd` | FST/(1-FST) against distance (log km), fitted line and 95% band. The seven DRIFT pairs sit high at every distance. |
| `pca` | Individuals on PC1 (11.8%) and PC2 (5.6%). Drift Creek frogs form a separate group on PC1. |
| `bic_k` | BIC for K = 1-16 clusters; K = 2 chosen by the diffNgroup rule (BIC is lowest at K = 4-5, and the curve is flat there). |
| `dapc_assignment` | DAPC membership: Drift Creek frogs in cluster 2 (one frog partly in cluster 1), all others in cluster 1; a few MARYS and YACHT frogs are partly assigned to cluster 2. |
| `sample_map` | Streams on the map, coloured by He and sized by sample size. |

All seven figures were opened and checked. No change to the `figures:` block
was needed.

**9. Gate (eco-acceptance-gate).** Planted truth PASS at this design (7/7);
convergence not applicable (closed form); goodness of fit OK (0/94 HWE);
independent route: estimators match hierfstat 0.5.11 to < 1e-10
(VALIDATION.md section 6); negative control: a naive estimator fails the
gate (VALIDATION.md section 4); reproducible commands below.
**Verdict: DONE**, with the limits in (f).

**10. Deliver.** Summary for the researcher:

> Drift Creek is the outlier. Its frogs carry much less genetic diversity
> (He 0.26, 95% CI 0.13-0.39, versus 0.70-0.75 in the other streams) and it
> is strongly differentiated from every other stream (pairwise FST
> 0.24-0.37). This is the signature of a small, isolated population. It is
> the candidate for monitoring and, if fitness problems appear, for genetic
> rescue from a nearby stream. The other seven streams are moderately
> differentiated (FST 0.01-0.08) and show no sign of inbreeding. Whether
> differentiation increases with distance is unclear: the slope is positive
> but the Mantel test is not significant, and Drift Creek drives the pattern.

## (e) Estimates versus the known truth

Three versions of "truth" exist for this dataset (`results/popgen-structure/si/true_values_with_design.csv`):
the census of all 1,040 simulated frogs; the census reweighted the way the
sample weights streams (the quantity the FST estimator targets; "design");
and the model expectation over many histories.

| quantity | estimate (95% CI) | truth: census | truth: design | truth: model |
|---|---|---|---|---|
| global FST | 0.094 (0.075 to 0.114) | 0.065 (outside) | 0.099 (inside) | 0.086 (inside) |
| global FIS | 0.003 (-0.028 to 0.034) | -0.004 | -0.004 | 0 |
| mean He | 0.678 (0.627 to 0.728) | 0.673 | 0.673 | 0.660 |
| He, Drift Creek | 0.256 (0.126 to 0.386) | 0.220 | 0.220 | 0.259 |
| He, North Fork Siuslaw | 0.755 (0.694 to 0.815) | 0.743 | 0.743 | 0.734 |
| FST, Drift vs Smith | 0.366 (0.276 to 0.456) | 0.334 | 0.388 | 0.287 |
| FST, Five Mile vs NF Siuslaw | 0.007 (-0.004 to 0.018) | 0.019 (outside) | 0.019 (outside) | 0.021 (outside) |
| IBD slope | 0.057 (0.022 to 0.092) | 0.058 | 0.069 | 0.047 |

Over all 57 quantities with intervals, the design truth is inside the
interval for 54 (95%). The 3 misses: FST Five Mile-NF Siuslaw and FST
Five Mile-Smith (both estimated lower than truth, intervals just short), and
FIS in Drift Creek (true census value -0.11, interval -0.099 to -0.050).

**Why global FST differs from the census value.** Drift Creek holds 3% of
all frogs but 9% of the genotyped frogs. Weir & Cockerham's FST weights each
stream by its sample size. Drift Creek is by far the most differentiated
stream, so giving it more weight raises FST. Against the matching
(design-weighted) truth the estimator is unbiased (+0.8%, coverage 0.970 in
200 replicate histories). Against the census-weighted truth it is +47% with
coverage 0.24. **Global FST describes the populations as sampled.** When
one very different population is over-sampled, report pairwise FST and
per-population diversity alongside it.

## (f) Limits from VALIDATION.md that matter here

* **Pairwise FST intervals are slightly narrow with about 10-12 loci.**
  Here 4 of 28 pairs had coverage 0.890-0.900 (bound 0.904) in the
  exact-process check. Treat pairwise intervals as approximate.
* **Global FST weights populations by sample size.** With unequal samples
  and one outlier population it is not the census-weighted value
  (finding of this example, section e).
* **Isolation by distance**: 8 populations is the minimum; one outlier
  (Drift Creek) can create or hide the pattern. The slope interval and
  Mantel test disagree here, and the report says so.
* **FST is not migration.** Nm (2.4) is an island-model order of magnitude only.
* **Intervals describe these populations**, not the evolutionary process
  (He coverage against the model expectation was 0.900).
* **DAPC/K are summaries**, not management units. K was chosen by the
  diffNgroup rule (K = 2), although BIC itself is lowest at K = 4-5.
  DAPC is not part of the recovery gate.
* **Neutral, unlinked loci assumed.** Not validated: null alleles, scoring
  error, sampling of siblings. HWE tests (0/94 deviations) show no sign of
  null alleles here.

## (g) Reproduce

From `.` (R 4.5 via `scripts/Rscript`; on the
cluster wrap each in `scripts/ejob`):

```bash
# 0. regenerate the example data (identical CSVs) and the true values
scripts/ejob -c 2 -m 8 -t 0:30:00 -l logs/popgen-structure_si_truth.log -- \
  scripts/Rscript dev/popgen-structure/si/si_truth.R
# 1. planted-truth gate at the example's design (about 3 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/popgen-structure_si_recover.log -- \
  scripts/Rscript core/run.R recover popgen-structure \
  --config modules/popgen-structure/example/recovery.yaml --nsim 200 --cores 16 \
  --out results/popgen-structure/si/recovery
# 2. exact generating process, design-weighted truth, all 57 quantities
scripts/ejob -c 16 -m 32 -t 1:00:00 -l logs/popgen-structure_si_design.log -- \
  scripts/Rscript dev/popgen-structure/si/si_design_truth.R 200 16
# 3. the real-data fit: report + figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/popgen-structure_si_fit.log -- \
  scripts/Rscript core/run.R fit popgen-structure \
  --config modules/popgen-structure/example/config.yaml --out results/popgen-structure/si/fit
```

Supporting diagnostics (census-weighted and rescaled-weight truths) are in
`dev/popgen-structure/si/si_exact_recovery.R` and `si_exact_weights.R`.
Machine-readable summary: `results/popgen-structure/si/si_summary.json`.
