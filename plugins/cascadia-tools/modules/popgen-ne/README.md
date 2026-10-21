---
id: popgen-ne
title: Effective population size and relatedness
when_to_use: Genotypes (microsatellites or SNPs) from one or more samples of a population; question is how large the genetically effective population is (Ne), how fast it is losing diversity, whether Ne is above conservation thresholds (50/100, 500/1000), and which sampled individuals are close relatives (full sibs, half sibs, parent-offspring).
not_for: Population structure, FST, isolation by distance, assignment (use popgen-structure); demographic history over many generations (coalescent methods such as GONE/SNeP); census abundance (use capture-recapture or close-kin mark-recapture models); pedigree reconstruction with parentage likelihoods (COLONY, sequoia).
data_needed: one row per individual with population, optional sampling time and optional flag to exclude juveniles/siblings from Ne; diploid genotypes as two allele columns per locus, "152/160" single columns, SNP 0/1/2 dosages, GENEPOP or VCF; optional list of known relative pairs; optional census size and generation time for the temporal method
engine: R implementation of the NeEstimator v2 formulas (LDNe with jackknife CIs; temporal Jorde-Ryman, Nei-Tajima, Pollak), Queller-Goodnight and Lynch-Ritland moment estimators; cross-checked against the NeEstimator C program
---

# Effective population size (Ne) and relatedness

## Question it answers

- **How big is the genetically effective population?** Ne is the size of an
  ideal population that would lose genetic diversity at the same rate as the
  real one. It sets the rate of inbreeding and of loss of heterozygosity,
  1/(2Ne) per generation, and is the quantity behind the 50/500 (now
  100/1000) conservation rules (Franklin 1980; Frankham et al. 2014).
- **Single sample (LD method).** Drift in a small population creates random
  associations (linkage disequilibrium) between unlinked loci; their size,
  after removing the part due to sampling a finite number of individuals,
  gives Ne (Waples 2006; Waples & Do 2008).
- **Two samples t generations apart (temporal method).** Allele frequencies
  drift faster in small populations; the standardised variance of change,
  corrected for sampling, gives the harmonic mean Ne over the interval
  (Nei & Tajima 1981; Pollak 1983; Waples 1989; Jorde & Ryman 2007).
- **Who is related to whom?** Pairwise relatedness r (0 unrelated, 0.25 half
  sibs, 0.5 full sibs / parent-offspring) from moment estimators
  (Queller & Goodnight 1989; Lynch & Ritland 1999), with the misclassification
  rate of kin calls at *your* markers estimated by simulation.
- Heterozygosity (Ho, He, F_IS) and, optionally, a heterozygosity-excess
  bottleneck test (Cornuet & Luikart 1996).

## Model

**LD method (LDNe).** For each pair of loci, Burrows' composite correlation
r between allele dosages (equivalently the Pearson correlation of 0/1/2 allele
counts), multiplied by S/(S-1), is squared and averaged over allele pairs
(alleles below the frequency cut-off Pcrit are excluded). Locus pairs are
averaged with weights equal to the number of independent allelic comparisons
(K1-1)(K2-1) (times S^2 when data are missing; Peel et al. 2013), giving R^2.
Sampling contributes E(R^2_sample) = 1/S + 3.19/S^2 (S >= 30) or
0.0018 + 0.907/S + 4.44/S^2 (S < 30), where S is the weighted harmonic mean
sample size. With R^2' = R^2 - E(R^2_sample) (the drift signal),

  Ne = [1/3 + sqrt(1/9 - 2.76 R^2')] / (2 R^2')   (S >= 30)
  Ne = [0.308 + sqrt(0.308^2 - 2.08 R^2')] / (2 R^2')   (S < 30)

(Waples 2006, random mating; a monogamy variant is available). When R^2' <= 0
the estimate is infinite. 95% CIs: jackknife over individuals, with the
empirical variance correction 0.84^2 and an effective number of degrees of
freedom 2/CV^2 for a chi-square interval on R^2 (Jones et al. 2016, as in
NeEstimator v2.1); a parametric chi-square CI (df = number of independent
comparisons) is also computed but is far too narrow with many loci.

**Temporal method.** With x and y the allele frequencies in the two samples
and z = (x+y)/2, Jorde & Ryman's Fs = sum (x-y)^2 / sum z(1-z) over alleles
and loci, corrected for sampling:

  Fs' = [Fs (1 - 1/(4S~) + 1/(4N)) - 1/S~ + 1/N] / [(1 + Fs/4)(1 - 1/(2 S_t))],

Ne = t / (2 Fs'), with S~ the harmonic mean sample size, S_t the size of the
second sample, N the census size at the first sample (sampling plan I; omit
N for plan II, where sampled individuals do not contribute offspring).
Nei-Tajima Fc and Pollak Fk (F' = F - 1/S~ [+ 1/N]) are reported for
comparison. CI for Jorde-Ryman: chi-square interval on the uncorrected Fs
with df = 2/CV^2 from a jackknife over loci, mapped to Fs' (NeEstimator's
normal jackknife interval is kept internally; it under-covered with 15
microsatellites); Fc/Fk: chi-square with df = number of independent alleles. **With hundreds of SNPs these CIs are too
narrow** because allele-frequency changes at unlinked loci are correlated
through linkage disequilibrium in a small population; the module inflates the
jackknife variance by a design effect estimated from the LD in the samples
(see VALIDATION.md).

**Relatedness.** Allele frequencies come from the sample being analysed.
Queller-Goodnight: r_xy = sum_l [ (Sac+Sad+Sbc+Sbd)/2 - pa - pb ] /
sum_l [ 1 + Sab - pa - pb ] (alleles a,b of x and c,d of y; S = 1 if
identical), averaged over the two reference directions. Lynch-Ritland: the
per-locus regression estimator [pa(Sbc+Sbd) + pb(Sac+Sad) - 4 pa pb] /
[(1+Sab)(pa+pb) - 4 pa pb] combined across loci with weights
[(1+Sab)(pa+pb) - 4 pa pb] / (2 pa pb), averaged over the two directions.
Pairs with r >= 0.354 are called first-degree and 0.177 <= r < 0.354
second-degree relatives (KING thresholds, Manichaikul et al. 2010).
Unrelated, half-sib, full-sib and parent-offspring pairs are simulated from the
sample's allele frequencies (unlinked loci, Hardy-Weinberg) to estimate how
often each class is called correctly. Mean r of known pairs is reported with a
jackknife-over-loci CI.

## Assumptions (and what breaks if they fail)

- **Closed, randomly mating population; one gene pool per sample.**
  Immigrants or a mixture of subpopulations create LD unrelated to drift and
  bias LD estimates (Waples & England 2011). Check structure first
  (popgen-structure). F_IS clearly above 0 is a warning sign.
- **Unlinked loci.** The LD formula assumes recombination 0.5. Tightly linked
  SNPs (same RAD contig, same chromosome arm) inflate R^2 and bias Ne
  downward; thin SNPs to one per contig and, if a map exists, use pairs on
  different chromosomes.
- **Discrete generations.** In iteroparous species a single cohort's LD
  estimate is closer to the number of breeders that produced it (Nb); a
  mixed-age sample gives roughly Ne per generation (Waples et al. 2014). The
  temporal method needs the number of generations between samples; give
  `generation_time` when times are in years.
- **Random sample of the population.** Over-representing families (many
  sibs from one nest) biases LD Ne downward; mark known juveniles/sibs with
  `use_for_ne = no` rather than deleting all putative sibs, which also causes
  bias (Waples & Anderson 2017).
- **Neutral markers, Hardy-Weinberg, no null alleles or genotyping error.**
  Null alleles and allelic dropout inflate F_IS and distort LD and relatedness.
- **Relatedness is relative to the sample.** Moment estimators measure
  relatedness relative to the sample's allele frequencies, so the average
  pair scores slightly below 0 and values are comparable only within the
  same reference. They assume no inbreeding.
- **Bottleneck test** assumes mutation-drift equilibrium before the
  bottleneck and the infinite-alleles model; the stepwise model is more
  conservative. Treat it as a weak, supporting signal.
  The P-value is a Monte Carlo test of the mean standardised difference
  (He - E[He|k]) / SD[He|k] (the standardised-differences statistic of
  BOTTLENECK, Piry et al. 1999) against its equilibrium distribution given each
  locus' allele number; the Wilcoxon test is also reported but was
  anti-conservative in our checks (VALIDATION.md).

## Data format

Genotype CSV, one row per individual (two-column format shown):

| sample_id | site | year | use_for_ne | Sco102_a | Sco102_b | Sco105_a | Sco105_b | ... |
|---|---|---|---|---|---|---|---|---|
| UEC-14-001 | UpperElkCreek | 2014 | yes | 154 | 162 | 0 | 0 | ... |
| UEC-24-101 | UpperElkCreek | 2024 | no  | 158 | 158 | 220 | 228 | ... |

- `0` (or any code in `missing`) = missing; both alleles must be missing together.
- Alleles can be fragment sizes, numbers or letters (A/C/G/T).
- `one_column`: one column per locus holding `154/162` (`allele_sep`).
- `dosage`: SNP allele counts 0/1/2 (NA missing).
- `genepop`: a standard .gen file (populations from Pop blocks or a `popmap`).
- `vcf`: biallelic SNP VCF plus `popmap` (columns id, pop[, time, use_for_ne]).

Known relatives (optional) `known_relationships.csv`: `id1, id2, relationship`
with relationship fullsib, halfsib, parent_offspring or unrelated.

## Config

```yaml
data:
  genotypes: bulltrout_microsats.csv
  format: two_column            # two_column | one_column | dosage | genepop | vcf
  id_column: sample_id          # default: first column
  pop_column: site              # optional; default one population
  time_column: year             # optional; two or more times -> temporal method
  sample_column: null           # optional; default pop_time
  use_for_ne_column: use_for_ne # optional yes/no
  first_locus_column: Sco102_a  # optional; otherwise all non-metadata columns
  missing: [0]
  allele_sep: "/"               # one_column only
  popmap: null                  # genepop/vcf: id,pop[,time,use_for_ne]
  known_kin: known_relationships.csv
  generation_time: 5            # time units per generation
  census_N: {UpperElkCreek: 120}  # optional: sampling plan I for the temporal method
  marker: msat                  # msat | snp (auto-detected if omitted)
settings:
  pcrit: auto        # or a number; auto = 0.02 if S > 25, else drop singletons
  jackknife: true    # jackknife CI for LD Ne
  mating: random     # random | monogamy
  max_loci: 3000     # LD method uses a random subset above this (speed)
  rel_primary: LR    # LR or QG for kin calls and figures
  kin_thresholds: {first: 0.354, second: 0.177}
  kin_sim_pairs: 300 # simulated pairs per class for misclassification rates (0 = skip)
  pcrit_grid: [0.01, 0.02, 0.05, 0.1]   # sensitivity of LD Ne to Pcrit
  bottleneck: false  # heterozygosity-excess test
  ci: 0.95
```

## Outputs

Estimates table (95% intervals):

- `Ne_LD_<sample>`: LD estimate of Ne (Inf = no drift signal detected).
- `het_loss_LD_<sample>` = 1/(2 Ne): expected proportion of heterozygosity lost
  per generation. This is the scale on which the estimator is nearly unbiased
  and finite even when Ne is infinite; Ne = 1/(2 x het_loss).
- `Ne_temporal[_<pop>]`, `het_loss_temporal[_<pop>]` (Jorde-Ryman) and the
  `_NT`, `_Pollak` variants.
- `relatedness_<class>_<LR|QG>`: mean r of known pairs of each class.
- `Ho_`, `He_`, `FIS_<sample>`.

Figures: see the Figures section below.

## Figures

Exported by `core/run.R fit` to `<out>/figures/` (PDF + PNG, `gallery.html`);
redraw without refitting with `core/run.R figures popgen-ne --fit <out>/popgen-ne_fit.rds --config <cfg>`.

| id | what it shows |
|---|---|
| `ne_estimates` | Ne with confidence intervals for every sample (LD) and sample pair (temporal: Jorde-Ryman, Nei-Tajima, Pollak), log scale; dashed reference lines at 50/500, dotted at 100/1,000; `>` marks an infinite upper limit. Colour and symbol = method. |
| `ne_by_maf` | LD estimate of Ne (with jackknife CI) against the allele-frequency cutoff Pcrit, one line per sample; shows whether rare alleles drive the estimate. |
| `relatedness_distribution` | Histogram of pairwise r (LR or QG, per `rel_primary`) per sample, pedigree expectations (0, 0.25, 0.5) and the kin-calling thresholds; known-relationship pairs as coloured rug ticks. |
| `kin_resolution` | Densities of r for unrelated, half-sib, full-sib and parent-offspring pairs simulated from the sample's allele frequencies, with the thresholds; overlap = misclassification. |
| `kin_network` | Individuals on a circle joined by putative first-degree (solid) and second-degree (dashed) kin links; line width proportional to r (first-degree only when links are too many). |
| `design_precision` | CV of 1/Ne against individuals sampled, by number of loci and true Ne, from the validation grid (`design_grid.csv`), with this study's sample size marked. |

Relabel anything through `figures: labels:` in the config: sample names
(e.g. `UpperElkCreek_2014`), population names, methods in the `ne_estimates` row labels (`short_LD`,
`short_JR`, `short_NT`, `short_Pollak`; the `method_*` keys name the colour
groups, whose legend is hidden because the row labels already say the method), relationship classes
(`unrelated`, `halfsib`, `fullsib`, `parent_offspring`, `first_degree`,
`second_degree`) and axes (`Ne`, `Ne_LD`, `pcrit`, `relatedness`,
`relatedness_sim`, `density`, `pairs`, `n_ind`, `cv_inv_ne`, `infinite`).
Relationship classes keep one colour and line type in every figure; greyscale
works (symbols and line types carry the groups; threshold lines turn black).
A global `figures: width:` other than `single` overrides the module's default
widths; `per_figure: {<id>: {width: ...}}` overrides both. Interval level is the one
the fit used (`settings: ci`, default 0.95); changing `figures: ci` alone
does not recompute intervals (a message says so). These are confidence
intervals, so there is no inner (50%) band.

## Simulation and power

`simulate()` runs a Wright-Fisher population of `Ne` individuals (random union
of gametes, so Ne = N), `n_loci` unlinked SNPs (`snp_maf_min`) or
microsatellites (`n_alleles`, Dirichlet `allele_dirichlet`), `burn_in`
generations, then samples `n_ind` individuals ("early"), and again after
`generations_between` generations ("late"), plus a kin panel of
`n_fullsib_pairs` and `n_halfsib_pairs` born in the late generation and
`n_unrelated_pairs` unrelated pairs from the late sample, with
`missing_rate`. Relatedness truth is the pedigree relatedness of each pair
relative to the sample (kinship tracked exactly). Design knobs for power:
`n_ind`, `n_loci`, `generations_between`, `Ne`. The `power` command reports
the probability that the CI for het_loss excludes 0, i.e. that the upper limit
of Ne is finite. `design_grid.csv` holds the CV of 1/Ne by S, loci and Ne.

## Pitfalls

- Infinite Ne is not an error: it means Ne is large relative to the sample.
  Report the lower limit. Sample more individuals rather than more loci when
  Ne >> S (see design_precision).
- Do not compare Ne from different Pcrit values or different marker sets.
- With fewer than ~25 individuals or ~10 microsatellites, LD estimates are
  very imprecise; the temporal method needs at least 2-3 generations.
- Relatedness calls from 10-15 microsatellites misclassify many unrelated
  pairs as half sibs: read the simulated misclassification table before
  acting on individual kin links.
- The simulated misclassification table is optimistic for half sibs: with
  1,000 SNPs in a population of Ne = 100, 10% of real half sibs fell below
  the second-degree threshold against 2% predicted (VALIDATION.md).

## Validation status (2026-10-19; details in VALIDATION.md)

- Default design (Ne 100, 1,000 SNPs, S = 60, 5 generations; nsim 200): PASS all 9 estimands
  (LD and temporal 1/(2Ne) bias -3.4% / -0.7%, coverage 0.91-0.97).
- Hard design (15 microsatellites, S = 30, Ne 200): PASS, but intervals for Ne nearly always reach infinity.
- Ne 500 with S = 50: LD bias in 1/Ne -12% for one sample (Ne overestimated ~10-15%): FAIL, documented.
- Negative controls fail as they should (wrong sampling plan; family over-representation).
- Matches the NeEstimator v2 C program exactly on its 8Ne50 worked example; Lynch-Ritland matches Demerelate.
- Wang (2002) relatedness estimator is not implemented (QG and LR only).
- Independent verification (separate sexes, no selfing): LD Ne biased upward 6-25% (gate FAIL for LD at 20 microsatellites and skewed sex ratios); temporal Ne and relatedness pass. Verdict VERIFIED_WITH_LIMITS.

## References

Waples 2006; Waples & Do 2008, 2010; Peel et al. 2013; Do et al. 2014 (NeEstimator v2);
Jones, Ovenden & Wang 2016; Nei & Tajima 1981; Pollak 1983; Waples 1989;
Jorde & Ryman 2007; Queller & Goodnight 1989; Lynch & Ritland 1999; Wang 2014;
Manichaikul et al. 2010; Cornuet & Luikart 1996; Piry et al. 1999; Franklin 1980; Frankham et al. 2014;
Franklin, Allendorf & Jamieson 2014;
Frankham 1995; Waples & England 2011; Waples, Antao & Luikart 2014;
Waples & Anderson 2017. Full citations in `module.R` (`references`).
