# Worked example: effective population size of a bull trout population (module `popgen-ne`)

This is the Appendix S1 worked example and the website's "try an example" walkthrough.
**The data are simulated, and the true values are known** (see `README_data.md`). So
every estimate below can be checked against the real answer.

---

## (a) The question, in the researcher's words

> "I monitor a bull trout spawning population in Upper Elk Creek. Redd counts say
> there are about 120 spawning adults. We fin-clipped 50 adults in 2014 and 50 more
> in 2024, and genotyped them at 16 microsatellites. Is the population genetically
> big enough, or is it losing diversity? How fast? In 2024 we also took two juveniles
> from each of six redds. Can the markers tell me which fish are related?"

## (b) The files the researcher has

`bulltrout_microsats.csv`: one row per fish, 112 rows. There are two allele columns per locus, and `0` means missing.
The first 11 of 37 columns are shown:
```
"sample_id","site","year","life_stage","use_for_ne","Sco102_a","Sco102_b","Sco105_a","Sco105_b","Sco106_a","Sco106_b"
"UEC-14-001","UpperElkCreek",2014,"adult","yes",226,226,180,184,282,282
"UEC-14-002","UpperElkCreek",2014,"adult","yes",226,238,196,224,262,282
"UEC-14-003","UpperElkCreek",2014,"adult","yes",202,210,196,212,286,286
"UEC-14-004","UpperElkCreek",2014,"adult","yes",210,226,196,196,286,286
```

`known_relationships.csv`: the juvenile nest-mates, 6 rows:
```
"id1","id2","relationship","note"
"UEC-24-101","UEC-24-102","fullsib","juveniles from redd R01"
"UEC-24-103","UEC-24-104","fullsib","juveniles from redd R02"
"UEC-24-105","UEC-24-106","fullsib","juveniles from redd R03"
"UEC-24-107","UEC-24-108","fullsib","juveniles from redd R04"
```

## (c) What the researcher types to Claude

> "I have microsatellite genotypes for bull trout from Upper Elk Creek in
> bulltrout_microsats.csv: 50 adults sampled in 2014 and 50 in 2024, 16 loci, two
> columns per locus, 0 = missing. There are also 12 juveniles from 2024, two from each of six
> redds; known_relationships.csv lists the nest-mate pairs. Generation time is about
> 5 years and redd counts suggest about 120 adults. Please estimate the effective
> population size and how fast the population is losing genetic diversity, and tell
> me whether the markers can identify relatives. I'd like a report with figures for a
> paper."

## (d) What the agent does, step by step (eco-analyze)

### 1. Project folder
The agent makes `projects/<name>/` with the raw data, a plan, a config and a log.
In this example, the same files live in `modules/popgen-ne/example/` (inputs) and `results/popgen-ne/si/` (outputs).

### 2. Data intake: what the agent reports back
- 112 fish: 50 adults in 2014, 50 adults in 2024, and 12 juveniles from 2024 (6 redds x 2).
- 16 microsatellite loci with 5 to 9 alleles each. 57 of 1,792 genotypes are missing (3.2%).
- The juveniles are full siblings. Siblings in an Ne sample make Ne look too small, so they are used only for relatedness (`use_for_ne = no`).
- 2014 to 2024 at 5 years per generation is 2 generations.
- Census N = 120 at the first sampling. This lets the temporal method use "sampling plan I" (fish are sampled before they spawn, and the census size is known).

### 3. Choosing the module
`popgen-ne` covers the question. It gives:
- single-sample linkage-disequilibrium Ne (LDNe), one estimate per year;
- temporal Ne from the change in allele frequencies between 2014 and 2024;
- pairwise relatedness and close-kin calls;
- heterozygosity and F_IS.

### 4. Plan (written before fitting)
- **Estimands:** the heterozygosity lost per generation, 1/(2Ne), for each LD sample and for the temporal method. Each is also reported as Ne.
  The module validates on the 1/(2Ne) scale because it stays finite when Ne is very large.
  Also: the mean relatedness of the known full sibs, and F_IS in each year.
- **Settings:** Pcrit = 0.02 (alleles rarer than 2% are dropped; Waples & Do 2010). Jackknife confidence intervals. Lynch-Ritland relatedness as the main estimator.
- **Adequate fit:** F_IS interval includes 0 (random mating). Known full sibs score near 0.5. LD estimates are stable across Pcrit.

### 5. Planted truth at this design (the gate, step 1)
The agent simulated 200 populations with Ne = 120 and the researcher's exact design.
The design was 2 samples of 50, 16 microsatellites, 2 generations apart, 6 full-sib pairs and 3% missing data.
It refitted each one (0.2 s per replicate, 200/200 ran). Gate rules: bias within 10% of the truth (or 3 Monte Carlo SE), coverage of 95% intervals at least 0.904.

| estimand | truth | mean estimate | relative bias | coverage | power* | gate |
|---|---|---|---|---|---|---|
| het. loss, LD 2014 sample | 0.00417 | 0.00362 | **-13.1%** (MC SE 1.8%) | 0.960 | 0.925 | **FAIL (bias)** |
| het. loss, LD 2024 sample | 0.00417 | 0.00369 | **-11.5%** (MC SE 1.8%) | 0.955 | 0.925 | **FAIL (bias)** |
| het. loss, temporal | 0.00417 | 0.00423 | +1.6% (MC SE 3.4%) | 0.920 | 0.835 | PASS |
| full-sib relatedness (Lynch-Ritland) | 0.469 | 0.459 | -2.1% | 0.975 | n/a | PASS |
| full-sib relatedness (Queller-Goodnight) | 0.469 | 0.487 | +3.9% | 0.960 | n/a | PASS |

*Power = share of replicates whose interval gives Ne a finite upper limit.

**In one sentence:** at this design, the temporal estimate and relatedness are recovered without bias, but each single-year LD estimate overstates Ne by about 15% on average.
That is a harmonic mean of about 137 when the truth is 120. Its intervals still cover the truth 96% of the time.
The cause is a documented limit: each sample of 50 fish is a large share (42%) of a population of 120.
The agent tells the researcher this before fitting. It names the temporal estimate as the main result and treats the LD estimates as supporting evidence that runs slightly high.

### 6. Fit
```
scripts/Rscript core/run.R fit popgen-ne --config modules/popgen-ne/example/config.yaml --out results/popgen-ne/si/fit
```
Moment estimators are used, so there is no MCMC. The full fit with figures takes under a minute.

### 7. Checks: convergence and goodness of fit
- **Numerical health:** the report says "Model converged". Every estimate and interval is finite, and nothing sits on a boundary.
  The 2024 LD upper limit is 500, which is wide but finite.
- **Random mating (Hardy-Weinberg):**
  - 2014: F_IS = 0.028 (95% CI -0.026 to 0.081).
  - 2024: F_IS = 0.004 (95% CI -0.025 to 0.032).
  - Both intervals include 0, so there is no sign of inbreeding, null alleles or mixed populations.
- **Positive control:** the 6 known full-sib pairs have mean r = 0.494 (95% CI 0.296 to 0.693). The pedigree expectation is 0.5.
- **Pcrit sensitivity:**
  - 2014 LD Ne ranges from 64 to 177 across Pcrit 0.01-0.10.
  - 2024 LD Ne ranges from 137 to 191.
  - The 2014 sample is the more sensitive one. The default Pcrit of 0.02 is the recommended value for samples above 25 fish.
- **Agreement among methods:** the three temporal estimators agree within 5%: Jorde-Ryman 122, Nei-Tajima 127, Pollak 124.

### 8. Interpretation: excerpts from the report (`popgen-ne_report.html`)
> "Temporal method (population UpperElkCreek, samples at times 2014 and 2024, 2.0 generations apart, sampling plan I, census N = 120): Jorde & Ryman's Fs' = 0.0082 gives Ne = 122 (95% CI 55.4 to 407; jackknife over loci (chi-square)). Nei-Tajima and Pollak estimators give 127 and 124."

> "Sample UpperElkCreek_2014 (50 individuals, 16 loci after excluding alleles rarer than 0.020): the linkage-disequilibrium estimate of Ne is 121 (95% CI 74.4 to 271; jackknife over individuals). An Ne of this size loses about 0.41% of its heterozygosity per generation (1/(2Ne)), about 4.1% over 10 generations, unless immigration or mutation replaces it."

> "The interval lies above 50 but reaches below 100: the population passes the classic short-term threshold (50) but not clearly the revised one (100; Frankham et al. 2014)."

> "Relatedness in sample UpperElkCreek_2024 (1891 pairs, 16 loci, Lynch-Ritland estimator): 15 pairs have r >= 0.354 (likely first-degree relatives: full sibs or parent-offspring) and 86 have 0.177 <= r < 0.354 (likely second-degree, e.g. half sibs). [...] about 7.0% of truly unrelated pairs would be called second-degree or closer (~132 false calls among 1891 pairs); 80% of full sibs and 46% of half sibs fall in their own class."

> "Known fullsib pairs (n = 6) have mean r = 0.494 (95% CI 0.296 to 0.693), against a pedigree expectation of 0.50."

> "Ne is not census size: Ne/N averages roughly 0.1-0.2 in wild populations (Frankham 1995) [...] In ecotools' independent validation with separate sexes (no selfing), LD estimates ran about 5-25% above the true Ne [...] treat a single LD estimate as slightly optimistic."

The agent checked every number in these sentences against `popgen-ne_estimates.csv`.

**What the agent tells the researcher (3-6 sentences):**
- The population's effective size is about 120. The best estimate is temporal Ne 122 (95% CI 55 to 407).
- At that size it loses about 0.4% of its heterozygosity per generation. That is about 4% over 10 generations (50 years) without immigration.
- The two single-year LD estimates agree (121 and 139). At this design, though, they run about 15% high, so they are supporting evidence only.
- The interval clears Ne = 50 but not clearly 100. So the data cannot yet say whether the population meets the revised short-term guideline.
- 16 microsatellites separate full sibs from unrelated fish well (80% correct class). Half sibs are only 46% correct, so treat single second-degree calls as tentative.

### 9. Figures (in `results/popgen-ne/si/fit/figures/`, PNG and PDF)

| id | file | caption (shortened) |
|---|---|---|
| `ne_estimates` | `ne_estimates.png` | Ne by sample and method with 95% CIs, log scale. Dashed lines at 50/500 (Franklin 1980), dotted at 100/1,000 (Frankham et al. 2014). **Key figure.** |
| `ne_by_maf` | `ne_by_maf.png` | Sensitivity of the LD estimate to the minimum allele frequency Pcrit, with 95% jackknife CIs, one line per year. |
| `relatedness_distribution` | `relatedness_distribution.png` | Pairwise Lynch-Ritland relatedness in each year. Grey lines: expectations 0, 0.25, 0.5. Red lines: kin-calling thresholds 0.177 and 0.354. Ticks: the known full-sib pairs. |
| `kin_resolution` | `kin_resolution.png` | Relatedness of simulated unrelated, half-sib, full-sib and parent-offspring pairs at these markers. Where the classes overlap, misclassification happens. |
| `kin_network` | `kin_network.png` | The 35 fish linked by first-degree relatedness (r >= 0.354). Line width is proportional to r. |
| `design_precision` | `design_precision.png` | Precision (CV of 1/Ne) of the LD method against the number of fish sampled, for 10 and 20 microsatellites and true Ne of 50, 100 and 500. The red line marks this study (50 fish, 16 loci). |

The agent opened every PNG. Labels, legends and scales read correctly, so the config's `figures:` block needed no changes.

### 10. The acceptance gate verdict: **OPEN**

| gate step | result |
|---|---|
| 0. plan declared before fitting | yes (section 4) |
| 1. planted truth at this design | temporal Ne and relatedness PASS; **LD Ne (both years) FAIL bias, -13% and -12% in 1/(2Ne)** |
| 2. numerical health | finite estimates, no boundaries |
| 3. model checks | F_IS includes 0 in both years; known full sibs r = 0.49 |
| 4. independent route | module-level: LD and temporal formulas identical to NeEstimator v2 (VALIDATION.md section 8) |
| 5. negative control | module-level: the gate fails when siblings are left in the Ne sample or the census is withheld (VALIDATION.md section 4) |
| 6. reproducible | commands below. The data generator reproduces the CSVs byte for byte. |

The temporal Ne, heterozygosity loss and relatedness results are fit to report.
The single-year LD estimates are not, because they fail the bias rule at this design.
The agent says so plainly: "With 50 fish from a population of about 120, the LD method overstates Ne by about 15% on average. Use the temporal estimate (122, 95% CI 55-407) as the answer. Read the LD estimates as upper-leaning."
Sampling more fish does not fix this here, because the sample is already a large share of the population.
The fix is to rely on the temporal method or to add sampling years.

## (e) Estimates versus the known truth (true Ne = 120)

| quantity | estimate (95% CI) | truth | interval covers truth? |
|---|---|---|---|
| temporal Ne (Jorde-Ryman) | 122 (55 to 407) | 120 | yes |
| temporal Ne (Nei-Tajima) | 127 (65 to 327) | 120 | yes |
| temporal Ne (Pollak) | 124 (64 to 311) | 120 | yes |
| LD Ne, 2014 | 121 (74 to 271) | 120 | yes |
| LD Ne, 2024 | 139 (76 to 500) | 120 | yes |
| het. loss per generation, temporal | 0.0041 (0.0012 to 0.0090) | 0.00417 | yes |
| het. loss per generation, LD 2014 | 0.0041 (0.0018 to 0.0067) | 0.00417 | yes |
| het. loss per generation, LD 2024 | 0.0036 (0.0010 to 0.0066) | 0.00417 | yes |
| full-sib relatedness, Lynch-Ritland | 0.494 (0.296 to 0.693) | 0.465 (realised) | yes |
| full-sib relatedness, Queller-Goodnight | 0.516 (0.379 to 0.652) | 0.465 (realised) | yes |

All intervals cover the truth in this one data set. The 2024 LD estimate is high, as the planted-truth runs predicted.
A single data set cannot show bias. The 200 replicates in step 5 can.

## (f) Limits from VALIDATION.md that matter at this design

- **The LD method overstates Ne when the sample is a large share of the population.** Here that share is 50/120. In planted truth, LD 1/Ne was biased -12% to -13% (Ne about +15%). This is why the verdict is OPEN for LD.
- **The LD method overstates Ne in species with separate sexes** by 6-25% (independent verification). The simulator used here has no separate sexes, but bull trout do, so the real LD bias could be larger. The temporal method is unaffected.
- **Overlapping generations have not been validated.** Bull trout spawn in several years, so an LD estimate from one sample is closer to the number of breeders (Nb) than to per-generation Ne (Waples et al. 2014).
- **The temporal method with fewer than 3 generations between samples** was outside the module validation. This example's own planted-truth run covers 2 generations at this design (bias +1.6%, coverage 0.92, PASS).
- **The temporal CI depends on a simulation-calibrated design effect.** With 15-20 microsatellites coverage was 0.90-0.92 (nominal 0.95). Here it was 0.92.
- **Missing data with microsatellites:** LD Ne can differ from NeEstimator by a median 1.7% (max 12%), because of how alleles are screened.
- **Half-sib calls are less reliable than the simulated-pair table suggests.** Read second-degree kin calls as tentative.
- **Not validated:** migration or population structure, null alleles or genotyping error, selection.
- The module estimates Ne, never census N. Wild Ne/N is typically 0.1-0.2. The simulated Ne/N here is 1.

## (g) Commands to reproduce

Run from `.`. Use `scripts/ejob` on this cluster.

```
# 1. regenerate the example data (byte-for-byte check) and write the truth table
scripts/ejob -c 1 -m 4 -t 0:20:00 -l logs/popgen-ne_si_truth.log -- scripts/Rscript dev/popgen-ne/si_truth.R

# 2. planted truth at the example's design (200 replicates)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/popgen-ne_si_recover.log -- scripts/Rscript core/run.R recover popgen-ne --config modules/popgen-ne/example/recovery.yaml --nsim 200 --cores 16 --out results/popgen-ne/si/recovery

# 3. fit the example data: report + figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/popgen-ne_si_fit.log -- scripts/Rscript core/run.R fit popgen-ne --config modules/popgen-ne/example/config.yaml --out results/popgen-ne/si/fit
```

Outputs:
- `results/popgen-ne/si/recovery/popgen-ne_recover_summary.csv` (gate table)
- `results/popgen-ne/si/fit/popgen-ne_report.html`
- `results/popgen-ne/si/fit/popgen-ne_estimates.csv`
- `results/popgen-ne/si/fit/figures/`
- `results/popgen-ne/si/si_summary.json`
