---
id: popgen-structure
title: Genetic diversity and population structure (FST, FIS, isolation by distance, DAPC)
when_to_use: Multilocus genotypes (microsatellites, or SNPs coded 0/1/2 or in a VCF) from individuals sampled in several populations or sites; questions about genetic diversity per population (Ho, He, allelic richness, private alleles), inbreeding / Hardy-Weinberg departures (FIS), how differentiated populations are (global and pairwise FST), whether differentiation grows with distance (isolation by distance), and how many genetic clusters there are (PCA, DAPC). Connectivity and genetic-rescue questions in conservation.
not_for: Effective population size or relatedness/kinship (use popgen-ne); assigning individuals to known source populations for forensics; landscape-resistance models (circuitscape/ResistanceGA); haploid, polyploid or dominant markers (AFLP); selection scans / outlier loci; datasets with one population (no structure to estimate) or fewer than about 3 typed individuals per population.
data_needed: genotype table (one row per individual with population and loci as 'a/b', two columns per locus, or SNPs 0/1/2) or a VCF plus a sample sheet; optional site coordinates (latitude/longitude or km) for isolation by distance and maps
engine: hand-written Weir & Cockerham (1984) / Nei (1987) estimators (agree with hierfstat to 1e-11), delete-one-locus jackknife + shared-individual covariance (delete-one-individual jackknife) intervals; adegenet/ade4 (PCA, DAPC), vegan (Mantel), poppr (AMOVA, optional)
---

# Genetic diversity and population structure

## Question it answers

* **How diverse is each population?** Observed heterozygosity (Ho), expected
  heterozygosity / gene diversity (He), allelic richness rarefied to a common
  sample size (AR) and private alleles. Low He and AR flag small or isolated
  populations that are candidates for monitoring and genetic rescue.
* **Are individuals inbred?** FIS per population and overall, plus
  Hardy-Weinberg exact tests per population and locus.
* **How different are populations?** Global FST and FST for every pair of
  populations, with 95% intervals. Island-model "effective migrants per
  generation" (Nm) is reported as an order-of-magnitude translation only.
* **Does differentiation grow with distance?** Regression of pairwise
  FST/(1-FST) on log distance (2-D habitat) or distance (1-D habitat such as
  a river), with an interval for the slope and a Mantel test.
* **How many genetic groups are there?** PCA of individuals, k-means
  clustering with BIC to choose K, DAPC membership probabilities.
* **Optional:** AMOVA (PhiST) with a permutation test.

## Model and estimators

Notation: populations i = 1..r, loci l, alleles u, n_i genotyped individuals,
p_iu allele frequency, h_iu proportion of individuals heterozygous for allele u.

* **Ho** = proportion of heterozygous individuals, averaged over loci.
* **He** (Nei 1987 unbiased gene diversity):
  He_il = n/(n-1) [1 - sum_u p_u^2 - Ho_il/(2n)], averaged over loci.
  It estimates the probability that two genes from *different* individuals
  of the population are different alleles.
* **AR** (El Mousadik & Petit 1996; Hurlbert 1971): expected number of
  alleles in a subsample of g genes, g = the smallest population sample
  (2 x min n), AR = sum_u [1 - C(2n - x_u, g) / C(2n, g)].
* **FST, FIS** (Weir & Cockerham 1984): variance components a (among
  populations), b (among individuals within populations), c (within
  individuals) per allele, summed over alleles and loci;
  FST = theta = sum a / sum (a + b + c), FIS = f = 1 - sum c / sum (b + c).
  Pairwise FST uses the same formulas on each pair of populations. These
  are ratio-of-sums estimators: unbiased for the ratio of expected identity
  differences, (Q_within - Q_between) / (1 - Q_between), the parameter that
  island and stepping-stone theory predicts. Results agree with
  `hierfstat::wc`, `pairwise.WCfst`, `basic.stats`, `allelic.richness`
  to < 1e-10 (VALIDATION.md, section 6).
* **Isolation by distance** (Rousset 1997): OLS slope of FST/(1-FST) on
  ln(distance) for 2-D habitats (`ibd_distance: log`) or distance for 1-D
  habitats (`linear`). Under Rousset's model the slope estimates 1/(4 pi D sigma^2)
  (2-D), i.e. 1 / neighbourhood size, or 1/(4 D sigma^2) per unit distance (1-D).
  Mantel test (vegan, Pearson, 999 permutations) for the association.
* **Uncertainty.** Loci are the independent replicates of the evolutionary
  process, so FST and He vary from locus to locus far more than sampling of
  individuals alone suggests (Weir 1996). But loci are *not* fully
  independent: they share the same individuals and pedigree (relatives in
  the sample, inbred individuals). The module therefore uses
  SE^2 = SE^2(delete-one-locus jackknife) + max(0, V_joint - V_perlocus),
  where V_joint is the delete-one-individual jackknife variance (stratified
  by population) when an individual is removed at all loci at once, and
  V_perlocus the same variance when the individual is removed at one locus
  at a time (summed over loci). Their difference is the cross-locus
  covariance created by sharing individuals (identity disequilibrium from
  inbreeding/selfing, kinship). Intervals are estimate +/- t(L-1, 0.975) x SE
  (L = number of loci). A bootstrap of individuals was tried first and
  rejected: the number of duplicated individuals is common to all loci and
  inflated FST intervals near 0. Simulation showed that the usual percentile
  bootstrap over loci alone undercovers (0.87-0.92, and 0.66-0.68 for FIS
  with selfing); see VALIDATION.md.
* **HWE** per population x locus: exact test for biallelic loci (Wigginton
  et al. 2005); Monte Carlo conditional exact test for multiallelic loci
  (Guo & Thompson 1992; permutation of alleles among individuals);
  Benjamini-Hochberg correction within each population.
* **PCA** on centred allele frequencies (adegenet `tab`, missing = mean;
  ade4 `dudi.pca`). **Clusters**: `adegenet::find.clusters` (k-means on all
  PCs, BIC for K = 1..max_k, K chosen by `k_criterion`, default
  "diffNgroup"); **DAPC** with K - 1 retained PCs (Thia 2023: only K - 1 PC
  axes carry among-group signal; more axes overfit) and n.da = K - 1.
* **AMOVA** (optional; Excoffier et al. 1992) via `poppr::poppr.amova`
  with an ade4 permutation test.

## Assumptions (and what breaks)

* **Diploid, codominant, neutral, unlinked loci.** Linked SNPs (RAD loci on
  the same contig, dense SNP chips) are not independent: thin to one SNP per
  contig / LD-prune first, or intervals are too narrow.
* **Populations are biological units sampled at random.** Sampling families
  (siblings) inflates FST and FIS; pooling two breeding groups as one
  "population" produces a heterozygote deficit (Wahlund effect).
* **Null alleles / allele dropout** (microsatellites) produce FIS > 0
  concentrated in a few loci and in HWE tests; the report flags loci that
  deviate in more than half the populations.
* **FST is not migration.** Converting FST to Nm assumes Wright's island
  model at equilibrium (equal sizes, equal migration, no selection or
  mutation); real systems rarely fit, and FST lags behind recent changes in
  connectivity (Whitlock & McCauley 1999). Report it only as an order of
  magnitude. A small, drifting population can have high FST to all others
  even with ongoing gene flow.
* **IBD** needs at least ~6-8 populations spread over a range of
  distances; one outlier population (very small or isolated) can create or
  hide a distance pattern. The report states when the slope interval and
  the Mantel test disagree.
* **DAPC/k-means** find discrete groups even under continuous IBD and can
  merge weakly differentiated populations; treat K as a summary, not as
  management units.

## Data format

**Genotype table** (CSV or tab-delimited), one row per individual.

Two columns per locus (GenAlEx-style; 0 = missing). The second column of
each pair can have any name; locus names are taken from the first column
(suffixes `_a/_b`, `.1/.2`, `_1/_2` are dropped):

```
individual,population,Atru01,Atru01.2,Atru02,Atru02.2
ATNFSIU-01,NFSIU,158,158,240,244
ATNFSIU-02,NFSIU,158,162,0,0
```

One column per locus with 'a/b' (also 'a|b'; missing '0/0', NA, '?', '-9'):

```
individual,population,Loc1,Loc2
F001,Upper,152/156,201/201
```

SNPs coded 0/1/2 (count of the alternative allele; set `format: snp012`,
missing NA or -9; 0 is a genotype here, not missing):

```
individual,population,snp0001,snp0002
S01,North,0,2
```

**VCF** (`vcf:`), read with vcfR; GT field only. Populations come from the
metadata sheet (columns `individual`, `population`).

**Metadata / coordinates** (`metadata:`), optional: one row per population
(population, latitude, longitude) or per individual (individual, latitude,
longitude; averaged per population). Projected coordinates must be in km
(`x_col`, `y_col`). Distances from lat/long are great-circle (haversine).

## Config

```yaml
data:
  genotypes: tailed_frog_microsats.csv   # or vcf: my.vcf
  format: two_column                     # auto | two_column | slash | snp012
  id_col: individual
  pop_col: population
  exclude: [sex, date]                   # non-genotype columns
  missing: ["000"]                       # extra missing codes
  marker: microsat                       # microsat | snp (default from format)
  metadata: sites.csv
  lat_col: latitude                      # or x_col / y_col in km
  lon_col: longitude
  ibd_distance: log                      # log (2-D habitat) | linear (1-D: rivers, coasts)
settings:
  individual_covariance: true   # add the shared-individual covariance to the loci jackknife SE
  ci: 0.95
  hwe: true          # Hardy-Weinberg tests (hwe_B Monte Carlo permutations, default 1000)
  pca: true
  dapc: true
  k_criterion: diffNgroup   # find.clusters rule: diffNgroup | min | goodfit | goodfitbic ...
  max_k: 20                 # default max(2 x number of populations, 10)
  mantel_perm: 999          # 0 = skip Mantel
  amova: false              # true = poppr AMOVA (slow for > 1000 individuals)
  seed: 1
figures:
  labels: {population: "Stream", NFSIU: "North Fork Siuslaw"}
```

## Outputs

Estimates table (`<id>_estimates.csv`): `Fst_global`, `Fis_global`,
`He_mean`, `Ho_mean`; per population `He_<pop>`, `Ho_<pop>`, `AR_<pop>`,
`Fis_<pop>`, `PrivateAlleles_<pop>` (count, no interval); pairwise
`Fst_<pop1>_vs_<pop2>`; `ibd_slope` (if coordinates); `Nm_island_model`
(derived, island-model translation, interval from the FST interval). The
report adds HWE results, Mantel r and P, AMOVA, K and plain-language
interpretation in management terms. Figures: see [Figures](#figures).

## Figures

`core/run.R fit` exports these to `<out>/figures/` (PDF + PNG and a
`gallery.html`); `core/run.R figures popgen-structure --fit <out>/popgen-structure_fit.rds
--config <cfg>` redraws them from the saved fit after editing the `figures:` block.

| id | What it shows |
|---|---|
| `diversity` | Ho, He (Nei unbiased) and rarefied allelic richness per population: point = mean over loci, bar = `ci` confidence interval (t(L-1) x jackknife SE). |
| `fst_matrix` | Lower-triangle heatmap of pairwise Weir & Cockerham FST (viridis, labelled cells; intervals in the estimates table). |
| `ibd` | FST/(1-FST) per population pair vs geographic distance (log axis when `ibd_distance: log`), least-squares line and `ci` confidence band; Mantel r and P in the caption. Needs coordinates. |
| `pca` | PC1 vs PC2 of individual genotypes, colour + symbol = sampled population, % of variance on the axes. Needs `pca: true`. |
| `bic_k` | BIC vs number of k-means clusters K; chosen K marked in orange. Needs `dapc: true`. |
| `dapc_assignment` | STRUCTURE-style bar plot of DAPC membership probabilities per individual, grouped by population. A cluster made up mainly of one population takes that population's PCA colour; clusters shared by several populations are grey. |
| `sample_map` | Sampling locations (lon/lat in degrees, or km), fill = He, size = individuals genotyped, labelled by population. Needs coordinates. |

Customisation: `labels:` keys are `population` (legend/axis title), each
population code (e.g. `NFSIU: "North Fork Siuslaw"`), `Ho`, `He`, `AR` (facet
titles), `diversity_x`, `fst`, `distance`, `ibd_y`, `pca_x`, `pca_y`, `bic_x`,
`bic_y`, `dapc_x`, `dapc_y`, `cluster` (legend title), `cluster_<k>` (one
cluster), `map_x`, `map_y`, `map_fill`, `map_size`. `ci:` redraws every interval
at that level from the stored SEs without refitting. `greyscale: true` keeps
populations distinguishable by symbol. With more than 8 populations the colours
fall back to a viridis ramp (use the `diversity` and `fst_matrix` figures, which
do not depend on colour). Geographic distance is in km (from coordinates), so
no `units:` back-transformation applies; nothing in this module is standardised.

## Simulation and power

`simulate()` is a Wright-Fisher simulator of r demes of Ne diploid
monoecious individuals (random union of gametes; optional partial selfing)
with backward migration matrix M: `island` (rate m spread over all other
demes), `stepping_stone_1d`, `stepping_stone_2d` (lattice, rate m spread
over 2 / 4 neighbours, reflecting edges), or `panmixia` (one population
split into sampling sites). Microsatellites: K-allele mutation (rate mu,
K = n_alleles), founders from Dirichlet(alpha) frequencies. SNPs: two
alleles, no mutation. After T_gen generations, n_ind individuals per deme
are sampled at n_loci loci.

Two definitions of truth are available (`design$truth`):

* `realised` (default): the genome-wide values of the simulated
  populations themselves, from a census of every individual at the sampled
  plus `n_shadow_loci` (300) extra loci on the same pedigree. This is what
  a manager wants to know about *these* populations.
* `expected`: the expectation over replicate histories, computed exactly
  from the identity-by-descent recursion of the same life cycle
  (F' = s sum_k M_ik (1+F_k)/2 + (1-s) Q'_ii, Q' = M H M^T, ...; see
  simulate.R). Used to show that the estimators are unbiased and that the
  simulator reproduces Wright's island-model FST.

Params: `Ne`, `m`, `selfing`, `mu`, `n_alleles`, `alpha`, `snp_alpha`,
`T_gen`. Design: `model`, `n_pops`, `spacing_km`, `n_ind`, `n_loci`,
`marker`, `missing_rate`, `ibd_distance`, `truth`, `n_shadow_loci`.

Power: `power_params` are `Fst_global` (power = interval excludes 0, i.e.
structure detected), `ibd_slope` and `Fis_global`. Typical knobs:
`--vary n_loci=...` or `n_ind=...` at the expected FST (set via `m`).
Rough guide (island model, 6 populations, 20 individuals each,
FST ~ 0.01): see VALIDATION.md section 7.

## Pitfalls

* Allelic richness and private-allele counts depend on sample size; AR is
  rarefied, private-allele counts are not.
* Per-population FIS intervals are tested separately: with many
  populations some will exclude 0 by chance.
* Negative FST estimates are possible when true FST is ~0; they mean "no
  detectable differentiation".
* Monomorphic loci add nothing to FST but count as He = 0; drop them first
  if He must be comparable with studies that excluded them.
* For SNP panels, remove loci with very low minor-allele count and
  LD-prune before analysis; the module does not filter loci for you.

## References

Benjamini & Hochberg (1995) J R Stat Soc B 57:289. El Mousadik & Petit
(1996) Theor Appl Genet 92:832. Excoffier, Smouse & Quattro (1992)
Genetics 131:479. Frankham et al. (2011) Conserv Biol 25:465. Goudet (2005)
Mol Ecol Notes 5:184. Guo & Thompson (1992) Biometrics 48:361. Hurlbert
(1971) Ecology 52:577. Jombart (2008) Bioinformatics 24:1403. Jombart,
Devillard & Balloux (2010) BMC Genetics 11:94. Kamvar et al. (2014) PeerJ
2:e281. Knaus & Grunwald (2017) Mol Ecol Resour 17:44. Mantel (1967)
Cancer Res 27:209. Mills & Allendorf (1996) Conserv Biol 10:1509. Nei
(1987) Molecular Evolutionary Genetics. Rousset (1996) Genetics 142:1357.
Rousset (1997) Genetics 145:1219. Thia (2023) Mol Ecol Resour 23:523.
Weir (1996) Genetic Data Analysis II. Weir & Cockerham (1984) Evolution
38:1358. Whitlock & McCauley (1999) Heredity 82:117. Wigginton, Cutler &
Abecasis (2005) Am J Hum Genet 76:887. Wright (1931) Genetics 16:97.
