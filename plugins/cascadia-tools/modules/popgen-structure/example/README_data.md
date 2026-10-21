# Example data: coastal tailed frogs in eight Oregon Coast Range streams

> **These genotypes are SIMULATED.** No real frogs were sampled. The data were
> generated with a known population history so that the right answers are known
> and the analysis can be checked against them. Stream names and coordinates
> are real places, used only to give the example a realistic geography.

## The system

The coastal tailed frog (*Ascaphus truei*) lives in cold, fast, forested
headwater streams of the Pacific Northwest. Adults rarely leave the stream
corridor, and dispersal between catchments is limited. Logging, road
building and warming water can isolate populations. Managers want to know
which streams are genetically diverse, which are isolated, and whether
any population is small enough that it might need genetic rescue.

## The (simulated) study

A field crew swabbed 18-32 frogs (adults and late-stage tadpoles) in each of
eight streams in the central Oregon Coast Range, 6-81 km apart. Every frog
was genotyped at 12 microsatellite loci (Atru01-Atru12). About 3% of
genotypes failed and are recorded as missing.

The researcher's questions:

1. How genetically diverse is each stream population?
2. Are any populations inbred?
3. How different are the streams from one another, and does that grow with distance?
4. How many distinct genetic groups are there?
5. Is any stream isolated and losing diversity?

## Design

| dimension | value |
|---|---|
| populations (streams) | 8 |
| frogs genotyped per stream | NFSIU 30, LOBST 24, FIVEM 28, DRIFT 18, MARYS 32, ALSEA 26, YACHT 22, SMITH 20 (200 total) |
| loci | 12 microsatellites (dinucleotide or tetranucleotide repeats) |
| missing genotypes | 3% per genotype at random (realised 2.6%) |
| distances between streams | 6.2-80.8 km (great-circle); nearest neighbour median 15.2 km |
| habitat geometry | 2-D (several catchments), so isolation by distance uses log distance |

## Files

### `tailed_frog_microsats.csv` (200 rows, one per frog)

| column | meaning |
|---|---|
| `individual` | frog ID, `AT<stream>-<number>` |
| `population` | stream code (matches `tailed_frog_sites.csv`) |
| `Atru01` ... `Atru12` | first allele at each locus (fragment size in base pairs) |
| `Atru01.2` ... `Atru12.2` | second allele at the same locus |

Two columns per locus (GenAlEx style). `0` in both columns = missing genotype.
Alleles are fragment sizes: each locus has its own base size and repeat
step (2 or 4 bp).

### `tailed_frog_sites.csv` (8 rows, one per stream)

| column | meaning |
|---|---|
| `population` | stream code |
| `stream` | stream name |
| `latitude`, `longitude` | sampling site, decimal degrees (WGS84) |
| `elevation_m` | site elevation (m); descriptive only, not used in the analysis |

### Other files

* `config.yaml`: the analysis settings used by `core/run.R fit`.
* `recovery.yaml`: planted-truth check at this design (`core/run.R recover`).
* `true_values.csv`: the true values (see below).

## How the data were generated (TRUE values, simulated)

Script: `dev/popgen-structure/gen_example.R`, `set.seed(2026)`.
Re-running `dev/popgen-structure/si/si_truth.R` regenerates both CSVs
byte for byte (checked: identical).

* Forward Wright-Fisher simulation, 150 generations, diploid, random mating
  (no selfing), the module's `ps_wf_forward()`.
* Effective (census) size per stream: NFSIU 150, LOBST 150, FIVEM 120,
  **DRIFT 30**, MARYS 150, ALSEA 150, YACHT 120, SMITH 100.
* Dispersal: each generation 5% of a stream's parents are immigrants, spread
  over the other streams with weight exp(-distance / 12 km).
  **Drift Creek is nearly isolated: 0.4% immigrants per generation.**
* Mutation: K-allele model, 14 alleles per locus, rate 5 x 10^-4 per generation;
  starting allele frequencies drawn from a flat Dirichlet.
* Sampling: the frogs listed above drawn at random from the final generation;
  3% of genotypes set to missing.

### True values (`true_values.csv`)

"census_12loci" = the module's estimators applied to every frog alive at
the end of the simulation (1,040 frogs) at the 12 typed loci. "expected_model"
= exact expectation over many replicate histories (identity-by-descent
recursion). Genetic quantities are themselves random, so the two differ.

| quantity | census (all frogs) | expected (model) |
|---|---|---|
| global FIS (inbreeding) | -0.004 | 0 |
| mean He (gene diversity) | 0.673 | 0.660 |
| He, Drift Creek | 0.220 | 0.259 |
| He, other seven streams | 0.713-0.759 | 0.666-0.734 |
| pairwise FST with Drift Creek | 0.239-0.334 | 0.193-0.287 |
| pairwise FST among the other seven | 0.019-0.082 | 0.021-0.060 |
| global FST | 0.065 | 0.086 |
| IBD slope (FST/(1-FST) per unit log km) | 0.058 | 0.047 |

**Important note on global FST.** The census value (0.065) weights each
stream by how many frogs live there; Drift Creek is 3% of all frogs. The
Weir & Cockerham estimator weights each stream by how many frogs were
*genotyped*; Drift Creek is 9% of the sample. Because Drift Creek is by far
the most differentiated stream, the estimator targets a larger value at
this design. That design-weighted true value (census allele frequencies,
streams weighted as sampled; `dev/popgen-structure/si/si_design_truth.R`)
is **0.099** for global FST and **0.069** for the IBD slope. All three
"truths" are listed in `results/popgen-structure/si/true_values_with_design.csv`.

What a correct analysis should find: Drift Creek has far lower diversity
than the others and is strongly differentiated from every other stream;
there is no inbreeding within streams (FIS about 0); the other seven
streams are moderately differentiated (FST 0.02-0.08) and connected.
