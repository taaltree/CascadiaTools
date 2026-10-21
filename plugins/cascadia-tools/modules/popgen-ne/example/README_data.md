# Example data: bull trout, Upper Elk Creek (SIMULATED)

**These data are simulated.** No real fish were sampled. The values were generated
with the module's own Wright-Fisher simulator so that the true answer is known. The
place, the people and the redd numbers are invented for the story.

## The (simulated) study

Bull trout (*Salvelinus confluentus*) are listed as threatened in the United States.
They are cold-water char that spawn in small headwater streams. A biologist monitors
one spawning population, Upper Elk Creek. Redd counts suggest about 120 spawning
adults.

The biologist fin-clipped 50 adults in 2014 and 50 adults in 2024. All 100 fish were
genotyped at 16 microsatellite loci. In 2024 the crew also dug six redds and took
two juveniles from each one. Juveniles from the same redd are nest-mates, so they are
treated as known full siblings.

Bull trout mature at about 5 years, so 2014 to 2024 is about 2 generations.

The biologist wants to know three things:
1. How large is the genetically effective population (Ne)?
2. How fast is the population losing genetic diversity?
3. Which sampled fish are close relatives, and can the markers tell kin apart?

## Design

| dimension | value |
|---|---|
| population | 1 (Upper Elk Creek) |
| samples | 2 (2014 and 2024) |
| adults per sample (used for Ne) | 50 and 50 |
| juveniles (2024, known full sibs, not used for Ne) | 12 (6 redds x 2) |
| total genotyped fish | 112 |
| loci | 16 microsatellites, 5-9 alleles each in the data |
| missing genotypes | 57 of 1,792 (3.2%) |
| generations between samples | 2 (10 years at 5 years per generation) |
| census size at first sampling | 120 adults (from redd counts) |

## Files

### `bulltrout_microsats.csv` (112 rows, one per fish)

The two-column layout used by GenAlEx and CERVUS: each locus has two columns, one per allele.

| column | meaning |
|---|---|
| `sample_id` | fish ID. `UEC-14-...` = 2014, `UEC-24-...` = 2024. `UEC-24-101` to `UEC-24-112` are the juveniles. |
| `site` | population name (`UpperElkCreek`) |
| `year` | sampling year (2014 or 2024) |
| `life_stage` | `adult` or `juvenile` |
| `use_for_ne` | `yes` = use for Ne; `no` = juvenile nest-mate, use for relatedness only |
| `Sco102_a`, `Sco102_b`, ... `Smm22_a`, `Smm22_b` | the two allele sizes (bp) at each of 16 loci: Sco102, Sco105, Sco106, Sco107, Sco109, Sco200, Sco202, Sco212, Sco215, Sco216, Sco218, Sco220, Omm1128, Omm1130, Sfo18, Smm22. `0` = missing. |

The locus names follow real salmonid microsatellite naming. The allele sizes are
relabelled simulated alleles (a different size offset per locus, 4 bp steps). The
labels do not change any estimate.

### `known_relationships.csv` (6 rows)

| column | meaning |
|---|---|
| `id1`, `id2` | the two juveniles from one redd |
| `relationship` | `fullsib` |
| `note` | which redd (R01-R06) |

### `config.yaml`

Tells ecotools which columns hold what: genotype format `two_column`, missing code
`0`, generation time 5 years, census N 120, known kin file, figure labels.

### Other files in this folder

`snp_panel_dosage.csv` and `config_snp.yaml` are a second, separate example (a RADseq
SNP panel, 80 frogs, one sample). They are not part of this walkthrough.
`recovery.yaml` is the planted-truth design used in the walkthrough.

## TRUE values (simulated)

| quantity | true value | note |
|---|---|---|
| Ne | 120 | Wright-Fisher population, N = Ne = 120 every generation |
| heterozygosity lost per generation, 1/(2Ne) | 0.00417 (0.42%) | the same for both LD samples and for the temporal method |
| census N at first sampling | 120 | supplied in `config.yaml` (sampling plan I) |
| generations between samples | 2 | |
| relatedness of the 6 known full-sib pairs | 0.465 | realised pedigree relatedness of these 6 pairs, relative to the 2024 sample. The pedigree expectation is 0.5. |

The simulated population has random union of gametes (no separate sexes, no
overlapping generations, no migration, unlinked loci). Real bull trout have separate
sexes and overlapping generations; see the walkthrough's limits section.

## How the data were made (reproducible)

- Generator: `dev/popgen-ne/make_example.R`.
  It calls `wf_simulate(list(Ne = 120), design, seed = 2024)` from the module's `simulate.R`.
  The design is 16 msats, 12 alleles, Dirichlet 0.8, 50 individuals per sample, 2 generations, 6 full-sib pairs and 3% missing.
  Allele-size offsets use `set.seed(7)`.
- Check: `dev/popgen-ne/si_truth.R` regenerates both CSVs into `dev/popgen-ne/si_check/`.
  The regenerated files are byte-for-byte identical to the files here (md5 match, log `logs/popgen-ne_si_truth.log`).
  The same script writes `results/popgen-ne/si/truth.csv`.
