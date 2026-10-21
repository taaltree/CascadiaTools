# popgen-structure: genetic diversity and population structure (conservation genetics).
#
# Engine: hand-written, vectorised Weir & Cockerham (1984) / Nei (1987)
# estimators (verified to match hierfstat exactly); intervals = estimate +/-
# t(L-1) x SE with SE from a delete-one-locus jackknife plus the shared-individual
# covariance (delete-one-individual jackknife). adegenet/ade4 for PCA and DAPC,
# vegan for Mantel tests, poppr for optional AMOVA. simulate() is a Wright-Fisher
# island / stepping-stone simulator; truth is the realised census value (default)
# or the exact identity-by-descent expectation (design$truth = "expected").
# See README.md (method) and VALIDATION.md (evidence and limits).

source(file.path(MODULE_DIR, "simulate.R"), local = TRUE)
source(file.path(MODULE_DIR, "engine.R"), local = TRUE)
source(file.path(MODULE_DIR, "io.R"), local = TRUE)
source(file.path(MODULE_DIR, "output.R"), local = TRUE)

module <- ecomodule(
  id = "popgen-structure",
  title = "Genetic diversity and population structure (FST, FIS, IBD, DAPC)",
  description = paste(
    "How genetically diverse is each population, how differentiated are populations from one another,",
    "are individuals inbred, does differentiation grow with geographic distance (isolation by distance),",
    "and how many genetic clusters are there? Microsatellite or SNP genotypes from several sampled populations."),
  defaults = list(
    params = list(Ne = 100, m = 0.05, selfing = 0, mu = 5e-4, n_alleles = 12, alpha = 1,
                  snp_alpha = 1, T_gen = 200),
    design = list(model = "stepping_stone_2d", n_pops = 9, spacing_km = 25, n_ind = 25,
                  n_loci = 15, marker = "microsat", missing_rate = 0.02, ibd_distance = "log")),
  simulate = ps_simulate,
  fit = ps_fit,
  settings = list(ci = 0.95, individual_covariance = TRUE, hwe = TRUE, hwe_B = 1000, pca = TRUE, dapc = TRUE,
                  k_criterion = "diffNgroup", mantel_perm = 999, amova = FALSE, amova_perm = 999,
                  estimator = "WC", seed = 1),
  fast_settings = list(ci = 0.95, individual_covariance = TRUE, hwe = FALSE, pca = FALSE, dapc = FALSE,
                       mantel_perm = 0, amova = FALSE, estimator = "WC", seed = 1),
  read_data = ps_read_data,
  check_data = ps_check_data,
  interpret = ps_interpret,
  figures = ps_figures,
  plots = function(fit, data) ps_figures(fit, data, fig_opts()),
  methods_text = ps_methods_text,
  power_params = c("Fst_global", "Fis_global", "ibd_slope"),
  null_values = c(Fst_global = 0, Fis_global = 0, ibd_slope = 0),
  engine = "R: Weir & Cockerham / Nei estimators (cross-checked with hierfstat), adegenet (PCA, DAPC), vegan (Mantel), poppr (AMOVA)",
  references = c(
    "Weir BS, Cockerham CC (1984) Estimating F-statistics for the analysis of population structure. Evolution 38:1358-1370. doi:10.1111/j.1558-5646.1984.tb05657.x",
    "Nei M (1987) Molecular Evolutionary Genetics. Columbia University Press, New York.",
    "El Mousadik A, Petit RJ (1996) High level of genetic differentiation for allelic richness among populations of the argan tree. Theoretical and Applied Genetics 92:832-839.",
    "Hurlbert SH (1971) The nonconcept of species diversity: a critique and alternative parameters. Ecology 52:577-586.",
    "Goudet J (2005) HIERFSTAT, a package for R to compute and test hierarchical F-statistics. Molecular Ecology Notes 5:184-186.",
    "Rousset F (1996) Equilibrium values of measures of population subdivision for stepwise mutation processes. Genetics 142:1357-1362.",
    "Rousset F (1997) Genetic differentiation and estimation of gene flow from F-statistics under isolation by distance. Genetics 145:1219-1228. doi:10.1093/genetics/145.4.1219",
    "Wright S (1931) Evolution in Mendelian populations. Genetics 16:97-159.",
    "Whitlock MC, McCauley DE (1999) Indirect measures of gene flow and migration: FST != 1/(4Nm+1). Heredity 82:117-125.",
    "Mills LS, Allendorf FW (1996) The one-migrant-per-generation rule in conservation and management. Conservation Biology 10:1509-1518.",
    "Frankham R, Ballou JD, Eldridge MDB, Lacy RC, Ralls K, Dudash MR, Fenster CB (2011) Predicting the probability of outbreeding depression. Conservation Biology 25:465-475. doi:10.1111/j.1523-1739.2011.01662.x",
    "Guo SW, Thompson EA (1992) Performing the exact test of Hardy-Weinberg proportion for multiple alleles. Biometrics 48:361-372.",
    "Wigginton JE, Cutler DJ, Abecasis GR (2005) A note on exact tests of Hardy-Weinberg equilibrium. American Journal of Human Genetics 76:887-893.",
    "Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate. Journal of the Royal Statistical Society B 57:289-300.",
    "Mantel N (1967) The detection of disease clustering and a generalized regression approach. Cancer Research 27:209-220.",
    "Jombart T (2008) adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics 24:1403-1405.",
    "Jombart T, Devillard S, Balloux F (2010) Discriminant analysis of principal components: a new method for the analysis of genetically structured populations. BMC Genetics 11:94. doi:10.1186/1471-2156-11-94",
    "Thia JA (2023) Guidelines for standardizing the application of discriminant analysis of principal components to genotype data. Molecular Ecology Resources 23:523-538. doi:10.1111/1755-0998.13706",
    "Excoffier L, Smouse PE, Quattro JM (1992) Analysis of molecular variance inferred from metric distances among DNA haplotypes. Genetics 131:479-491.",
    "Kamvar ZN, Tabima JF, Grunwald NJ (2014) Poppr: an R package for genetic analysis of populations with clonal, partially clonal, and/or sexual reproduction. PeerJ 2:e281.",
    "Knaus BJ, Grunwald NJ (2017) vcfR: a package to manipulate and visualize variant call format data in R. Molecular Ecology Resources 17:44-53.")
)
