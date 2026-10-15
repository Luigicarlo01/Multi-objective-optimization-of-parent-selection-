#análise pairwise#

packages <- c("AlphaSimR", "openxlsx")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(char = packages)

strategy_codes <- c("Strategy 1", "Strategy 2", "EA", "Local", "SA", "GA", "GRASP", "TABU", "EDA", "PSO")
strategies_to_run <- strategy_codes

classic_strategies <- c("Strategy 1", "Strategy 2")
cross_plan_strategies <- setdiff(strategy_codes, classic_strategies)

target_inbreeding_values <- c(0, 0.25, 0.50, 0.75)
scenario_tolerance <- 0.025
max_strategy1_cycles_to_find <- 100

nQTL <- 550
nSNP <- 5500
gen_length <- 13 / 11
nis <- 11
n_progeny_per_fam <- 40
n_crosses <- 50
nThreads <- 1

lambda_fixed <- 0.3
top_prop_optimization <- 1
scenario_candidate_n <- 300
max_ind_per_family <- 4
max_parent_use_cross_plan <- 4
allow_duplicate_crosses <- FALSE
min_unique_parents_cross_plan <- ceiling((2 * n_crosses) / max_parent_use_cross_plan)
n_recurrent_cycles <- 5
CurrentRecurrentCycle <- NA_integer_
pasta_saida_base <- "Resultados_cenarios_endogamia_5_ciclos_candidatos300_pairwise_reduzido5x_sem_repeticoes"

Sys.setenv(
  OMP_NUM_THREADS = "4",
  OPENBLAS_NUM_THREADS = "4",
  MKL_NUM_THREADS = "4",
  VECLIB_MAXIMUM_THREADS = "4",
  NUMEXPR_NUM_THREADS = "4"
)

set.seed(3532357)

# Parametros reduzidos em aproximadamente 5x no orçamento de avaliações por ciclo.
# A estrutura biologica da analise foi mantida: 4 cenarios, 5 ciclos recorrentes,
# 300 candidatos para caracterizar os cenarios, matriz de cruzamentos pairwise e mesmas metricas.
eval_budget <- 60000
ea_pop_size <- 200
ea_n_gen <- 300
ea_elite_n <- 20
ea_tournament_k <- 3
ea_mutation_swaps <- 1
ea_random_injection_rate <- 0.05
local_n_swaps <- 5

sa_n_iter <- eval_budget
sa_temp_ini <- 1
sa_cooling <- 0.95

ga_pop_size <- 200
ga_n_gen <- 300
ga_elite_n <- 20
ga_mutation_rate <- 0.10

grasp_alpha <- 0.20

tabu_n_iter <- 12000
tabu_n_neighbors <- 5
tabu_tenure_value <- 10

eda_pop_size <- 200
eda_n_iter <- 300
eda_elite_prop <- 0.20

pso_n_particles <- 300
pso_n_iter <- 200
pso_w <- 0.72
pso_c1 <- 1.49
pso_c2 <- 1.49
pso_vmax <- 5

convergence_record_every <- 100
cross_construction_pair_candidates <- 60

founderPop <- runMacs2(
  nInd = 100,
  nChr = 11,
  segSites = (nQTL + nSNP) / 11,
  Ne = 100,
  bp = 6.9e+08,
  genLen = gen_length,
  mutRate = 4.8e-08,
  ploidy = 2L,
  returnCommand = FALSE,
  nThreads = nThreads
)

SP <- SimParam$new(founderPop)

SP$addTraitAD(
  nQtlPerChr = nQTL / 11,
  gamma = TRUE,
  mean = 0,
  var = 0.5,
  meanDD = 0.5,
  varDD = 0.5
)

SP$addSnpChip(50)
SP$setVarE(H2 = 0.5)

Egra_Parents <- newPop(founderPop)
Prog <- randCross(Egra_Parents, 100, ignoreSexes = TRUE)

for (p in 1:4) {
  Prog <- randCross(Prog, 100, ignoreSexes = TRUE)
}

base_snp_freqs <- colMeans(pullSnpGeno(Prog, simParam = SP), na.rm = TRUE) / 2

maf_threshold <- 0.05
initial_qtl_freqs <- colMeans(pullQtlGeno(Prog, simParam = SP), na.rm = TRUE) / 2

rare_alt_indices <- which(initial_qtl_freqs > 0 & initial_qtl_freqs < maf_threshold)
rare_ref_indices <- which(initial_qtl_freqs < 1 & initial_qtl_freqs > (1 - maf_threshold))
rare_loci_indices <- c(rare_alt_indices, rare_ref_indices)
total_initial_rare_alleles <- length(rare_loci_indices)

minmax_01 <- function(x) {
  x <- as.numeric(x)
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  if (!is.finite(xmin) || !is.finite(xmax) || xmax == xmin) return(rep(0.5, length(x)))
  out <- (x - xmin) / (xmax - xmin)
  out[!is.finite(out)] <- 0.5
  out
}

safe_div <- function(a, b) {
  if (!is.finite(a) || !is.finite(b) || b == 0) return(NA_real_)
  a / b
}

get_family_ids <- function(pop) {
  n <- length(pop@iid)
  dam <- as.character(pop@mother)
  sire <- as.character(pop@father)
  iid <- as.character(pop@iid)
  if (length(dam) != n) dam <- rep(NA_character_, n)
  if (length(sire) != n) sire <- rep(NA_character_, n)
  no_family <- is.na(dam) | is.na(sire) | dam == "" | sire == "" | dam == "0" | sire == "0"
  parent_1 <- pmin(dam, sire)
  parent_2 <- pmax(dam, sire)
  family_id <- paste(parent_1, parent_2, sep = "__")
  family_id[no_family] <- paste0("IndSemFamilia__", iid[no_family])
  family_id
}

violates_family_limit <- function(indices, family_ids, max_ind_per_family = 4) {
  indices <- unique(indices)
  if (is.null(family_ids) || is.infinite(max_ind_per_family) || length(indices) == 0) return(FALSE)
  fam <- family_ids[indices]
  any(table(fam) > max_ind_per_family)
}


effective_family_capacity <- function(indices, family_ids, max_ind_per_family = 4) {
  if (length(indices) == 0) return(0L)
  fam_tab <- table(family_ids[indices])
  sum(pmin(as.integer(fam_tab), max_ind_per_family))
}

sample_valid_subset <- function(candidate_pool, n_select, family_ids, max_ind_per_family = 4, max_attempts = 1000) {
  if (length(candidate_pool) < n_select) stop("candidate_pool possui menos individuos do que o necessario.")
  for (attempt in seq_len(max_attempts)) {
    sel <- sample(candidate_pool, n_select)
    if (!violates_family_limit(sel, family_ids, max_ind_per_family)) return(sel)
  }
  ordered_pool <- sample(candidate_pool, length(candidate_pool))
  sel <- integer(0)
  fam_count <- list()
  for (idx in ordered_pool) {
    fam <- family_ids[idx]
    current_count <- ifelse(is.null(fam_count[[fam]]), 0, fam_count[[fam]])
    if (current_count < max_ind_per_family) {
      sel <- c(sel, idx)
      fam_count[[fam]] <- current_count + 1
    }
    if (length(sel) == n_select) break
  }
  if (length(sel) < n_select) stop("Nao foi possivel montar um subconjunto valido com o limite de individuos por familia.")
  sel
}

repair_family_limit <- function(indices, candidate_pool, n_select, family_ids, max_ind_per_family = 4) {
  indices <- unique(indices)
  sel <- integer(0)
  fam_count <- list()
  for (idx in indices) {
    fam <- family_ids[idx]
    current_count <- ifelse(is.null(fam_count[[fam]]), 0, fam_count[[fam]])
    if (current_count < max_ind_per_family) {
      sel <- c(sel, idx)
      fam_count[[fam]] <- current_count + 1
    }
    if (length(sel) == n_select) break
  }
  if (length(sel) < n_select) {
    remaining_pool <- sample(setdiff(candidate_pool, sel))
    for (idx in remaining_pool) {
      fam <- family_ids[idx]
      current_count <- ifelse(is.null(fam_count[[fam]]), 0, fam_count[[fam]])
      if (current_count < max_ind_per_family) {
        sel <- c(sel, idx)
        fam_count[[fam]] <- current_count + 1
      }
      if (length(sel) == n_select) break
    }
  }
  if (length(sel) < n_select) stop("Nao foi possivel reparar a solucao respeitando o limite de individuos por familia.")
  sel
}

calc_genomic_relationship <- function(snp_geno, base_freqs = NULL) {
  M <- as.matrix(snp_geno)
  if (nrow(M) < 1 || ncol(M) < 2) return(diag(nrow(M)))
  if (is.null(base_freqs)) {
    p <- colMeans(M, na.rm = TRUE) / 2
  } else {
    p <- base_freqs
    if (length(p) != ncol(M)) stop("base_freqs deve ter o mesmo numero de marcadores que snp_geno.")
  }
  keep <- which(is.finite(p) & p > 0 & p < 1)
  if (length(keep) < 2) return(diag(nrow(M)))
  M_keep <- M[, keep, drop = FALSE]
  p_keep <- p[keep]
  for (j in seq_len(ncol(M_keep))) {
    miss_j <- is.na(M_keep[, j])
    if (any(miss_j)) M_keep[miss_j, j] <- 2 * p_keep[j]
  }
  Z <- sweep(M_keep, 2, 2 * p_keep, FUN = "-")
  denom <- 2 * sum(p_keep * (1 - p_keep))
  if (!is.finite(denom) || denom <= 0) return(diag(nrow(M)))
  tcrossprod(Z) / denom
}

calc_genomic_inbreeding <- function(pop, SP, base_snp_freqs) {
  snp_geno <- pullSnpGeno(pop, simParam = SP)
  G <- calc_genomic_relationship(snp_geno = snp_geno, base_freqs = base_snp_freqs)
  diag_G <- diag(G)
  individual_F <- diag_G - 1
  mean_diag_G <- mean(diag_G, na.rm = TRUE)
  genomic_F <- mean(individual_F, na.rm = TRUE)
  genomic_F_variance <- if (length(individual_F) > 1) var(individual_F, na.rm = TRUE) else NA_real_
  genomic_F_sd <- if (length(individual_F) > 1) sd(individual_F, na.rm = TRUE) else NA_real_
  genomic_F_min <- min(individual_F, na.rm = TRUE)
  genomic_F_max <- max(individual_F, na.rm = TRUE)
  
  if (!is.finite(mean_diag_G)) mean_diag_G <- NA_real_
  if (!is.finite(genomic_F)) genomic_F <- NA_real_
  if (!is.finite(genomic_F_variance)) genomic_F_variance <- NA_real_
  if (!is.finite(genomic_F_sd)) genomic_F_sd <- NA_real_
  if (!is.finite(genomic_F_min)) genomic_F_min <- NA_real_
  if (!is.finite(genomic_F_max)) genomic_F_max <- NA_real_
  
  list(
    mean_diag_G = mean_diag_G,
    genomic_F = genomic_F,
    genomic_F_variance = genomic_F_variance,
    genomic_F_sd = genomic_F_sd,
    genomic_F_min = genomic_F_min,
    genomic_F_max = genomic_F_max,
    individual_F = individual_F
  )
}


calc_observed_heterozygosity <- function(snp_geno) {
  M <- as.matrix(snp_geno)
  if (nrow(M) < 1 || ncol(M) < 1) return(NA_real_)
  ho <- mean(M == 1, na.rm = TRUE)
  if (!is.finite(ho)) ho <- NA_real_
  ho
}

calc_mean_hamming_distance <- function(snp_geno) {
  M <- as.matrix(snp_geno)
  n_ind <- nrow(M)
  n_mark <- ncol(M)
  if (n_ind < 2 || n_mark < 1) return(NA_real_)
  n_nonmiss <- colSums(!is.na(M))
  total_pairs_marker <- n_nonmiss * (n_nonmiss - 1) / 2
  n0 <- colSums(M == 0, na.rm = TRUE)
  n1 <- colSums(M == 1, na.rm = TRUE)
  n2 <- colSums(M == 2, na.rm = TRUE)
  same_pairs_marker <- n0 * (n0 - 1) / 2 + n1 * (n1 - 1) / 2 + n2 * (n2 - 1) / 2
  different_pairs_marker <- total_pairs_marker - same_pairs_marker
  valid_markers <- which(total_pairs_marker > 0)
  if (length(valid_markers) < 1) return(NA_real_)
  hamming_marker <- different_pairs_marker[valid_markers] / total_pairs_marker[valid_markers]
  hamming_mean <- mean(hamming_marker, na.rm = TRUE)
  if (!is.finite(hamming_mean)) hamming_mean <- NA_real_
  hamming_mean
}

get_top_candidate_indices <- function(pop, SP, top_prop = 0.30, min_required_candidates = 1, merit_source = "gv", base_snp_freqs = NULL) {
  if (merit_source == "gv") {
    merit <- as.numeric(gv(pop))
  } else if (merit_source == "pheno") {
    merit <- as.numeric(pheno(pop))
  } else {
    stop("merit_source deve ser 'gv' ou 'pheno'.")
  }
  merit_01 <- minmax_01(merit)
  snp_geno <- pullSnpGeno(pop, simParam = SP)
  G_for_index <- calc_genomic_relationship(snp_geno = snp_geno, base_freqs = base_snp_freqs)
  individual_F <- diag(G_for_index) - 1
  individual_F[!is.finite(individual_F)] <- 0
  selection_index <- merit_01 - individual_F
  selection_index[!is.finite(selection_index)] <- -Inf
  ordered_indices <- order(selection_index, decreasing = TRUE)
  n_top <- max(min_required_candidates, ceiling(length(merit) * top_prop))
  n_top <- min(n_top, length(merit))
  family_ids <- get_family_ids(pop)
  while (n_top < length(ordered_indices) && effective_family_capacity(ordered_indices[seq_len(n_top)], family_ids, max_ind_per_family) < min_required_candidates) {
    n_top <- n_top + 1
  }
  ordered_indices[seq_len(n_top)]
}

get_top_n_candidate_indices <- function(pop, SP, n_candidates = 300, nis, merit_source = "gv", base_snp_freqs = NULL) {
  if (merit_source == "gv") {
    merit <- as.numeric(gv(pop))
  } else if (merit_source == "pheno") {
    merit <- as.numeric(pheno(pop))
  } else {
    stop("merit_source deve ser 'gv' ou 'pheno'.")
  }
  
  merit_01 <- minmax_01(merit)
  snp_geno <- pullSnpGeno(pop, simParam = SP)
  G_for_index <- calc_genomic_relationship(snp_geno = snp_geno, base_freqs = base_snp_freqs)
  individual_F <- diag(G_for_index) - 1
  individual_F[!is.finite(individual_F)] <- 0
  
  selection_index <- merit_01 - individual_F
  selection_index[!is.finite(selection_index)] <- -Inf
  
  n_top <- min(length(merit), max(nis, n_candidates))
  order(selection_index, decreasing = TRUE)[seq_len(n_top)]
}

calc_scenario_candidate_inbreeding <- function(pop, SP, base_snp_freqs, n_candidates = 300) {
  candidate_indices <- get_top_n_candidate_indices(
    pop = pop,
    SP = SP,
    n_candidates = n_candidates,
    nis = nis,
    merit_source = "gv",
    base_snp_freqs = base_snp_freqs
  )
  
  safe_numeric <- function(x) {
    out <- tryCatch(as.numeric(x), error = function(e) NA_real_)
    out
  }
  
  safe_stat_mean <- function(x) {
    x <- as.numeric(x)
    out <- mean(x, na.rm = TRUE)
    if (!is.finite(out)) NA_real_ else out
  }
  
  safe_stat_var <- function(x) {
    x <- as.numeric(x)
    if (sum(is.finite(x)) <= 1) return(NA_real_)
    out <- var(x, na.rm = TRUE)
    if (!is.finite(out)) NA_real_ else out
  }
  
  safe_stat_sd <- function(x) {
    x <- as.numeric(x)
    if (sum(is.finite(x)) <= 1) return(NA_real_)
    out <- sd(x, na.rm = TRUE)
    if (!is.finite(out)) NA_real_ else out
  }
  
  safe_stat_min <- function(x) {
    x <- as.numeric(x)
    if (!any(is.finite(x))) return(NA_real_)
    out <- min(x, na.rm = TRUE)
    if (!is.finite(out)) NA_real_ else out
  }
  
  safe_stat_max <- function(x) {
    x <- as.numeric(x)
    if (!any(is.finite(x))) return(NA_real_)
    out <- max(x, na.rm = TRUE)
    if (!is.finite(out)) NA_real_ else out
  }
  
  candidate_pop <- pop[candidate_indices]
  candidate_inb <- calc_genomic_inbreeding(candidate_pop, SP, base_snp_freqs)
  full_inb <- calc_genomic_inbreeding(pop, SP, base_snp_freqs)
  
  candidate_gv <- safe_numeric(gv(candidate_pop))
  full_gv <- safe_numeric(gv(pop))
  
  candidate_pheno <- tryCatch(
    as.numeric(pheno(candidate_pop)),
    error = function(e) rep(NA_real_, length(candidate_indices))
  )
  
  full_pheno <- tryCatch(
    as.numeric(pheno(pop)),
    error = function(e) rep(NA_real_, length(pop@iid))
  )
  
  list(
    candidate_indices = candidate_indices,
    candidate_n = length(candidate_indices),
    
    candidate_mean_F = candidate_inb$genomic_F,
    candidate_variance_F = candidate_inb$genomic_F_variance,
    candidate_sd_F = candidate_inb$genomic_F_sd,
    candidate_min_F = candidate_inb$genomic_F_min,
    candidate_max_F = candidate_inb$genomic_F_max,
    
    candidate_mean_GV = safe_stat_mean(candidate_gv),
    candidate_variance_GV = safe_stat_var(candidate_gv),
    candidate_sd_GV = safe_stat_sd(candidate_gv),
    candidate_min_GV = safe_stat_min(candidate_gv),
    candidate_max_GV = safe_stat_max(candidate_gv),
    
    candidate_mean_pheno = safe_stat_mean(candidate_pheno),
    candidate_variance_pheno = safe_stat_var(candidate_pheno),
    candidate_sd_pheno = safe_stat_sd(candidate_pheno),
    candidate_min_pheno = safe_stat_min(candidate_pheno),
    candidate_max_pheno = safe_stat_max(candidate_pheno),
    
    full_mean_F = full_inb$genomic_F,
    full_variance_F = full_inb$genomic_F_variance,
    full_sd_F = full_inb$genomic_F_sd,
    full_min_F = full_inb$genomic_F_min,
    full_max_F = full_inb$genomic_F_max,
    
    full_mean_GV = safe_stat_mean(full_gv),
    full_variance_GV = safe_stat_var(full_gv),
    full_sd_GV = safe_stat_sd(full_gv),
    full_min_GV = safe_stat_min(full_gv),
    full_max_GV = safe_stat_max(full_gv),
    
    full_mean_pheno = safe_stat_mean(full_pheno),
    full_variance_pheno = safe_stat_var(full_pheno),
    full_sd_pheno = safe_stat_sd(full_pheno),
    full_min_pheno = safe_stat_min(full_pheno),
    full_max_pheno = safe_stat_max(full_pheno)
  )
}


build_selection_components <- function(pop, SP, merit_source = "gv", candidate_indices = NULL, base_snp_freqs = NULL) {
  snp_geno <- pullSnpGeno(pop, simParam = SP)
  K <- calc_genomic_relationship(snp_geno = snp_geno, base_freqs = base_snp_freqs)
  if (is.null(candidate_indices)) candidate_indices <- seq_len(nrow(K))
  K_candidate <- K[candidate_indices, candidate_indices, drop = FALSE]
  off_diag_vals <- K_candidate[lower.tri(K_candidate)]
  kin_mean <- mean(off_diag_vals, na.rm = TRUE)
  kin_sd <- sd(off_diag_vals, na.rm = TRUE)
  if (!is.finite(kin_mean)) kin_mean <- 0
  if (!is.finite(kin_sd) || kin_sd <= 0) kin_sd <- 1
  if (merit_source == "gv") {
    merit <- as.numeric(gv(pop))
  } else if (merit_source == "pheno") {
    merit <- as.numeric(pheno(pop))
  } else {
    stop("merit_source deve ser 'gv' ou 'pheno'.")
  }
  merit_mean <- mean(merit[candidate_indices], na.rm = TRUE)
  merit_sd <- sd(merit[candidate_indices], na.rm = TRUE)
  if (!is.finite(merit_mean)) merit_mean <- 0
  if (!is.finite(merit_sd) || merit_sd <= 0) merit_sd <- 1
  merit_z <- (merit - merit_mean) / merit_sd
  merit_z[is.na(merit_z)] <- 0
  list(K = K, merit = merit, merit_z = merit_z, kin_mean = kin_mean, kin_sd = kin_sd, candidate_indices = candidate_indices)
}

calc_fitness <- function(indices, merit_z, K, lambda, kin_mean, kin_sd, family_ids = NULL, max_ind_per_family = Inf) {
  if (violates_family_limit(indices, family_ids, max_ind_per_family)) return(-Inf)
  mean_merit_z <- mean(merit_z[indices], na.rm = TRUE)
  if (length(indices) < 2) return(mean_merit_z)
  subK <- K[indices, indices, drop = FALSE]
  kin_m <- mean(subK[lower.tri(subK)], na.rm = TRUE)
  kin_z <- (kin_m - kin_mean) / kin_sd
  if (!is.finite(kin_z)) kin_z <- 0
  mean_merit_z - lambda * kin_z
}

calc_solution_stats <- function(indices, merit_z, K, lambda, kin_mean, kin_sd, family_ids = NULL, max_ind_per_family = Inf) {
  objective <- calc_fitness(indices, merit_z, K, lambda, kin_mean, kin_sd, family_ids, max_ind_per_family)
  mean_merit_z <- mean(merit_z[indices], na.rm = TRUE)
  subK <- K[indices, indices, drop = FALSE]
  mean_kinship <- if (length(indices) < 2) NA_real_ else mean(subK[lower.tri(subK)], na.rm = TRUE)
  kin_z <- (mean_kinship - kin_mean) / kin_sd
  if (!is.finite(mean_merit_z)) mean_merit_z <- NA_real_
  if (!is.finite(mean_kinship)) mean_kinship <- NA_real_
  if (!is.finite(kin_z)) kin_z <- NA_real_
  if (!is.finite(objective)) objective <- NA_real_
  list(objective = objective, mean_merit_z = mean_merit_z, mean_kinship = mean_kinship, kin_z = kin_z)
}

new_convergence_row <- function(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy, iteration, evaluation, current_objective, best_objective) {
  data.frame(
    Replication = 1L,
    TargetInbreeding = target_inbreeding,
    ObservedInitialInbreeding = observed_initial_inbreeding,
    ScenarioCycleFound = scenario_cycle_found,
    Generation = CurrentRecurrentCycle,
    Strategy = strategy,
    Iteration = iteration,
    Evaluation = evaluation,
    CurrentObjective = current_objective,
    BestObjective = best_objective,
    stringsAsFactors = FALSE
  )
}


record_generation_metrics_one_cycle <- function(new_pop, initial_pop, SP, qtl_effects, total_qtl, base_snp_freqs) {
  snp_geno <- pullSnpGeno(new_pop, simParam = SP)
  qtl_freqs <- colMeans(pullQtlGeno(new_pop, simParam = SP), na.rm = TRUE) / 2
  fav_f <- ifelse(qtl_effects > 0, qtl_freqs, 1 - qtl_freqs)
  gen_inb <- calc_genomic_inbreeding(pop = new_pop, SP = SP, base_snp_freqs = base_snp_freqs)
  initial_inb <- calc_genomic_inbreeding(pop = initial_pop, SP = SP, base_snp_freqs = base_snp_freqs)
  F_current <- gen_inb$genomic_F
  F_previous <- initial_inb$genomic_F
  delta_F <- (F_current - F_previous) / (1 - F_previous)
  ne_val <- ifelse(is.finite(delta_F) && delta_F > 0, 1 / (2 * delta_F), NA_real_)
  if (total_initial_rare_alleles > 0) {
    lost_alt <- if (length(rare_alt_indices) > 0) sum(qtl_freqs[rare_alt_indices] <= 0) else 0
    lost_ref <- if (length(rare_ref_indices) > 0) sum(qtl_freqs[rare_ref_indices] >= 1) else 0
    rare_allele_loss_rate <- (lost_alt + lost_ref) / total_initial_rare_alleles
  } else {
    rare_allele_loss_rate <- NA_real_
  }
  list(
    MeanPhenotypicValue = mean(new_pop@gv),
    GeneticGain = mean(new_pop@gv) - mean(initial_pop@gv),
    AdditiveVariance = varA(new_pop),
    PhenotypicVariance = varP(new_pop),
    DominanceVariance = varD(new_pop),
    TotalGeneticVariance = varG(new_pop),
    MeanDiagonalG = gen_inb$mean_diag_G,
    GenomicEndogamyF = gen_inb$genomic_F,
    EffectiveCrossings = ne_val,
    MeanHammingDistance = calc_mean_hamming_distance(snp_geno),
    MeanFavorableAlleleFreq = mean(fav_f),
    MeanUnfavorableAlleleFreq = 1 - mean(fav_f),
    FixedFavorableRate = sum((qtl_effects > 0 & qtl_freqs == 1) | (qtl_effects < 0 & qtl_freqs == 0)) / total_qtl,
    FixedUnfavorableRate = sum((qtl_effects < 0 & qtl_freqs == 1) | (qtl_effects > 0 & qtl_freqs == 0)) / total_qtl,
    RareAlleleLossRate = rare_allele_loss_rate,
    ObservedHeterozygosity = calc_observed_heterozygosity(snp_geno)
  )
}

record_generation_metrics_recurrent <- function(new_pop, previous_pop, baseline_pop, SP, qtl_effects, total_qtl, base_snp_freqs) {
  snp_geno <- pullSnpGeno(new_pop, simParam = SP)
  
  qtl_freqs <- colMeans(pullQtlGeno(new_pop, simParam = SP), na.rm = TRUE) / 2
  fav_f <- ifelse(qtl_effects > 0, qtl_freqs, 1 - qtl_freqs)
  
  gen_inb <- calc_genomic_inbreeding(pop = new_pop, SP = SP, base_snp_freqs = base_snp_freqs)
  prev_inb <- calc_genomic_inbreeding(pop = previous_pop, SP = SP, base_snp_freqs = base_snp_freqs)
  base_inb <- calc_genomic_inbreeding(pop = baseline_pop, SP = SP, base_snp_freqs = base_snp_freqs)
  
  F_current <- gen_inb$genomic_F
  F_previous <- prev_inb$genomic_F
  F_baseline <- base_inb$genomic_F
  
  delta_F_per_cycle <- (F_current - F_previous) / (1 - F_previous)
  delta_F_cumulative <- (F_current - F_baseline) / (1 - F_baseline)
  
  ne_per_cycle <- ifelse(
    is.finite(delta_F_per_cycle) && delta_F_per_cycle > 0,
    1 / (2 * delta_F_per_cycle),
    NA_real_
  )
  
  ne_cumulative <- ifelse(
    is.finite(delta_F_cumulative) && delta_F_cumulative > 0,
    1 / (2 * delta_F_cumulative),
    NA_real_
  )
  
  if (total_initial_rare_alleles > 0) {
    lost_alt <- if (length(rare_alt_indices) > 0) sum(qtl_freqs[rare_alt_indices] <= 0) else 0
    lost_ref <- if (length(rare_ref_indices) > 0) sum(qtl_freqs[rare_ref_indices] >= 1) else 0
    rare_allele_loss_rate <- (lost_alt + lost_ref) / total_initial_rare_alleles
  } else {
    rare_allele_loss_rate <- NA_real_
  }
  
  list(
    MeanPhenotypicValue = mean(new_pop@gv),
    GeneticGain = mean(new_pop@gv) - mean(baseline_pop@gv),
    GeneticGainPerCycle = mean(new_pop@gv) - mean(previous_pop@gv),
    AdditiveVariance = varA(new_pop),
    PhenotypicVariance = varP(new_pop),
    DominanceVariance = varD(new_pop),
    TotalGeneticVariance = varG(new_pop),
    MeanDiagonalG = gen_inb$mean_diag_G,
    GenomicEndogamyF = gen_inb$genomic_F,
    GenomicEndogamyVariance = gen_inb$genomic_F_variance,
    GenomicEndogamySD = gen_inb$genomic_F_sd,
    GenomicEndogamyMin = gen_inb$genomic_F_min,
    GenomicEndogamyMax = gen_inb$genomic_F_max,
    DeltaFPerCycle = delta_F_per_cycle,
    DeltaFCumulative = delta_F_cumulative,
    EffectiveCrossings = ne_per_cycle,
    EffectiveCrossingsCumulative = ne_cumulative,
    MeanHammingDistance = calc_mean_hamming_distance(snp_geno),
    MeanFavorableAlleleFreq = mean(fav_f),
    MeanUnfavorableAlleleFreq = 1 - mean(fav_f),
    FixedFavorableRate = sum((qtl_effects > 0 & qtl_freqs == 1) | (qtl_effects < 0 & qtl_freqs == 0)) / total_qtl,
    FixedUnfavorableRate = sum((qtl_effects < 0 & qtl_freqs == 1) | (qtl_effects > 0 & qtl_freqs == 0)) / total_qtl,
    RareAlleleLossRate = rare_allele_loss_rate,
    ObservedHeterozygosity = calc_observed_heterozygosity(snp_geno)
  )
}


prepare_selection_env <- function(pop) {
  candidate_indices <- get_top_candidate_indices(
    pop = pop,
    SP = SP,
    top_prop = top_prop_optimization,
    min_required_candidates = max(nis, min_unique_parents_cross_plan),
    merit_source = "gv",
    base_snp_freqs = base_snp_freqs
  )
  comps <- build_selection_components(
    pop = pop,
    SP = SP,
    merit_source = "gv",
    candidate_indices = candidate_indices,
    base_snp_freqs = base_snp_freqs
  )
  list(
    K = comps$K,
    merit_z = comps$merit_z,
    kin_mean = comps$kin_mean,
    kin_sd = comps$kin_sd,
    candidate_pool = comps$candidate_indices,
    family_ids = get_family_ids(pop)
  )
}

normalize_cross_plan <- function(cross_plan) {
  if (is.null(cross_plan) || length(cross_plan) == 0) return(matrix(integer(0), ncol = 2))
  out <- as.matrix(cross_plan)
  storage.mode(out) <- "integer"
  if (ncol(out) != 2) stop("cross_plan deve ter duas colunas.")
  out
}

canonical_pair_key <- function(cross_plan) {
  cross_plan <- normalize_cross_plan(cross_plan)
  if (nrow(cross_plan) == 0) return(character(0))
  paste(pmin(cross_plan[, 1], cross_plan[, 2]), pmax(cross_plan[, 1], cross_plan[, 2]), sep = "__")
}

parent_use_table <- function(cross_plan) {
  cross_plan <- normalize_cross_plan(cross_plan)
  if (nrow(cross_plan) == 0) return(integer(0))
  table(as.vector(cross_plan))
}

violates_parent_use_limit <- function(cross_plan, max_parent_use = 4) {
  tab <- parent_use_table(cross_plan)
  if (length(tab) == 0) return(FALSE)
  any(tab > max_parent_use)
}

violates_selfing_cross_plan <- function(cross_plan) {
  cross_plan <- normalize_cross_plan(cross_plan)
  if (nrow(cross_plan) == 0) return(FALSE)
  any(cross_plan[, 1] == cross_plan[, 2])
}

violates_duplicate_crosses <- function(cross_plan, allow_duplicates = FALSE) {
  if (allow_duplicates) return(FALSE)
  keys <- canonical_pair_key(cross_plan)
  any(duplicated(keys))
}

validate_cross_plan <- function(cross_plan, family_ids, max_ind_per_family = 4, max_parent_use = 4, allow_duplicates = FALSE, required_n_crosses = n_crosses) {
  cross_plan <- normalize_cross_plan(cross_plan)
  if (nrow(cross_plan) != required_n_crosses) return(FALSE)
  if (ncol(cross_plan) != 2) return(FALSE)
  if (any(!is.finite(cross_plan))) return(FALSE)
  if (violates_selfing_cross_plan(cross_plan)) return(FALSE)
  if (violates_parent_use_limit(cross_plan, max_parent_use)) return(FALSE)
  if (violates_duplicate_crosses(cross_plan, allow_duplicates)) return(FALSE)
  unique_parents <- unique(as.vector(cross_plan))
  if (violates_family_limit(unique_parents, family_ids, max_ind_per_family)) return(FALSE)
  TRUE
}

get_cross_plan_unique_parents <- function(cross_plan) {
  cross_plan <- normalize_cross_plan(cross_plan)
  unique(as.vector(cross_plan))
}

can_use_parent_in_partial_plan <- function(parent, partial_plan, family_ids, max_ind_per_family = 4, max_parent_use = 4) {
  partial_plan <- normalize_cross_plan(partial_plan)
  current_uses <- if (nrow(partial_plan) == 0) 0 else sum(as.vector(partial_plan) == parent)
  if (current_uses >= max_parent_use) return(FALSE)
  current_unique <- if (nrow(partial_plan) == 0) integer(0) else unique(as.vector(partial_plan))
  if (!(parent %in% current_unique)) {
    fam_current <- family_ids[current_unique]
    fam_parent <- family_ids[parent]
    if (sum(fam_current == fam_parent) >= max_ind_per_family) return(FALSE)
  }
  TRUE
}

can_add_pair_to_partial_plan <- function(pair, partial_plan, family_ids, max_ind_per_family = 4, max_parent_use = 4, allow_duplicates = FALSE) {
  if (length(pair) != 2) return(FALSE)
  if (pair[1] == pair[2]) return(FALSE)
  partial_plan <- normalize_cross_plan(partial_plan)
  candidate <- rbind(partial_plan, as.integer(pair))
  if (violates_parent_use_limit(candidate, max_parent_use)) return(FALSE)
  if (violates_duplicate_crosses(candidate, allow_duplicates)) return(FALSE)
  unique_parents <- unique(as.vector(candidate))
  if (violates_family_limit(unique_parents, family_ids, max_ind_per_family)) return(FALSE)
  TRUE
}

sample_valid_pair <- function(partial_plan, candidate_pool, family_ids, max_ind_per_family = 4, max_parent_use = 4, allow_duplicates = FALSE, weights = NULL, max_attempts = 5000) {
  candidate_pool <- as.integer(candidate_pool)
  if (is.null(weights)) {
    weights <- rep(1, length(candidate_pool))
  } else {
    weights <- as.numeric(weights)
    weights[!is.finite(weights) | weights < 0] <- 0
    if (sum(weights) <= 0) weights <- rep(1, length(candidate_pool))
  }
  for (attempt in seq_len(max_attempts)) {
    p1 <- sample(candidate_pool, 1, prob = weights)
    p2 <- sample(candidate_pool, 1, prob = weights)
    if (can_add_pair_to_partial_plan(c(p1, p2), partial_plan, family_ids, max_ind_per_family, max_parent_use, allow_duplicates)) {
      return(as.integer(c(p1, p2)))
    }
  }
  ordered_pool <- sample(candidate_pool, length(candidate_pool), prob = weights)
  for (p1 in ordered_pool) {
    for (p2 in sample(setdiff(candidate_pool, p1), length(setdiff(candidate_pool, p1)))) {
      if (can_add_pair_to_partial_plan(c(p1, p2), partial_plan, family_ids, max_ind_per_family, max_parent_use, allow_duplicates)) {
        return(as.integer(c(p1, p2)))
      }
    }
  }
  stop("Nao foi possivel construir um par valido para o plano de cruzamentos.")
}

sample_valid_cross_plan <- function(candidate_pool, family_ids, n_crosses, max_ind_per_family = 4, max_parent_use = 4, allow_duplicates = FALSE, weights = NULL, max_plan_attempts = 100) {
  candidate_pool <- as.integer(candidate_pool)
  min_needed <- ceiling((2 * n_crosses) / max_parent_use)
  if (length(candidate_pool) < min_needed) {
    stop("candidate_pool possui menos individuos que o minimo necessario para o limite de uso parental.")
  }
  if (effective_family_capacity(candidate_pool, family_ids, max_ind_per_family) < min_needed) {
    stop("candidate_pool nao possui diversidade de familias suficiente para o limite de individuos por familia.")
  }
  for (attempt_plan in seq_len(max_plan_attempts)) {
    plan <- matrix(integer(0), ncol = 2)
    ok <- TRUE
    for (i in seq_len(n_crosses)) {
      pair <- tryCatch(
        sample_valid_pair(plan, candidate_pool, family_ids, max_ind_per_family, max_parent_use, allow_duplicates, weights),
        error = function(e) NULL
      )
      if (is.null(pair)) {
        ok <- FALSE
        break
      }
      plan <- rbind(plan, pair)
    }
    if (!ok) next
    storage.mode(plan) <- "integer"
    if (validate_cross_plan(plan, family_ids, max_ind_per_family, max_parent_use, allow_duplicates, n_crosses)) return(plan)
  }
  stop("Nao foi possivel gerar um plano de cruzamentos valido apos varias tentativas.")
}

calc_cross_plan_fitness <- function(cross_plan, merit_z, K, lambda, kin_mean, kin_sd, family_ids = NULL, max_ind_per_family = 4, max_parent_use = 4, allow_duplicates = FALSE, required_n_crosses = n_crosses) {
  cross_plan <- normalize_cross_plan(cross_plan)
  if (nrow(cross_plan) != required_n_crosses) return(-Inf)
  if (!validate_cross_plan(cross_plan, family_ids, max_ind_per_family, max_parent_use, allow_duplicates, required_n_crosses)) return(-Inf)
  parent_uses <- as.vector(cross_plan)
  mean_merit_z <- mean(merit_z[parent_uses], na.rm = TRUE)
  pair_kinships <- K[cbind(cross_plan[, 1], cross_plan[, 2])]
  mean_pair_kinship <- mean(pair_kinships, na.rm = TRUE)
  kin_z <- (mean_pair_kinship - kin_mean) / kin_sd
  if (!is.finite(mean_merit_z)) mean_merit_z <- -Inf
  if (!is.finite(kin_z)) kin_z <- 0
  mean_merit_z - lambda * kin_z
}

calc_cross_plan_solution_stats <- function(cross_plan, merit_z, K, lambda, kin_mean, kin_sd, family_ids = NULL, max_ind_per_family = 4, max_parent_use = 4, allow_duplicates = FALSE, required_n_crosses = n_crosses) {
  objective <- calc_cross_plan_fitness(cross_plan, merit_z, K, lambda, kin_mean, kin_sd, family_ids, max_ind_per_family, max_parent_use, allow_duplicates, required_n_crosses)
  cross_plan <- normalize_cross_plan(cross_plan)
  parent_uses <- as.vector(cross_plan)
  unique_parents <- unique(parent_uses)
  use_counts <- as.integer(parent_use_table(cross_plan))
  mean_merit_z <- mean(merit_z[parent_uses], na.rm = TRUE)
  pair_kinships <- K[cbind(cross_plan[, 1], cross_plan[, 2])]
  mean_pair_kinship <- mean(pair_kinships, na.rm = TRUE)
  kin_z <- (mean_pair_kinship - kin_mean) / kin_sd
  unique_subK <- K[unique_parents, unique_parents, drop = FALSE]
  mean_unique_kinship <- if (length(unique_parents) < 2) NA_real_ else mean(unique_subK[lower.tri(unique_subK)], na.rm = TRUE)
  if (!is.finite(mean_merit_z)) mean_merit_z <- NA_real_
  if (!is.finite(mean_pair_kinship)) mean_pair_kinship <- NA_real_
  if (!is.finite(mean_unique_kinship)) mean_unique_kinship <- NA_real_
  if (!is.finite(kin_z)) kin_z <- NA_real_
  if (!is.finite(objective)) objective <- NA_real_
  list(
    objective = objective,
    mean_merit_z = mean_merit_z,
    mean_kinship = mean_pair_kinship,
    kin_z = kin_z,
    mean_unique_kinship = mean_unique_kinship,
    num_unique_parents = length(unique_parents),
    max_parent_use_observed = ifelse(length(use_counts) > 0, max(use_counts), NA_integer_)
  )
}

add_valid_rows_from_source <- function(start_plan, source_rows, env) {
  plan <- normalize_cross_plan(start_plan)
  source_rows <- normalize_cross_plan(source_rows)
  if (nrow(source_rows) == 0) return(plan)
  ord <- sample(seq_len(nrow(source_rows)))
  for (i in ord) {
    if (nrow(plan) >= n_crosses) break
    pair <- source_rows[i, ]
    if (can_add_pair_to_partial_plan(pair, plan, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses)) {
      plan <- rbind(plan, pair)
    }
  }
  plan
}

complete_cross_plan <- function(partial_plan, env, weights = NULL) {
  plan <- normalize_cross_plan(partial_plan)
  while (nrow(plan) < n_crosses) {
    pair <- sample_valid_pair(
      partial_plan = plan,
      candidate_pool = env$candidate_pool,
      family_ids = env$family_ids,
      max_ind_per_family = max_ind_per_family,
      max_parent_use = max_parent_use_cross_plan,
      allow_duplicates = allow_duplicate_crosses,
      weights = weights
    )
    plan <- rbind(plan, pair)
  }
  storage.mode(plan) <- "integer"
  if (!validate_cross_plan(plan, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, n_crosses)) {
    plan <- sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, weights)
  }
  plan
}

mutate_cross_plan <- function(plan, env, n_changes = 1, weights = NULL) {
  plan <- normalize_cross_plan(plan)
  if (nrow(plan) == 0) return(sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, weights))
  for (change in seq_len(n_changes)) {
    remove_row <- sample(seq_len(nrow(plan)), 1)
    partial <- plan[-remove_row, , drop = FALSE]
    new_pair <- sample_valid_pair(
      partial_plan = partial,
      candidate_pool = env$candidate_pool,
      family_ids = env$family_ids,
      max_ind_per_family = max_ind_per_family,
      max_parent_use = max_parent_use_cross_plan,
      allow_duplicates = allow_duplicate_crosses,
      weights = weights
    )
    plan <- rbind(partial, new_pair)
  }
  plan <- plan[sample(seq_len(nrow(plan))), , drop = FALSE]
  storage.mode(plan) <- "integer"
  if (!validate_cross_plan(plan, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, n_crosses)) {
    plan <- sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, weights)
  }
  plan
}

crossover_cross_plans <- function(parent1, parent2, env) {
  p1 <- normalize_cross_plan(parent1)
  p2 <- normalize_cross_plan(parent2)
  rows <- rbind(p1[sample(seq_len(nrow(p1))), , drop = FALSE], p2[sample(seq_len(nrow(p2))), , drop = FALSE])
  child <- matrix(integer(0), ncol = 2)
  child <- add_valid_rows_from_source(child, rows, env)
  complete_cross_plan(child, env)
}

construct_greedy_cross_plan <- function(env, pair_candidates_per_step = cross_construction_pair_candidates) {
  plan <- matrix(integer(0), ncol = 2)
  while (nrow(plan) < n_crosses) {
    candidates <- vector("list", pair_candidates_per_step)
    fits <- rep(-Inf, pair_candidates_per_step)
    for (i in seq_len(pair_candidates_per_step)) {
      pair <- tryCatch(
        sample_valid_pair(plan, env$candidate_pool, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses),
        error = function(e) NULL
      )
      if (is.null(pair)) next
      tmp <- rbind(plan, pair)
      if (nrow(tmp) < n_crosses) {
        parent_uses <- as.vector(tmp)
        mean_merit_z <- mean(env$merit_z[parent_uses], na.rm = TRUE)
        pair_kinships <- env$K[cbind(tmp[, 1], tmp[, 2])]
        mean_pair_kinship <- mean(pair_kinships, na.rm = TRUE)
        kin_z <- (mean_pair_kinship - env$kin_mean) / env$kin_sd
        if (!is.finite(kin_z)) kin_z <- 0
        fits[i] <- mean_merit_z - lambda_fixed * kin_z
      } else {
        fits[i] <- calc_cross_plan_fitness(tmp, env$merit_z, env$K, lambda_fixed, env$kin_mean, env$kin_sd, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, n_crosses)
      }
      candidates[[i]] <- pair
    }
    if (!any(is.finite(fits))) {
      return(sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses))
    }
    chosen <- candidates[[which.max(fits)]]
    plan <- rbind(plan, chosen)
  }
  storage.mode(plan) <- "integer"
  if (!validate_cross_plan(plan, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, n_crosses)) {
    plan <- sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses)
  }
  plan
}

make_cross_from_plan <- function(pop, cross_plan, n_progeny, SP) {
  cross_plan <- normalize_cross_plan(cross_plan)
  out <- tryCatch(
    makeCross(pop = pop, crossPlan = cross_plan, nProgeny = n_progeny, simParam = SP),
    error = function(e1) {
      tryCatch(
        makeCross2(females = pop, males = pop, crossPlan = cross_plan, nProgeny = n_progeny, simParam = SP),
        error = function(e2) {
          stop(paste0("Falha ao executar makeCross/makeCross2: ", conditionMessage(e1), " | ", conditionMessage(e2)))
        }
      )
    }
  )
  out
}

select_ea_cross_plan <- function(pop, target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name) {
  t0 <- proc.time()[3]
  env <- prepare_selection_env(pop)
  evals <- 0
  convergence <- list()
  eval_fit <- function(x) {
    evals <<- evals + 1
    calc_cross_plan_fitness(x, env$merit_z, env$K, lambda_fixed, env$kin_mean, env$kin_sd, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, n_crosses)
  }
  tournament_parent <- function(pop_list, fits, k) {
    k <- min(k, length(pop_list))
    idx <- sample(seq_along(pop_list), k)
    pop_list[[idx[which.max(fits[idx])]]]
  }
  ea_pop <- vector("list", ea_pop_size)
  ea_pop[[1]] <- construct_greedy_cross_plan(env)
  if (ea_pop_size > 1) {
    for (i in 2:ea_pop_size) {
      ea_pop[[i]] <- sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses)
    }
  }
  best_global <- NULL
  best_global_fit <- -Inf
  child_n <- ea_pop_size - ea_elite_n
  for (gen_ea in seq_len(ea_n_gen)) {
    fits_ea <- sapply(ea_pop, eval_fit)
    if (max(fits_ea, na.rm = TRUE) > best_global_fit) {
      best_global_fit <- max(fits_ea, na.rm = TRUE)
      best_global <- ea_pop[[which.max(fits_ea)]]
    }
    elite_order <- order(fits_ea, decreasing = TRUE)[seq_len(ea_elite_n)]
    elites <- ea_pop[elite_order]
    children <- vector("list", child_n)
    for (c_idx in seq_len(child_n)) {
      if (runif(1) < ea_random_injection_rate) {
        children[[c_idx]] <- sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses)
      } else {
        parent <- tournament_parent(ea_pop, fits_ea, ea_tournament_k)
        children[[c_idx]] <- mutate_cross_plan(parent, env, ea_mutation_swaps)
      }
    }
    ea_pop <- c(elites, children)
    if (gen_ea == 1 || gen_ea == ea_n_gen || gen_ea %% convergence_record_every == 0) {
      convergence[[length(convergence) + 1]] <- new_convergence_row(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name, gen_ea, evals, max(fits_ea, na.rm = TRUE), best_global_fit)
    }
  }
  fits_final <- sapply(ea_pop, eval_fit)
  if (max(fits_final, na.rm = TRUE) > best_global_fit) {
    best_global_fit <- max(fits_final, na.rm = TRUE)
    best_global <- ea_pop[[which.max(fits_final)]]
  }
  stats <- calc_cross_plan_solution_stats(best_global, env$merit_z, env$K, lambda_fixed, env$kin_mean, env$kin_sd, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, n_crosses)
  list(indices = get_cross_plan_unique_parents(best_global), cross_plan = best_global, stats = stats, elapsed = proc.time()[3] - t0, evals = evals, convergence = do.call(rbind, convergence))
}

select_metaheuristic_cross_plan <- function(pop, method, target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name) {
  t0 <- proc.time()[3]
  env <- prepare_selection_env(pop)
  evals <- 0
  convergence <- list()
  eval_fit <- function(x) {
    evals <<- evals + 1
    calc_cross_plan_fitness(x, env$merit_z, env$K, lambda_fixed, env$kin_mean, env$kin_sd, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, n_crosses)
  }
  if (method == "Local") {
    n_restarts_local <- max(1, ceiling(eval_budget / (n_crosses * max(1, local_n_swaps) * 10)))
    best_global <- NULL
    best_global_fit <- -Inf
    for (restart in seq_len(n_restarts_local)) {
      if (restart == 1) {
        sol <- construct_greedy_cross_plan(env)
      } else {
        sol <- sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses)
      }
      improved <- TRUE
      while (improved && evals < eval_budget) {
        improved <- FALSE
        fit_old <- eval_fit(sol)
        for (j in seq_len(local_n_swaps)) {
          new_sol <- mutate_cross_plan(sol, env, 1)
          fit_new <- eval_fit(new_sol)
          if (fit_new > fit_old) {
            sol <- new_sol
            fit_old <- fit_new
            improved <- TRUE
          }
        }
      }
      fit_sol <- eval_fit(sol)
      if (fit_sol > best_global_fit) {
        best_global_fit <- fit_sol
        best_global <- sol
      }
      if (restart == 1 || restart == n_restarts_local || restart %% convergence_record_every == 0) {
        convergence[[length(convergence) + 1]] <- new_convergence_row(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name, restart, evals, fit_sol, best_global_fit)
      }
      if (evals >= eval_budget) break
    }
    selected_cross_plan <- best_global
  } else if (method == "SA") {
    current <- sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses)
    best <- current
    Temp <- sa_temp_ini
    fit_current <- eval_fit(current)
    best_fit <- fit_current
    convergence[[length(convergence) + 1]] <- new_convergence_row(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name, 0, evals, fit_current, best_fit)
    for (iter in seq_len(sa_n_iter)) {
      new_cand <- mutate_cross_plan(current, env, 1)
      fit_new <- eval_fit(new_cand)
      diff <- fit_new - fit_current
      if (diff > 0 || runif(1) < exp(diff / (Temp + 1e-10))) {
        current <- new_cand
        fit_current <- fit_new
      }
      if (fit_current > best_fit) {
        best <- current
        best_fit <- fit_current
      }
      Temp <- Temp * sa_cooling
      if (iter == 1 || iter == sa_n_iter || iter %% convergence_record_every == 0) {
        convergence[[length(convergence) + 1]] <- new_convergence_row(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name, iter, evals, fit_current, best_fit)
      }
    }
    selected_cross_plan <- best
  } else if (method == "GA") {
    ga_pop <- replicate(
      ga_pop_size,
      sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses),
      simplify = FALSE
    )
    best_global <- NULL
    best_global_fit <- -Inf
    child_n <- ga_pop_size - ga_elite_n
    for (gen_ga in seq_len(ga_n_gen)) {
      fits_ga <- sapply(ga_pop, eval_fit)
      if (max(fits_ga, na.rm = TRUE) > best_global_fit) {
        best_global_fit <- max(fits_ga, na.rm = TRUE)
        best_global <- ga_pop[[which.max(fits_ga)]]
      }
      parents_ga <- ga_pop[order(fits_ga, decreasing = TRUE)[seq_len(ga_elite_n)]]
      children <- vector("list", child_n)
      for (c_idx in seq_len(child_n)) {
        p1 <- parents_ga[[sample(seq_len(ga_elite_n), 1)]]
        p2 <- parents_ga[[sample(seq_len(ga_elite_n), 1)]]
        child <- crossover_cross_plans(p1, p2, env)
        if (runif(1) < ga_mutation_rate) child <- mutate_cross_plan(child, env, 1)
        children[[c_idx]] <- child
      }
      ga_pop <- c(parents_ga, children)
      if (gen_ga == 1 || gen_ga == ga_n_gen || gen_ga %% convergence_record_every == 0) {
        convergence[[length(convergence) + 1]] <- new_convergence_row(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name, gen_ga, evals, max(fits_ga, na.rm = TRUE), best_global_fit)
      }
    }
    fits_final <- sapply(ga_pop, eval_fit)
    if (max(fits_final, na.rm = TRUE) > best_global_fit) {
      best_global_fit <- max(fits_final, na.rm = TRUE)
      best_global <- ga_pop[[which.max(fits_final)]]
    }
    selected_cross_plan <- best_global
  } else if (method == "GRASP") {
    grasp_n_iter <- max(1, ceiling(eval_budget / (n_crosses * cross_construction_pair_candidates)))
    best_global <- NULL
    best_global_fit <- -Inf
    for (iter in seq_len(grasp_n_iter)) {
      sol <- matrix(integer(0), ncol = 2)
      while (nrow(sol) < n_crosses) {
        sampled_pairs <- vector("list", cross_construction_pair_candidates)
        fits <- rep(-Inf, cross_construction_pair_candidates)
        for (i in seq_len(cross_construction_pair_candidates)) {
          pair <- tryCatch(sample_valid_pair(sol, env$candidate_pool, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses), error = function(e) NULL)
          if (is.null(pair)) next
          tmp <- rbind(sol, pair)
          if (nrow(tmp) == n_crosses) {
            fits[i] <- eval_fit(tmp)
          } else {
            parent_uses <- as.vector(tmp)
            mean_merit_z <- mean(env$merit_z[parent_uses], na.rm = TRUE)
            pair_kinships <- env$K[cbind(tmp[, 1], tmp[, 2])]
            mean_pair_kinship <- mean(pair_kinships, na.rm = TRUE)
            kin_z <- (mean_pair_kinship - env$kin_mean) / env$kin_sd
            if (!is.finite(kin_z)) kin_z <- 0
            fits[i] <- mean_merit_z - lambda_fixed * kin_z
          }
          sampled_pairs[[i]] <- pair
        }
        if (!any(is.finite(fits))) {
          sol <- sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses)
          break
        }
        threshold <- max(fits, na.rm = TRUE) - grasp_alpha * (max(fits, na.rm = TRUE) - min(fits[is.finite(fits)], na.rm = TRUE))
        rcl <- which(fits >= threshold & is.finite(fits))
        chosen <- sampled_pairs[[sample(rcl, 1)]]
        sol <- rbind(sol, chosen)
      }
      improved <- TRUE
      while (improved && evals < eval_budget) {
        improved <- FALSE
        fit_old <- eval_fit(sol)
        for (j in seq_len(local_n_swaps)) {
          new_sol <- mutate_cross_plan(sol, env, 1)
          fit_new <- eval_fit(new_sol)
          if (fit_new > fit_old) {
            sol <- new_sol
            fit_old <- fit_new
            improved <- TRUE
          }
        }
      }
      fit_sol <- eval_fit(sol)
      if (fit_sol > best_global_fit) {
        best_global_fit <- fit_sol
        best_global <- sol
      }
      if (iter == 1 || iter == grasp_n_iter || iter %% convergence_record_every == 0) {
        convergence[[length(convergence) + 1]] <- new_convergence_row(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name, iter, evals, fit_sol, best_global_fit)
      }
      if (evals >= eval_budget) break
    }
    selected_cross_plan <- best_global
  } else if (method == "TABU") {
    cur <- sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses)
    best_sol <- cur
    best_fit <- eval_fit(cur)
    tabu_list <- character(0)
    convergence[[length(convergence) + 1]] <- new_convergence_row(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name, 0, evals, best_fit, best_fit)
    for (iter in seq_len(tabu_n_iter)) {
      neighbors <- vector("list", tabu_n_neighbors)
      for (v in seq_len(tabu_n_neighbors)) {
        neighbors[[v]] <- mutate_cross_plan(cur, env, 1)
      }
      fits <- sapply(neighbors, eval_fit)
      ord <- order(fits, decreasing = TRUE)
      moved <- FALSE
      for (idx in ord) {
        cand_sol <- neighbors[[idx]]
        cand_hash <- paste(sort(canonical_pair_key(cand_sol)), collapse = "_")
        is_tabu <- cand_hash %in% tabu_list
        if (!is_tabu || fits[idx] > best_fit) {
          cur <- cand_sol
          if (fits[idx] > best_fit) {
            best_sol <- cur
            best_fit <- fits[idx]
          }
          tabu_list <- c(tabu_list, cand_hash)
          if (length(tabu_list) > tabu_tenure_value) tabu_list <- tabu_list[-1]
          moved <- TRUE
          break
        }
      }
      if (!moved) cur <- neighbors[[ord[1]]]
      if (iter == 1 || iter == tabu_n_iter || iter %% convergence_record_every == 0) {
        convergence[[length(convergence) + 1]] <- new_convergence_row(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name, iter, evals, max(fits, na.rm = TRUE), best_fit)
      }
    }
    selected_cross_plan <- best_sol
  } else if (method == "EDA") {
    probs <- rep(1 / length(env$candidate_pool), length(env$candidate_pool))
    names(probs) <- as.character(env$candidate_pool)
    eda_pop <- NULL
    best_global <- NULL
    best_global_fit <- -Inf
    for (iter in seq_len(eda_n_iter)) {
      eda_pop <- replicate(
        eda_pop_size,
        sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, weights = probs),
        simplify = FALSE
      )
      fits <- sapply(eda_pop, eval_fit)
      if (max(fits, na.rm = TRUE) > best_global_fit) {
        best_global_fit <- max(fits, na.rm = TRUE)
        best_global <- eda_pop[[which.max(fits)]]
      }
      elite_n <- max(1, round(eda_pop_size * eda_elite_prop))
      elite <- eda_pop[order(fits, decreasing = TRUE)[seq_len(elite_n)]]
      freq <- table(unlist(lapply(elite, as.vector)))
      probs <- rep(1e-4, length(env$candidate_pool))
      names(probs) <- as.character(env$candidate_pool)
      idx <- match(names(freq), names(probs))
      probs[idx] <- as.numeric(freq)
      probs <- probs / sum(probs)
      if (iter == 1 || iter == eda_n_iter || iter %% convergence_record_every == 0) {
        convergence[[length(convergence) + 1]] <- new_convergence_row(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name, iter, evals, max(fits, na.rm = TRUE), best_global_fit)
      }
    }
    selected_cross_plan <- best_global
  } else if (method == "PSO") {
    n_particles <- pso_n_particles
    n_cand <- length(env$candidate_pool)
    decode_particle <- function(position_vector) {
      score <- as.numeric(position_vector)
      names(score) <- as.character(env$candidate_pool)
      candidate_order <- env$candidate_pool[order(score, decreasing = TRUE)]
      plan <- matrix(integer(0), ncol = 2)
      for (cross_i in seq_len(n_crosses)) {
        added <- FALSE
        for (p1 in candidate_order) {
          if (!can_use_parent_in_partial_plan(p1, plan, env$family_ids, max_ind_per_family, max_parent_use_cross_plan)) next
          for (p2 in candidate_order) {
            if (p1 == p2) next
            if (!can_add_pair_to_partial_plan(c(p1, p2), plan, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses)) next
            plan <- rbind(plan, as.integer(c(p1, p2)))
            added <- TRUE
            break
          }
          if (added) break
        }
        if (!added) {
          weights <- score[match(as.character(env$candidate_pool), names(score))] - min(score, na.rm = TRUE) + 1e-6
          plan <- complete_cross_plan(plan, env, weights = weights)
          break
        }
      }
      storage.mode(plan) <- "integer"
      if (!validate_cross_plan(plan, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, n_crosses)) {
        weights <- score[match(as.character(env$candidate_pool), names(score))] - min(score, na.rm = TRUE) + 1e-6
        plan <- sample_valid_cross_plan(env$candidate_pool, env$family_ids, n_crosses, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, weights = weights)
      }
      plan
    }
    positions <- matrix(runif(n_particles * n_cand, min = -1, max = 1), nrow = n_particles, ncol = n_cand)
    velocities <- matrix(rnorm(n_particles * n_cand, mean = 0, sd = 0.1), nrow = n_particles, ncol = n_cand)
    pbest_positions <- positions
    pbest_fitness <- rep(-Inf, n_particles)
    pbest_solutions <- vector("list", n_particles)
    for (particle in seq_len(n_particles)) {
      sol_particle <- decode_particle(positions[particle, ])
      fit_particle <- eval_fit(sol_particle)
      pbest_fitness[particle] <- fit_particle
      pbest_solutions[[particle]] <- sol_particle
    }
    gbest_index <- which.max(pbest_fitness)
    gbest_position <- pbest_positions[gbest_index, ]
    gbest_fitness <- pbest_fitness[gbest_index]
    gbest_solution <- pbest_solutions[[gbest_index]]
    convergence[[length(convergence) + 1]] <- new_convergence_row(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name, 0, evals, gbest_fitness, gbest_fitness)
    for (iter_pso in seq_len(pso_n_iter)) {
      gbest_matrix <- matrix(gbest_position, nrow = n_particles, ncol = n_cand, byrow = TRUE)
      r1 <- matrix(runif(n_particles * n_cand), nrow = n_particles, ncol = n_cand)
      r2 <- matrix(runif(n_particles * n_cand), nrow = n_particles, ncol = n_cand)
      velocities <- pso_w * velocities + pso_c1 * r1 * (pbest_positions - positions) + pso_c2 * r2 * (gbest_matrix - positions)
      velocities[velocities > pso_vmax] <- pso_vmax
      velocities[velocities < -pso_vmax] <- -pso_vmax
      positions <- positions + velocities
      positions[positions > 10] <- 10
      positions[positions < -10] <- -10
      current_best_iter <- -Inf
      for (particle in seq_len(n_particles)) {
        sol_particle <- decode_particle(positions[particle, ])
        fit_particle <- eval_fit(sol_particle)
        if (fit_particle > current_best_iter) current_best_iter <- fit_particle
        if (fit_particle > pbest_fitness[particle]) {
          pbest_fitness[particle] <- fit_particle
          pbest_positions[particle, ] <- positions[particle, ]
          pbest_solutions[[particle]] <- sol_particle
          if (fit_particle > gbest_fitness) {
            gbest_fitness <- fit_particle
            gbest_position <- positions[particle, ]
            gbest_solution <- sol_particle
          }
        }
      }
      if (iter_pso == 1 || iter_pso == pso_n_iter || iter_pso %% convergence_record_every == 0) {
        convergence[[length(convergence) + 1]] <- new_convergence_row(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name, iter_pso, evals, current_best_iter, gbest_fitness)
      }
    }
    selected_cross_plan <- gbest_solution
  } else {
    stop("Metodo nao reconhecido.")
  }
  stats <- calc_cross_plan_solution_stats(selected_cross_plan, env$merit_z, env$K, lambda_fixed, env$kin_mean, env$kin_sd, env$family_ids, max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, n_crosses)
  list(
    indices = get_cross_plan_unique_parents(selected_cross_plan),
    cross_plan = selected_cross_plan,
    stats = stats,
    elapsed = proc.time()[3] - t0,
    evals = evals,
    convergence = if (length(convergence) > 0) do.call(rbind, convergence) else NULL
  )
}

select_classic_one_cycle <- function(pop, type, target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name) {
  t0 <- proc.time()[3]
  if (type == "Strategy 1") {
    parents_selected <- selectInd(pop, nis, use = "gv", simParam = SP)
  } else if (type == "Strategy 2") {
    fam <- selectFam(pop, nis, use = "gv", famType = "B", simParam = SP)
    parents_selected <- selectWithinFam(fam, 1, use = "gv", famType = "B", simParam = SP)
  } else {
    stop("Estrategia classica nao reconhecida.")
  }
  selected_ids <- as.character(parents_selected@iid)
  pop_ids <- as.character(pop@iid)
  selected_indices <- match(selected_ids, pop_ids)
  env <- prepare_selection_env(pop)
  stats <- calc_solution_stats(selected_indices, env$merit_z, env$K, lambda_fixed, env$kin_mean, env$kin_sd, env$family_ids, Inf)
  convergence <- new_convergence_row(target_inbreeding, observed_initial_inbreeding, scenario_cycle_found, strategy_name, 1, 1, stats$objective, stats$objective)
  list(indices = selected_indices, cross_plan = NULL, stats = stats, elapsed = proc.time()[3] - t0, evals = 1, convergence = convergence)
}

advance_strategy1_to_inbreeding_scenarios <- function(start_pop, SP, base_snp_freqs, targets, tolerance, max_cycles) {
  scenario_pops <- list()
  scenario_info <- list()
  
  start_summary <- calc_scenario_candidate_inbreeding(
    pop = start_pop,
    SP = SP,
    base_snp_freqs = base_snp_freqs,
    n_candidates = scenario_candidate_n
  )
  
  F0 <- start_summary$candidate_mean_F
  
  make_scenario_info_row <- function(target, summary_obj, cycle, status) {
    data.frame(
      TargetInbreeding = target,
      ObservedInitialInbreeding = summary_obj$candidate_mean_F,
      ScenarioCycleFound = cycle,
      AbsoluteDifference = abs(summary_obj$candidate_mean_F - target),
      Status = status,
      ScenarioCandidateN = summary_obj$candidate_n,
      ScenarioCandidateMeanF = summary_obj$candidate_mean_F,
      ScenarioCandidateVarianceF = summary_obj$candidate_variance_F,
      ScenarioCandidateSDF = summary_obj$candidate_sd_F,
      ScenarioCandidateMinF = summary_obj$candidate_min_F,
      ScenarioCandidateMaxF = summary_obj$candidate_max_F,
      ScenarioCandidateMeanGV = summary_obj$candidate_mean_GV,
      ScenarioCandidateVarianceGV = summary_obj$candidate_variance_GV,
      ScenarioCandidateSDGV = summary_obj$candidate_sd_GV,
      ScenarioCandidateMinGV = summary_obj$candidate_min_GV,
      ScenarioCandidateMaxGV = summary_obj$candidate_max_GV,
      ScenarioCandidateMeanPheno = summary_obj$candidate_mean_pheno,
      ScenarioCandidateVariancePheno = summary_obj$candidate_variance_pheno,
      ScenarioCandidateSDPheno = summary_obj$candidate_sd_pheno,
      ScenarioCandidateMinPheno = summary_obj$candidate_min_pheno,
      ScenarioCandidateMaxPheno = summary_obj$candidate_max_pheno,
      FullPopulationGenomicEndogamyF = summary_obj$full_mean_F,
      FullPopulationGenomicEndogamyVariance = summary_obj$full_variance_F,
      FullPopulationGenomicEndogamySD = summary_obj$full_sd_F,
      FullPopulationGenomicEndogamyMin = summary_obj$full_min_F,
      FullPopulationGenomicEndogamyMax = summary_obj$full_max_F,
      FullPopulationMeanGV = summary_obj$full_mean_GV,
      FullPopulationVarianceGV = summary_obj$full_variance_GV,
      FullPopulationSDGV = summary_obj$full_sd_GV,
      FullPopulationMinGV = summary_obj$full_min_GV,
      FullPopulationMaxGV = summary_obj$full_max_GV,
      FullPopulationMeanPheno = summary_obj$full_mean_pheno,
      FullPopulationVariancePheno = summary_obj$full_variance_pheno,
      FullPopulationSDPheno = summary_obj$full_sd_pheno,
      FullPopulationMinPheno = summary_obj$full_min_pheno,
      FullPopulationMaxPheno = summary_obj$full_max_pheno,
      stringsAsFactors = FALSE
    )
  }
  
  scenario_pops[[as.character(0)]] <- start_pop
  scenario_info[[as.character(0)]] <- make_scenario_info_row(
    target = 0,
    summary_obj = start_summary,
    cycle = 0,
    status = "InitialPopulation"
  )
  
  remaining_targets <- setdiff(targets, 0)
  best_pop <- list()
  best_info <- list()
  
  for (tg in remaining_targets) {
    key <- as.character(tg)
    best_pop[[key]] <- start_pop
    best_info[[key]] <- make_scenario_info_row(
      target = tg,
      summary_obj = start_summary,
      cycle = 0,
      status = "ClosestPopulation"
    )
  }
  
  pop <- start_pop
  found_targets <- character(0)
  
  history <- data.frame(
    Cycle = 0,
    ScenarioCandidateN = start_summary$candidate_n,
    CandidateGenomicEndogamyF = start_summary$candidate_mean_F,
    CandidateGenomicEndogamyVariance = start_summary$candidate_variance_F,
    CandidateGenomicEndogamySD = start_summary$candidate_sd_F,
    CandidateGenomicEndogamyMin = start_summary$candidate_min_F,
    CandidateGenomicEndogamyMax = start_summary$candidate_max_F,
    CandidateMeanGV = start_summary$candidate_mean_GV,
    CandidateVarianceGV = start_summary$candidate_variance_GV,
    CandidateSDGV = start_summary$candidate_sd_GV,
    CandidateMinGV = start_summary$candidate_min_GV,
    CandidateMaxGV = start_summary$candidate_max_GV,
    CandidateMeanPheno = start_summary$candidate_mean_pheno,
    CandidateVariancePheno = start_summary$candidate_variance_pheno,
    CandidateSDPheno = start_summary$candidate_sd_pheno,
    CandidateMinPheno = start_summary$candidate_min_pheno,
    CandidateMaxPheno = start_summary$candidate_max_pheno,
    FullPopulationGenomicEndogamyF = start_summary$full_mean_F,
    FullPopulationGenomicEndogamyVariance = start_summary$full_variance_F,
    FullPopulationGenomicEndogamySD = start_summary$full_sd_F,
    FullPopulationGenomicEndogamyMin = start_summary$full_min_F,
    FullPopulationGenomicEndogamyMax = start_summary$full_max_F,
    FullPopulationMeanGV = start_summary$full_mean_GV,
    FullPopulationVarianceGV = start_summary$full_variance_GV,
    FullPopulationSDGV = start_summary$full_sd_GV,
    FullPopulationMinGV = start_summary$full_min_GV,
    FullPopulationMaxGV = start_summary$full_max_GV,
    FullPopulationMeanPheno = start_summary$full_mean_pheno,
    FullPopulationVariancePheno = start_summary$full_variance_pheno,
    FullPopulationSDPheno = start_summary$full_sd_pheno,
    FullPopulationMinPheno = start_summary$full_min_pheno,
    FullPopulationMaxPheno = start_summary$full_max_pheno,
    MeanGv = mean(pop@gv),
    stringsAsFactors = FALSE
  )
  
  for (cycle in seq_len(max_cycles)) {
    parents_selected <- selectInd(pop, nis, use = "gv", simParam = SP)
    
    pop <- randCross(
      parents_selected,
      n_crosses,
      n_progeny_per_fam,
      ignoreSexes = TRUE
    )
    
    current_summary <- calc_scenario_candidate_inbreeding(
      pop = pop,
      SP = SP,
      base_snp_freqs = base_snp_freqs,
      n_candidates = scenario_candidate_n
    )
    
    F_current <- current_summary$candidate_mean_F
    
    history <- rbind(
      history,
      data.frame(
        Cycle = cycle,
        ScenarioCandidateN = current_summary$candidate_n,
        CandidateGenomicEndogamyF = current_summary$candidate_mean_F,
        CandidateGenomicEndogamyVariance = current_summary$candidate_variance_F,
        CandidateGenomicEndogamySD = current_summary$candidate_sd_F,
        CandidateGenomicEndogamyMin = current_summary$candidate_min_F,
        CandidateGenomicEndogamyMax = current_summary$candidate_max_F,
        CandidateMeanGV = current_summary$candidate_mean_GV,
        CandidateVarianceGV = current_summary$candidate_variance_GV,
        CandidateSDGV = current_summary$candidate_sd_GV,
        CandidateMinGV = current_summary$candidate_min_GV,
        CandidateMaxGV = current_summary$candidate_max_GV,
        CandidateMeanPheno = current_summary$candidate_mean_pheno,
        CandidateVariancePheno = current_summary$candidate_variance_pheno,
        CandidateSDPheno = current_summary$candidate_sd_pheno,
        CandidateMinPheno = current_summary$candidate_min_pheno,
        CandidateMaxPheno = current_summary$candidate_max_pheno,
        FullPopulationGenomicEndogamyF = current_summary$full_mean_F,
        FullPopulationGenomicEndogamyVariance = current_summary$full_variance_F,
        FullPopulationGenomicEndogamySD = current_summary$full_sd_F,
        FullPopulationGenomicEndogamyMin = current_summary$full_min_F,
        FullPopulationGenomicEndogamyMax = current_summary$full_max_F,
        FullPopulationMeanGV = current_summary$full_mean_GV,
        FullPopulationVarianceGV = current_summary$full_variance_GV,
        FullPopulationSDGV = current_summary$full_sd_GV,
        FullPopulationMinGV = current_summary$full_min_GV,
        FullPopulationMaxGV = current_summary$full_max_GV,
        FullPopulationMeanPheno = current_summary$full_mean_pheno,
        FullPopulationVariancePheno = current_summary$full_variance_pheno,
        FullPopulationSDPheno = current_summary$full_sd_pheno,
        FullPopulationMinPheno = current_summary$full_min_pheno,
        FullPopulationMaxPheno = current_summary$full_max_pheno,
        MeanGv = mean(pop@gv),
        stringsAsFactors = FALSE
      )
    )
    
    for (tg in remaining_targets) {
      key <- as.character(tg)
      diff_current <- abs(F_current - tg)
      diff_best <- best_info[[key]]$AbsoluteDifference[1]
      
      if (is.finite(diff_current) && diff_current < diff_best) {
        best_pop[[key]] <- pop
        best_info[[key]] <- make_scenario_info_row(
          target = tg,
          summary_obj = current_summary,
          cycle = cycle,
          status = "ClosestPopulation"
        )
      }
      
      if (!(key %in% found_targets) && is.finite(diff_current) && diff_current <= tolerance) {
        scenario_pops[[key]] <- pop
        scenario_info[[key]] <- make_scenario_info_row(
          target = tg,
          summary_obj = current_summary,
          cycle = cycle,
          status = "WithinTolerance"
        )
        
        found_targets <- c(found_targets, key)
      }
    }
    
    if (length(found_targets) == length(remaining_targets)) break
  }
  
  for (tg in remaining_targets) {
    key <- as.character(tg)
    
    if (!(key %in% names(scenario_pops))) {
      scenario_pops[[key]] <- best_pop[[key]]
      scenario_info[[key]] <- best_info[[key]]
    }
  }
  
  scenario_order <- as.character(targets)
  scenario_pops <- scenario_pops[scenario_order]
  scenario_info_df <- do.call(rbind, scenario_info[scenario_order])
  rownames(scenario_info_df) <- NULL
  
  list(
    pops = scenario_pops,
    info = scenario_info_df,
    history = history
  )
}

make_strategy_label <- function(strategy_to_run) {
  filter_label <- paste0(
    "Top", round(top_prop_optimization * 100),
    "IndexFam", max_ind_per_family,
    "_CrossPlanUse", max_parent_use_cross_plan,
    "_L", lambda_fixed,
    "_Balanced"
  )
  
  if (strategy_to_run %in% classic_strategies) {
    strategy_to_run
  } else if (strategy_to_run == "EA") {
    paste0("EA_", filter_label)
  } else {
    paste0(strategy_to_run, "_", filter_label)
  }
}


collapse_cross_plan <- function(cross_plan) {
  if (is.null(cross_plan)) return(NA_character_)
  cross_plan <- normalize_cross_plan(cross_plan)
  paste(apply(cross_plan, 1, function(x) paste(x, collapse = "x")), collapse = ";")
}

collapse_cross_plan_ids <- function(cross_plan, pop) {
  if (is.null(cross_plan)) return(NA_character_)
  cross_plan <- normalize_cross_plan(cross_plan)
  ids <- as.character(pop@iid)
  paste(apply(cross_plan, 1, function(x) paste(ids[x], collapse = "x")), collapse = ";")
}

collapse_parent_use_counts <- function(cross_plan, pop) {
  if (is.null(cross_plan)) return(NA_character_)
  ids <- as.character(pop@iid)
  tab <- parent_use_table(cross_plan)
  if (length(tab) == 0) return(NA_character_)
  parent_indices <- as.integer(names(tab))
  paste(paste(ids[parent_indices], as.integer(tab), sep = ":"), collapse = ";")
}

run_one_model_recurrent_cycles <- function(scenario_pop, strategy_to_run, target_inbreeding, observed_initial_inbreeding, scenario_cycle_found) {
  strategy_label <- make_strategy_label(strategy_to_run)
  
  cat("\n----------------------------------------\n")
  cat(sprintf("F alvo: %.2f | F observado: %.4f | Modelo: %s\n", target_inbreeding, observed_initial_inbreeding, strategy_label))
  cat(sprintf("Ciclos de selecao recorrente: %d\n", n_recurrent_cycles))
  cat("----------------------------------------\n")
  
  baseline_pop <- scenario_pop
  current_pop <- scenario_pop
  
  result_list <- list()
  selection_list <- list()
  convergence_list <- list()
  
  qtl_effects <- SP$traits[[1]]@addEff
  total_qtl <- length(qtl_effects)
  
  for (cycle in seq_len(n_recurrent_cycles)) {
    CurrentRecurrentCycle <<- cycle
    
    cat(sprintf(
      "\nCenario F alvo %.2f | Modelo %s | Ciclo recorrente %d/%d\n",
      target_inbreeding,
      strategy_label,
      cycle,
      n_recurrent_cycles
    ))
    
    if (strategy_to_run %in% classic_strategies) {
      sel_out <- select_classic_one_cycle(
        pop = current_pop,
        type = strategy_to_run,
        target_inbreeding = target_inbreeding,
        observed_initial_inbreeding = observed_initial_inbreeding,
        scenario_cycle_found = scenario_cycle_found,
        strategy_name = strategy_label
      )
    } else if (strategy_to_run == "EA") {
      sel_out <- select_ea_cross_plan(
        pop = current_pop,
        target_inbreeding = target_inbreeding,
        observed_initial_inbreeding = observed_initial_inbreeding,
        scenario_cycle_found = scenario_cycle_found,
        strategy_name = strategy_label
      )
    } else {
      sel_out <- select_metaheuristic_cross_plan(
        pop = current_pop,
        method = strategy_to_run,
        target_inbreeding = target_inbreeding,
        observed_initial_inbreeding = observed_initial_inbreeding,
        scenario_cycle_found = scenario_cycle_found,
        strategy_name = strategy_label
      )
    }
    
    selected_indices <- unique(sel_out$indices)
    
    if (strategy_to_run %in% cross_plan_strategies) {
      if (!validate_cross_plan(sel_out$cross_plan, get_family_ids(current_pop), max_ind_per_family, max_parent_use_cross_plan, allow_duplicate_crosses, n_crosses)) {
        stop("O plano de cruzamentos final violou alguma restricao.")
      }
      parents_selected <- current_pop[selected_indices]
      new_pop <- make_cross_from_plan(
        pop = current_pop,
        cross_plan = sel_out$cross_plan,
        n_progeny = n_progeny_per_fam,
        SP = SP
      )
    } else {
      parents_selected <- current_pop[selected_indices]
      new_pop <- randCross(
        parents_selected,
        n_crosses,
        n_progeny_per_fam,
        ignoreSexes = TRUE
      )
    }
    
    pop_ids <- as.character(current_pop@iid)
    selected_ids <- pop_ids[selected_indices]
    selected_families <- get_family_ids(current_pop)[selected_indices]
    number_selected_families <- length(unique(selected_families))
    selected_parents_inb <- calc_genomic_inbreeding(parents_selected, SP, base_snp_freqs)
    
    metric_out <- record_generation_metrics_recurrent(
      new_pop = new_pop,
      previous_pop = current_pop,
      baseline_pop = baseline_pop,
      SP = SP,
      qtl_effects = qtl_effects,
      total_qtl = total_qtl,
      base_snp_freqs = base_snp_freqs
    )
    
    result_row <- data.frame(
      Replication = 1L,
      TargetInbreeding = target_inbreeding,
      ObservedInitialInbreeding = observed_initial_inbreeding,
      ScenarioCycleFound = scenario_cycle_found,
      Generation = cycle,
      Strategy = strategy_label,
      SelectionMode = ifelse(strategy_to_run %in% classic_strategies, "nis_randCross", "cross_plan_makeCross"),
      MeanPhenotypicValue = metric_out$MeanPhenotypicValue,
      GeneticGain = metric_out$GeneticGain,
      GeneticGainPerCycle = metric_out$GeneticGainPerCycle,
      AdditiveVariance = metric_out$AdditiveVariance,
      PhenotypicVariance = metric_out$PhenotypicVariance,
      DominanceVariance = metric_out$DominanceVariance,
      TotalGeneticVariance = metric_out$TotalGeneticVariance,
      MeanDiagonalG = metric_out$MeanDiagonalG,
      GenomicEndogamyF = metric_out$GenomicEndogamyF,
      GenomicEndogamyVariance = metric_out$GenomicEndogamyVariance,
      GenomicEndogamySD = metric_out$GenomicEndogamySD,
      GenomicEndogamyMin = metric_out$GenomicEndogamyMin,
      GenomicEndogamyMax = metric_out$GenomicEndogamyMax,
      SelectedParentsMeanF = selected_parents_inb$genomic_F,
      SelectedParentsVarianceF = selected_parents_inb$genomic_F_variance,
      SelectedParentsSDF = selected_parents_inb$genomic_F_sd,
      SelectedParentsMinF = selected_parents_inb$genomic_F_min,
      SelectedParentsMaxF = selected_parents_inb$genomic_F_max,
      NumberSelectedFamilies = number_selected_families,
      DeltaFPerCycle = metric_out$DeltaFPerCycle,
      DeltaFCumulative = metric_out$DeltaFCumulative,
      EffectiveCrossings = metric_out$EffectiveCrossings,
      EffectiveCrossingsCumulative = metric_out$EffectiveCrossingsCumulative,
      MeanHammingDistance = metric_out$MeanHammingDistance,
      MeanFavorableAlleleFreq = metric_out$MeanFavorableAlleleFreq,
      MeanUnfavorableAlleleFreq = metric_out$MeanUnfavorableAlleleFreq,
      FixedFavorableRate = metric_out$FixedFavorableRate,
      FixedUnfavorableRate = metric_out$FixedUnfavorableRate,
      RareAlleleLossRate = metric_out$RareAlleleLossRate,
      ObservedHeterozygosity = metric_out$ObservedHeterozygosity,
      BestObjective = sel_out$stats$objective,
      BestSolutionMeanMeritZ = sel_out$stats$mean_merit_z,
      BestSolutionMeanKinship = sel_out$stats$mean_kinship,
      BestSolutionKinshipZ = sel_out$stats$kin_z,
      BestSolutionMeanUniqueKinship = ifelse(is.null(sel_out$stats$mean_unique_kinship), NA_real_, sel_out$stats$mean_unique_kinship),
      NumUniqueParents = ifelse(is.null(sel_out$stats$num_unique_parents), length(selected_indices), sel_out$stats$num_unique_parents),
      MaxParentUseObserved = ifelse(is.null(sel_out$stats$max_parent_use_observed), NA_integer_, sel_out$stats$max_parent_use_observed),
      ElapsedTime = sel_out$elapsed,
      NumEvaluations = sel_out$evals,
      ObjectivePerSecond = safe_div(sel_out$stats$objective, sel_out$elapsed),
      ObjectivePerEvaluation = safe_div(sel_out$stats$objective, sel_out$evals),
      stringsAsFactors = FALSE
    )
    
    selection_row <- data.frame(
      Replication = 1L,
      TargetInbreeding = target_inbreeding,
      ObservedInitialInbreeding = observed_initial_inbreeding,
      ScenarioCycleFound = scenario_cycle_found,
      Generation = cycle,
      Strategy = strategy_label,
      SelectionMode = ifelse(strategy_to_run %in% classic_strategies, "nis_randCross", "cross_plan_makeCross"),
      SelectedIndices = paste(selected_indices, collapse = ";"),
      SelectedIDs = paste(selected_ids, collapse = ";"),
      SelectedFamilies = paste(selected_families, collapse = ";"),
      NumberSelectedFamilies = number_selected_families,
      SelectedParentsMeanF = selected_parents_inb$genomic_F,
      SelectedParentsVarianceF = selected_parents_inb$genomic_F_variance,
      SelectedParentsSDF = selected_parents_inb$genomic_F_sd,
      SelectedParentsMinF = selected_parents_inb$genomic_F_min,
      SelectedParentsMaxF = selected_parents_inb$genomic_F_max,
      CrossPlanIndices = collapse_cross_plan(sel_out$cross_plan),
      CrossPlanIDs = collapse_cross_plan_ids(sel_out$cross_plan, current_pop),
      ParentUseCounts = collapse_parent_use_counts(sel_out$cross_plan, current_pop),
      BestObjective = sel_out$stats$objective,
      BestSolutionMeanMeritZ = sel_out$stats$mean_merit_z,
      BestSolutionMeanKinship = sel_out$stats$mean_kinship,
      BestSolutionMeanUniqueKinship = ifelse(is.null(sel_out$stats$mean_unique_kinship), NA_real_, sel_out$stats$mean_unique_kinship),
      NumUniqueParents = ifelse(is.null(sel_out$stats$num_unique_parents), length(selected_indices), sel_out$stats$num_unique_parents),
      MaxParentUseObserved = ifelse(is.null(sel_out$stats$max_parent_use_observed), NA_integer_, sel_out$stats$max_parent_use_observed),
      stringsAsFactors = FALSE
    )
    
    result_list[[cycle]] <- result_row
    selection_list[[cycle]] <- selection_row
    
    if (!is.null(sel_out$convergence) && nrow(sel_out$convergence) > 0) {
      convergence_list[[length(convergence_list) + 1]] <- sel_out$convergence
    }
    
    current_pop <- new_pop
  }
  
  CurrentRecurrentCycle <<- NA_integer_
  
  list(
    result = do.call(rbind, result_list),
    selection = do.call(rbind, selection_list),
    convergence = if (length(convergence_list) > 0) do.call(rbind, convergence_list) else data.frame()
  )
}


calc_jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  inter <- length(intersect(a, b))
  union_size <- length(union(a, b))
  if (union_size == 0) return(NA_real_)
  inter / union_size
}

build_diversity_tables <- function(selection_log) {
  pairwise_list <- list()
  group_cols <- intersect(c("Replication", "TargetInbreeding", "ObservedInitialInbreeding", "ScenarioCycleFound", "Generation"), names(selection_log))
  groups <- unique(selection_log[, group_cols, drop = FALSE])
  counter <- 1
  for (g in seq_len(nrow(groups))) {
    sub <- selection_log
    for (gc in group_cols) sub <- sub[sub[[gc]] == groups[[gc]][g], , drop = FALSE]
    if (nrow(sub) < 2) next
    for (i in seq_len(nrow(sub) - 1)) {
      for (j in (i + 1):nrow(sub)) {
        set_i <- unlist(strsplit(sub$SelectedIDs[i], ";", fixed = TRUE))
        set_j <- unlist(strsplit(sub$SelectedIDs[j], ";", fixed = TRUE))
        jac <- calc_jaccard(set_i, set_j)
        pairwise_list[[counter]] <- cbind(
          groups[g, , drop = FALSE],
          data.frame(StrategyA = sub$Strategy[i], StrategyB = sub$Strategy[j], JaccardSimilarity = jac, JaccardDistance = 1 - jac, SharedSelected = length(intersect(set_i, set_j)), UnionSelected = length(union(set_i, set_j)), stringsAsFactors = FALSE)
        )
        counter <- counter + 1
      }
    }
  }
  pairwise <- if (length(pairwise_list) > 0) do.call(rbind, pairwise_list) else data.frame()
  if (nrow(pairwise) > 0) {
    numeric_cols <- c("JaccardSimilarity", "JaccardDistance", "SharedSelected", "UnionSelected")
    summary_group_cols <- intersect(c("Replication", "TargetInbreeding", "Generation"), names(pairwise))
    summary_out <- aggregate(pairwise[numeric_cols], by = pairwise[summary_group_cols], FUN = mean, na.rm = TRUE)
  } else {
    summary_out <- data.frame()
  }
  list(pairwise = pairwise, summary = summary_out)
}

build_selected_individuals_long <- function(selection_log) {
  out_list <- vector("list", nrow(selection_log))
  for (i in seq_len(nrow(selection_log))) {
    selected_indices <- unlist(strsplit(selection_log$SelectedIndices[i], ";", fixed = TRUE))
    selected_ids <- unlist(strsplit(selection_log$SelectedIDs[i], ";", fixed = TRUE))
    selected_families <- unlist(strsplit(selection_log$SelectedFamilies[i], ";", fixed = TRUE))
    n_sel <- length(selected_ids)
    out_list[[i]] <- data.frame(
      Replication = if ("Replication" %in% names(selection_log)) selection_log$Replication[i] else 1L,
      TargetInbreeding = selection_log$TargetInbreeding[i],
      ObservedInitialInbreeding = selection_log$ObservedInitialInbreeding[i],
      ScenarioCycleFound = selection_log$ScenarioCycleFound[i],
      Generation = if ("Generation" %in% names(selection_log)) selection_log$Generation[i] else NA_integer_,
      Strategy = selection_log$Strategy[i],
      SelectionMode = selection_log$SelectionMode[i],
      SelectedRank = seq_len(n_sel),
      SelectedIndex = selected_indices,
      SelectedID = selected_ids,
      SelectedFamily = selected_families,
      NumberSelectedFamilies = if ("NumberSelectedFamilies" %in% names(selection_log)) selection_log$NumberSelectedFamilies[i] else NA_integer_,
      SelectedParentsMeanF = if ("SelectedParentsMeanF" %in% names(selection_log)) selection_log$SelectedParentsMeanF[i] else NA_real_,
      SelectedParentsVarianceF = if ("SelectedParentsVarianceF" %in% names(selection_log)) selection_log$SelectedParentsVarianceF[i] else NA_real_,
      SelectedParentsSDF = if ("SelectedParentsSDF" %in% names(selection_log)) selection_log$SelectedParentsSDF[i] else NA_real_,
      BestObjective = selection_log$BestObjective[i],
      BestSolutionMeanMeritZ = selection_log$BestSolutionMeanMeritZ[i],
      BestSolutionMeanKinship = selection_log$BestSolutionMeanKinship[i],
      stringsAsFactors = FALSE
    )
  }
  if (length(out_list) > 0) return(do.call(rbind, out_list))
  data.frame()
}

build_cross_plan_long <- function(selection_log) {
  out <- list()
  counter <- 1
  for (i in seq_len(nrow(selection_log))) {
    if (is.na(selection_log$CrossPlanIDs[i]) || selection_log$CrossPlanIDs[i] == "") next
    pairs_id <- unlist(strsplit(selection_log$CrossPlanIDs[i], ";", fixed = TRUE))
    pairs_idx <- unlist(strsplit(selection_log$CrossPlanIndices[i], ";", fixed = TRUE))
    for (k in seq_along(pairs_id)) {
      id_pair <- unlist(strsplit(pairs_id[k], "x", fixed = TRUE))
      idx_pair <- unlist(strsplit(pairs_idx[k], "x", fixed = TRUE))
      out[[counter]] <- data.frame(
        Replication = selection_log$Replication[i],
        TargetInbreeding = selection_log$TargetInbreeding[i],
        ObservedInitialInbreeding = selection_log$ObservedInitialInbreeding[i],
        ScenarioCycleFound = selection_log$ScenarioCycleFound[i],
        Generation = if ("Generation" %in% names(selection_log)) selection_log$Generation[i] else NA_integer_,
        Strategy = selection_log$Strategy[i],
        CrossNumber = k,
        Parent1Index = idx_pair[1],
        Parent2Index = idx_pair[2],
        Parent1ID = id_pair[1],
        Parent2ID = id_pair[2],
        stringsAsFactors = FALSE
      )
      counter <- counter + 1
    }
  }
  if (length(out) > 0) return(do.call(rbind, out))
  data.frame()
}

build_jaccard_matrix_long <- function(selection_log) {
  matrix_list <- list()
  group_cols <- intersect(c("Replication", "TargetInbreeding", "ObservedInitialInbreeding", "ScenarioCycleFound", "Generation"), names(selection_log))
  groups <- unique(selection_log[, group_cols, drop = FALSE])
  counter <- 1
  for (g in seq_len(nrow(groups))) {
    sub <- selection_log
    for (gc in group_cols) sub <- sub[sub[[gc]] == groups[[gc]][g], , drop = FALSE]
    if (nrow(sub) < 1) next
    strategies <- sub$Strategy
    for (i in seq_along(strategies)) {
      for (j in seq_along(strategies)) {
        set_i <- unlist(strsplit(sub$SelectedIDs[i], ";", fixed = TRUE))
        set_j <- unlist(strsplit(sub$SelectedIDs[j], ";", fixed = TRUE))
        jac <- calc_jaccard(set_i, set_j)
        matrix_list[[counter]] <- cbind(
          groups[g, , drop = FALSE],
          data.frame(StrategyA = strategies[i], StrategyB = strategies[j], JaccardSimilarity = jac, JaccardDistance = 1 - jac, SharedSelected = length(intersect(set_i, set_j)), UnionSelected = length(union(set_i, set_j)), stringsAsFactors = FALSE)
        )
        counter <- counter + 1
      }
    }
  }
  if (length(matrix_list) > 0) return(do.call(rbind, matrix_list))
  data.frame()
}

build_upset_input <- function(selection_log) {
  upset_list <- list()
  group_cols <- intersect(c("Replication", "TargetInbreeding", "ObservedInitialInbreeding", "ScenarioCycleFound", "Generation"), names(selection_log))
  groups <- unique(selection_log[, group_cols, drop = FALSE])
  counter <- 1
  for (g in seq_len(nrow(groups))) {
    sub <- selection_log
    for (gc in group_cols) sub <- sub[sub[[gc]] == groups[[gc]][g], , drop = FALSE]
    all_ids <- unique(unlist(strsplit(sub$SelectedIDs, ";", fixed = TRUE)))
    if (length(all_ids) == 0) next
    out <- cbind(groups[rep(g, length(all_ids)), , drop = FALSE], data.frame(IndividualID = all_ids, stringsAsFactors = FALSE))
    for (s in sub$Strategy) {
      ids_s <- unlist(strsplit(sub$SelectedIDs[sub$Strategy == s][1], ";", fixed = TRUE))
      safe_s <- gsub("[^A-Za-z0-9_]", "_", s)
      out[[safe_s]] <- as.integer(out$IndividualID %in% ids_s)
    }
    upset_list[[counter]] <- out
    counter <- counter + 1
  }
  if (length(upset_list) > 0) return(do.call(rbind, upset_list))
  data.frame()
}

sanitize_filename <- function(x) {
  x <- gsub(" ", "_", x)
  x <- gsub("\\.", "_", x)
  x <- gsub("[^A-Za-z0-9_-]", "_", x)
  x <- gsub("_+", "_", x)
  x
}

make_unique_sheet_names <- function(x) {
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  x <- gsub("\\.", "_", x)
  x <- gsub("[^A-Za-z0-9_-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- ifelse(nchar(x) == 0, "Sheet", x)
  out <- character(length(x))
  used <- character(0)
  for (i in seq_along(x)) {
    base <- substr(x[i], 1, 31)
    candidate <- base
    k <- 1
    while (candidate %in% used) {
      suffix <- paste0("_", k)
      candidate <- paste0(substr(base, 1, 31 - nchar(suffix)), suffix)
      k <- k + 1
    }
    out[i] <- candidate
    used <- c(used, candidate)
  }
  out
}

cat("\n============================================================\n")
cat("EXECUCAO PAIRWISE COM CENARIOS DE ENDOGAMIA INICIAL E SELECAO RECORRENTE\n")
cat("============================================================\n")
cat("Modelos:", paste(strategies_to_run, collapse = ", "), "\n")
cat("Endogamias alvo:", paste(target_inbreeding_values, collapse = ", "), "\n")
cat("Ciclos recorrentes por modelo:", n_recurrent_cycles, "\n")
cat("Cenarios definidos pela media de F dos", scenario_candidate_n, "candidatos melhores pelo indice GV padronizado - F individual.\n")
cat("Strategy 1 e Strategy 2 mantem nis e randCross.\n")
cat("EA, Local, SA, GA, GRASP, TABU, EDA e PSO usam matriz 50 x 2 de cruzamentos.\n")
cat("Cada genitor pode aparecer no maximo", max_parent_use_cross_plan, "vezes na matriz.\n")
cat("============================================================\n\n")

all_results_list <- list()
all_selection_list <- list()
all_convergence_list <- list()

result_counter <- 1
selection_counter <- 1
convergence_counter <- 1

set.seed(3532357)

scenario_out <- advance_strategy1_to_inbreeding_scenarios(
  start_pop = Prog,
  SP = SP,
  base_snp_freqs = base_snp_freqs,
  targets = target_inbreeding_values,
  tolerance = scenario_tolerance,
  max_cycles = max_strategy1_cycles_to_find
)

scenario_info <- scenario_out$info
scenario_history <- scenario_out$history

for (tg in target_inbreeding_values) {
  key <- as.character(tg)
  scenario_pop <- scenario_out$pops[[key]]
  info_row <- scenario_info[scenario_info$TargetInbreeding == tg, ]
  observed_initial_inbreeding <- info_row$ObservedInitialInbreeding[1]
  scenario_cycle_found <- info_row$ScenarioCycleFound[1]
  cat("\n============================================================\n")
  cat(sprintf(
    "Cenario alvo F = %.2f | F observado nos %d candidatos = %.4f | ciclo Strategy 1 = %d | status = %s\n",
    tg,
    scenario_candidate_n,
    observed_initial_inbreeding,
    scenario_cycle_found,
    info_row$Status[1]
  ))
  cat("============================================================\n")
  for (strategy_to_run in strategies_to_run) {
    strategy_id <- match(strategy_to_run, strategy_codes)
    set.seed(3532357 + 100000 * strategy_id + round(tg * 10000))
    out <- run_one_model_recurrent_cycles(
      scenario_pop = scenario_pop,
      strategy_to_run = strategy_to_run,
      target_inbreeding = tg,
      observed_initial_inbreeding = observed_initial_inbreeding,
      scenario_cycle_found = scenario_cycle_found
    )
    all_results_list[[result_counter]] <- out$result
    result_counter <- result_counter + 1
    all_selection_list[[selection_counter]] <- out$selection
    selection_counter <- selection_counter + 1
    if (!is.null(out$convergence) && nrow(out$convergence) > 0) {
      all_convergence_list[[convergence_counter]] <- out$convergence
      convergence_counter <- convergence_counter + 1
    }
    cat(sprintf("Concluido: F alvo %.2f | %s | %d ciclos recorrentes\n", tg, unique(out$result$Strategy)[1], n_recurrent_cycles))
  }
}

combined_results <- do.call(rbind, all_results_list)
rownames(combined_results) <- NULL

selection_log <- do.call(rbind, all_selection_list)
rownames(selection_log) <- NULL

convergence_log <- if (length(all_convergence_list) > 0) do.call(rbind, all_convergence_list) else data.frame()

diversity_tables <- build_diversity_tables(selection_log)
diversity_pairwise <- diversity_tables$pairwise
diversity_summary <- diversity_tables$summary
selected_individuals_long <- build_selected_individuals_long(selection_log)
cross_plan_long <- build_cross_plan_long(selection_log)
jaccard_matrix_long <- build_jaccard_matrix_long(selection_log)
upset_input <- build_upset_input(selection_log)

run_parameters <- data.frame(
  Parameter = c(
    "target_inbreeding_values",
    "scenario_tolerance",
    "max_strategy1_cycles_to_find",
    "n_recurrent_cycles",
    "scenario_candidate_n",
    "nQTL",
    "nSNP",
    "gen_length",
    "nis_classic_strategies",
    "n_progeny_per_fam",
    "n_crosses",
    "nThreads",
    "lambda_fixed",
    "top_prop_optimization",
    "max_ind_per_family",
    "max_parent_use_cross_plan",
    "allow_duplicate_crosses",
    "min_unique_parents_cross_plan",
    "eval_budget",
    "ea_pop_size",
    "ea_n_gen",
    "ea_elite_n",
    "ea_tournament_k",
    "ea_mutation_swaps",
    "ea_random_injection_rate",
    "local_n_swaps",
    "sa_n_iter",
    "sa_temp_ini",
    "sa_cooling",
    "ga_pop_size",
    "ga_n_gen",
    "ga_elite_n",
    "ga_mutation_rate",
    "grasp_alpha",
    "tabu_n_iter",
    "tabu_n_neighbors",
    "tabu_tenure_value",
    "eda_pop_size",
    "eda_n_iter",
    "eda_elite_prop",
    "pso_n_particles",
    "pso_n_iter",
    "pso_w",
    "pso_c1",
    "pso_c2",
    "pso_vmax",
    "convergence_record_every",
    "cross_construction_pair_candidates"
  ),
  Value = c(
    paste(target_inbreeding_values, collapse = ", "),
    scenario_tolerance,
    max_strategy1_cycles_to_find,
    n_recurrent_cycles,
    scenario_candidate_n,
    nQTL,
    nSNP,
    gen_length,
    nis,
    n_progeny_per_fam,
    n_crosses,
    nThreads,
    lambda_fixed,
    top_prop_optimization,
    max_ind_per_family,
    max_parent_use_cross_plan,
    allow_duplicate_crosses,
    min_unique_parents_cross_plan,
    eval_budget,
    ea_pop_size,
    ea_n_gen,
    ea_elite_n,
    ea_tournament_k,
    ea_mutation_swaps,
    ea_random_injection_rate,
    local_n_swaps,
    sa_n_iter,
    sa_temp_ini,
    sa_cooling,
    ga_pop_size,
    ga_n_gen,
    ga_elite_n,
    ga_mutation_rate,
    grasp_alpha,
    tabu_n_iter,
    tabu_n_neighbors,
    tabu_tenure_value,
    eda_pop_size,
    eda_n_iter,
    eda_elite_prop,
    pso_n_particles,
    pso_n_iter,
    pso_w,
    pso_c1,
    pso_c2,
    pso_vmax,
    convergence_record_every,
    cross_construction_pair_candidates
  ),
  stringsAsFactors = FALSE
)

parameter_balance <- data.frame(
  Model = c("EA", "Local", "SA", "GA", "GRASP", "TABU", "EDA", "PSO"),
  Main_Parameter_1 = c(ea_pop_size, eval_budget, sa_n_iter, ga_pop_size, eval_budget, tabu_n_iter, eda_pop_size, pso_n_particles),
  Main_Parameter_2 = c(ea_n_gen, NA, NA, ga_n_gen, NA, tabu_n_neighbors, eda_n_iter, pso_n_iter),
  Approx_Main_Evaluations_Per_Cycle = c(
    ea_pop_size * ea_n_gen,
    eval_budget,
    sa_n_iter,
    ga_pop_size * ga_n_gen,
    eval_budget,
    tabu_n_iter * tabu_n_neighbors,
    eda_pop_size * eda_n_iter,
    pso_n_particles * pso_n_iter
  ),
  Approx_Main_Evaluations_5_Cycles = c(
    ea_pop_size * ea_n_gen,
    eval_budget,
    sa_n_iter,
    ga_pop_size * ga_n_gen,
    eval_budget,
    tabu_n_iter * tabu_n_neighbors,
    eda_pop_size * eda_n_iter,
    pso_n_particles * pso_n_iter
  ) * n_recurrent_cycles,
  stringsAsFactors = FALSE
)

pasta_saida_geral <- file.path(pasta_saida_base, "Todos_os_modelos")
if (!dir.exists(pasta_saida_geral)) dir.create(pasta_saida_geral, recursive = TRUE)

arquivo_excel <- file.path(
  pasta_saida_geral,
  "resultados_cenarios_endogamia_5_ciclos_candidatos300_pairwise_reduzido5x_sem_repeticoes.xlsx"
)

strategy_tables <- split(combined_results, combined_results$Strategy)
names(strategy_tables) <- make_unique_sheet_names(names(strategy_tables))

sheets_to_write <- c(
  list(
    Combined_Results = combined_results,
    Scenario_Info = scenario_info,
    Scenario_History = scenario_history,
    Selection_Log = selection_log,
    Selected_Individuals_Long = selected_individuals_long,
    CrossPlan_Long = cross_plan_long,
    Convergence_Log = convergence_log,
    Diversity_Pairwise = diversity_pairwise,
    Diversity_Summary = diversity_summary,
    Jaccard_Matrix_Long = jaccard_matrix_long,
    UpSet_Input = upset_input,
    Run_Parameters = run_parameters,
    Parameter_Balance = parameter_balance
  ),
  strategy_tables
)

openxlsx::write.xlsx(sheets_to_write, file = arquivo_excel, overwrite = TRUE)

cat("\n============================================================\n")
cat("EXECUCAO FINALIZADA\n")
cat("============================================================\n")
cat("Arquivo Excel:", normalizePath(arquivo_excel), "\n")
cat("============================================================\n\n")
