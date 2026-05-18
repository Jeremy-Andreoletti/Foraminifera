suppressMessages({
  library(ape); library(hisse); library(phangorn); library(paleotree); library(dplyr)
})

functional_df <- read.csv("../4-OBDD/Functional_data_newTaxo.csv", row.names = 1)
load("../../3-Data_processed/Morphospecies_phylogenies/NonSpinoseOnly_sampledAncestors.RData")
load("../../3-Data_processed/Morphospecies_phylogenies/SpinoseOnly_sampledAncestors.RData")
load("../../3-Data_processed/Morphospecies_phylogenies/MixedSpinose_sampledAncestors.RData")

eco_to_obs <- function(x) {
  ifelse(x %in% c("Symbiotic SML", "Asymbiotic SML"), 1L,
  ifelse(x == "Thermocline",                          2L,
  ifelse(x == "Sub-thermocline",                      3L, NA_integer_)))
}
encode_3_to_2 <- function(s) {
  M <- matrix(NA_integer_, nrow = length(s), ncol = 2)
  i <- !is.na(s) & s == 1; M[i, 1] <- 0L; M[i, 2] <- 1L
  i <- !is.na(s) & s == 2; M[i, 1] <- 0L; M[i, 2] <- 0L
  i <- !is.na(s) & s == 3; M[i, 1] <- 1L; M[i, 2] <- 0L
  M
}

extract_k_and_collapse <- function(tree, functional_df) {
  keep <- !is.na(eco_to_obs(functional_df$eco[match(tree$tip.species, functional_df$Morphospecies)]))
  tree <- drop.tip(tree, tree$tip.label[!keep], collapse.singles = FALSE)
  tree$tip.species  <- tree$tip.species[tree$tip.label]
  tree$node.species <- tree$node.species[tree$node.label]
  tree_labels <- c(tree$tip.label, tree$node.label)
  k_nodes <- grep("F", tree$node.label, value = TRUE)
  k_descendants <- sapply(k_nodes, function(k_node) {
    k_children <- phangorn::Children(tree, k_node)
    while (any(grepl("F|N", tree_labels[k_children]))) {
      if (length(k_children) == 1) k_children <- phangorn::Children(tree, k_children)
      else if (grepl("F|N", tree_labels[k_children])[1]) k_children[1] <- phangorn::Children(tree, k_children[1])[[1]]
      else k_children[2] <- phangorn::Children(tree, k_children[2])[[1]]
    }
    if (length(k_children) == 1) c(tree_labels[k_children], tree_labels[k_children]) else tree_labels[k_children]
  })
  k_samples <- data.frame(k_nodes, t(k_descendants))
  names(k_samples) <- c("node_label", "taxon1", "taxon2")
  all_ages <- paleotree::dateNodes(tree = tree); names(all_ages) <- tree_labels
  k_samples$timefrompresent <- all_ages[k_samples$node_label]
  k_samples$species <- tree$node.species[k_samples$node_label]
  k_samples$eco_state <- eco_to_obs(functional_df$eco[match(k_samples$species, functional_df$Morphospecies)])
  k_samples$state1 <- encode_3_to_2(k_samples$eco_state)[, 1]
  k_samples$state2 <- encode_3_to_2(k_samples$eco_state)[, 2]
  k_samples <- k_samples[!is.na(k_samples$state1), ]
  list(phy = ape::collapse.singles(tree), k.samples = k_samples)
}

clade_inputs <- list(
  Mixed      = extract_k_and_collapse(Phylo_MixedSpinose_sampledAncestors,   functional_df),
  Spinose    = extract_k_and_collapse(Phylo_SpinoseOnly_sampledAncestors,    functional_df),
  NonSpinose = extract_k_and_collapse(Phylo_NonSpinoseOnly_sampledAncestors, functional_df)
)

# Three model variants
# - null_shared: hidden=1, make.null=T, rates shared between hidden states (baseline)
# - null_free  : hidden=1, make.null=T, rates FREE between hidden states A and B
# - no_hidden  : hidden=0 (plain MuSSE), 3 observed states only
fit_model <- function(tree, k.samples, model) {
  obs <- eco_to_obs(functional_df$eco[match(tree$tip.species, functional_df$Morphospecies)])
  X   <- encode_3_to_2(obs)
  dat <- data.frame(taxon = tree$tip.label, stateA = X[, 1], stateB = X[, 2],
                    stringsAsFactors = FALSE, row.names = tree$tip.label)
  ks  <- k.samples[c("taxon1", "taxon2", "timefrompresent", "state1", "state2")]

  if (model == "no_hidden") {
    trans     <- TransMatMakerMuHiSSE(hidden.traits = 0)
    trans.mod <- ParDrop(trans, c(4, 6, 7, 8))
    turnover  <- c(1, 2, 3, 0)
    eps       <- c(1, 2, 3, 0)
    root.p    <- c(0.5, 0, 0, 0)  # Surface or Thermocline at root (50/50 of valid two)
    hidden    <- FALSE
  } else if (model == "null_shared") {
    trans     <- TransMatMakerMuHiSSE(hidden.traits = 1, make.null = TRUE)
    trans.mod <- ParDrop(trans, c(4, 6, 7, 8))
    turnover  <- c(1, 2, 3, 0, 1, 2, 3, 0)
    eps       <- c(1, 2, 3, 0, 1, 2, 3, 0)
    root.p    <- c(0.5, 0, 0, 0, 0.5, 0, 0, 0)
    hidden    <- TRUE
  } else if (model == "null_free") {
    trans     <- TransMatMakerMuHiSSE(hidden.traits = 1, make.null = TRUE)
    trans.mod <- ParDrop(trans, c(4, 6, 7, 8))
    turnover  <- c(1, 2, 3, 0, 4, 5, 6, 0)
    eps       <- c(1, 2, 3, 0, 4, 5, 6, 0)
    root.p    <- c(0.5, 0, 0, 0, 0.5, 0, 0, 0)
    hidden    <- TRUE
  }

  fit <- MuHiSSE(phy = tree, data = dat, f = c(1, 1, 1, 0),
                 turnover = turnover, turnover.upper = 10,
                 eps = eps, eps.upper = 2,
                 trans.rate = trans.mod, k.samples = ks,
                 hidden.states = hidden, includes.fossils = TRUE,
                 root.p = root.p,
                 sann = FALSE, starting.vals = c(0.1, 0.5, 0.01),
                 ode.eps = 0)
  fit$tree <- tree; fit$dat <- dat; fit$k.samples <- ks; fit$model <- model
  fit
}

summary_tab <- data.frame()
for (model in c("no_hidden", "null_shared", "null_free")) {
  for (clade in names(clade_inputs)) {
    out_file <- sprintf("FMuHiSSE_fit_%s_%s.RData", model, clade)
    var_name <- sprintf("fit_%s_%s", model, clade)
    if (file.exists(out_file)) {
      message("Skip ", out_file)
      load(out_file)
    } else {
      message(sprintf("[%s] Fitting %s / %s ...", format(Sys.time(), "%H:%M:%S"), model, clade))
      t0  <- Sys.time()
      fit <- fit_model(clade_inputs[[clade]]$phy, clade_inputs[[clade]]$k.samples, model)
      elapsed <- format(Sys.time() - t0, digits = 3)
      message(sprintf("[%s]   done in %s. LogLik=%.2f AIC=%.2f",
                      format(Sys.time(), "%H:%M:%S"), elapsed, fit$loglik, fit$AIC))
      assign(var_name, fit)
      save(list = var_name, file = out_file)
    }
    fit <- get(var_name)
    summary_tab <- rbind(summary_tab, data.frame(
      model = model, clade = clade, loglik = fit$loglik, AIC = fit$AIC,
      n_params = sum(fit$solution != 0)
    ))
  }
}

cat("\n=== Summary ===\n")
print(summary_tab, row.names = FALSE)

cat("\n=== Delta AIC per clade (vs best model) ===\n")
for (cl in unique(summary_tab$clade)) {
  s <- summary_tab[summary_tab$clade == cl, ]
  s$dAIC <- s$AIC - min(s$AIC)
  cat(sprintf("--- %s ---\n", cl)); print(s[order(s$AIC), c("model","loglik","AIC","dAIC","n_params")], row.names = FALSE)
}
