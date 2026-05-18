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
      if (length(k_children) == 1) {
        k_children <- phangorn::Children(tree, k_children)
      } else if (grepl("F|N", tree_labels[k_children])[1]) {
        k_children[1] <- phangorn::Children(tree, k_children[1])[[1]]
      } else {
        k_children[2] <- phangorn::Children(tree, k_children[2])[[1]]
      }
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
  NonSpinose = extract_k_and_collapse(Phylo_NonSpinoseOnly_sampledAncestors, functional_df),
  Spinose    = extract_k_and_collapse(Phylo_SpinoseOnly_sampledAncestors,    functional_df)
)

run_simple <- function(tree, k.samples) {
  obs <- eco_to_obs(functional_df$eco[match(tree$tip.species, functional_df$Morphospecies)])
  X   <- encode_3_to_2(obs)
  dat <- data.frame(taxon = tree$tip.label, stateA = X[, 1], stateB = X[, 2],
                    stringsAsFactors = FALSE, row.names = tree$tip.label)

  trans     <- TransMatMakerMuHiSSE(hidden.traits = 1, make.null = TRUE)
  trans.mod <- ParDrop(trans, c(4, 6, 7, 8, 12, 14, 15, 16))

  # Shared across hidden states: 3 free turnover, 3 free eps (vs 6 + 6 before)
  turnover <- c(1, 2, 3, 0, 1, 2, 3, 0)
  eps      <- c(1, 2, 3, 0, 1, 2, 3, 0)

  fit <- MuHiSSE(phy = tree, data = dat, f = c(1, 1, 1, 0),
                 turnover = turnover, turnover.upper = 10,
                 eps = eps, eps.upper = 2,
                 trans.rate = trans.mod,
                 k.samples = k.samples[c("taxon1", "taxon2", "timefrompresent", "state1", "state2")],
                 hidden.states = TRUE, includes.fossils = TRUE,
                 root.p = c(0.5, 0, 0, 0, 0.5, 0, 0, 0),
                 sann = FALSE, starting.vals = c(0.1, 0.5, 0.01),
                 ode.eps = 0)
  fit$tree <- tree; fit$dat <- dat
  fit$k.samples <- k.samples[c("taxon1", "taxon2", "timefrompresent", "state1", "state2")]
  fit
}

for (clade in names(clade_inputs)) {
  out_file <- sprintf("FMuHiSSE_fit_simple_%s.RData", clade)
  if (file.exists(out_file)) { message("Skip ", clade, " (exists)"); next }
  message(sprintf("[%s] Fitting %s (tips=%d, k_samples=%d)...",
                  format(Sys.time(), "%H:%M:%S"), clade,
                  Ntip(clade_inputs[[clade]]$phy), nrow(clade_inputs[[clade]]$k.samples)))
  t0 <- Sys.time()
  fit_obj <- run_simple(clade_inputs[[clade]]$phy, clade_inputs[[clade]]$k.samples)
  message(sprintf("[%s]   done in %s. LogLik=%.3f AIC=%.3f",
                  format(Sys.time(), "%H:%M:%S"),
                  format(Sys.time() - t0, digits = 3),
                  fit_obj$loglik, fit_obj$AIC))
  assign(paste0("fit_", clade), fit_obj)
  save(list = paste0("fit_", clade), file = out_file)
}

message("All done.")
