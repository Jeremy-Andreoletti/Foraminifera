# Set working directory

include(homedir()*"/Nextcloud/Recherche/1_Methods/INSANE/Source_INSANE.jl");

using Distributions
using DataFrames
using DelimitedFiles
using StatsBase
using Random: seed!, randperm
using Distributed
using Plots

# using Profile
# using PProf

seed_nb = isempty(ARGS) ? 0 : parse(Int64, ARGS[1])
seed!(seed_nb)

## Read in data
# Trees and occurrences
treesDLM_0 = readdlm("../../4-Phylogenetic_reconstruction/output_pf95_pt95_nodeDating/forams_Spinose_seed0.trees", header=true)
treesDLM_1 = readdlm("../../4-Phylogenetic_reconstruction/output_pf95_pt95_nodeDating/forams_Spinose_seed1.trees", header=true)
treesDLM_2 = readdlm("../../4-Phylogenetic_reconstruction/output_pf95_pt95_nodeDating/forams_Spinose_seed2.trees", header=true)
trees = DataFrame(vcat(treesDLM_0[1], treesDLM_1[1], treesDLM_2[1]), vec(treesDLM_0[2]))
occurrence_ages = readdlm("../../3-Data_processed/Triton_occurrences/TritonDB_subsampled_STT_43817occurrences_Spinose.csv", ';')[:]

# Origin times
logsDLM_0 = readdlm("../../4-Phylogenetic_reconstruction/output_pf95_pt95_nodeDating/forams_Spinose_seed0.log", header=true)
logsDLM_1 = readdlm("../../4-Phylogenetic_reconstruction/output_pf95_pt95_nodeDating/forams_Spinose_seed1.log", header=true)
logsDLM_2 = readdlm("../../4-Phylogenetic_reconstruction/output_pf95_pt95_nodeDating/forams_Spinose_seed2.log", header=true)
logs_df = DataFrame(vcat(logsDLM_0[1], logsDLM_1[1], logsDLM_2[1]), vec(logsDLM_0[2]))

# # Sample a tree and occurrences
idx = rand(findall(trees.Iteration .>= 5000))
tree = _parse_newick(String(trees.obd_tree[idx]), accerr, false)
n_samples = 1000
ωtimes = occurrence_ages[randperm(length(occurrence_ages))[1:n_samples]]
    
# Adjust the stem branch length
sete!(tree, logs_df.origin_time[idx] - treeheight(tree))


### OBDD inference ###

veryShortMCMC = false
shortMCMC = false

λa_prior = (1.5, 1.0)
μa_prior = (1.5, 1.0)
αλ_prior = (0.0, 0.1)
αμ_prior = (0.0, 0.1)
σλ_prior = (5.0, 0.5)
σμ_prior = (5.0, 0.5)
ψ_prior  = (1.0, 10000.0)
ω_prior  = (1.0, 0.5)
# ψω_epoch = Float64[33.9, 28.1, 23.03, 20.44, 7.246, 5.333]
ψω_epoch = Float64[]
f_epoch  = Int64[0]
niter    = veryShortMCMC ? 100 : (shortMCMC ? 50_000 : 50_000_000)
nthin    = veryShortMCMC ? 1 : (shortMCMC ? niter+1 : niter+1)
nburn    = veryShortMCMC ? 10 : (shortMCMC ? 1_000 : 50_000)
nflushθ  = Int64(ceil(niter/20_000))
nflushΞ  = Int64(ceil(niter/100))
tune_int = 100
ϵi       = 0.2
λi       = prod(λa_prior)
μi       = prod(μa_prior)
ψi       = 0.0
ωi       = ω_prior[1]/ω_prior[2]
αλi      = 0.0
αμi      = 0.0
σλi      = 0.1
σμi      = 0.1
pupdp    = (0.02, 0.02, 0.1, 0.0, 0.01, 0.02, 0.1, 0.2)
δt       = 1e-3
survival = true
mxthf    = 0.1
prints   = 5
stnλ     = 0.5
stnμ     = 0.5
tρ       = Dict("" => 1.0)

ofile    = "OBDD_Forams_Spinose_pf95_pt95_nodeDating_$(niter)iter_mxthf$(join(split(string(mxthf), ".")))_seed$(seed_nb)" ; isdir("outputs/") || mkdir("outputs/")

writedlm("outputs/$(ofile)_ωtimes.txt", ωtimes)
write_newick(tree, "outputs/$(ofile)")

seed!(seed_nb); insane_gbmobd(tree::sTf_label,
                              ωtimes,
                              λa_prior = λa_prior,
                              μa_prior = μa_prior,
                              αλ_prior = αλ_prior,
                              αμ_prior = αμ_prior,
                              σλ_prior = σλ_prior,
                              σμ_prior = σμ_prior,
                              ψ_prior  = ψ_prior,
                              ω_prior  = ω_prior,
                              ψω_epoch = ψω_epoch,
                              f_epoch  = f_epoch,
                              niter    = niter,
                              nthin    = nthin,
                              nburn    = nburn,
                              nflushθ  = nflushθ,
                              nflushΞ  = nflushΞ,
                              ofile    = "outputs/"*ofile,
                              tune_int = tune_int,
                              ϵi       = ϵi,
                              λi       = λi,
                              μi       = μi,
                              ψi       = ψi,
                              ωi       = ωi,
                              αλi      = αλi,
                              αμi      = αμi,
                              σλi      = σλi,
                              σμi      = σμi,
                              pupdp    = pupdp,
                              δt       = δt,
                              survival = survival,
                              mxthf    = mxthf,
                              prints   = prints,
                              stnλ     = stnλ,
                              stnμ     = stnμ,
                              tρ       = tρ)

include("Forams_plot.jl")
Forams_plot(ofile, tree, ωtimes)


