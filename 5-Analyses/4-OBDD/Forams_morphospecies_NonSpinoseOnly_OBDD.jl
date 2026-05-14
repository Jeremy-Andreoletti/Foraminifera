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
tree = read_newick("../../3-Data_processed/Morphospecies_phylogenies/NonSpinoseOnly_sampledAncestors.tree", true)
occurrence_ages = readdlm("../../3-Data_processed/Triton_occurrences/TritonDB_subsampled_STT_34846occurrences_NonSpinoseOnly.csv", ';')[:]

# # Sample a tree and occurrences
n_samples = 5000
ωtimes = occurrence_ages[randperm(length(occurrence_ages))[1:n_samples]]

### OBDD inference ###

veryShortMCMC = false
shortMCMC = false

λa_prior = (1.5, 1.0)
μa_prior = (1.5, 1.0)
αλ_prior = (0.0, 0.1)
αμ_prior = (0.0, 0.1)
σλ_prior = (5.0, 0.5)
σμ_prior = (5.0, 0.5)
ψ_prior  = (1.0, 0.5)
ω_prior  = (1.0, 0.1)
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
ψi       = ψ_prior[1]/ψ_prior[2]
ωi       = ω_prior[1]/ω_prior[2]
αλi      = 0.0
αμi      = 0.0
σλi      = σλ_prior[2]/(σλ_prior[1]-1)
σμi      = σμ_prior[2]/(σμ_prior[1]-1)
pupdp    = (0.02, 0.02, 0.1, 0.01, 0.01, 0.02, 0.1, 0.2)
δt       = 1e-3
survival = true
mxthf    = 0.05
prints   = 5
stnλ     = 0.5
stnμ     = 0.5
tρ       = Dict("" => 1.0)

ofile    = "OBDD_Forams_NonSpinoseOnly_sampledAncestors_$(niter)iter_seed$(seed_nb)" ; isdir("outputs/") || mkdir("outputs/")

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


