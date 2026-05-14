include(homedir()*"/Nextcloud/Recherche/1_Methods/INSANE/Source_INSANE.jl");

# using Distributions
# using DataFrames
using DelimitedFiles
# using StatsBase
# using Random: seed!, randperm
# using Distributed
using Plots

include("Forams_plot.jl")

# Automatically find all output files in the "outputs" directory with a specific pattern
output_files = readdir("outputs"; join=true)  # Get full paths
output_files = filter(x -> occursin(r".txt", x), output_files)
output_files = filter(x -> !occursin(r"\.\_", x), output_files)
output_files = filter(x -> !occursin(r"ωtimes", x), output_files)

# Relevant parameters
δt       = 1e-3

# Group files by their common model name, ignoring the seed part
model_groups = Dict{String, Vector{String}}()
for file in output_files
    if occursin(r"_seed", file)
        # Extract the part before "_seed" as the model name
        model_name = match(r"(.*)_seed", basename(file)).captures[1]
        if haskey(model_groups, model_name)
            push!(model_groups[model_name], file)
        else
            model_groups[model_name] = [file]
        end
    end
end
model_groups

# for (model_name, files) in Dict("OBDD_Forams_Spinose_50000000" => model_groups["OBDD_Forams_Spinose_50000000"])
for (model_name, files) in model_groups
  println("Processing model: $model_name")
  
  for file in files
    println("File: $file")

    ofile = match(r"OBDD(.*)[0-9]", basename(file)).match
    println("ofile: $ofile")
    
    out_trees = iread("outputs/"*ofile*".txt")
    tree = read_newick("outputs/$(ofile).tre", true)
    ωtimes = readdlm("outputs/$(ofile)_ωtimes.txt")[:,1]

    Forams_plot(ofile, out_trees, tree, ωtimes)

  end
end

# Merge outputs, trees, and ωtimes from seeds for each model group before plotting merged results.
for (model_name, files) in model_groups
    println("Processing model (merged seeds): $model_name")
    
    # Prepare arrays to store data per seed.
    out_trees_by_seed = [];
    trees_by_seed = [];
    ωtimes_by_seed = [];
    
    for file in files
        println("Reading file: $file")
        trees = iread(file);
        push!(out_trees_by_seed, trees)
        
        # Assume each seed has its own ωtimes file.
        seed_ω_file = replace(file, ".txt" => "_ωtimes.txt");
        seed_ω = readdlm(seed_ω_file)[:,1];
        push!(ωtimes_by_seed, seed_ω)
        
        # Extract the common identifier and read the tree for this seed.
        ofile = match(r"OBDD(.*)[0-9]", basename(file)).match;
        tree_i = read_newick("outputs/$(ofile).tre", true);
        push!(trees_by_seed, tree_i)
    end

    # Merge out_trees and ωtimes across seeds.
    out_trees_merged = reduce(vcat, out_trees_by_seed);
    ωtimes_merged = reduce(vcat, ωtimes_by_seed);
    
    # Use the first file in the group for the common identifier.
    ofile = match(r"(.*)_seed", basename(files[1])).captures[1];
    println("Using ofile: $ofile")
    
    Forams_plot(ofile, out_trees_merged, trees_by_seed, ωtimes_merged)
end

