# -------------------------------------------------------------------------
# Purpose: For each seed replicate, read the chain of trees, pick
#          the median tree, extract birth and death rates for each tip, 
#          and save them along with the tip labels.
# -------------------------------------------------------------------------

using DelimitedFiles
using DataFrames
using CSV

include(homedir()*"/Nextcloud/Recherche/1_Methods/INSANE/Source_INSANE.jl")

# 1) Identify all relevant output files in the "outputs" directory.
output_files = readdir("outputs"; join=true)
output_files = filter(x -> occursin(r".txt", x), output_files)
output_files = filter(x -> !occursin(r"\.\_", x), output_files)
output_files = filter(x -> !occursin(r"_ωtimes.txt", x), output_files)

# 2) Group files by their common model name
model_groups = Dict{String, Vector{String}}()
for file in output_files
    if occursin(r"_seed", file)
        # Extract everything before “_seed”
        model_name = match(r"(.*)_seed", basename(file)).captures[1]
        if haskey(model_groups, model_name)
            push!(model_groups[model_name], file)
        else
            model_groups[model_name] = [file]
        end
    end
end

# 3) For each model group and each seed replicate
for (model_name, files) in model_groups
    println("Processing model: $model_name")

    for file in files
        println("  Processing replicate: $file")

        # Load the chain of trees, get tip rates
        out_trees = iread(file);
        tip_birth_rates = reduce(hcat, tip_rates.(remove_unsampled.(out_trees), b))
        tip_death_rates = reduce(hcat, tip_rates.(remove_unsampled.(out_trees), d))

        tree = read_newick(replace(file, "txt"=>"tre"), true);
        tip_labels = tiplabels(tree);
        
        ## Compute median tip rates ##
        df = DataFrame(
            "tip_label"  => tip_labels,
            "birth_rate" => [mapslices(median, tip_birth_rates, dims=2)...],
            "death_rate" => [mapslices(median, tip_death_rates, dims=2)...]
        )

        # Save
        out_csv = replace(file, ".txt" => "_tiprates.csv")
        CSV.write(out_csv, df)

        ## Compute average rates through time ##
        δt = 0.2
        average_b_rates = aggr_rates.(out_trees, b, δt, af = mean)
        average_d_rates = aggr_rates.(out_trees, d, δt, af = mean)

        # Merge across samples
        b_rates_matrix = reduce(hcat, [t[2] for t in average_b_rates])
        d_rates_matrix = reduce(hcat, [t[2] for t in average_d_rates])
        # Add tip rates
        b_rates_matrix = vcat(b_rates_matrix, mapslices(mean, tip_birth_rates, dims=1))
        d_rates_matrix = vcat(d_rates_matrix, mapslices(mean, tip_death_rates, dims=1))
        times = vcat(average_b_rates[1][1], 0)
        rates_df = hcat(DataFrame(time = times), 
                        DataFrame(b_rates_matrix, ["b_$i" for i in 1:size(b_rates_matrix, 2)]), 
                        DataFrame(d_rates_matrix, ["d_$i" for i in 1:size(d_rates_matrix, 2)]))

        # Save
        CSV.write(replace(file, ".txt" => "_averageRatesThroughTime.csv"), rates_df)

        ## Compute lineages‑through‑time (LTT) ##
        ltt_curves = ltt.(out_trees)
        for c in ltt_curves
            popfirst!(c.t)
            popfirst!(c.n)
        end

        times_ltt  = round.(average_b_rates[1][1], digits=1)

        # For one curve, pick the lineage count valid at each grid time.
        counts_at_grid(curve, grid) = [curve.n[searchsortedfirst(round.(curve.t, digits=7), τ, rev=true)] for τ in grid]

        ltt_matrix  = reduce(hcat, counts_at_grid.(ltt_curves, Ref(times_ltt)))


        ltt_df = hcat(
            DataFrame(time = times_ltt),
            DataFrame(ltt_matrix, ["ltt_$i" for i in 1:size(ltt_matrix, 2)]),
        )
        CSV.write(replace(file, ".txt" => "_LTT.csv"), ltt_df)

    end
end
