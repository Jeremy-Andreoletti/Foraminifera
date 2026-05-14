ENV["GKSwstype"] = "nul"
gr(dpi=400, size=(500,300))

function Forams_plot(ofile, out_trees, tree, ωtimes)
    if isa(tree, AbstractVector)
        reorder!.(out_trees);
    else
        median_tree = iquantile(remove_unsampled.(out_trees), 0.5);
        reorder!.(out_trees);
        reorder!(median_tree);
    end
    # anim_tree = @animate for tree_i in out_trees
    #   plot(tree_i, shownodes=(true, true, true), showda=true, shsizes=[1.0, 1.0, 1.0])
    # end
    # mp4(anim_tree, "Animations/$(ofile)_anim_tree.mp4", fps=5)

    # λmax = ceil(maximum(b.(out_trees))[1]*2)
    # anim_tree = @animate for tree_i in out_trees
    #   plot(tree_i, b, shownodes=(false, false, false), clim=(0,λmax))
    # end
    # mp4(anim_tree, "Animations/$(ofile)_anim_tree_λ.mp4", fps=5)

    if !isa(tree, AbstractVector)
        plotω(out_trees[lastindex(out_trees)], ωtimes, shownodes=(false, false, true), showda=true)
        png("Images/$(ofile)_tree.png")

        plot(median_tree, b, shownodes=(false, false, false), clim=(0,1.5), size = (450, 300))
        png("Images/$(ofile)_tree_λ.png")
        plot(median_tree, d, shownodes=(false, false, false), clim=(0,1.5), size = (450, 300))
        png("Images/$(ofile)_tree_μ.png")
        plot(median_tree, nd, shownodes=(false, false, false), clim=(-0.5,0.5), size = (450, 300))
        png("Images/$(ofile)_tree_nd.png")
    end
    
    plot(out_trees, b, δt, fillcolor="#005AB5", linecolor=:darkblue)
    plot!(out_trees, d, δt, fillcolor="#DC3220", linecolor=:darkred, xlab = "Time before present", ylab = "Speciation and extinction rates")
    png("Images/$(ofile)_λTT_μTT.png")

    plotω(ltt.(out_trees), ωtimes, 0.05)
    if isa(tree, AbstractVector)
        for t in tree
            plot!(ltt(t), scale=:identity)
        end
    else
        plot!(ltt(tree), scale=:identity)
    end
    png("Images/$(ofile)_LTT.png")
end


