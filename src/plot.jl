# Haven't gotten solid/dashed lines to work yet
function findsequentialresidues(biochain::BioStructures.Chain, inds)
    residues = BioStructures.collectresidues(biochain, MatchRBDs.bioselector)

    return Bool[BioStructures.sequentialresidues(residues[fr], residues[to]) for (fr, to) in zip(inds[1:end-1], inds[2:end])]
end

function plottingdata(biochain::BioStructures.Chain, chain, inds)
    issequential = findsequentialresidues(biochain, inds)

    xs = chain.coordinates[1, inds]
    ys = chain.coordinates[2, inds]
    zs = chain.coordinates[3, inds]

    return xs, ys, zs, issequential
end

function plotchain!(pl, biochain, chain, args...; inds = collect(eachindex(chain.sequence)), complexname = nothing, kwargs...)
    xs, ys, zs, issequential = plottingdata(biochain, chain, inds)

    label = isnothing(complexname) ? chain.name : "$(chain.name) of $complexname"
    plot!(pl, xs, ys, zs, linestyle = :solid, label = label; kwargs...)
end

plotchain!(pl, biochain, chain, matchedchain::ProteinChain, args...; kwargs...) = plotchain!(pl, biochain, chain, args...; inds = first(matchedindices(chain, matchedchain)), kwargs...)

plotchain(args...; kwargs...) = plotchain!(plot(), args...; kwargs...)

function plotcomplex!(pl, biostruc, complex, matchedcomplex = nothing, args...; trimtomatchedcomplex = false, kwargs...)
    for (i, ch) in enumerate(chains(complex))
        bioch = biostruc[ch.name]

        inds = ch == complex.specialchain && trimtomatchedcomplex ? first(matchedindices(ch, matchedcomplex.specialchain)) : collect(eachindex(ch.sequence))

        if isnothing(matchedcomplex)
            plotchain!(pl, bioch, ch, args...; inds = inds, complexname = complex.name, kwargs...)
        else
            plotchain!(pl, bioch, ch, chains(matchedcomplex)[i], args...; inds = inds, complexname = complex.name, kwargs...)
        end
    end
    pl
end

plotcomplex(args...; kwargs...) = plotcomplex!(plot(), args...; kwargs...)

function plotcomplexes!(pl, complextuples, args...; kwargs...)
    for complextuple in complextuples
        plotcomplex!(pl, complextuple..., args...; kwargs...)
    end
    pl
end

plotcomplexes(complextuples, args...; kwargs...) = plotcomplexes!(plot(), complextuples, args...; kwargs...)

function plotmatch!(pl, q_struc, query, t_struc, target, args...; trimtoquery = true, kwargs...)
    plotcomplex!(pl, q_struc, query, target, trimtomatchedcomplex = false, args...; kwargs...)
    plotcomplex!(pl, t_struc, target, query, trimtomatchedcomplex = trimtoquery, args...; kwargs...)
end

plotmatch(args...; kwargs...) = plotmatch!(plot(), args...; kwargs...)

function plotepitope!(pl, biochain::BioStructures.Chain, complex::ProcessedComplex)
    epitopeinds = epitope_indices(complex)
    
    nonepitopeinds = setdiff(eachindex(complex.specialchain.sequence), epitopeinds)
    #@show epitopeinds nonepitopeinds
    plotchain!(pl, biochain, complex.specialchain; inds = epitopeinds, alpha = 1.0)
    plotchain!(pl, biochain, complex.specialchain; inds = nonepitopeinds, alpha = 1.0)
end

plotepitope!(pl, biostruc, complex) = plotepitope!(pl, biostruc[complex.specialchain.name], complex)