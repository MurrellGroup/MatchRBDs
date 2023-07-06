module MatchRBDs
    import BioStructures
    import NextGenSeqUtils
    using Distances: euclidean
    using Plots
    
    include("proteincomplex.jl")
    include("pdb.jl")
    include("alignment.jl")
    include("interface.jl")
    include("spatial.jl")
    include("plot.jl")

    export 
        # types
        ProcessedComplex, 
        UnprocessedComplex,
        ProteinChain,
        # typeutils
        showids,
        # io, pdb
        loadentirepdb,
        writermsdanswer,
        alignable_targets_in_pdb,
        # spatial
        alignabletargetinstances,
        minimum_rmsd,
        rmsd,
        epitope_overlap_ratio,
        # interface
        findmatches,
        findbestmatch,
        # plot 
        plotchain!,
        plotchain,
        plotcomplex!,
        plotcomplex,
        plotcomplexes!,
        plotcomplexes,
        plotmatch!,
        plotmatch,
        plotepitope!,
        plotepitope

    
end