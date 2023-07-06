findmatches(struc::BioStructures.ProteinStructure, chainidofspecialchain::String; specialchain_minmatch::Real = 0.3, otherchains_minmatch::Real = 0.3, maxrmsd = 50.0) = 
    threshold_rmsd(ProcessedComplex(struc, chainidofspecialchain), specialchain_minmatch = specialchain_minmatch, otherchains_minmatch = otherchains_minmatch, maxrmsd = maxrmsd)

findmatches(pathtopdbfile::String, chainidofspecialchain::String; specialchain_minmatch::Real = 0.3, otherchains_minmatch::Real = 0.3, maxrmsd = 50.0) = 
    findmatches(BioStructures.read(pathtopdbfile, BioStructures.PDB, structure_name = split(splitpath(pathtopdbfile)[end], ".")[1]), chainidofspecialchain, specialchain_minmatch = specialchain_minmatch, otherchains_minmatch = otherchains_minmatch, maxrmsd = maxrmsd)

function findbestmatch(struc::BioStructures.ProteinStructure, chainidofspecialchain::String; specialchain_minmatch::Real = 0.3, otherchains_minmatch::Real = 0.3, maxrmsd::Real = Inf, outputfile::Union{String, Nothing} = nothing)
    answer = minimum_rmsd(ProcessedComplex(struc, chainidofspecialchain), specialchain_minmatch = specialchain_minmatch, otherchains_minmatch = otherchains_minmatch, maxrmsd = maxrmsd)
    
    isnothing(outputfile) || open(io->writermsdanswer(io, answer), outputfile, "w")

    return answer
end

findbestmatch(pathtopdbfile::String, chainidofspecialchain::String; specialchain_minmatch::Real = 0.3, otherchains_minmatch::Real = 0.3, maxrmsd = Inf, outputfile = nothing) =
    findbestmatch(BioStructures.read(pathtopdbfile, BioStructures.PDB, structure_name = first(split(splitpath(pathtopdbfile)[end], "."))), chainidofspecialchain, specialchain_minmatch = specialchain_minmatch, otherchains_minmatch = otherchains_minmatch, maxrmsd = maxrmsd, outputfile = outputfile)
