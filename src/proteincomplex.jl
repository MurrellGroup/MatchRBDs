struct ProteinChain
    name::String
    sequence::String
    coordinates::Matrix{Float64}
    resids::Vector{String}
end

# Do we want to keep disordered atoms or not?
bioselector(at::BioStructures.AbstractAtom) = BioStructures.calphaselector(at) && !BioStructures.isdisorderedatom(at)

function bioselector(res::BioStructures.AbstractResidue)
    ats = BioStructures.collectatoms(res, BioStructures.calphaselector)
    !isempty(ats) && bioselector(only(ats))
end

function ProteinChain(chain::BioStructures.Chain)
    name = BioStructures.chainid(chain)
    sequence = string(BioStructures.LongAA(chain, bioselector, gaps = false))
    coordinates = BioStructures.coordarray(chain, bioselector)
    ids = BioStructures.resid.(BioStructures.collectresidues(chain, bioselector))

    @assert length(sequence) == size(coordinates, 2) "Sequence and coordinates are not the same length, they are $(length(sequence)) and $(size(coordinates, 2)) respectively."

    return ProteinChain(name, sequence, coordinates, ids)
end

abstract type AbstractProteinComplex end

struct UnprocessedComplex
    name::String
    chains::Vector{ProteinChain}
end

UnprocessedComplex(struc::BioStructures.ProteinStructure) = UnprocessedComplex(struc.name, ProteinChain.(values(BioStructures.chains(struc))))

struct ProcessedComplex
    name::String
    specialchain::ProteinChain
    otherchains::Vector{ProteinChain}
end

function ProcessedComplex(struc::BioStructures.ProteinStructure, idofspecialchain::String)
    otherchains = ProteinChain[ProteinChain(ch) for (chid, ch) in BioStructures.chains(struc) if chid != idofspecialchain]
    specialchain = ProteinChain(struc[idofspecialchain])
    return ProcessedComplex(struc.name, specialchain, otherchains)
end

function showids(complex::ProcessedComplex)
    println("Name: $(complex.name)")
    println("Special chain: $(complex.specialchain.name)")
    println("Other chains: $(join([chain.name for chain in complex.otherchains], ", "))")
end

chains(complex::ProcessedComplex) = vcat(complex.specialchain, complex.otherchains...)

joincoords(chains::Vector{ProteinChain}) = hcat(map(chain -> chain.coordinates, chains)...)

coordarray(complex::ProcessedComplex) = joincoords(chains(complex))