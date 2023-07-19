struct ProteinChain
    id::String
    sequence::String
    coords::Matrix{Float64}
end

struct ProteinComplex
    name::String
    chains::Vector{ProteinChain}
end

struct ComplexInstance
    proteincomplex::ProteinComplex # shallow copy
    specialchain_index::Int
end

specialchain(instance::ComplexInstance) = instance.proteincomplex.chains[instance.specialchain_index]

otherchains(instance::ComplexInstance) = ProteinChain[chain for (i, chain) in enumerate(instance.proteincomplex.chains) if i != instance.specialchain_index]

# Matches

struct ChainMatch
    query::ProteinChain # shallow copy
    target::ProteinChain # shallow copy
    queryindices::Vector{Int}
    targetindices::Vector{Int}
end

struct ComplexInstanceMatch
    query::ComplexInstance
    target::ComplexInstance
    specialchainmatch::ChainMatch
    otherchainmatches::Vector{ChainMatch}
end

matchedcoords_query(chainmatch::ChainMatch) = chainmatch.query.coords[:, chainmatch.queryindices]
matchedcoords_target(chainmatch::ChainMatch) = chainmatch.target.coords[:, chainmatch.targetindices]

matchedcoords_query_otherchains(instancematch::ComplexInstanceMatch) = hcat([matchedcoords_query(matchedchain) for matchedchain in instancematch.otherchainmatches]...)
matchedcoords_target_otherchains(instancematch::ComplexInstanceMatch) = hcat([matchedcoords_target(matchedchain) for matchedchain in instancematch.otherchainmatches]...)

function rmsd(instancematch::ComplexInstanceMatch)
    
    transformation = BioStructures.Transformation(matchedcoords_target(instancematch.specialchainmatch), matchedcoords_query(instancematch.specialchainmatch))

    query_coords_otherchains = matchedcoords_query_otherchains(instancematch)
    transformed_target_coords_otherchains = BioStructures.applytransform(matchedcoords_target_otherchains(instancematch), transformation)

    return BioStructures.rmsd(query_coords_otherchains, transformed_target_coords_otherchains)
end

using Base: show

function Base.show(io::IO, x::ComplexInstanceMatch)
    print(io, "ComplexMatch( ")
    print(io, "target name: ", x.target.proteincomplex.name, ", ")
    print(io, "specialchain: ", x.specialchainmatch.target.id, ", ")
    print(io, "otherchains: ", [" "*m.target.id for m in x.otherchainmatches]..., " )")
end