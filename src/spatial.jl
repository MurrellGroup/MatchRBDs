rmsd(target::ProteinChain, query::ProteinChain) = BioStructures.rmsd(matchedcoords(query, target)...)

rmsd(target::ProcessedComplex, query::ProcessedComplex) = BioStructures.rmsd(matchedcoords(query.otherchains, target.otherchains)...)

function applytransform!(complex::ProcessedComplex, transform::BioStructures.Transformation)
    complex.specialchain.coordinates .= BioStructures.applytransform(complex.specialchain.coordinates, transform)
    for i in eachindex(complex.otherchains)
        complex.otherchains[i].coordinates .= BioStructures.applytransform(complex.otherchains[i].coordinates, transform)
    end
    return complex
end

function superimpose!(target::ProcessedComplex, query::ProcessedComplex)
    transform = BioStructures.Transformation(matchedcoords(target.specialchain, query.specialchain)...) # FIXME: Check if applying one transformation to a separate set of coordinates works as expected

    applytransform!(target, transform)

    return transform
end

# The answer target_instance will have been superimposed onto the query
function threshold_rmsd(targetinstances::Vector{ProcessedComplex}, query::ProcessedComplex; maxrmsd::T) where T <: Real
    ans = Vector{NamedTuple{(:rmsd, :targetinstance, :appliedtransformation), Tuple{T, ProcessedComplex, BioStructures.Transformation}}}()

    for targetinstance in targetinstances
        transform = superimpose!(targetinstance, query)
        rmsd_ = rmsd(targetinstance, query)
        if rmsd_ < maxrmsd
            push!(ans, (rmsd = rmsd_, targetinstance = targetinstance, appliedtransformation = transform))
        end
    end

    return ans
end

threshold_rmsd(targets::Vector{UnprocessedComplex}, query::ProcessedComplex; specialchain_minmatch::Real = 0.3, otherchains_minmatch::Real = 0.3, maxrmsd::Real) = 
    threshold_rmsd(alignabletargetinstances(targets, query, specialchain_minmatch = specialchain_minmatch, otherchains_minmatch = otherchains_minmatch), query, maxrmsd = maxrmsd)

threshold_rmsd(query::ProcessedComplex; specialchain_minmatch::Real = 0.3, otherchains_minmatch::Real = 0.3, maxrmsd::Real) = 
    threshold_rmsd(alignable_targets_in_pdb(query, specialchain_minmatch = specialchain_minmatch), query, specialchain_minmatch = specialchain_minmatch, otherchains_minmatch = otherchains_minmatch, maxrmsd = maxrmsd)

# mainly this function will be used
# should have separate minmatch for specialchain and otherchains
# could superimpose otherchain and then filter by rmsd of otherchains
function minimum_rmsd(args...; specialchain_minmatch = 0.3, otherchains_minmatch = 0.3, maxrmsd = Inf)
    answers = threshold_rmsd(args..., specialchain_minmatch = specialchain_minmatch, otherchains_minmatch = otherchains_minmatch, maxrmsd = maxrmsd)
    isempty(answers) ? nothing : argmin(x -> x.rmsd, answers)
end

# Perhaps this could be optimized a lot by using a k-d tree or something. Currently it's O(n * m), where n is the number of residues in the query and m is the number of residues in the target.
# We pass specialchaincoords because it can be multiple different things or have been filtered in some way.
epitope_indices(specialchaincoords::Matrix{Float64}, complex::ProcessedComplex; epitope_distance_threshold::Real = 10.0) = findall( eachcol(specialchaincoords) ) do p1
    any( p2 -> euclidean(p1, p2) < epitope_distance_threshold, eachcol(joincoords(complex.otherchains)) )
end

epitope_indices(complex::ProcessedComplex; epitope_distance_threshold::Real = 10.0) = epitope_indices(complex.specialchain.coordinates, complex; epitope_distance_threshold = epitope_distance_threshold)

# TODO: Figure out what to do when there are two antibodies on a single antigen.
# TODO: Use matchedcoords. Already implemented that. Perhaps use it in epitope_indices, meaning that you would need a new implementation for joincoords.
function matched_epitope_indices(query::ProcessedComplex, target::ProcessedComplex)
    query_matchedspecialchaincoords, target_matchedspecialchaincoords = matchedcoords(query.specialchain, target.specialchain)
    return epitope_indices(query_matchedspecialchaincoords, query), epitope_indices(target_matchedspecialchaincoords, target)
end

function epitope_overlap_ratio(query::ProcessedComplex, target::ProcessedComplex)
    queryepitopeindices, targetepitopeindices = matched_epitope_indices(query, target)

    isempty(queryepitopeindices) && isempty(targetepitopeindices) && (println("Warning: both complexes have empty epitope indices");  return 0.0)

    return length(intersect(queryepitopeindices, targetepitopeindices)) / length(union(queryepitopeindices, targetepitopeindices))
end