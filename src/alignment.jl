
# FIXME: Try choosing the min subset. Perhaps that will be better.
function yeet_edges(short_seq::String, long_seq::String)
    left_edge = 1 + (length(short_seq) - length(lstrip(short_seq, '-')))
    right_edge = length(rstrip(short_seq, '-'))

    return short_seq[left_edge:right_edge], long_seq[left_edge:right_edge]
end

function yeet_all_edges(seq1, seq2) # currently not used
    seq1, seq2 = yeet_edges(seq1, seq2)
    seq2, seq1 = yeet_edges(seq2, seq1)

    return seq1, seq2
end

struct Alignment
    originalquery::String
    originaltarget::String
    alignedquery::String
    alignedtarget::String    
end

Alignment(; query::String, target::String) = Alignment(query, target, yeet_all_edges(NextGenSeqUtils.affine_nw_align(query, target)...)...)

nummatches(alignment::Alignment) = count(x -> x[1] == x[2], zip(alignment.alignedquery, alignment.alignedtarget))

alignmentscore(alignment::Alignment) = nummatches(alignment) / length(alignment.alignedquery)

alignmentscore(target::ProteinChain, query::ProteinChain) = alignmentscore(Alignment(query = query.sequence, target = target.sequence))

isalignable(alignment::Alignment; minalignmentscore = 0.3, minalignmentsizefactor = 0.45) = alignmentscore(alignment) >= minalignmentscore && length(alignment.alignedquery) > length(alignment.originalquery) * minalignmentsizefactor

isalignable(target::ProteinChain, query::ProteinChain; minalignmentscore = 0.3) = isalignable(Alignment(query = query.sequence, target = target.sequence); minalignmentscore = minalignmentscore)

# This method is sketchy. It doesn't check alignability for the entire complex, just for the special chain
isalignable(target::UnprocessedComplex, query::ProcessedComplex; minalignmentscore = 0.3) = any(isalignable(target_chain, query.specialchain; minalignmentscore = minalignmentscore) for target_chain in target.chains)

alignable_targets_in_pdb(query::ProcessedComplex; specialchain_minmatch = 0.3) = filter(target -> isalignable(target, query; minalignmentscore = specialchain_minmatch), loadentirepdb())

function alignabletargetinstances(choseninds::Vector{Int}, specialchain::ProteinChain, target::UnprocessedComplex, query::ProcessedComplex; otherchains_minmatch::Real = 0.3)
    length(choseninds) == length(query.otherchains) && return [deepcopy(ProcessedComplex(target.name, specialchain, target.chains[choseninds]))]

    targetinstances = ProcessedComplex[]
    for (i, targetchain) in enumerate(target.chains)
        i ∉ choseninds || continue # could optimize this by using a Set instead of a Vector
        targetchain != specialchain || continue
        isalignable(targetchain, query.otherchains[length(choseninds) + 1], minalignmentscore = otherchains_minmatch) || continue

        append!(targetinstances, alignabletargetinstances(vcat(choseninds, [i]), specialchain, target, query, otherchains_minmatch = otherchains_minmatch)) # Could optimize this by making a tree instead of using vcat, but that's a lot of work
    end
    return targetinstances
end

function alignabletargetinstances(target::UnprocessedComplex, query::ProcessedComplex; specialchain_minmatch::Real = 0.3, otherchains_minmatch::Real = 0.3) # TODO
    targetinstances = ProcessedComplex[]
    for specialchain in target.chains
        isalignable(specialchain, query.specialchain, minalignmentscore = specialchain_minmatch) || continue

        append!(targetinstances, alignabletargetinstances(Int[], specialchain, target, query, otherchains_minmatch = otherchains_minmatch))
    end
    return targetinstances
end

alignabletargetinstances(targets::Vector{UnprocessedComplex}, query::ProcessedComplex; specialchain_minmatch::Real = 0.3, otherchains_minmatch::Real = 0.3) = vcat(Vector{ProcessedComplex}[alignabletargetinstances(target, query, specialchain_minmatch = specialchain_minmatch, otherchains_minmatch = otherchains_minmatch) for target in targets]...)

function matchedindices(queryseq::String, targetseq::String)
    query, target = NextGenSeqUtils.affine_nw_align(queryseq, targetseq) # We don't want to use the yeet_edges function here because we want to index into the original sequences, not the aligned ones

    queryindices = Int[]
    targetindices = Int[]

    query_index = 1
    target_index = 1

    for (q, t) in zip(query, target)
        if q == t
            push!(queryindices, query_index)
            push!(targetindices, target_index)
        end

        query_index += q != '-'
        target_index += t != '-'
    end

    return queryindices, targetindices
end

matchedindices(query::ProteinChain, target::ProteinChain) = matchedindices(query.sequence, target.sequence)

function matchedindices(query::Vector{ProteinChain}, target::Vector{ProteinChain})
    queryindices = Int[]
    targetindices = Int[]

    q_offset = 0
    t_offset = 0
    for (querychain, targetchain) in zip(query, target)
        queryindices_, targetindices_ = matchedindices(querychain, targetchain)
        append!(queryindices, queryindices_ .+ q_offset)
        append!(targetindices, targetindices_ .+ t_offset)

        q_offset += length(querychain.sequence)
        t_offset += length(targetchain.sequence)
    end

    return queryindices, targetindices
end

# We are not just concatenating the sequences and then calling matchedindices, because we want to make sure to align the chains in the query to the their corresponding chains in the target
# Perhaps we shouldn't call it query and target since it doesn't really matter which is which
function matchedcoords(query::ProteinChain, target::ProteinChain)
    queryindices, targetindices = matchedindices(query, target)
    return query.coordinates[:, queryindices], target.coordinates[:, targetindices]
end

# Perhaps we shouldn't call it query and target since it doesn't really matter which is which
function matchedcoords(query::Vector{ProteinChain}, target::Vector{ProteinChain})
    queryindices, targetindices = matchedindices(query, target)

    querycoords = joincoords(query)[:, queryindices]
    targetcoords = joincoords(target)[:, targetindices]

    return querycoords, targetcoords
end