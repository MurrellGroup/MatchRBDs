num_alignment_queries = 0

function alignabletargetinstances!(chainmatches, chaininds, query, target, minimum_otherchain_alignmentscore)
    # Base case
    if length(chaininds) == length(query.proteincomplex.chains)
        targetinstance = ComplexInstance(target, chaininds[1])
        complexinstancematch = ComplexInstanceMatch(query, targetinstance, chainmatches[(specialchain(query).id, specialchain(targetinstance).id)].chainmatch, Vector{ChainMatch}())

        for (querychain, targetchain) in zip(otherchains(query), target.chains[chaininds[2:end]])
            push!(complexinstancematch.otherchainmatches, chainmatches[(querychain.id, targetchain.id)].chainmatch)
        end
        return [complexinstancematch]
    end
        
    # DFS
    complexinstancematches = Vector{ComplexInstanceMatch}()
    for i in setdiff(eachindex(target.chains), chaininds)

        targetchain = target.chains[i]
        querychain = otherchains(query)[length(chaininds)]

        key = (querychain.id, targetchain.id)
        val = get!(chainmatches, key, align(querychain, targetchain, minimum_otherchain_alignmentscore))
        
        val.isalignable || continue

        append!(complexinstancematches, alignabletargetinstances!(chainmatches, vcat(chaininds, i), query, target, minimum_otherchain_alignmentscore))
    end
    return complexinstancematches
end

function alignabletargetinstances(query::ComplexInstance, target::ProteinComplex; minimum_specialchain_alignmentscore = 0.3, minimum_otherchain_alignmentscore = 0.3)
    complexinstancematches = Vector{ComplexInstanceMatch}()

    chainmatches = Dict{Tuple{String, String}, NamedTuple{(:isalignable, :chainmatch), Tuple{Bool, ChainMatch}}}()
    
    for (specialchain_index, targetspecialchain) in enumerate(target.chains)
        key = (specialchain(query).id, targetspecialchain.id)
        val = get!(chainmatches, key, align(specialchain(query), targetspecialchain, minimum_specialchain_alignmentscore))
        
        val.isalignable || continue

        append!(complexinstancematches, alignabletargetinstances!(chainmatches, [specialchain_index], query, target, minimum_otherchain_alignmentscore))
    end

    # This is ugly. Shouldn't be here
    #best_complexinstancematch = argmin(rmsd, complexinstancematches)

    @show collect(keys(chainmatches))
    #return [best_complexinstancematch]
    global num_alignment_queries += length(chainmatches)

    return complexinstancematches
end


alignabletargetinstances(query::ComplexInstance, targets::Vector{ProteinComplex}; minimum_specialchain_alignmentscore = 0.3, minimum_otherchain_alignmentscore = 0.3) = vcat(Vector{ComplexInstanceMatch}[alignabletargetinstances(query, target, minimum_specialchain_alignmentscore = minimum_specialchain_alignmentscore, minimum_otherchain_alignmentscore = minimum_otherchain_alignmentscore) for target in targets]...)

function matchedindices(aligned_queryseq::String, aligned_targetseq::String)
    queryindices = Int[]
    targetindices = Int[]

    query_index = 1
    target_index = 1

    for (q, t) in zip(aligned_queryseq, aligned_targetseq)

        q == t && (push!(queryindices, query_index); push!(targetindices, target_index))

        query_index += q != '-'
        target_index += t != '-'
    end

    return queryindices, targetindices
end

function alignmentscore(aligned_query::String, aligned_target::String)

    num_matches = count(i -> aligned_query[i] == aligned_target[i], eachindex(aligned_query))
    
    firstnongap = maximum((seq -> findfirst(!=('-'), seq)).([aligned_query, aligned_target]))
    lastnongap = minimum((seq -> findprev(!=('-'), seq, lastindex(seq))).([aligned_query, aligned_target]))

    alignmentsize_noedgegaps = 1 + (lastnongap - firstnongap)

    return num_matches / alignmentsize_noedgegaps
end

function align(query::ProteinChain, target::ProteinChain, minalignmentscore::Float64)
    aligned_query, aligned_target = NextGenSeqUtils_affine_nw_align(query.sequence, target.sequence)

    query_indices, target_indices = matchedindices(aligned_query, aligned_target)

    return (isalignable = alignmentscore(aligned_query, aligned_target) >= minalignmentscore, chainmatch = ChainMatch(query, target, query_indices, target_indices))
end

function NextGenSeqUtils_affine_nw_align(s1::String, s2::String;
    gap_open = -2.0,
    gap_extend = -0.2,
    match_cost = 1.0,
    mismatch_cost = -1.0)
   
    s1arr = collect(s1)  # vertical
    s1len = length(s1arr)
    s2arr = collect(s2)  # horizontal
    s2len = length(s2arr)
  
    M = zeros(s1len+1, s2len+1)
    IX = zeros(s1len+1, s2len+1)
    IY = zeros(s1len+1, s2len+1)

    IX[:,1] .= gap_open .+ gap_extend * (0:s1len)
    IX[1,:] .= -Inf * [1:s2len+1;]
    IY[1,:] .= gap_open .+ gap_extend * (0:s2len)
    IY[:,1] .= -Inf * [1:s1len+1;]
    M[1:end,1] .= gap_open .+ gap_extend * (0:s1len) #Not sure about this.
    M[1,1:end] .= gap_open .+ gap_extend * (0:s2len) #Not sure about this.
    
    traceM = zeros(Int, s1len+1, s2len+1)
    traceIX = zeros(Int, s1len+1, s2len+1)
    traceIY = zeros(Int, s1len+1, s2len+1)
    
    #ALSO UNSURE
    traceM[1,:] .+= 3 
    traceM[2:end,1] .+= 2
    
    traceIX[1,:] .+= 3
    traceIY[2:end,1] .+= 2
    
    traceIY[1,:] .+= 3 
    traceIX[2:end,1] .+= 2
    
    for i in 2:s1len+1
        for j in 2:s2len+1
            diag_cost = 0.0
            if s1arr[i-1] == s2arr[j-1]
                diag_cost = match_cost
            else
                diag_cost = mismatch_cost
            end
            diagM = M[i-1, j-1] + diag_cost
            IX2M = IX[i-1, j-1] + diag_cost
            IY2M = IY[i-1, j-1] + diag_cost
            
            M[i,j],traceM[i,j] = findmax((diagM,IX2M,IY2M)) #1: stay in M. 2: come from IX. 3: come from IY
            
            M2IX = M[i-1, j] + gap_open
            IXextend = IX[i-1, j] + gap_extend
            IX[i,j],traceIX[i,j] = findmax((M2IX,IXextend)) #1: gap open. 2: gap extend
            
            M2IY = M[i, j-1] + gap_open
            IYextend = IY[i, j-1] + gap_extend
            IY[i,j],traceIY[i,j] = findmax((M2IY,-Inf,IYextend)) #1: gap open. 3: gap extend
              
        end
    end
    
    rev_arr1 = Char[]
    rev_arr2 = Char[]
    sizehint!(rev_arr1, s1len+s2len)
    sizehint!(rev_arr2, s1len+s2len)
    
    mats = [traceM,traceIX,traceIY]
    x_i = s1len+1
    y_i = s2len+1
    m_i = argmax([M[x_i,y_i],IX[x_i,y_i],IY[x_i,y_i]])
    while x_i-1 > 0 &&  y_i-1 > 0
        next_m_i = mats[m_i][x_i,y_i]
        if m_i == 1
            push!(rev_arr1,s1arr[x_i-1])
            push!(rev_arr2,s2arr[y_i-1])
            
            x_i -= 1
            y_i -= 1
        elseif m_i == 2
            push!(rev_arr1,s1arr[x_i-1])
            push!(rev_arr2,'-')
            x_i -= 1
        elseif m_i == 3
            push!(rev_arr1,'-')
            push!(rev_arr2,s2arr[y_i-1])
            y_i -= 1
        end
        m_i = next_m_i
    end
    while x_i-1 > 0
        push!(rev_arr1,s1arr[x_i-1])
        push!(rev_arr2,'-')
        x_i -= 1
    end
    while y_i-1 > 0
        push!(rev_arr1,'-')
        push!(rev_arr2,s2arr[y_i-1])
        y_i -= 1
    end
    return join(reverse(rev_arr1)),join(reverse(rev_arr2))
end