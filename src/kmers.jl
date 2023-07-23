function kmers(seq::String, k::Integer)
    kmerset = Set{String}()

    first_index_range = 1 : 1 + length(seq) - k
    last_index_range = k : length(seq) 

    @inbounds for (firstindex, lastindex) in zip(first_index_range, last_index_range)
        push!(kmerset, seq[firstindex:lastindex])
    end

    return kmerset
end

function hasanyof(seq::String, kmerset::Set{String})
    k = length(first(kmerset))

    first_index_range = 1 : 1 + length(seq) - k
    last_index_range = k : length(seq) 

    @inbounds for (firstindex, lastindex) in zip(first_index_range, last_index_range)
        view(seq, firstindex:lastindex) ∈ kmerset && return true
    end

    return false
end