function kmers(seq::String, k::Integer)
    kmerset = Set{String}()

    first_index_range = 1 : 1 + length(seq) - k
    last_index_range = k : length(seq) 

    for (firstindex, lastindex) in zip(first_index_range, last_index_range)
        push!(kmerset, seq[firstindex:lastindex])
    end

    return kmerset
end