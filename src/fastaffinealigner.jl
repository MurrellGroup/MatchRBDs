# using StaticArrays

# const GAP, A, C, D, E, F, G = -1

function fastaffinenw(s1::String, s2::String,
    gapopen::Int = -20,
    gapextend::Int = -2,
    matchcost::Int = 10,
    mismatchcost::Int = -10
    )

    inf = Int(1e6)

    s1len = length(s1) # vertical
    s2len = length(s2) # horizontal

    M = zeros(Int, s1len + 1, s2len + 1)
    IX = zeros(Int, s1len + 1, s2len + 1)
    IY = zeros(Int, s1len + 1, s2len + 1)
    
    @. IX[:, 1] = gapopen + gapextend * (0:s1len)
    IX[1, :] .= -inf
    @. IY[1, :] = gapopen + gapextend * (0:s2len)
    IY[:, 1] .= -inf

    @. M[:, 1] = gapopen + gapextend * (0:s1len)
    @. M[1, :] = gapopen + gapextend * (0:s2len)

    @inbounds for c in 2:s2len + 1, r in 2:s1len + 1

        M[r, c] = max(M[r-1, c-1], M[r-1, c-1], M[r-1, c-1]) + ifelse(s1[r-1] == s2[c-1], matchcost, mismatchcost)

        IX[r, c] = max(M[r-1, c] + gapopen, IX[r-1, c] + gapextend)

        IY[r, c] = max(M[r, c-1] + gapopen, IY[r, c-1] + gapextend)

    end

    rev_arr1 = Char[]
    rev_arr2 = Char[]
    sizehint!(rev_arr1, s1len + s2len)
    sizehint!(rev_arr2, s1len + s2len)

    r = s1len + 1
    c = s2len + 1
    cur_dp = argmax((M[r, c], IX[r, c], IY[r, c]))

    @inbounds while r-1 > 0 && c-1 > 0
        if cur_dp == 1
            r -= 1
            c -= 1
            
            push!(rev_arr1, s1[r])
            push!(rev_arr2, s2[c])

            cur_dp = argmax((M[r, c], IX[r, c], IY[r, c]))
        elseif cur_dp == 2
            r -= 1

            push!(rev_arr1, s1[r])
            push!(rev_arr2, '-')

            cur_dp = argmax((M[r, c] + gapopen, IX[r, c] + gapextend))
        elseif cur_dp == 3
            c -= 1
            
            push!(rev_arr1, '-')
            push!(rev_arr2, s2[c])

            cur_dp = argmax((M[r, c] + gapopen, -inf, IY[r, c] + gapextend))
        end
    end

    @inbounds while (r -= 1) > 0
        push!(rev_arr1, s1[r])
        push!(rev_arr2, '-')
    end

    @inbounds while (c -= 1) > 0
        push!(rev_arr1, '-')
        push!(rev_arr2, s2[c])
    end

    reverse!(rev_arr1)
    reverse!(rev_arr2)

    return join(rev_arr1), join(rev_arr2)
end

        








