function findmatches(query::ComplexInstance; num_rand_kmers = 20, k = 6)
    @time "I/O" begin
        df = DataFrames.DataFrame(Arrow.Table(path_to_dataset))
        BAs_in_pdb = load_proteincomplexes(df)
        @info "Number of BAs in dataset: " length(BAs_in_pdb)
    end

    @time "BA chaincount filter" begin
        filter!(BAs_in_pdb) do BA
            length(BA.chains) >= length(query.proteincomplex.chains)
        end
        @info "Number of BAs after size filter: " length(BAs_in_pdb)
    end

    @time "Kmer filter" begin
        kmerset = kmers(specialchain(query).sequence, k)

        rand_kmers = Set([pop!(kmerset, rand(kmerset)) for _ in 1:min(num_rand_kmers, length(kmerset))])

        filter!(BAs_in_pdb) do BA
            any(BA.chains) do ch
                any(kmer ∈ rand_kmers for kmer in kmers(ch.sequence, k))
            end
        end
        @info "Number of BAs after kmer filter: " length(BAs_in_pdb)
    end

    @time "Alignment" begin
        matchedtargetinstances = alignabletargetinstances(query, BAs_in_pdb)
        @info "Number of possible matches: " length(matchedtargetinstances)
    end

    @time "Superimposition and rmsd calculation" begin
        matches = [(rmsd = rmsd(x), match = x) for x in matchedtargetinstances]
    end

    @time "Result sorting" begin
        sort!(matches, by = x -> x.rmsd)
    end

    global num_alignment_queries
    @info "Number of alignment_queries" num_alignment_queries

    return matches
end