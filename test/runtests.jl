using MatchRBDs
using Test
using BioStructures
using Rotations
using Dates

const querydir = joinpath("..", "pdb", "queries")
const targetdir = joinpath("..", "pdb", "targets")
const answerdir = joinpath("..", "pdb", "answers")

@testset verbose = true "MatchRBDs.jl" begin
    # Write your tests here.
    @testset "Alignment filter" begin
        struc = read(joinpath(querydir, "7MMO_ba1.pdb"), PDB, structure_name = "7MMO_ba1")
        @test (loadentirepdb(); true) skip=true
        @test (alignable_targets_in_pdb(ProcessedComplex(struc, "C")); true) skip=true
    end
    
    @testset "Target instances"  begin
        struc = read(joinpath(querydir, "7MMO.pdb"), PDB, structure_name = "7MMO")
    
        # make all chains be identical
        for (chid, ch) in chains(struc)
            struc[chid] = struc["C"]
            struc[chid].id = chid
        end
    
        q = ProcessedComplex(struc, "C")
        t = UnprocessedComplex(struc)

        # should give all possible permutations of the chains. n! / (n - k)!, where n is the number of chains in the target and k is the number of chains in the query
        @test length(MatchRBDs.alignabletargetinstances(t, q)) == factorial(length(t.chains)) / factorial(length(t.chains) - length(MatchRBDs.chains(q))) skip=true

    end
    
    @testset "BioStructures transformations" begin
        coords1 = randn(3, 10) .* 100
        coords2 = rand(QuatRotation) * coords1 .+ randn(3) .* 1000
    
        xform = Transformation(coords1, coords2)
    
        @test applytransform(coords1, xform) ≈ coords2
        @test !(applytransform(coords2, xform) ≈ coords1)
    end
    
    @testset "Epitope intersection" begin
        pth = joinpath(querydir, "7MMO_ba1.pdb")
        struc = read(pth, PDB, structure_name = "7MMO_ba1")
        q = ProcessedComplex(struc, "C")

        t = deepcopy(q)
    
        @test epitope_overlap_ratio(q, t) == 1.0
    
        # This test works, but it's not a good test.
        # TODO: set up a test that compares two different complexes
        
        t.specialchain.coordinates .+= 1e9
        
        @test epitope_overlap_ratio(q, t) == 0.0
    end

    @testset "superimpose!" begin
        pth = joinpath(querydir, "7MMO_ba1.pdb")
        struc = read(pth, PDB, structure_name = "7MMO_ba1")

        rot = rand(QuatRotation)
        tr = randn(3) .* 1000

        query = ProcessedComplex(struc, "C")
        target = deepcopy(query)

        target.specialchain.coordinates .= rot * target.specialchain.coordinates .+ tr
        for i in eachindex(target.otherchains)
            target.otherchains[i].coordinates .= rot * target.otherchains[i].coordinates .+ tr
        end
        rotatedtargetcopy = deepcopy(target)

        @test !(target.specialchain.coordinates ≈ query.specialchain.coordinates)
        @test !(all(target.otherchains[i].coordinates ≈ query.otherchains[i].coordinates for i in eachindex(target.otherchains)))

        MatchRBDs.superimpose!(target, query)

        @test target.specialchain.coordinates ≈ query.specialchain.coordinates
        @test all(target.otherchains[i].coordinates ≈ query.otherchains[i].coordinates for i in eachindex(target.otherchains))
    end

    @testset "RMSD" begin
        pth = joinpath(querydir, "7MMO_ba1.pdb")
        outputfile = joinpath(answerdir, "7MMO_ba1_answer_"*string(Dates.format(now(), "yyyymmddTHHMMSSsss"))*".pdb")

        struc = read(pth, PDB, structure_name = "7MMO_ba1")
        q = ProcessedComplex(struc, "C")
        
        @testset "Simple RMSD tests" begin
            t = deepcopy(q)

            t.specialchain.coordinates .+= [1e6, 0, 0]

            MatchRBDs.superimpose!(t, q)

            @test MatchRBDs.rmsd(q, q) ≈ 0.0

            @test MatchRBDs.rmsd(t, q) ≈ 1e6
        end

        @testset "I/O for `findbestmatch`" begin
            # export answer
            ans = findbestmatch(pth, "C", otherchains_minmatch = 0.0, outputfile = outputfile)

            # import answer
            biostruc = read(joinpath(MatchRBDs.pdb_dir, "$(ans.targetinstance.name).pdb"), PDB, structure_name = ans.targetinstance.name)

            #MatchRBDs.updatestruc!(biostruc, ans.targetinstance, ans.appliedtransformation)
            
            chidsinanswer = Set{String}( [ch.name for ch in MatchRBDs.chains(ans.targetinstance)] )

            for chid in keys(BioStructures.chains(biostruc))
                if chid ∉ chidsinanswer
                    pop!(BioStructures.chains(biostruc), chid)
                end
            end
            
            applytransform!(biostruc, ans.appliedtransformation)

            t = ProcessedComplex(biostruc, ans.targetinstance.specialchain.name)
            
            sort!(t.otherchains, by = ch -> findfirst(x -> ch.name == x.name, ans.targetinstance.otherchains))

            @assert all(t.otherchains[i].name == ans.targetinstance.otherchains[i].name for i in eachindex(t.otherchains, ans.targetinstance.otherchains)) "Chains in answer and targetinstance are not in the same order"
        
            # check that the two answers are the same
            @test MatchRBDs.rmsd(q, ans.targetinstance) ≈ MatchRBDs.rmsd(q, t)

            @test MatchRBDs.rmsd(q, ans.targetinstance) ≈ ans.rmsd
        end
    end
end
