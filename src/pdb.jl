# TODO: Function to write answer to pdb file

const pdb_dir = "C:\\Users\\arons\\.julia\\dev\\MatchRBDs\\pdb\\targets"

function loadentirepdb()
    pdb_files = readdir(pdb_dir)

    targets = UnprocessedComplex[]
    for f in pdb_files
        name = split(f, ".")[1]
        biostruc = BioStructures.read(joinpath(pdb_dir, f), BioStructures.PDB, structure_name = name)::BioStructures.ProteinStructure
        struc = UnprocessedComplex(biostruc)
        push!(targets, struc)
    end
    return targets
end

# Assert that multiple chains don't have the same name. This is a requirement for the rest of the code to work.
function writermsdanswer(io, answer)
    biostruc = BioStructures.read(joinpath(pdb_dir, "$(answer.targetinstance.name).pdb"), BioStructures.PDB, structure_name = answer.targetinstance.name)

    updatestruc!(biostruc, answer.targetinstance, answer.appliedtransformation)

    BioStructures.writepdb(io, biostruc)
end

function updatechain!(biochain::BioStructures.Chain, chain)
    for at in BioStructures.collectatoms(biochain, bioselector)
        ind = findfirst(resid -> resid == BioStructures.resid(at), chain.resids)
        BioStructures.coords!(at, chain.coordinates[:, ind])
    end
end

function updatestruc!(biostruc::BioStructures.ProteinStructure, complex, transformation::BioStructures.Transformation)
    
    @assert length(unique([ch.name for ch in chains(complex)])) == length(chains(complex)) "Multiple chains have the same name. This is not allowed."

    @assert length(unique([chid for chid in keys(BioStructures.chains(biostruc))])) == length(keys(BioStructures.chains(biostruc))) "Multiple chains have the same name. This is not allowed."

    chidsinanswer = Set{String}( [ch.name for ch in chains(complex)] )

    
    for (chid, ch) in BioStructures.chains(biostruc)
        if chid ∈ chidsinanswer
            #updatechain!(ch, only(filter(x -> x.name == chid, chains(complex))))
        else
            pop!(BioStructures.chains(biostruc), chid)
        end
    end
    

    BioStructures.applytransform!(biostruc, transformation)


    return biostruc
end