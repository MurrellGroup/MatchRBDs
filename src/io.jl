path_to_dataset = joinpath("pdb", "targets", "BAs.arrow")

path_to_pdb_dir = joinpath("pdb", "targets")

# Temp selector?
# Do we want to keep disordered atoms or not?
bioselector(at::BioStructures.AbstractAtom) = BioStructures.calphaselector(at) && !BioStructures.isdisorderedatom(at)

function bioselector(res::BioStructures.AbstractResidue)
    ats = BioStructures.collectatoms(res, BioStructures.calphaselector)
    !isempty(ats) && bioselector(only(ats))
end

function ProteinChain(chain::BioStructures.Chain)
    name = BioStructures.chainid(chain)
    sequence = string(BioStructures.LongAA(chain, bioselector, gaps = false))

    coordinates = BioStructures.coordarray(chain, bioselector)

    @assert length(sequence) == size(coordinates, 2) "Sequence and coordinates are not the same length, they are $(length(sequence)) and $(size(coordinates, 2)) respectively."

    mask = collect(sequence) .!= '-'

    return ProteinChain(name, join(collect(sequence)[mask]), coordinates[:, mask])
end

ProteinComplex(struc::BioStructures.ProteinStructure) = ProteinComplex(struc.name, ProteinChain[ProteinChain(ch) for ch in values(BioStructures.chains(struc))])

function ComplexInstance(struc::BioStructures.ProteinStructure, specialchain_id::String)
    proteincomplex = ProteinComplex(struc)
    specialchain_index = findfirst(ch -> ch.id == specialchain_id, proteincomplex.chains)
    ComplexInstance(proteincomplex, specialchain_index)
end

function DataFrames.DataFrame(BAs::Vector{BioStructures.ProteinStructure})
    df = DataFrames.DataFrame(:BAid => String[], :chainid => String[], :sequence => String[], :x => Vector{Float64}[], :y => Vector{Float64}[], :z => Vector{Float64}[], :resids => Vector{String}[])
    for ba in BAs
        for (chid, ch) in BioStructures.chains(ba)
            sequence = string(BioStructures.LongAA(ch, bioselector, gaps = false))
            coordinates = BioStructures.coordarray(ch, bioselector)
            resids = BioStructures.resid.(BioStructures.collectresidues(ch, bioselector))

            @assert length(sequence) == size(coordinates, 2) "Sequence and coordinates are not the same length, they are $(length(sequence)) and $(size(coordinates, 2)), respectively."

            DataFrames.push!(df, [ba.name, chid, sequence, coordinates[1, :], coordinates[2, :], coordinates[3, :], resids])
        end
    end
    return df
end

function load_proteincomplex(df::DataFrames.AbstractDataFrame)
    ba_id = only(unique(df[!, :BAid]))

    chains = ProteinChain[]#ProteinChain(ch.chainid, ch.sequence, ) ]
    for ch in eachrow(df)
        coords = zeros(3, length(ch.sequence))
        coords[1, :] .= ch.x
        coords[2, :] .= ch.y
        coords[3, :] .= ch.z

        push!(chains, ProteinChain(ch.chainid, ch.sequence, coords))
    end
    return ProteinComplex(ba_id, chains)
end

load_proteincomplexes(df::DataFrames.DataFrame) = ProteinComplex[load_proteincomplex(g) for g in DataFrames.groupby(df, :BAid)]
    
function generate_dataset(PDB_entries::Vector{Tuple{String, Vector{Int}}}; outfile = path_to_dataset, pdbentrydir = path_to_pdb_dir, overwrite = true)
    df = nothing

    for (pdbid, ba_numbers) in PDB_entries, ba_number in ba_numbers
        try
            BA = BioStructures.retrievepdb(pdbid, ba_number = ba_number, structure_name = "$(pdbid)_$(ba_number)", dir = pdbentrydir, overwrite = overwrite, read_het_atoms = false, remove_disorder = true)
            if isnothing(df)
                df = DataFrames.DataFrame([BA])
            else
                df = vcat(df, DataFrames.DataFrame([BA]))
            end
        catch
            @info "Not able to retrieve $(pdbid)_$(ba_number)"
        end
        GC.gc()
    end
    Arrow.write(outfile, df)
end