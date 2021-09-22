using DelimitedFiles
using LinearAlgebra
using Clustering
using StatsBase
using Distributed
using DataFrames
using CSV
import ClusterManagers
using LinearAlgebra
using Distances
using Random
using Base
using ArgParse

function ParseCommands()

    s = ArgParseSettings()

    @add_arg_table! s begin
        "--mutation_file", "-f"
            help = "File with the matrix count"
            required = true
        "--min_signatures_to_extract", "-s"
            help = "Minimum number of signatures to extract"
            arg_type = Int
            required = true
        "--max_signatures_to_extract", "-m"
            help = "Maximum number of signatures to extract"
            arg_type = Int
            required = true
        "--outpath", "-o"
            help = "Folder where the results will be saved"
            required = true
        "--iterations", "-i"
            help = "Total NMF iterations performed"
            arg_type = Int
            default = 1024
        "--workers", "-n"
            help = "Number of slurm parallel process or local cores"
            arg_type = Int
            default = 1
        "--weak", "-w"
            help = "Whether to remove weak mutation types. If unsure, use default"
            arg_type = Int
            default = 0
        "--slurm"
            help = "Whether to run the pipeline using SLURM"
            action = :store_true

    end
    return parse_args(s)
end

parsed_args = ParseCommands()
numberWorkers = parsed_args["workers"]
println("Adding workers...")
if parsed_args["slurm"]
    ClusterManagers.addprocs_slurm(numberWorkers;exeflags="--project")
else
    addprocs(numberWorkers;exeflags="--project")
end

# get the path where the script is located
@eval @everywhere path_root=@__DIR__

# add all the necessary scripts
include(string(path_root, "/src/removeWeak.jl"))
include(string(path_root,"/src/addWeak.jl"))
include(string(path_root,"/src/filterOutIterations.jl"))
include(string(path_root,"/src/evaluateStability.jl"))
include(string(path_root,"/src/deconvoluteMutCounts.jl"))
include(string(path_root,"/src/extractSignaturesParallelNoArray.jl"))
include(string(path_root,"/src/removeAllSingleSignatures.jl"))
include(string(path_root,"/src/removeOneSig.jl"))
include(string(path_root,"/src/evalSingleSample.jl"))
include(string(path_root,"/src/stabilities.jl"))

# mutation file
mutation_file = parsed_args["mutation_file"]

# number of signatures to extract
min_numberProcessesToExtract = parsed_args["min_signatures_to_extract"]
max_numberProcessesToExtract = parsed_args["max_signatures_to_extract"]

# total number of iterations
totalIterations = parsed_args["iterations"]

# outpath where they will be saved
outpath = parsed_args["outpath"]

# define whether is a Pan whole-genome, so weak types won't be removed
isWeak = parsed_args["weak"]
remove_weak_muttypes=0.001
if isWeak==0
    remove_weak_muttypes=0.00
end

# add seed so that the extraction is reproducible
Random.seed!(123456)

totalReplicates = 100
removeLastPercentage = 0.07

# Load input
genomesOriginal, sample_names = readdlm(mutation_file, '\t', Int, '\n'; header=true)
sample_header = vec(Symbol.(sample_names))

# Deconvolute mutation counts
genomes, sort_idx, n_types_to_remove = removeWeak(genomesOriginal, remove_weak_muttypes)

# start processing
name_ttype = splitext(basename(mutation_file))[1]

@info "Processing $name_ttype"
@info "Outpath $outpath"

if parsed_args["slurm"]
    @info "Slurm mode"
end

# Bootstrap variables
@eval @everywhere parallelGenomes=$genomes
@eval @everywhere originalGenomes=$genomesOriginal
@eval @everywhere parallelTotalIterations=$totalIterations
@eval @everywhere OutPath=$outpath

@eval @everywhere mutfile = $name_ttype
@eval @everywhere initial_mutfile = $mutation_file
@eval @everywhere SavePath = joinpath(OutPath, mutfile)
mkpath(SavePath)

@everywhere parallelGenomesColsSum = sum(parallelGenomes, dims=1)
@everywhere parallelNormGenomes = parallelGenomes ./ repeat(parallelGenomesColsSum, size(parallelGenomes, 1), 1 )

for numberProcessesToExtract = min_numberProcessesToExtract : max_numberProcessesToExtract

    @eval @everywhere parallelTotalProcesses=$numberProcessesToExtract

    @info "Extracting $numberProcessesToExtract mutational signatures for $totalIterations iterations using $numberWorkers cores"

    #signature extraction
    @time processes, exposures, processStab, processStabAvg, exposures_fitting, avgReconstructionError, avgReconstructionErrorPercentage, avgStability, reconError, reconErrorPercentage, simError = deconvoluteMutCounts(numberProcessesToExtract, totalIterations, sort_idx, n_types_to_remove)

    # save outfiles for each signature extracted
    file = string(SavePath, "/exposures_", parallelTotalProcesses)
    CSV.write(file, DataFrame(exposures), header = sample_header, delim = "\t")

    file = string(SavePath, "/processes_",  parallelTotalProcesses)
    CSV.write(file, DataFrame(processes), delim = "\t")

    file = string(SavePath, "/processesStabAvg_", parallelTotalProcesses)
    processes_stab = convert(DataFrame, processStabAvg)
    CSV.write(file, processes_stab, delim = "\t")

    file = string(SavePath, "/exposures_fitting_", parallelTotalProcesses)
    CSV.write(file, DataFrame(exposures_fitting), header = sample_header, delim = "\t")

    file = string(SavePath, "/reconError_", parallelTotalProcesses)
    CSV.write(file, DataFrame(reconError'), header = sample_header, delim="\t")

    file = string(SavePath, "/reconErrorPercentage_", parallelTotalProcesses)
    CSV.write(file, DataFrame(reconErrorPercentage'), header = sample_header, delim="\t")

    file = string(SavePath, "/simError_", parallelTotalProcesses)
    CSV.write(file, DataFrame(simError'), header = sample_header, delim="\t")

    file = string(SavePath, "/avgReconstructionError_", parallelTotalProcesses)
    writedlm(file, avgReconstructionError,  "\t")

    file = string(SavePath, "/avgReconstructionErrorPercentage_", parallelTotalProcesses)
    writedlm(file, avgReconstructionErrorPercentage,  "\t")

    file = string(SavePath, "/avgStability_", parallelTotalProcesses)
    writedlm(file, avgStability,  "\t")
end

println("Removing workers...")
if numberWorkers>1
    for i in workers()
        rmprocs(i)
    end
end

