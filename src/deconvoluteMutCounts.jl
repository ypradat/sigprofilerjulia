@everywhere function remove_sigs(j)

    exposuresOutput, accr, frob_rel_div, norm_one_dif = removeAllSingleSignatures(parallelExposure[:,j], parallelProcesses , originalGenomes[:, j])

    return j, exposuresOutput, accr, frob_rel_div, norm_one_dif
end

function deconvoluteMutCounts(
        processes,
        iterations,
        sort_idx,
        n_types_to_remove,
    )

    totalReplicates = 100
    removeLastPercentage = 0.07

    Wall, Hall, genomeErrors, genomesReconstructed = extractSignatures(genomes, processes, iterations)
    Wall, Hall, genomeErrors, genomesReconstructed = filterOutIterations( Wall, Hall, genomeErrors, processes, genomesReconstructed, removeLastPercentage )

    processes, processesStd, exposure, exposureStd, idx, idxS, processStab, processStabAvg, clusterCompactness = evaluateStability(
        Wall, Hall, processes, totalReplicates)

    processes, processesStd, Wall, genomeErrors, genomesReconstructed = addWeak(
        sort_idx, n_types_to_remove, processes, processesStd, Wall, genomeErrors, genomesReconstructed)


    # signature fitting
    exposure_copied = deepcopy(exposure)

    @eval @everywhere parallelExposure=$exposure_copied
    @eval @everywhere parallelProcesses=$processes

    for (j, exposuresOutput, accr, frob_rel_div, norm_one_dif) in pmap(remove_sigs, 1:size(parallelExposure,2))
        exposure[:, j] = exposuresOutput
    end

    # calculate stabilities and errors in reconstruction
    avgReconstructionError, avgReconstructionErrorPercentage, avgStability, reconError, reconErrorPercentage, simError = calculate_stabilities(processes, exposure, processStabAvg)

    return processes, exposure_copied, processStab, processStabAvg, exposure, avgReconstructionError, avgReconstructionErrorPercentage, avgStability, reconError, reconErrorPercentage, simError
end