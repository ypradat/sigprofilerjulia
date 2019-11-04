function addWeak(sort_idx, n_types_to_remove, processes_I, processesStd_I, Wall_I, genomeErrors_I, genomesReconstructed_I)
    
    totalMutTypes = size(Wall_I, 1) + n_types_to_remove
    processes = zeros(totalMutTypes, size(processes_I, 2))
    processesStd = zeros(totalMutTypes, size(processesStd_I, 2))
    Wall = zeros(totalMutTypes, size(Wall_I, 2))
    genomeErrors = zeros(totalMutTypes, size(genomeErrors_I, 2), size(genomeErrors_I, 3))
    genomesReconstructed = zeros(totalMutTypes, size(genomesReconstructed_I, 2), size(genomesReconstructed_I, 3))
   
    for i = 1 : size(Wall_I, 1)
        idx = sort_idx[i+n_types_to_remove]
        processes[idx, :] = processes_I[i, :]
        processesStd[idx, :] = processesStd_I[i, :]
        Wall[idx, :] = Wall_I[i, :]
        genomeErrors[idx, :, :] = genomeErrors_I[i, :, :]
        genomesReconstructed[idx, :, :] = genomesReconstructed_I[i, :, :]
    end

    return processes, processesStd, Wall, genomeErrors, genomesReconstructed
end