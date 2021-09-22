using Random
# Defining function constants
const BIG_NUMBER = 100
const CONVERG_ITER = 10
const CONVERG_CUTOFF = 0.005  # cosine distance
const TOTAL_INIT_CONDITIONS = 5

function evaluateStability(Wall, Hall, totalProcesses, totalReplicates)
    
    Wall = map(x->max(x, EPS), Wall)
    Hall = map(x->max(x, EPS), Hall)

    Random.seed!(123456)

    # Clustering mutational processes using custom clustering procedure
    minClusterDist = BIG_NUMBER
    totalIter = div(size(Wall, 2), totalProcesses)
    idx = zeros(size(Hall, 1))
    clusterCompactness = zeros(totalProcesses, totalIter)
    iStartigDataSet = 1 : totalProcesses : size(Wall,2)
    iStartingDataSet = iStartigDataSet[randperm(totalIter)]

    local centroids, centroidsFinal, clusterCompactnessFinal, idxFinal, centroidStd
    local exposure, exposureStd, idxS, processStab, processStabAvg
    
    for iInitData = 1 : min(TOTAL_INIT_CONDITIONS, totalIter) % size(Wall, 2)
        iStartingData = iStartingDataSet[iInitData]
        iEnd = iStartingData + totalProcesses - 1
        centroids = Wall[:, iStartingData:iEnd]

        # this is random process
        centroidsTest = rand(Float64, size(centroids))
        countIRep = 0
    
        for iRep = 1 : totalReplicates
            allDist = Distances.pairwise(CosineDist(), hcat(centroids, Wall), dims=2)
            centroidDist = allDist[ 1:size(centroids, 2), (size(centroids, 2)+1):size(allDist, 2) ]'

            jRange = randperm(totalProcesses)
            for jIndex = 1 : totalProcesses
                j = jRange[jIndex]
                for i = 1 : totalProcesses : size(Wall, 2)
                    iRange = i: (i + totalProcesses - 1)
                    Ind = argmin( centroidDist[iRange, j] )[1]
                    centroidDist[iRange[Ind], :] .= BIG_NUMBER
                    idx[iRange[Ind]] = j
                end
            end

            maxDistToNewCentroids = 0
            for i = 1 : size(centroids, 2)
                centroids[:, i] = mean( Wall[:, idx .== i], dims=2 )
                maxDistToNewCentroids = max(maxDistToNewCentroids, maximum(Distances.pairwise(CosineDist(), hcat(centroids[:, i], centroidsTest[:, i]), dims=2)))
            end

            if maxDistToNewCentroids < CONVERG_CUTOFF
                countIRep =  countIRep + 1
            else
                countIRep = 0
                centroidsTest = centroids
            end

            if countIRep == CONVERG_ITER
                break
            end
        end

        for i = 1 : size(centroids, 2)
            clusterDist = Distances.pairwise(CosineDist(), hcat(centroids[:,i], Wall[:, idx.==i]), dims=2)
            clusterCompactness[i, :] = clusterDist[1, 2:size(clusterDist, 2)]
        end

        meanCompactness = mean(clusterCompactness)
        if minClusterDist > meanCompactness
            centroidsFinal = centroids
            idxFinal = copy(idx)
            clusterCompactnessFinal = clusterCompactness
            minClusterDist = meanCompactness
        end
    end

    centroids = centroidsFinal'
    idx = copy(idxFinal)
    clusterCompactness = clusterCompactnessFinal
        
    centDist = mean( clusterCompactness, dims=2 )[:,1]
    centDistInd = sortperm(centDist)  # centDistSorted not used
    clusterCompactness = clusterCompactness[centDistInd,:]
    
    
    centroids = centroids[centDistInd, :]
    idxNew = zeros(size(idx));

    for i = 1 : totalProcesses
        idxNew[idx .== centDistInd[i]] .= i
    end
    
    idx = copy(idxNew)
        
    
    if ( totalProcesses > 1)
        idxCount = countmap(idx)
        countVector = [get(idxCount, i, 0) for i in 1:maximum(idx)]
        processStab = silhouettes(trunc.(Int, idx), countVector, Distances.pairwise(CosineDist(), Wall, dims=2))
        processStabAvg = zeros(1, totalProcesses)
        for i = 1 : totalProcesses
            processStabAvg[i] = mean(processStab[idx.==i]) 
        end
    else 
        allDist = Distances.pairwise(CosineDist(), hcat(centroids', Wall), dims=2)
        processStab = 1 - allDist[1:size(centroids', 2), (size(centroids', 2)+1): size(allDist, 2) ]'
        processStabAvg = mean(processStab)
    end
        
    centroidStd = zeros( size(centroids) )
    for i = 1 : totalProcesses
        centroidStd[i,:] = std(Wall[:, idx.==i]; dims=2 )
    end

    centroids = centroids'
    centroidStd = centroidStd'

    idxS = zeros(size(idx))
    for i = 1 : totalProcesses : size(Wall, 2)
        iEnd = i + totalProcesses - 1
        idxG = idx[i:iEnd]
        
        for j = 1 : totalProcesses
            idxS[i+j-1] = findfirst(idxG .== j)
        end
    end

    exposure = zeros(convert(Int, maximum(idxS)), size(Hall, 2))
    exposureStd = zeros(size(exposure))
        
    for i = 1 : convert(Int, maximum(idxS))
        exposure[i, :] = mean(Hall[idx.==i, :], dims=1)
        exposureStd[i, :] = std(Hall[idx.==i, :], dims=1)
    end
        
    return centroids, centroidStd, exposure, exposureStd, idx, idxS, processStab, processStabAvg, clusterCompactness
end
