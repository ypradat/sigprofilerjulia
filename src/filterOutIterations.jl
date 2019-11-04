function filterOutIterations(Wall, Hall, genomeErrors, numberProcessesToExtract, genomesReconstructed, removeLastPercentage)

    totalIterations = div(size(Wall, 2), numberProcessesToExtract)
    totalRemoveIter = convert(Int, round(removeLastPercentage * totalIterations ))
   
    closenessGenomes = zeros(totalIterations)
    for i = 1 : totalIterations
        closenessGenomes[i] = norm(genomeErrors[:, :, i])
    end

    indexClosenessGenomes = sortperm(closenessGenomes, rev=true)
    removeIterations = indexClosenessGenomes[1:totalRemoveIter]
  
    removeIterationSets = zeros(Int, numberProcessesToExtract * totalRemoveIter)
   
    for i = 1 : totalRemoveIter
        iStart = convert(Int, numberProcessesToExtract * ( removeIterations[i] - 1) + 1)
        iEnd = convert(Int, numberProcessesToExtract * removeIterations[i])
        removeIterationSets[(numberProcessesToExtract*(i-1)+1):(numberProcessesToExtract*i)] = iStart:iEnd
    end
   
    Wall = Wall[:, setdiff(1:end, removeIterationSets)] 
    Hall = Hall[setdiff(1:end, removeIterationSets), :]
    genomeErrors = genomeErrors[:, :, setdiff(1:end, removeIterations)]
    genomesReconstructed = genomesReconstructed[:, :, setdiff(1:end, removeIterations)]
    
    return Wall, Hall, genomeErrors, genomesReconstructed
end
