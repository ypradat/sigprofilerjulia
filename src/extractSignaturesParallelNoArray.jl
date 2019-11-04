
@everywhere using SharedArrays
@everywhere using Distributions
@everywhere using Random
@everywhere const EPS = 2.2204e-16

@everywhere function nmf(v,r)

    n, m = size(v)
    stopconv = 10000        # stopping criterion (can be adjusted)
    niter = 1000000         # maximum number of iterations (can be adjusted)

    cons = fill(false, (m,m))
    consold = fill(false, (m,m))
    inc = 0
    j = 0

    # initialize random w and h
    w = Random.rand(n, r)
    h = Random.rand(r, m)

    for i = 1:niter

        # divergence-reducing NMF iterations
        x1 = repeat(sum(w, dims=1)', 1, m)
        h = h.*(w'*(v./(w*h)))./x1
        x2 = repeat(sum(h, dims=2)', n, 1)
        w = w.*((v./(w*h))*h')./x2

        # test convergence every 10 iterations

        if i%10 == 0 
            j = j + 1

            # adjust small values to avoid undeflow
            h = map(x->max(x, EPS), h)
            w = map(x->max(x, EPS), w)

            # construct connectivity matrix
            index = [i[1] for i in argmax(h, dims=1)]    # find largest factor
            mat1 = repeat(index, m, 1)                   # spread index down
            mat2 = repeat(index', 1, m)                  # spread index right
            cons = mat1==mat2

            if cons == consold             # connectivity matrix has not changed
                inc = inc + 1              # accumulate count 
            else
                inc = 0                    # else restart count
            end
            
            if inc > stopconv 
                break                      # assume convergence if connectivity stops changing
            end 

            consold = cons

        end
    end

    return w, h
end


@everywhere function to_parallelize(i)
    
    if i%20 == 0
            @info "Iteration $i."
    end

    # Create a random seed for each iteration
    Random.seed!(i+123456)

    # Generating boostrapped genomes
    bootstrapGenomes = Array{Float64,2}(undef, size(parallelNormGenomes))
    for j = eachindex(parallelGenomesColsSum)
        bootstrapGenomes[:,j] = Random.rand(Multinomial(parallelGenomesColsSum[j], parallelNormGenomes[:,j]))
    end

    bootstrapGenomes = map(x->max(x, EPS), bootstrapGenomes)

    # Solving NMF for these boostrapped genomes
    W, H = nmf(bootstrapGenomes, parallelTotalProcesses)

    for j = 1 : parallelTotalProcesses
        total = sum( W[:, j] )
        W[:, j] = W[:, j] / total
        H[j, :] = H[j, :] * total
    end
    
    return i, bootstrapGenomes, W, H
end


function extractSignatures(genomes, processes, iterations)

    totalMutationTypes = size(genomes, 1)
    totalGenomes = size(genomes, 2)
    
    Wall = zeros( totalMutationTypes, processes * iterations )
    Hall = zeros( processes * iterations, totalGenomes )
    genomeErrors = zeros(totalMutationTypes, totalGenomes, iterations)
    genomesReconstructed = zeros(totalMutationTypes, totalGenomes, iterations)

    for (i, bootstrapGenomes, W, H) = pmap(to_parallelize, 1:iterations)
        processCount = (i-1)*processes + 1
        genomeErrors[:, :, i] = bootstrapGenomes -  W * H
        genomesReconstructed[:, :, i] = W * H
        Wall[ :, processCount : (processCount + processes - 1) ] = W
        Hall[ processCount : (processCount + processes - 1), : ] = H

    end

    return Wall, Hall, genomeErrors, genomesReconstructed
end
