@everywhere using JuMP
@everywhere using Ipopt
@everywhere using MathOptInterface
@everywhere using Optim
@everywhere using LinearAlgebra
@everywhere using DelimitedFiles
@everywhere using Distances

@everywhere function evalSingleSample(exposures, allSignatures, genome)

    maxMutations = sum(genome);

    numberOfSignatures = count(!iszero, exposures);

    x0 = maxMutations * rand(numberOfSignatures);
    x0 = x0 / sum(x0) * maxMutations;
    A = ones(1, numberOfSignatures);
    max_channels = size(genome, 1);

    # select only the signatures with exposure
    subSignatures = allSignatures[:,findall(!iszero, exposures)];

    # JuMP model instance
    model = Model(with_optimizer(Ipopt.Optimizer, print_level=0, sb="yes"));

    # set variables with bounds
    @variable(model, 0 <= x[i=1:numberOfSignatures] <= maxMutations, start=x0[i]);

    # set objective
    @NLobjective(model, Min,  sum( (sum(subSignatures[i,j]*x[j] for j=1:numberOfSignatures) - genome[i])^2 for i=1:max_channels));

    # set additional constraints
    @constraint(model, sum(A[i]*x[i] for i=1:numberOfSignatures) == maxMutations);

    # solve it
    JuMP.optimize!(model);
    x = JuMP.value.(x);

    x = round.(x);
    maximum_x_index = argmax(x);

    if sum(x) != maxMutations
         x[maximum_x_index] = x[maximum_x_index] + maxMutations - sum(x);
    end

    exposuresSample = zeros(size(allSignatures, 2));
    exposuresSample[findall(!iszero, exposures)] = x ;
    recon = allSignatures*exposuresSample;
    accr = 1-cosine_dist(recon, genome);
    error_fro = norm(recon-genome);
    norm_genome_fro = norm(genome);
    frob_rel_div = error_fro/norm_genome_fro;
    norm_one_dif = norm(recon - genome, 1);

    return accr, frob_rel_div, norm_one_dif, exposuresSample
end