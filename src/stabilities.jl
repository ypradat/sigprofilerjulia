using Distances

function calculate_stabilities(allSignatures, exposures, processStabAvg)

    recon = allSignatures*exposures;
    sz = size(exposures);
    totalSamples = sz[2];

    avgReconstructionError = norm(recon-originalGenomes);
    norm_genome_fro = norm(originalGenomes);
    avgReconstructionErrorPercentage = avgReconstructionError/norm_genome_fro;
    avgStability = mean(processStabAvg);

    reconError = zeros(totalSamples);
    reconErrorPercentage = zeros(totalSamples);
    simError = zeros(totalSamples);

    for j = 1 : totalSamples
        error_sample = norm(originalGenomes[:,j] - recon[:,j]);
        norm_sample_fro = norm(originalGenomes[:,j]);
        reconError[j] = error_sample;
        reconErrorPercentage[j] = error_sample/norm_sample_fro;
        accr = 1-cosine_dist(recon[:, j], originalGenomes[:, j]);
        simError[j] = accr;
    end

    return avgReconstructionError, avgReconstructionErrorPercentage, avgStability, reconError, reconErrorPercentage, simError
end
