@everywhere using StatsBase

@everywhere function removeOneSig(exposures, allSignatures, genome)

    numberOfSignatures = count(!iszero, exposures)
    sig_IDs = findall(!iszero, exposures)
    totalSignatures = size(allSignatures, 2)

    if numberOfSignatures > 1

        exposuresSample_wo = zeros(totalSignatures, numberOfSignatures)
        accr_wo = zeros(numberOfSignatures, 1)
        frob_rel_div_wo = zeros(numberOfSignatures, 1)
        norm_one_dif_wo = zeros(numberOfSignatures, 1)

        for j = 1 : numberOfSignatures

            exposuresRemoveOne = deepcopy(exposures)
            exposuresRemoveOne[sig_IDs[j]] = 0;

            accr1, frob1, norm_one_dif1, exposuresSample1 =  evalSingleSample(exposuresRemoveOne, allSignatures, genome)

            accr_wo[j] = accr1
            frob_rel_div_wo[j] = frob1
            norm_one_dif_wo[j] = norm_one_dif1
            exposuresSample_wo[:, j] = exposuresSample1

        end

        accr_wo_min = 1 .- accr_wo

        fID = argmin(mapslices(harmmean, hcat(accr_wo_min, frob_rel_div_wo), dims = 2)[:, ])
        accr = accr_wo[fID];
        frob_rel_div = frob_rel_div_wo[fID];
        norm_one_dif = norm_one_dif_wo[fID];
        exposuresSample = exposuresSample_wo[:,fID];
        fID = sig_IDs[fID];

    else
        accr, frob_rel_div, norm_one_dif, exposuresSample =  evalSingleSample(exposures, allSignatures, genome)
        fID = 0;

    end

    return accr, frob_rel_div, norm_one_dif, exposuresSample, fID
end