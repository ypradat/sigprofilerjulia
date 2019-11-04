@everywhere function removeAllSingleSignatures(exposures_init, allSignatures, genome,)

    accr_first, frob_rel_div_first, norm_one_dif_first, exposures = evalSingleSample(exposures_init, allSignatures, genome)

    exposuresOutput = exposures;
    accr = accr_first;
    frob_rel_div = frob_rel_div_first;
    norm_one_dif =  norm_one_dif_first;

    numSig = count(!iszero, exposures);

    for j = 1 : ( numSig - 1)

        accr_temp, frob_rel_div_temp, norm_one_dif_temp, exposuresSample_temp, fID = removeOneSig(exposuresOutput, allSignatures, genome);

        if ((accr_first - accr_temp ) < 0.01)

            exposuresOutput = exposuresSample_temp;
            accr = accr_temp;
            frob_rel_div = frob_rel_div_temp;
            norm_one_dif =  norm_one_dif_temp;
        else
            break
        end
    end

    return exposuresOutput, accr, frob_rel_div, norm_one_dif
end
