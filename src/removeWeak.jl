
function removeWeak(mut_counts, remove_weak_muttypes)

    sum_mut_counts = reshape(sum(mut_counts, dims=2), size(mut_counts, 1))
    sort_idx = sortperm(sum_mut_counts)
    sort_muts = sum_mut_counts[sort_idx]

    cumsum_sort_muts = cumsum(sort_muts) / sum(mut_counts)
    n_types_to_remove = length(cumsum_sort_muts[cumsum_sort_muts .< remove_weak_muttypes])

    return mut_counts[sort_idx[n_types_to_remove+1:end], :], sort_idx, n_types_to_remove
end
