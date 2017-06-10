module MCMCDiagnostics

using StatsBase

export chain_statistics, effective_sample_size, potential_scale_reduction

"""
Summary statistics from (part of) a chain for convergence diagnostics.
"""
struct ChainStatistics{T <: AbstractFloat}
    "Number of draws."
    sample_length::Int
    "Sample mean."
    sample_mean::T
    "Sample variance."
    sample_var::T
    "Scale factor for effective sample size, ∈ [0,1]."
    ess_factor::T
    """Last lag used for calculating scaling factor for effective
    sample size. For diagnostic purposes only."""
    last_lag::Int
end

"""
Calculate summary chain statistics (as a `ChainStatistics` object) for
vector of scalar draws.

# Recommended usage

1. Use multiple (3-5) chains, started from overdispersed initial
points.

2. After discarding burn-in, split the remaining sample in half.

3. Calculate the chain statistics for each.

4. Calculate [`potential_chain_reduction`](@ref) and
[`effective_sample_size`](@ref) using all the chains.
"""
function chain_statistics{T}(x::AbstractVector{T})
    N = length(x)
    m, v = mean_and_var(x, corrected = false)
    z = (x-m)/√v # normalize by variance for better numerical stability
    autocov(k) = dot(@view(z[1:(N-k)]), @view(z[(1+k):N])) / (N-k)
    invscale = 1 + 2*autocov(1)
    last_lag = 2
    while last_lag < N-2
        increment = autocov(last_lag) + autocov(last_lag + 1)
        if increment < 0
            break
        else
            invscale += 2*increment
            last_lag += 2
        end
    end
    ChainStatistics(N, m, v, 1 / invscale, last_lag)
end

"Effective sample size.

Estimated from autocorrelations. See Gelman et al (2013), section 11.4."
effective_sample_size(stat::ChainStatistics) = stat.sample_length * stat.ess_factor

effective_sample_size(stats::ChainStatistics...) = sum(effective_sample_size.(stats))

"""
Potential scale reduction factor (for possibly ragged chains).

Always ≥ 1 by construction, but values much larger than 1 (say 1.05)
indicate poor mixing.

Uses formula from Stan Development Team (2017), section 28.3.
"""
function potential_scale_reduction(stats::ChainStatistics...)
    W = mean(stat.sample_var for stat in stats)
    B = var(stat.sample_mean for stat in stats)
    √(1+B/W)
end

end # module
