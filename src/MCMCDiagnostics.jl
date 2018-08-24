module MCMCDiagnostics

using Statistics: mean, var
using StatsBase: mean_and_var

export ess_factor_estimate, effective_sample_size, potential_scale_reduction

"""
    autocorrelation(x, k, v = var(x))

Estimate of lag-`k` autocorrelation of `x` from a variogram. `v` is the variance
of `x`, used when supplied.

See Gelman et al (2013), section 11.4.
"""
function autocorrelation(x::AbstractVector, k::Integer, v = var(x))
    x1 = @view(x[1:(end-k)])
    x2 = @view(x[(1+k):end])
    V = sum((x1 .- x2).^2) / length(x1)
    1 - V / (2*v)
end

"""
    ess_factor_estimate(x, v = var(x))

Estimate for effective sample size factor.

Return `τ, K` where `τ` is estimated effective sample size / sample size, and
`K` is the last lag used for autocorrelation estimation.

# Notes

See Gelman et al (2013), section 11.4.

`τ` is capped at 1, this is relevant when the sample has large negative
autocorrelation (happens with HMC/NUTS).

Some implementations (eg Stan) use FFT for autocorrelations, which yields the
whole spectrum. In practice, a <50-100 lags are usually sufficient for
reasonable samplers, so the “naive” version may be more efficient.
"""
function ess_factor_estimate(x::AbstractVector, v = var(x))
    N = length(x)
    τ_inv = 1 + 2 * autocorrelation(x, 1, v)
    K = 2
    while K < N - 2
        Δ = autocorrelation(x, K, v) + autocorrelation(x, K + 1, v)
        if Δ < 0
            break
        else
            τ_inv += 2*Δ
            K += 2
        end
    end
    min(1 / τ_inv, one(τ_inv)), K
end

"""
    effective_sample_size(x, v = var(x))

Effective sample size of vector `x`.

Estimated from autocorrelations. See Gelman et al (2013), section 11.4.

When the variance `v` is supplied, it saves some calculation time.
"""
function effective_sample_size(x::AbstractVector, v = var(x))
    τ, _ = ess_factor_estimate(x, v)
    τ * length(x)
end

"""
    potential_scale_reduction(chains...)

Potential scale reduction factor (for possibly ragged chains).

Also known as R̂. Always ≥ 1 by construction, but values much larger than 1 (say
1.05) indicate poor mixing.

Uses formula from Stan Development Team (2017), section 28.3.
"""
function potential_scale_reduction(chains::AbstractVector...)
    mvs = mean_and_var.(chains)
    W = mean(last.(mvs))
    B = var(first.(mvs))
    √(1 + B / W)
end

end # module
