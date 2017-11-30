# MCMCDiagnostics

Markov Chain Monte Carlo convergence diagnostics in Julia.

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/tpapp/MCMCDiagnostics.jl.svg?branch=master)](https://travis-ci.org/tpapp/MCMCDiagnostics.jl)
[![Coverage Status](https://coveralls.io/repos/tpapp/MCMCDiagnostics.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/tpapp/MCMCDiagnostics.jl?branch=master)
[![codecov.io](http://codecov.io/github/tpapp/MCMCDiagnostics.jl/coverage.svg?branch=master)](http://codecov.io/github/tpapp/MCMCDiagnostics.jl?branch=master)

## Overview

This package contains two very useful diagnostics for Markov Chain Monte Carlo:

1. `potential_scale_reduction(chains...)`, which estimates the potential scale reduction factor, also known as R̂, for multiple scalar chains,
2. `effective_sample_size(chain)`, which calculates the effective sample size for scalar chains.

These are intended as *building blocks*, to be used by other libraries, and were organized into a separate library for testing and DRY.

## Installation

The package is not (yet) registered. If you find it useful and want me to register it, please open an issue. In the meantime, use

```julia
Pkg.clone("https://github.com/tpapp/MCMCDiagnostics.jl.git")
```

to install it.

## Related

You may find my other packages for MCMC interesting. See the documentation of [DynamicHMC.jl](https://github.com/tpapp/DynamicHMC.jl) for details.

## Bibliography

Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. Statistical science, 457-472.

Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B. (2013). Bayesian data analysis (3rd edition). Chapman & Hall/CRC.

Stan Development Team. (2017). Stan Modeling Language Users Guide and Reference Manual, Version 2.15.0. http://mc-stan.org
