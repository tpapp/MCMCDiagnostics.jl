# MCMCDiagnostics.jl

**MCMCDiagnostics.jl has been deprecated** in favor of [MCMCDiagnosticTools.jl](https://github.com/TuringLang/MCMCDiagnosticTools.jl/) and is no longer maintained.

Markov Chain Monte Carlo convergence diagnostics in Julia.
<!-- Tidyverse lifecycle badges, see https://www.tidyverse.org/lifecycle/ Uncomment or delete as needed. -->
![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![build](https://github.com/tpapp/MCMCDiagnostics.jl/workflows/CI/badge.svg)](https://github.com/tpapp/MCMCDiagnostics.jl/actions?query=workflow%3ACI)
[![codecov.io](http://codecov.io/github/tpapp/MCMCDiagnostics.jl/coverage.svg?branch=master)](http://codecov.io/github/tpapp/MCMCDiagnostics.jl?branch=master)

## Overview

This package contains two very useful diagnostics for Markov Chain Monte Carlo:

1. `potential_scale_reduction(chains...)`, which estimates the potential scale reduction factor, also known as `Rhat`, for multiple scalar chains,

2. `effective_sample_size(chain)`, which calculates the effective sample size for scalar chains.

These are intended as *building blocks*, to be used by other libraries, and were organized into a separate library for testing and DRY.

## Installation

The package is registered. You can install it with

```julia
Pkg.add("MCMCDiagnostics")
```

## Related

You may find my other packages for MCMC interesting. See the documentation of [DynamicHMC.jl](https://github.com/tpapp/DynamicHMC.jl) for details.

## Bibliography

Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. Statistical science, 457-472.

Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B. (2013). Bayesian data analysis (3rd edition). Chapman & Hall/CRC.

Stan Development Team. (2017). Stan Modeling Language Users Guide and Reference Manual, Version 2.15.0. http://mc-stan.org
