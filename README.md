# Parametric Bayesian Instrumental Variable (PBIV) Estimation

This `R` package extends [Li-Lu's PBIV estimation](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.6369), which was orignally for right-censored data only, to arbitrary censoring time-to-event outcome. The *PBIV* package deals with causal effect estimation using instrumental variables. 

**Installation:**

    install_github("https://github.com/ElvisCuiHan/PBIV/")

**Instrumental variables regression:**

    library("PBIV")
    ivreg(Q ~ P + D | D + F + A, data = Kmenta)
