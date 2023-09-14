# Parametric Bayesian Instrumental Variable (PBIV) Estimation

This `R` package extends [Li-Lu's PBIV estimation](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.6369), which was orignally for right-censored data only, to arbitrary censoring time-to-event outcome. The *PBIV* package deals with causal effect estimation using instrumental variables.

**Installation:**

    devtools::install_github("https://github.com/ElvisCuiHan/PBIV/")

**Instrumental variables regression:**

    library("PBIV")
    ivreg(L, R, d, X, G, U = NULL, m, wid, init = NULL, prior_1, prior_2)

where
- L: A $n\times1$ vector refers to left-observed time to event or censoring
- R: A $n\times1$ vector refers to right-observed time to event or censoring
- d: A $n\times1$ vector refers to censoring status, *4=event, 3=right-censored, 2=interval-censored, 1=left-censored*
- X: A $n\times1$ vector refers to the covariate of interest, i.e., the exposure using the language of randomized clinical trials (RCT)
- G: A $n\times kG$ matrix if multiple instruments are used
