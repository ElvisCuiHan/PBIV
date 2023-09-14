# Parametric Bayesian Instrumental Variable (PBIV)

This `R` package extends [Li-Lu's PBIV estimation](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.6369), which was orignally for right-censored data only, to arbitrary censoring time-to-event outcome. The *PBIV* package deals with causal effect estimation using instrumental variables.

**Installation:**

    devtools::install_github("https://github.com/ElvisCuiHan/PBIV/")

There is one main function `IV_MH_IC` in the package. It deals with bivariate normal random errors with multiple instruments and multiple observed potential confounders within the [two-stage linear model framework](https://www.tandfonline.com/doi/pdf/10.1080/01621459.1995.10476535?casa_token=eOGWmrHSCSQAAAAA:azlG4TqGnptppixI5srRZ7R47z_pv5bbAsDU-I8_oRFzrEC_slXF9MJ7e-YlpFRYX5UoyoW_DSuk).

**Instrumental variables regression:**

    library("PBIV")
    IV_MH_IC(L, R, d, X, G, U = NULL, m, wid, init = NULL, prior_1, prior_2)

where
- L: $n\times1$ vector refers to left-observed time to event or censoring
- R: $n\times1$ vector refers to right-observed time to event or censoring
- d: $n\times1$ vector refers to censoring status, *4=event, 3=right-censored, 2=interval-censored, 1=left-censored*
- X: $n\times1$ vector refers to the covariate of interest, i.e., the exposure using the language of randomized clinical trials (RCT)
- G: $n\times kG$ matrix if multiple instruments are used
- U: $n \times kU$ matrix if multiple observed confounders are used; the default is NULL where n is sample size, $kG$ is number of instruments, $kU$ is number of observed confounders A total of $(6+kG+2*kU)$ parameters to estimate: $a_0$, $a_1$ (length=$kG$), $a_2$ (length=$kU$), $\sigma_1^2$, $b_0$, $b_1$, $b_2$ (length=$kU$), $\sigma_2$, and $\rho$
- m: A scalr refers to the number of iterations of the MCMC algorithm
- wid: vector of the random walk width for $(a_0,a_1,a_2,\sigma_1^2,b_0,b_1,b_2,\sigma_2^2,\rho)$ in the MCMC algorithm
- init:	A vector of the initial values for $(a_0,a_1,a_2,\sigma_1^2,b_0,b_1,b_2,\sigma_2^2,\rho)$, default value is NULL
- prior_1: A vector of the first parameter of the priors for $(a_0,a_1,a_2,\sigma_1^2,b_0,b_1,b_2,\sigma^2_2)$: mean of the normal priors for $a_0,a_1,a_2,b_0,b_1,b_2$; shape parameter of the inverse-gamma priors for $\sigma_1^2, \sigma_2^2$
- prior_2: A vector of the second parameter of the priors for $(a_0,a_1,a_2,\sigma_1^2,b_0,b_1,b_2,\sigma_2^2)$: variance of the normal priors for $a_0,a_1,a_2,b_0,b_1,b_2$; scale parameter of the inverse-gamma priors for $\sigma_1^2, \sigma_2^2$

The ``
