# ImplantEFV

Function to calculate the incidence rate ratio of pregnancy for efavirenz-containing antiretroviral therapy compared to nevirapine via:

1. Poisson regression, analyzing each dataset seperately (unweighted analysis);
2. Inverse-probability weighting (IPW);
3. Generalized raking.

For both IPW and Generalized raking, we did a 2-phase analysis, using the EMR and chart review dataset as well as a three-phase analysis, using the EMR, chart review, and telephone interveiew datasets. For the 3-phase analysis, the sampling probability is the product of the two sampling probabilities (EMR to chart review and from chart review to telephone interview).
