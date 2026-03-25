# m6A Bayesian Mixture model

By Lisa Hülsmann for *Höhn and Pittroff et al.*

This repository contains the code to prepare the data and fit the Bayesian mixture model for the effect of m6A position on differential gene expression upon m6A loss.

It requires a dataset of site level methylation rates and log2 fold changes in transcript abundance that should be added to the folder [data](/data).

Step 1: Aggregate transcript level data with the script [1_site_level_prep.R](/1_site_level_prep.R).

Step 2: Fit and analyze the Bayesian mixture model with the script [2_mixture_of_experts.R](#0) that calls the STAN model [model.stan](#0).

Results are stored in the folder [output](/output).
