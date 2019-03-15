Statistical Methods for Replicability Assessment
========
This is the repository accompanying "Statistical Methods for Replicability Assessment"

Contributors
--------
Kenneth Hung, William Fithian

LaTeX Files
--------
The files related to LaTeX are `paper.tex` (LaTeX source for paper) `supplement.tex` (LaTeX source for supplement), `papers.bib` (bibliography) and `additional.bib` (additional bibliography). Furthermore, there are additional files related to submission to journals and slides for presentation.

R Files
--------
`shared.R` contains almost all functions, helper or not. `ecdf.R` plots the original and replication empirical cumulative distributions of the *p*-values. `effect-size-k.R` contains the real constants *k* used in `effect-size.R` that performs the selective test, and produces predictive interval and confidence interval. `fdp.R` computes the directional FDP and plots the histograms of *p*-values. `sel-bias-sim.R` simulates and shows how the RP:P metrics can be lower despite exact replications. `t-approx.R` investigates the precision of our normal approximations.

PDF Files
--------
`paper.pdf` is the paper and `supplement.pdf` is the supplement, the rest are all figures located in the `fig/` folder.
