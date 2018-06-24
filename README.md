# survbayes2
survbayes2

### Installing
The R package can be installed either via git clone, see tutorial here (https://bit.ly/2MgUPIu) or by using devtools in Rstudio:
```
if(!require("devtools")) install.packages("devtools")
devtools::install_github("csetraynor/survbayes2")
```

## Getting Started

We fltered a datset of iC-2 (n = 72) patients from the METABRIC cohort [1], which included
clinical and genomic covariates. In addition, treatment efects were also considered such as, chemotherapy,
radio-therapy, hormone-therapy or surgery (mastectomy or breast conservation). We propose a Poisson
generalised additive model with log link function that relates the hazard ratio to a linear combination of the
log-hazard ratios, or β parameters, and X the n*p matrix of covariates; the logartithm of the diferential time
τ i as an of-set variable; and a low-rank thin-plate splines function f(·), where the fxed knots k k are shrinked
towards a frst degree polynomial to avoid over-ftting.
[HormoneSurgeryInteraction.pdf](https://github.com/csetraynor/survbayes2/files/2130697/HormoneSurgeryInteraction.pdf)

Model performance was measured via the BS, which is defned as the squared distance
between the predicted and observed survival outcomes. Three covariate models were compared: C, CG and
CGwT. A Bayesian hierarchical model was built on the Monte Carlo cross-validation results by assuming the BS
measures, x i , to be jointly multivariate normal with average diference, μ 0 , the quantity of interest.

[MC_Results.pdf](https://github.com/csetraynor/survbayes2/files/2130698/MC_Results.pdf)

## Acknowledgment

I would like to thank my supervisors Prof Michael Chappell,
Dr Neil Evans, Dr Tarj Sahota and Ms Helen Tomkinson for giving me the opportunity to study for this PhD at University of Warwick. 

In addition, many thanks to the original authors of the study METABRIC (Pereira et al) and the creators of cBioPortal (Gao et al) for making easier to share knowledge in biology and promote the development of science that may find cures for the difficult cancerous diseases.

## References

Pereira, Bernard, et al. "The somatic mutation profiles of 2,433 breast
 cancers refine their genomic and transcriptomic landscapes." 
Nature communications 7 (2016): 11479.

Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012). Evaluating
  Random Forests for Survival Analysis Using Prediction Error Curves.
  Journal of Statistical Software, 50(11), 1-23. URL
  http://www.jstatsoft.org/v50/i11/.

Max Kuhn and Hadley Wickham (2017). rsample: General Resampling
  Infrastructure. R package version 0.0.2.
  https://CRAN.R-project.org/package=rsample
