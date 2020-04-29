A SEIR model with time-varying coefficients for analysing the SARS-CoV-2 epidemic
=====================================================================================

[Paolo Girardi](mailto://paolo.girardi@unipd.it)(1), [Carlo Gaetan](mailto://gaetan@unive.it)(2)

  1. Department of Developmental and Social Psychology  
     Università di Padova
     Via Venezia 8, 
     I-35131 Padua, Italy  
  2. Dipartimento di Scienze Ambientali, Informatica e Statistica,
     Università "Ca' Foscari" di Venezia,
     Campus Scientifico, Via Torino 155,
     I-30172 Mestre-Venezia, Italy
  

Corresponding author: [Paolo Girardi, paolo.girardi@unipd.it](mailto://paolo.girardi@unipd.it)
[Personal page](https://paolin83.github.io)

(analysis updated on 29 April 2020, 22.00)  

This repository contained the functions and the statistical tools employed in the manuscript.
The analysis can be found here  
[a SEIR analysis Italy](main_analysis.md)  
(html R-Markdown can be downloaded [here](main_analysis.Rmd))  

Abstract
--------

In this study, we propose a time-dependent Susceptible-Infected-Exposed-Recovered (SEIR)
differential model for the analysis of the SARS-CoV-2 epidemic outbreak in three different countries,
the United States of America, Italy and Iceland using public data inherent the numbers of the
epidemic wave. Several types and grades of actions were adopted by the governments, including
travel restrictions, social distancing, or limitation of movement. We want to investigate how these
measures can affect the epidemic curve of the infectious population. The parameters of interest for
the SEIR model were estimated employing a composite likelihood approach. Moreover, standard
errors have been corrected for the temporal dependence. The adoption of restrictive measures
results in flatten epidemic curves, and the future evolution indicated a decrease in the number of
cases.

License
-------

The package is distributed on CRAN at [this web page](https://cran.r-project.org/web/packages/spMC/index.html) under the terms of the [GENERAL PUBLIC LICENSE (>= 2)](https://cran.r-project.org/web/licenses/GPL-2). The package spMC version 0.3.6.3153 was edited in order to satisfy the requirement for the manuscript submission to *Computer and Geosciences* (Elsevier). This specific version will never appear on the CRAN servers, but it can be available on the Elsevier servers with the same terms and conditions specified by the [GENERAL PUBLIC LICENSE (>= 2)](https://cran.r-project.org/web/licenses/GPL-2).



