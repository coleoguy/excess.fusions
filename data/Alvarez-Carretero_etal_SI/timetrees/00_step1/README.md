This directory contains the timetrees (72 taxa) estimated with `MCMCtree` for
the seven tree hypotheses (see corresponding calibrated tree files in `trees/00_step1`).

Note that file `00_main_tree_T2-updated-geochronology-time.tree` is the main tree
used as a backbone to generate the large 4,705-taxon tree. 

Files `00_main_tree_T2-narrowLcalibs-updated-geochronology-time.tree` and
`00_main_tree_T2-wideLcalibs-updated-geochronology-time.tree` 
are the timetrees estimated as part of the sensitivity analysis to the fossil
calibrations we carried out. With this analysis, we assessed the impact of fossil
calibration strategies on node age estimates.

File `FigTree-72s-updated-geochronology-prior.tree` is the timetree estimated 
with `MCMCtree` when sampling from the prior (i.e., no data) and using the calibrated 
main tree topology (file `00_main_tree_T2_72sp-updated-geochronology.tree` in `trees/00_step1`).

File `00_main_tree_T2-time.tree` is the timetree estimated with `MCMCtree` when using 
the non-updated version of the main tree topology (file `00_main_tree_T2_72sp.tree` in 
`trees/00_step1`). Note that this file is **not** used in the subsequent steps of the sequential 
Bayesian dating approach, but we provide it here to keep track of the differences 
in time estimates when using the non-updated and the updated set of calibrations as 
of September 2021.

## **NOTE** 
We have created a GitHub repository,
[mammals_dating](https://github.com/sabifo4/mammals_dating), where we explain
all the steps involved in the Bayesian sequential-subtree analysis described
in this study.

[Here](https://github.com/sabifo4/mammals_dating/tree/main/01_SeqBayes_S1/02_MCMCtree)
you will find more information about how we ran `MCMCtree` with this dataset
under each tree hypothesis.
