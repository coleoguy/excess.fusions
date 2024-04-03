This directory contains the thirteen timetrees (Nexus format) estimated with
`MCMCtree` in the second step of the sequential Bayesian dating approach, as well as
the final stitched tree of 4,705 taxa (`4705sp_mammal-time.tree`, Nexus format).

File `4705sp_colours_mammal-time.tree` (Nexus format) has the instructions to plot the
stitched tree with `FigTree` in radial format, with 95% CIs of node ages, and with
the main clades coloured; all as in the main figure in the paper.

The `*.nwk` files contain the stitched tree in Newick format. There are three
versions depending on the information plotted for each node on the mammal phylogeny:
posterior mean of node ages (`4705sp_mean.nwk`), and the 2.5% and 97.5% quantiles of
the HPD-Credibility intervals (`4705sp_ci025.nwk` and `4705sp_ci975.nwk`,
respectively). You can use these trees in macroevolutionary analysis
(e.g., trait regression).

---
## **NOTE**
We have created a GitHub repository,
[mammals_dating](https://github.com/sabifo4/mammals_dating), where we explain
all the steps involved in the Bayesian sequential-subtree analysis described in
this study.

[Here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2), in
the last section, you will find more information about how we ran `MCMCtree`
with each data subset.
[Here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/03_Generate_final_mammal_tree)
you will find more details about how we generated the 4,705-taxon mammal tree of
life.
