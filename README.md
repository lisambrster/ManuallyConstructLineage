# ManuallyConstructLineage
Hayden's Matlab Code for manually constructing lineage trees from segmentation output with lots of nice visualizations

start by installing CPD for matlab https://sites.google.com/site/myronenko/research/cpd You need to fill out the form, download the stuff, install Xcode on your mac (this can take hours) -

Then run: track_nuclei_based_on_CPD_HardFromMasha_make_graph.m

update the paths to your installation of CPD and where you downloaded the data

This will enable you to do registration and initial matching of nuclei between two frames. You can set breakpoints in matlab for each iteration and update the graph as needed.
