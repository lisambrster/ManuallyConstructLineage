# ManuallyConstructLineage
Hayden's Matlab Code for manually constructing lineage trees from segmentation output with lots of nice visualizations

start by installing CPD for matlab https://sites.google.com/site/myronenko/research/cpd You need to fill out the form, download the stuff, install Xcode on your mac (this can take hours) -

First run: PrecomputeRegistrationTransforms.m
update the paths to your installation of CPD, where you downloaded the data and where you want to store the registration transforms
you may need to update which time steps you want to process
This will determine the registration transform for every consecutive pair of images and can take a long
time depending on the number of time steps and their difficulty.

Then run: MakeLineage.m
Again - update the paths.
You can set breakpoints in matlab for each iteration and update the graph as needed.
An example of adding a new node and link is given as a comment at the end of the program
