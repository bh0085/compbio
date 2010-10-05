This directory has two directories:
inputnetworks: has the 4 networks, bRN (chip-binding), mRN (motif), kRN (REDfly), fRN (functional network) used as input for our predictive model
Each line is a an edge: first column is the TF and the second column is the target.

regressionwts: has directories for each network.
For each network we have 10 learned networks: nw_0-9.sif and bias_0-9.sif
The format of the nw_0-9.sif file is as follows:
Target TF regression_wt/coefficient

To make a prediction for a target using a TF's expression levels, take the
linear combination of the expression levels using the coeffients as the
weights, and add the bias term.
