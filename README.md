# Preprocessing
data preprocessing scripts for jet images construction

This is the adaptation of the previous pre-processing code.
It has been altered to be able to run with the ROOT 6 framework yet still be backwards compatible with ROOT 5.

Code usage. 

Requires an input root data file, either W_signal.root or QCD_background.root. 

Currently the scripts are split up into 2 forms.

1. for the QCD background images
2. for the W Signal images

The plan is to merge them two in the future for better scalability and maintenance

- previous bugs. 

1. redefinition of variables, forbidden in ROOT 6
2. Histogram axes where made from -1 to -1, ROOT 5 ignores this, but ROOT 6 defaults to the full scale histogram
3. name on data file was not correct, i.e. typo.
