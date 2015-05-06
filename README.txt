The supplied MALTAB code and files allow for the execution of
a. Procedure for estimation of Beta parameters for the three expression classes which are part of the probabilistic model used in the paper,
b. Filtering and imprinting analysis on the Geuvadis LCLs dataset.

Below we describe how to execute each of these two procedures.
Execution requires MATLAB.

Additionally, the downstream processing steps and R code in the downstream_analysis directory documents the analysis for the results presented in the paper. 

Parameter estimation:
====================

1. Replace the variable “refdir” with the path to the imprinting_code directory in your computer.
2. Cd into the “params4model” directory and execute the script param4model.m -
either type
/Applications/MATLAB_R2014a.app/bin/matlab -nodesktop -nosplash -r "param4model()"
or start MATLAB and type “param4model()” in the command window.
The procedure should begin running and output intermediate and final parameter estimates.

Filtering and imprinting analysis:
=================================

1. Replace the variable “refdir” with the path to the imprinting_code directory in your computer.
2. Choose between a sample run on two genes to a run on the complete Geuvadis LCLs dataset by setting the “datafile” variable.
3. Cd into the “analysis” directory and execute the script filter_model.m -
either type
/Applications/MATLAB_R2014a.app/bin/matlab -nodesktop -nosplash -r "filter_model()"
or start MATLAB and type “filter_model()” in the command window.

Analysis results are written to the “results” directory. The output includes a list of genes ranked according to the impglr statistic, as well as additional statistics and relevant information per gene. All output fields are explain in the supplementary text.
In addition, three figures will be generated per each of the top twenty ranked genes (while considering only genes for which at least one SNPs with both allele expressed monoallelically, see supplementary text). In all figures each dot depicts one informative sites, and the x,y axes give the ref,alt counts for the site, respectively. The colors are individual-specific and are different for each of the three figures; colors in the *BAL.jpg,*IMB.jpg and *IMP.jpg denote the probability (per-gene, per-individual Z, see supplementary text) of this individual belonging to the balanced, imbalanced and imprinted classes, respectively.
