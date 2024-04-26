# Preclinical_Spiral_Correction

This repository contains instructions and source code to reproduce the results presented in 

> On the correction of spiral trajectories on a preclinical MRI scanner with a high-performance gradient insert
> Scholten H, Wech T, Köhler S, Smart S S, Boyle J H, Teh I, Köstler H, Schneider J E.
> (currently under review)

Please cite this work if you use the content of this repository in your project.

The code is written in MATLAB (R2023b).

## Preparation
In order for all scripts to run faultless, the *Michigan image reconstruction toolbox (MIRT)* has to be downloaded from [https://web.eecs.umich.edu/~fessler/irt/fessler.tgz](https://web.eecs.umich.edu/~fessler/irt/fessler.tgz) or [http://web.eecs.umich.edu/~fessler/irt/irt](http://web.eecs.umich.edu/~fessler/irt/irt). Information about the toolbox can be found [here](https://web.eecs.umich.edu/~fessler/code/). The content of the toolbox (the folder named *irt*) needs to be placed inside the folder named *MIRT_toolbox* of this repository.

Furthermore, the data used for this publication need to be downloaded from [zenodo](https://zenodo.org/). (Link will be updated once the zenodo repository is created.)

The repository should finally contain the following folders:
* GSTF_code
* GSTF_data
* invivo_code
* invivo_data
* MIRT_toolbox/irt/...
* phantom_code
* phantom_data
* Spiral_Recon_NUFFT

## Scripts to reproduce the figures from the paper

### GSTF_code/compare_GSTFs.m

This script reproduces Figure 1 of the paper, i.e. it plots the magnitude and phase of the GSTFs of the x-, y-, and z-axis of the gradient system, it compares the forward calculation of the triangle waveforms used for the GSTF determination to the measured ones, it plots one interleaf each of the three spiral trajectories examined in the paper, and it calculates the spectra of the associated gradient waveforms and displays them im conjunction with the GSTFs.

### phantom_code/delay_optimization.m

This script reproduces Figure 2 of the paper. It takes the nRMSE values calculated for the isotropic delay correction and the GSTF+delay correction for different delays, fits a polynomial of 4th degree to them and finds the respective minima. The nRMSE values are stored in the folders *phantom_data* and *invivo_data* and have been calculated by the scripts *phantom_code/reconstruct_singlePhantom_multiDelay.m* and *invivo_code/reconstruct_mouse_multiDelay.m*, respectively.

### phantom_code/reconstruct_phantoms_loop.m

This script reproduces Figures 3, 4, 5, and 6 of the paper by reconstructing the phantom images acquired with the different spiral trajectories with the different trajectory correction methods. The required raw data and trajectory information are stored in the folder *phantom_data*. The script then plots the phantom images, line profiles, gradient waveforms and trajecotries.

### phantom_code/reconstruct_numericalphantom_loop.m

This script reproduces Figure 7 of the paper by repeating the reconstructions from *reconstruct_phantoms_loop.m* on a modified Shepp-Logan phantom.

### invivo_code/reconstruct_mouse.m

This script reproduces Figures 8 and 9 of the paper. It reconstructs the spiral and Cartesian in vivo images, and displays them as well as the gradient waveforms and k-space trajectories. The necessary data files are stored in the folder *invivo_data*.
