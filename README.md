# Spatial subdomains in the Optic Tectum for the encoding of visual information

This repository contains the code used in the analysis of [Shallcross et al 2023](https://www.biorxiv.org/content/10.1101/2023.05.15.540762v1)

## Data Analysis

The analysis pipeline takes the output of 2photon imaging data which has been processed using [suite2p](https://suite2p.readthedocs.io/en/latest/#). There is also the option to input a binary mask to select subsets of the neurons. The mask should be the same dimensions as was used during 2photon imaging and saved as a tiff, with one mask per imaging plane.

## Data structure

/path/to/data/experiment_1/suite2p/combined/suite2p_data
/path/to/data/experiment_1/experiment_1_slice_i_mask.tiff

/path/to/data/experiment_2/suite2p/combined/suite2p_data
/path/to/data/experiment_2/experiment_2_slice_i_mask.tiff

prefix_list <- c("experiment_1", "experiment_2")

prefix <- prefix_list[1]

## Pipeline

0_run_analysis.R - runs the analysis pipeline

1_preprocess.R - takes the output from suite2p and extracts tectal neuron calcium transients and spatial coordinates

2_binarise_traces.R - takes the timeseries and binarises the activity, then finds those neurons whose firing correlates with the presentation of the stimuli

3_calculate_MI.R - calculates the mutual information between neuronal response and the presentation of dot stimuli and grating stimuli

4_ant_post_bias.R - estimates the distribution of visually responsive neurons across the anterior-posterior axis of the tectum 

5_gaussian_process.R - estimates how the average amount of information each neuron encodes changes across the anterior-posterior and medial-lateral axis of the tectum

6_neuronal_subtypes.R - splits up visually responsive neurons into feature selective subtypes

7_ant_post_bias_subtypes.R - estimates the distribution of neuronal subtypes across the anterior-posterior axis of the tectum

8_linear_regression.R - calculates how the amount of information per subtype changes as a function of its location across the anterior-posterior axis of the tectum

