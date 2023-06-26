## Spatial subdomains in the Optic Tectum for the encoding of visual information

This repository contains the code used in the analysis of [Shallcross et al 2023](https://www.biorxiv.org/content/10.1101/2023.05.15.540762v1)

1_preprocess.R - takes the output from suite2p and extracts tectal neuron calcium transients and spatial coordinates

2_binarise_traces.R - takes the timeseries and binarises the activity, then finds those neurons whose firing correlates with the presentation of the stimuli

3_calculate_MI.R - calculates the mutual information between neuronal response and the presentation of dot stimuli and grating stimuli

