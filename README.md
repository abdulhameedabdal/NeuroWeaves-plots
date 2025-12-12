# NeuroWeaves MATLAB Code

MATLAB scripts for processing, analyzing, and visualizing electrophysiological data used in the **NeuroWeaves** project. The code supports recordings acquired with **Intan systems** and **custom wireless EEG hardware**, and was developed for acute and chronic neural recordings from minimally invasive, suturable microthread electrodes.

## Overview
This repository includes scripts for:
- Loading and appending multi-file recordings
- Signal preprocessing (downsampling, filtering, artifact handling)
- Trial-based evoked response analysis
- Time–frequency analysis (PSD, spectrograms, scalograms)
- Latency and amplitude quantification
- Figure generation for main and supplementary figures

## Data
Raw data are not included. Scripts assume access to:
- Intan `.rhs` or exported `.csv` files
- Wireless EEG `.csv` recordings

Update file paths, sampling rates, and channel indices to match your dataset.

## Usage
Scripts are modular and can be run independently once required variables are loaded. A typical workflow is:
1. Load and preprocess raw data
2. Segment trials and compute metrics
3. Generate figures

## Requirements
- MATLAB (R2021a or later)
- Signal Processing Toolbox

## Notes
This code is provided to support transparency and reproducibility of the NeuroWeaves study. It is intended for research use and may require minor adaptation for different hardware or experimental protocols.

## Citation
If you use this code, please cite the associated NeuroWeaves publication.

