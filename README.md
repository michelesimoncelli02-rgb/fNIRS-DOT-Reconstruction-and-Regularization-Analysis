# fNIRS DOT Reconstruction and Regularization Analysis

MATLAB implementation of a **Diffuse Optical Tomography (DOT)** pipeline for functional Near-Infrared Spectroscopy (**fNIRS**) data, including channel quality assessment, signal preprocessing, image reconstruction, and comparison of different regularization strategies for solving the optical inverse problem.

The project was developed as part of the *Imaging for Neuroscience* course at the **University of Padova**. The implemented pipeline reproduces the complete analysis requested in the assignment, from raw fNIRS measurements to cortical hemoglobin reconstructions, with particular emphasis on evaluating the effect of the regularization parameter on image quality.

---

## Project Overview

Diffuse Optical Tomography (DOT) extends conventional fNIRS by reconstructing three-dimensional maps of cortical hemodynamic activity from optical measurements acquired at the scalp. Since the inverse problem is highly ill-posed, image quality strongly depends on the choice of the regularization parameter.

The objective of this project is to investigate how different regularization strengths affect DOT image reconstruction by:

- assessing channel quality,
- preprocessing raw optical signals,
- reconstructing cortical **HbO** and **HbR** concentration changes,
- comparing reconstructed activation maps obtained with different regularization parameters.

The analysis is performed on data acquired from one adult participant performing a **color-naming task**, during which the subject verbally responded according to the colour of visual stimuli.

---

## Analysis Pipeline

The MATLAB script implements the complete workflow described in the assignment.

### 1. Array Configuration

The three-dimensional geometry of the acquisition system is visualized by plotting:

- optical sources,
- optical detectors,
- measurement channels connecting each source-detector pair.

This allows verification of the probe layout before processing the measurements.

---

### 2. Source-Detector Distance Analysis

The Euclidean distance between every source-detector pair is computed and summarized with a histogram.

This provides an overview of the spatial sampling characteristics of the acquisition array and allows inspection of channel-length variability.

---

### 3. Channel Quality Assessment

Signal quality is evaluated for every channel.

Channels are classified as **bad** if they satisfy either of the following criteria:

- mean intensity outside the range **[0.03, 3]**
- signal-to-noise ratio (SNR) below **7**

The resulting binary vector is stored in the `SD.MeasListAct` field and subsequently used throughout the processing pipeline.

The array configuration is then displayed:

- highlighting bad channels,
- showing only the retained channels after quality control.

---

### 4. fNIRS Signal Preprocessing

The optical measurements undergo a standard preprocessing pipeline composed of several sequential steps.

#### Optical Density Conversion

Raw light intensity measurements are converted into optical density changes using the logarithmic Beer-Lambert relationship.

#### Motion Artifact Correction

Motion artifacts are corrected using the **Wavelet Motion Correction** algorithm implemented in Homer2 with:

- **IQR = 0.5**

#### Band-Pass Filtering

The optical density signals are filtered using a band-pass filter with cut-off frequencies:

- **0.01 Hz**
- **0.5 Hz**

to suppress slow baseline drifts and high-frequency physiological noise.

#### Block Averaging

Stimulus-locked hemodynamic responses are computed using a block-average approach over the interval:

- **−2 s to 15 s**

Baseline correction is applied prior to averaging.

---

### 5. Head Mesh Visualization

The anatomical head volume mesh is displayed together with the optical probe configuration.

The visualization overlays:

- head volume mesh,
- sources,
- detectors,

allowing inspection of probe placement relative to the anatomical model.

---

### 6. Whole-Array Sensitivity Analysis

The Jacobian matrix provided with the dataset is used to visualize the sensitivity distribution of the optical array on the grey matter volumetric mesh.

Sensitivity maps are generated for:

- all channels,
- only the channels retained after quality assessment.

This comparison illustrates the impact of removing noisy measurements on spatial sensitivity.

---

### 7. DOT Image Reconstruction

The inverse problem is solved using a regularized inversion of the Jacobian matrix.

Images of:

- oxygenated hemoglobin (**HbO**)
- deoxygenated hemoglobin (**HbR**)

are reconstructed and mapped onto the cortical grey matter surface.

Five regularization parameters are investigated:

| Regularization | λ |
|---------------:|---:|
| Very weak | 0.0001 |
| Weak | 0.01 |
| Moderate | 0.1 |
| Strong | 1 |
| Very strong | 10 |

For each reconstruction, cortical activation maps are generated at:

- **0 s**
- **7 s**
- **15 s**

for both HbO and HbR.

---

### 8. Regularization Comparison

The reconstructed images obtained with different regularization strengths are qualitatively compared.

The analysis highlights the classical trade-off between:

- **under-regularization**, which increases spatial detail but amplifies noise,
- **over-regularization**, which suppresses noise at the cost of spatial resolution.

The implementation suggests that regularization values approximately between **0.01** and **0.1** provide the best compromise between localization accuracy and image smoothness.

---

## Repository Structure

```text
fNIRS-DOT-Reconstruction-and-Regularization-Analysis/
│
├── README.md
├── main.m                          # Main analysis pipeline
│
├── DATASET/                        # Private dataset (not included)
│
├── homer2/                         # Homer2 toolbox
├── iso2mesh/                       # Iso2Mesh toolbox
│
├── *.mat                           # Supporting MATLAB data files
│
└── figures/                        # Optional output figures
```

---

## Dataset

The analysis relies on an experimental fNIRS dataset acquired during a color-naming task, together with anatomical models, Jacobian matrices, mesh files, and additional reconstruction resources.

**The dataset is not included in this repository because it is not publicly distributable.**

It can be shared **upon request**, subject to the original data-sharing policies and permissions.

All supplementary MATLAB functions and `.mat` files required to reproduce the complete analysis are included in the repository. Only the experimental dataset must be obtained separately.

---

## Requirements

The project was developed in **MATLAB** and requires:

- MATLAB
- **Homer2**
- **Iso2Mesh**
- Statistics and Signal Processing Toolboxes (recommended)

The repository includes all auxiliary MATLAB scripts and supporting `.mat` files used throughout the analysis. Only the external toolboxes (Homer2 and Iso2Mesh) and the experimental dataset need to be provided separately.

---

## Outputs

The pipeline produces:

- 3D probe configuration visualizations
- Source-detector distance histogram
- Channel quality assessment
- Optical density signals
- Motion-corrected and filtered fNIRS data
- Block-averaged hemodynamic responses
- Whole-array sensitivity maps
- HbO and HbR cortical reconstructions
- DOT images reconstructed with multiple regularization parameters
- Qualitative comparison of reconstruction quality across regularization strengths

---

## Author

**Michele Simoncelli**  
Department of Information Engineering  
**University of Padova**
