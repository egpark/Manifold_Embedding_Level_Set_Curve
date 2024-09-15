# README

## Associated Code for "Bridging Field Observations and Manifold Theory: A Novel Approach to Non-stationary Subsurface Characterization Using Level-Set Curves" by Park et al. (2024)

This repository contains the MATLAB code associated with the manuscript titled "Bridging Field Observations and Manifold Theory: A Novel Approach to Non-stationary Subsurface Characterization Using Level-Set Curves" submitted to *Computers & Geosciences*.

### Repository Contents

1. **Est_z_geodist.m**: Main script for executing the level-set curve analysis.
2. **D_M_T.m**: Script for computing geodesic distances.
3. **GPR_est_ok_seis.m**: Script for computing Gaussian Process Regression (GPR) estimation.
4. **GradientToManifold.m**: Function for interpolating the gradient field.
5. **ImageToGradient.m**: Function for reading images and obtaining features (strikes) for computing gradients.
6. **lgwt.m**: Function for computing Gauss-Legendre abscissae and weights.
7. **Sigma.m**: Script for computing matrices for GPR computation.
8. **MyFigureFormat.m**: Script for formatting figures.
9. **geological_map.png**: Example geological map used in the analysis.

### Requirements

- **MATLAB**: The code is designed to run in MATLAB.
- **Toolboxes Required**:
  1. **Image Processing Toolbox**: Required for working with images in `ImageToGradient.m`.
  2. **Statistics and Machine Learning Toolbox**: Required for working with images in `ImageToGradient.m`.

### Instructions

1. **Directory Setup**: 
   - Ensure that **all files** in the repository are placed in the **same directory**. The main script, `Est_z_geodist.m`, depends on the other files being available in the same directory for successful execution.
   
2. **Execution**:
   - Run the main script `Est_z_geodist.m` to generate the results. This script coordinates the overall analysis, including reading the geological map, computing the geodesic distances, and performing the Gaussian Process Regression (GPR) estimation.
   - Make sure the image file `geological_map.png` is in the same directory, as it will be processed by `ImageToGradient.m` for feature extraction.

3. **Output**:
   - The main script generates various results including figures and data visualizations of the level-set curve analysis and GPR-based estimates.

### Notes

- Before running the code, ensure you have the **required toolboxes** installed in your MATLAB environment. Failure to do so will result in errors during the execution of specific scripts.
- All scripts and functions are interdependent, so make sure no files are missing or renamed to avoid execution errors.

