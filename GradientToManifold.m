function [dfdx_grid,dfdy_grid]=GradientToManifold(Io,type,sig_max,epsil,len)
% GradientToManifold computes interpolated gradient fields (df/dx and df/dy)
% from an input image, using geodesic kernel-based Gaussian Process Regression.

% Inputs:
% - Io: Input image matrix (geological map)
% - type: Type of GPR (0 = conventional, 1 = geodesic-based)
% - sig_max: Maximum correlation scale for GPR
% - epsil: Error variance for GPR
% - len: Scaling factor for gradients

nx=size(Io,2); ny=size(Io,1); % Get dimensions of the image (nx: width, ny: height)

% Step 1: Extract gradients from the input image
imgout_sel=ImageToGradient(Io); 
% Compute gradients (df/dx, df/dy) from the input image using ImageToGradient function

% Step 2: Assign computed gradients to appropriate variables
data_points(:,1) = imgout_sel(:,1); % X-coordinates of gradient points
data_points(:,2) = ny+1-imgout_sel(:,2); % Y-coordinates (inverted for image coordinates)
dfdx_values=imgout_sel(:,3); % Gradient in x-direction (df/dx)
dfdy_values=-imgout_sel(:,4); % Gradient in y-direction (df/dy) (negated to match orientation)

% Step 3: Multiply gradient values by a scaling factor
dfdx_values=dfdx_values*len; % Scale df/dx values
dfdy_values=dfdy_values*len; % Scale df/dy values

%% Step 4: Interpolation of the gradient field
% Interpolate both df/dx and df/dy values using Gaussian Process Regression
% Perform interpolation using the GPR_est_ok_seis function for scattered data
dfdx_grid=GPR_est_ok_seis(type,nx,ny,[data_points dfdx_values],sig_max,epsil,[]); 
dfdy_grid=GPR_est_ok_seis(type,nx,ny,[data_points dfdy_values],sig_max,epsil,[]); 
