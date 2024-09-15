function [z_est,z_unc]=GPR_est_ok_seis(type,nx,ny,dat_trn,sig_max,err,df_total)
% GPR_est_ok_seis performs geodesic kernel-based Gaussian Process Regression (GPR)
% to estimate values and uncertainty at given locations based on training data.

% Inputs:
% - type: Type of GPR (0 = conventional, 1 = geodesic-based)
% - nx, ny: Number of grid points along x and y (dimensions of the domain)
% - dat_trn: Training data (x, y coordinates and observed values)
% - sig_max: Maximum correlation scale (hyperparameter for GPR kernel)
% - err: Error variance (used in regularization during matrix inversion)
% - df_total: Interpolated gradient data used for geodesic distance calculation (optional)

% Outputs:
% - z_est: Estimated values at each grid point
% - z_unc: Estimated uncertainty (variance) at each grid point

xo=dat_trn(:,1); % Extract x-coordinates of training points
yo=dat_trn(:,2); % Extract y-coordinates of training points
val=dat_trn(:,3); % Extract observed values at training points
ntrn=size(dat_trn,1); % Number of training data points
idx0=(xo-1)*ny+yo; % Index mapping from 2D coordinates to 1D

%% Step 1: Initialize kernel matrices
% Kernel matrix holds the covariance between points. 
% Sig_ab stores covariances between each grid point and the training data.
Sig_ab=ones(nx*ny,ntrn+1); % Initialize covariance matrix between all points and training data
Sig_bb=zeros(ntrn+1); % Initialize covariance matrix for the training data

mm=0;

%% Step 2: Compute Kernel Matrix for Geodesic GPR
% Fill Sig_ab with kernel values (similarity between training points and grid points)
for ii=1:ntrn
    mm=mm+1;
    Sig_ab(:,mm)=Sigma(type,idx0(ii),nx,ny,sig_max,df_total)'; 
    % Compute kernel function (Sigma) between the ii-th training point and all grid points
end

%% Step 3: Build Full Covariance Matrix for Training Data
Sig_bb(1:end-1,:)=Sig_ab(idx0,:); % Covariance between training points
Sig_bb(end,1:end-1)=1; % Add bias term to ensure positive definiteness in the covariance matrix

%% Step 4: Matrix Inversion with Regularization
% Invert the covariance matrix (with regularization for numerical stability)
iSig_bb=pinv(Sig_bb,err); % Compute pseudo-inverse with error tolerance `err`

%% Step 5: Compute Estimated Values
z_est=Sig_ab*(iSig_bb*[val;0]); % Compute the estimated values using the GPR formula

%% Step 6: Compute Uncertainty (Variance)
% Calculate the relative uncertainty of the estimates
z_unc=max(0,1-sum(Sig_ab*iSig_bb.*Sig_ab,2)); 
% Uncertainty is derived from the diagonal elements of the covariance matrix
