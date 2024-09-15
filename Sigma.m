function [Sig,dist]=Sigma(type,idx,nx,ny,sig_max,df_total)
% Sigma computes the kernel values (covariance) based on either Euclidean or geodesic distance.

% Inputs:
% - type: Type of kernel (0 = Euclidean, 1 = Geodesic)
% - idx: Index of the reference point (1D index)
% - nx, ny: Grid dimensions (number of grid points in x and y directions)
% - sig_max: Maximum correlation scale (kernel parameter)
% - df_total: Gradient data for geodesic distance calculation (optional, used if type == 1)

% Outputs:
% - Sig: Computed kernel values (covariances) between the reference point and all grid points
% - dist: Distances between the reference point and all grid points

% Step 1: Convert 1D index to 2D coordinates (x0, y0)
y0=double(mod(idx,ny)); % Compute y-coordinate from 1D index
x0=double(ceil(idx/ny)); % Compute x-coordinate from 1D index
if y0==0, y0=ny; end % Correct if y-coordinate is zero (edge case)

% Step 2: Generate meshgrid of all grid points
[xx,yy]=meshgrid(1:nx,1:ny); % Create a grid of all (x, y) coordinates
xx=xx(:); % Flatten the grid's x-coordinates into a column vector
yy=yy(:); % Flatten the grid's y-coordinates into a column vector

%% Step 3: Calculate distances (Euclidean or geodesic)
if type==0
    % Euclidean distance between reference point (x0, y0) and all grid points
    dist=sqrt((xx-x0).^2+(yy-y0).^2); 
else
    % Geodesic distance using precomputed gradients from df_total
    dist=D_M_T(nx,ny,x0,y0,xx,yy,df_total(:,:,1),df_total(:,:,2)); 
    % Call D_M_T to calculate geodesic distances based on the gradient data
end

%% Step 4: Compute the kernel values (covariance matrix)
Sig=exp(-1/2*dist.^2/sig_max^2); % Compute the kernel using the Gaussian RBF function
