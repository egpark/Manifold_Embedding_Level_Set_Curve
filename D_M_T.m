function dist=D_M_T(nx,ny,u0,v0,u1,v1,df_du_tot,df_dv_tot)
% D_M_T computes the geodesic distance between two points in a spatial domain
% using Legendre-Gauss quadrature for integration and interpolated gradient fields.

% Inputs:
% - nx, ny: Dimensions of the grid (number of points along x and y directions)
% - u0, v0: Coordinates of the starting point
% - u1, v1: Coordinates of the end point
% - df_du_tot: Interpolated gradient field in u direction (df/du)
% - df_dv_tot: Interpolated gradient field in v direction (df/dv)

% Outputs:
% - dist: Computed geodesic distance between points (u0, v0) and (u1, v1)

nd=5; % Number of segments for the Legendre-Gauss quadrature

%% Step 1: Legendre-Gauss Quadrature Nodes and Weights
% This method is used to perform numerical integration over a specified range
[absc,wght]=lgwt(nd,-1,1); % Get nodes (abscissae) and weights for Legendre-Gauss quadrature

%% Step 2: Parametrization of the curve
% Define lambda functions that interpolate between the start (u0, v0) and end (u1, v1) points
ulam=@(t) u0+(u1-u0)*t'; % Linear interpolation for u based on parameter t
vlam=@(t) v0+(v1-v0)*t'; % Linear interpolation for v based on parameter t
t=0.5*(absc+1); % Shift abscissae to range [0, 1]
ult=ulam(t); % Parametric u values along the path (discretized)
vlt=vlam(t); % Parametric v values along the path (discretized)

%% Step 3: Interpolation of gradient fields at intermediate points
% Interpolate the gradient fields (df/du and df/dv) at the parametric points (ult, vlt)
df_du=interp2(reshape(u1,ny,nx),reshape(v1,ny,nx),reshape(df_du_tot,ny,nx),ult(:),vlt(:)); 
% Interpolate df/du at the parametric points

df_dv=interp2(reshape(u1,ny,nx),reshape(v1,ny,nx),reshape(df_dv_tot,ny,nx),ult(:),vlt(:)); 
% Interpolate df/dv at the parametric points

%% Step 4: Metric tensor components
% Compute the metric tensor components (g11, g12, g22) for the geodesic distance calculation
g11=1+reshape(df_du,nx*ny,nd).^2; % First metric tensor component (g11)
g12=reshape(df_du,nx*ny,nd).*reshape(df_dv,nx*ny,nd); % Cross-term (g12)
g22=1+reshape(df_dv,nx*ny,nd).^2; % Second metric tensor component (g22)

%% Step 5: Geodesic distance calculation
% Compute the geodesic distance using the metric tensor components and the Legendre-Gauss weights
dist=sqrt((ult-u0).^2.*g11 + 2*(ult-u0).*(vlt-v0).*g12 + (vlt-v0).^2.*g22)*wght; 
% Apply quadrature to compute the total distance based on the metric tensor and path parametrization

end
