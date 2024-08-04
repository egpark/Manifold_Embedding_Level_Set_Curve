function [dfdx_grid,dfdy_grid]=GradientToManifold(Io,type,sig_max,epsil,len)


nx=size(Io,2); ny=size(Io,1);

imgout_sel=ImageToGradient(Io);
% imgout_sel=ImageToGradient_seis2(inputIMG,nx,ny,win_size);

% assign the gradient dx dy calculated from ImageTOGradient to the dfdx dfdy values
data_points(:,1) = imgout_sel(:,1);
data_points(:,2) = ny+1-imgout_sel(:,2);
dfdx_values=imgout_sel(:,3);
dfdy_values=-imgout_sel(:,4);

% multiply the length
dfdx_values=dfdx_values*len;
dfdy_values=dfdy_values*len;

%% 2. Gradient Field Interpolation
% Perform data interpolation for both dfdx and dfdy using scatteredInterpolant
dfdx_grid=GPR_est_ok_seis(type,nx,ny,[data_points dfdx_values],sig_max,epsil,[]);
dfdy_grid=GPR_est_ok_seis(type,nx,ny,[data_points dfdy_values],sig_max,epsil,[]);




 
