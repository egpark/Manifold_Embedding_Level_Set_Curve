clear
close all

%% Parameters
inputIMG='geological_map.png'; % Input image file, used for geological map
seed=456; % Random seed for reproducibility
Io=imread(inputIMG); % Load the input image

% Image dimensions
nx=size(Io,2); 
ny=size(Io,1);

type1=0; % =0: Conventional Gaussian Process Regression (GPR)
sig_max1=40; % Correlation scale for gradient interpolation
err1=1e-1; % Error variance for gradient interpolation
discon=3; % Discontinuity parameter (0: ignore, 1: consider, 3: no discontinuity)

ntrn=50; % Number of training data points
sig_max2=100; % Correlation scale for geodesic-based GPR
len=10; % Scaling factor for gradients from image
err2=1e-1; % Error variance for GPR

%% Software options
type2=1; % =0: Conventional GPR, =1: Geodesic kernel-based GPR
cond=0; % =0: Unconditional estimation, =1: Conditional estimation

%% Loading gradient data from input image (building manifold)
[dfdx_grid,dfdy_grid]=GradientToManifold(Io,type1,sig_max1,err1,len); 
% Computes gradients (df/dx, df/dy) from the image to construct the manifold

df_grid=cat(3,dfdx_grid, dfdy_grid); 
% Combines gradients into a single 3D array

%% Conditional or unconditional estimation
if cond==0 % Unconditional estimation
    fn="uncon"; % Output file name indicator
    rng(seed) % Set random seed for reproducibility

    aa=ceil(rand(ntrn,2).*[nx ny]); % Randomly generate training point coordinates
    bb=randn(ntrn,1); % Generate random training data
    dat_trn=[aa(:,1) aa(:,2) bb]; % Training data: x, y, z (random values)
else % Conditional estimation
    case_name=inputIMG(1:6); % Extract case name from input image
    fn="con"; % Output file name indicator
    z_est=load(case_name+"_z_est_uncon.txt"); % Load pre-existing estimation data
    res=25;
    ns=ny/res;
    yw=(res:res:ny)'; % Sample points along y-direction
    dat_trn=nan(ns*3,3); % Initialize training data
    nn=0;
    for kk=50:(nx-100)/2:nx-50 % Loop to fill training data with conditional values
        nn=nn+1;
        xw=ones(ns,1)*kk; % X-coordinates
        d=z_est(yw,kk); % Z-values from estimation
        dat_trn((nn-1)*ns+1:nn*ns,:)=[xw yw d]; % Store in training data
    end
    ntrn=size(dat_trn,1); % Update number of training points
end

%% Running geodesic kernel-based GPR
[z_est,z_unc]=GPR_est_ok_seis(type2,nx,ny,dat_trn,sig_max2,err2,df_grid); 
% Estimate the field (z_est) and its uncertainty (z_unc) using geodesic-based GPR

z_est=reshape(z_est,ny,nx); % Reshape the estimation to match image dimensions
z_unc=reshape(z_unc,ny,nx); % Reshape uncertainty estimation

%% Drawing Figure for estimated Z
x = 1:1:nx;
y = 1:1:ny;
figure('position',[250 250 800 700],'color','w') % Create a figure for visualization
contourf(x,y,z_est,20,'LineColor','none') % Plot the estimated field as a filled contour plot
hold on
plot(dat_trn(:,1),dat_trn(:,2),'ko','markersize',8,'markerfacecolor','w','linewidth',2) 
% Overlay training points on the plot
hc=colorbar;
axis equal
axis tight
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')
xlabel('\itu \rm\bf(pixel)','fontweight','bold','fontsize',32)
ylabel('\itv \rm\bf(pixel)','fontweight','bold','fontsize',32)
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000])
set(hc,'linewidth',2)
title('Hydraulic Conductivity Field (1)','fontsize',32)

figure('position',[250 250 800 700],'color','w')
surf(x, y, z_est,'FaceAlpha',0.8) % 3D surface plot of the estimated field
axis equal
axis tight
box on
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')
hold on
plot3(dat_trn(:,1),dat_trn(:,2),dat_trn(:,3),'ko','markersize',8,'markerfacecolor','w','linewidth',2) 
% Overlay training points in 3D
set(gca,'CameraPosition', [154.5000  109.5000  194.5109])
xlabel('\itu \rm\bf(pixel)','fontweight','bold','fontsize',32)
ylabel('\itv \rm\bf(pixel)','fontweight','bold','fontsize',32)
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000])
set(gca,'dataaspectratio',[1 1 0.1])
shading flat
camlight
title('Hydraulic Conductivity Field (2)','fontsize',32)

%% Plotting uncertainty (variance)
figure('position',[250 250 800 700],'color','w')
contourf(x,y,z_unc,10,'LineColor','none') % Plot uncertainty as filled contour
hold on
plot(dat_trn(:,1),dat_trn(:,2),'ko','markersize',8,'markerfacecolor','w','linewidth',2)
hc=colorbar;
axis equal
axis tight
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')
xlabel('\itu \rm\bf(pixel)','fontweight','bold','fontsize',32)
ylabel('\itv \rm\bf(pixel)','fontweight','bold','fontsize',32)
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000])
set(hc,'linewidth',2)
title('Variance (Uncertainty)','fontsize',32)

%% Plotting interpolated gradients (df/dx)
figure('position',[250 250 800 700],'color','w')
contourf(x,y,reshape(dfdx_grid,ny,nx),20,'LineColor','none') % Plot interpolated df/dx
hc=colorbar;
axis equal
axis tight
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')
xlabel('\itu \rm\bf(pixel)','fontweight','bold','fontsize',32)
ylabel('\itv \rm\bf(pixel)','fontweight','bold','fontsize',32)
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000])
set(hc,'linewidth',2)
title('Interpolated Gradients (df/du)_e_s_t','fontsize',32)

%% Plotting interpolated gradients (df/dy)
figure('position',[250 250 800 700],'color','w')
contourf(x,y,reshape(dfdy_grid,ny,nx),20,'LineColor','none') % Plot interpolated df/dy
hc=colorbar;
axis equal
axis tight
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')
xlabel('\itu \rm\bf(pixel)','fontweight','bold','fontsize',32)
ylabel('\itv \rm\bf(pixel)','fontweight','bold','fontsize',32)
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000])
set(hc,'linewidth',2)
title('Interpolated Gradients (df/dv)_e_s_t','fontsize',32)
