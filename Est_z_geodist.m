clear
close all

%% Parameters
inputIMG='geological_map.png';% Case 1
seed=456; % random seed number
Io=imread(inputIMG);

nx=size(Io,2); ny=size(Io,1);
type1=0; % =0 is conventional GPR
sig_max1=40; % correlation scale for interpolating df/dx df/dy
err1=1e-1; % error variance for interpolating df/dx df/dy
discon=3; % =0 ignore discontinuity, =1 considering discon, =3 no discontinuity image

ntrn=50; % number of training 
sig_max2=100; % 80 opt- case04opti 135.23
len=10; % the value multiply to the dfdx_values from image 
err2=1e-1; % 1e-1 case04opti 0.01

%% Software options
type2=1; % =0 is conventional GPR ~0 is geodesic based GPR
cond=0; % =0 is unconditional est. =1 is conditional est.

%% Loading gradient data from input image (building manifold)
[dfdx_grid,dfdy_grid]=GradientToManifold(Io,type1,sig_max1,err1,len);
df_grid=cat(3,dfdx_grid, dfdy_grid);

%% Conditional or unconditional?
if cond==0 % unconditional 
    fn="uncon";
    rng(seed)

    aa=ceil(rand(ntrn,2).*[nx ny]);
    bb=randn(ntrn,1);
    dat_trn=[aa(:,1) aa(:,2) bb];
else %conditioning on input data
    case_name=inputIMG(1:6);
    fn="con";
    z_est=load(case_name+"_z_est_uncon.txt");
    res=25;
    ns=ny/res;
    yw=(res:res:ny)';
    dat_trn=nan(ns*3,3);
    nn=0;
    for kk=50:(nx-100)/2:nx-50
        nn=nn+1;
        xw=ones(ns,1)*kk;
        d=z_est(yw,kk);
        dat_trn((nn-1)*ns+1:nn*ns,:)=[xw yw d];
    end
    ntrn=size(dat_trn,1);
end

%% Running geodesic kernel-based GPR
[z_est,z_unc]=GPR_est_ok_seis(type2,nx,ny,dat_trn,sig_max2,err2,df_grid);

z_est=reshape(z_est,ny,nx);
z_unc=reshape(z_unc,ny,nx);

%% Drawing Figure for estimated Z
x = 1:1:nx;
y = 1:1:ny;
figure('position',[250 250 800 700],'color','w')
contourf(x,y,z_est,20,'LineColor','none')
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
title('Hydraulic Conductivity Field (1)','fontsize',32)

figure('position',[250 250 800 700],'color','w')
surf(x, y, z_est,'FaceAlpha',0.8)
axis equal
axis tight
box on
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')
hold on
plot3(dat_trn(:,1),dat_trn(:,2),dat_trn(:,3),'ko','markersize',8,'markerfacecolor','w','linewidth',2)
set(gca,'CameraPosition', [154.5000  109.5000  194.5109])
xlabel('\itu \rm\bf(pixel)','fontweight','bold','fontsize',32)
ylabel('\itv \rm\bf(pixel)','fontweight','bold','fontsize',32)
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000])
set(gca,'dataaspectratio',[1 1 0.1])
shading flat
camlight
title('Hydraulic Conductivity Field (2)','fontsize',32)

figure('position',[250 250 800 700],'color','w')
contourf(x,y,z_unc,10,'LineColor','none')
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

figure('position',[250 250 800 700],'color','w')
contourf(x,y,reshape(dfdx_grid,ny,nx),20,'LineColor','none')
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

figure('position',[250 250 800 700],'color','w')
contourf(x,y,reshape(dfdy_grid,ny,nx),20,'LineColor','none')
% hold on
% plot(dat_trn(:,1),dat_trn(:,2),'ko','markersize',8,'markerfacecolor','w','linewidth',2)
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

