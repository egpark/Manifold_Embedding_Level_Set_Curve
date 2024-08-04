function [z_est,z_unc]=GPR_est_ok_seis(type,nx,ny,dat_trn,sig_max,err,df_total)
xo=dat_trn(:,1);% x-coordinate
yo=dat_trn(:,2);% y-coordinate
val=dat_trn(:,3);% value
ntrn=size(dat_trn,1);
idx0=(xo-1)*ny+yo;
%% Computations
% Kernel matrix computations (kernel trick)
Sig_ab=ones(nx*ny,ntrn+1);% Sigma_ab
Sig_bb=zeros(ntrn+1);% Sigma_bb; note that the dimension is ndat '+1'
mm=0;

for ii=1:ntrn
    mm=mm+1;
    Sig_ab(:,mm)=Sigma(type,idx0(ii),nx,ny,sig_max,df_total)';
end

Sig_bb(1:end-1,:)=Sig_ab(idx0,:);
Sig_bb(end,1:end-1)=1;
iSig_bb=pinv(Sig_bb,err);
z_est=Sig_ab*(iSig_bb*[val;0]);% value

% Relative uncertainty assessment
z_unc=max(0,1-sum(Sig_ab*iSig_bb.*Sig_ab,2));






