function [Sig,dist]=Sigma(type,idx,nx,ny,sig_max,df_total)

y0=double(mod(idx,ny));
x0=double(ceil(idx/ny));
if y0==0,y0=ny;end
[xx,yy]=meshgrid(1:nx,1:ny);
xx=xx(:);
yy=yy(:);

if type==0
    dist=sqrt((xx-x0).^2+(yy-y0).^2);
else
    dist=D_M_T(nx,ny,x0,y0,xx,yy,df_total(:,:,1),df_total(:,:,2));
end

Sig=exp(-1/2*dist.^2/sig_max^2);
