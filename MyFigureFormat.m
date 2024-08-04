function MyFigureFormat(ttl,XLabel,YLabel,ZLabel,xmin,xmax,ymin,ymax,zmin,zmax)

shading flat
camlight
grid on
axis([xmin xmax ymin ymax])
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')

xlabel(XLabel, 'FontSize', 32,'FontWeight','bold')
ylabel(YLabel, 'FontSize', 32,'FontWeight','bold')
zlabel(ZLabel, 'FontSize', 32,'FontWeight','bold')
set(gca,'fontname','times new roman','fontsize',18,'linewidth',2)
ylim([ymin ymax]); xlim([xmin xmax]); %zlim([zmin zmax])
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000])
axis equal
axis tight
set(gca,'dataaspectratio',[1 1 0.1])
% clim([zmin zmax])
% set(gca,'ztick',[zmin:(zmax-zmin)/2:zmax])
% clim([zmin-0.1 zmax+0.1])
% set(gca, 'DataAspectRatio')
% title(ttl, 'FontSize', 32,'FontWeight', 'bold');