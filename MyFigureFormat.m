function MyFigureFormat(ttl,XLabel,YLabel,ZLabel,xmin,xmax,ymin,ymax,zmin,zmax)
% MyFigureFormat is used to format a MATLAB figure with predefined settings
% for axes limits, labels, and overall appearance.

% Inputs:
% - ttl: Title of the figure (currently not used)
% - XLabel, YLabel, ZLabel: Axis labels for x, y, and z
% - xmin, xmax: Limits for the x-axis
% - ymin, ymax: Limits for the y-axis
% - zmin, zmax: Limits for the z-axis (z-axis limits currently commented out)

%% Step 1: General figure settings
shading flat % Apply flat shading (useful for surface plots)
camlight % Add lighting to the figure for 3D plots
grid on % Display the grid on the figure
axis([xmin xmax ymin ymax]) % Set the axis limits for x and y

% Step 2: Set the font size, line width, and font for the axes
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')

%% Step 3: Axis labels
xlabel(XLabel, 'FontSize', 32, 'FontWeight', 'bold') % X-axis label with bold font
ylabel(YLabel, 'FontSize', 32, 'FontWeight', 'bold') % Y-axis label with bold font
zlabel(ZLabel, 'FontSize', 32, 'FontWeight', 'bold') % Z-axis label with bold font (for 3D plots)

% Step 4: Additional axis formatting
set(gca,'fontname','times new roman','fontsize',18,'linewidth',2) % Set font and size for tick labels
ylim([ymin ymax]); % Set y-axis limits
xlim([xmin xmax]); % Set x-axis limits
% zlim([zmin zmax]); % Set z-axis limits (commented out for now)

%% Step 5: Tick labels and axis scaling
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000]) % Set x-axis tick labels and intervals
set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000]) % Set y-axis tick labels and intervals

%% Step 6: Ensure equal axis scaling and adjust aspect ratio
axis equal % Ensure equal scaling for x and y axes
axis tight % Adjust the axis limits to fit the data tightly
set(gca,'dataaspectratio',[1 1 0.1]) % Set aspect ratio (x:y:z) to 1:1:0.1
