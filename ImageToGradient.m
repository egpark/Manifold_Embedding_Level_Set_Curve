function grad_dat=ImageToGradient(Io)
% ImageToGradient extracts gradient directions (strikes) from an input image 
% using edge detection and Hough transformation.

% Inputs:
% - Io: Input image matrix (geological map)
% Outputs:
% - grad_dat: Gradient data, containing x, y coordinates and gradient directions

% Step 1: Convert input image to grayscale
I=rgb2gray(Io); % Convert RGB image to grayscale for processing

[imH,imW]=size(I); % Get the image dimensions (height and width)

%% Step 2: Set window size for image sampling
% A sampling window is used to analyze local features in the image
min_dim=min(imH,imW); % Identify the smaller dimension of the image
win_siz=ceil(0.025*min_dim); % Set window size as 2.5% of the smaller dimension
num_win=floor(imH/win_siz)*floor(imW / win_siz)*2; % Calculate the number of windows

%% Step 3: Generate sample points using Sobol sequence
% Sobol sequence is a quasi-random sampling method to select sample points
sobol_seq=sobolset(2,'leap',1000); % Initialize the Sobol sequence generator
samples=net(sobol_seq,num_win); % Generate quasi-random sample points
x_sam=round(win_siz/4+samples(:,2)*(imW-win_siz/2)); % Scale x-coordinates of the samples
y_sam=round(win_siz/4+samples(:,1)*(imH-win_siz/2)); % Scale y-coordinates of the samples

%% Step 4: Perform edge detection and Hough transformation for line detection
grad_dat=nan(num_win,4); % Initialize output gradient data
mm=0;
for ii=1:num_win
    xx=x_sam(ii); % x-coordinate of the sampling window center
    yy=y_sam(ii); % y-coordinate of the sampling window center
    
    % Define the boundaries of the sampling window
    y_fm=round(max(1,yy-win_siz/2));
    y_to=round(min(imH,yy+win_siz/2));
    x_fm=round(max(1,xx-win_siz/2));
    x_to=round(min(imW,xx+win_siz/2));
    
    I_win=I(y_fm:y_to,x_fm:x_to); % Extract the image within the sampling window

    % Apply Canny edge detection on the sampled window
    edges=edge(I_win,'canny'); % Detect edges using the Canny method

    % Compute the Hough transform to detect lines in the edge-detected image
    [H,theta,rho]=hough(edges,'rhoresolution',1); % Perform Hough transform

    % Find peaks in the Hough transform (potential line candidates)
    peaks=houghpeaks(H,100,'threshold',ceil(0.25 * max(H(:)))); % Identify Hough peaks

    % Extract line segments based on the detected peaks
    min_len=(x_to-x_fm)/2; % Minimum line length for valid detection
    lines=houghlines(edges,theta,rho,peaks,'fillgap',2,'minlength',min_len); 

    % Step 5: Compute the gradient directions (slopes) from the detected lines
    slp=[];
    for jj=1:length(lines)
        tmp=lines(jj).point1-lines(jj).point2; % Calculate line direction
        slp=[slp;tmp(2)/tmp(1)]; % Compute slope of the line
    end
    
    % If at least one valid line is detected, calculate median slope
    if ~isempty(lines) && length(slp)>1
        mm=mm+length(lines); % Debugging info: counts number of valid lines
        slp=median(slp); % Use median slope as a robust estimate
        angles = atan(slp); % Calculate the angle from the slope

        % Convert angles to directional vectors (dx, dy)
        dx = cos(angles);
        dy = sin(angles);

        grad_dat(ii,:)= [xx yy -dy dx]; % Store the coordinates and gradients (dy, dx)
    end
end
grad_dat=grad_dat(~isnan(grad_dat(:,4)),:); % Remove empty rows

%% Step 6: Plot detected lines and their directions (quiver plot)
figure('color','w','position',[200 200 800 700])
imagesc(imcomplement(Io)) % Display the inverted original image
hold on
quiver(grad_dat(:,1), grad_dat(:,2), grad_dat(:,4), -grad_dat(:,3), 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0); 
% Overlay gradient directions (strikes) using quiver plot
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')
xlabel('x (pixel)','fontweight','bold','fontsize',32)
ylabel('y (pixel)','fontweight','bold','fontsize',32)
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
ny=floor(imH/100)*100;
set(gca,'ytick',[imH-ny:50:imH],'yticklabel',fliplr([0 50:50:ny]))
grid on
axis equal
axis tight
title('Boundary Identification Results (Strikes)','fontsize',32)
