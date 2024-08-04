function grad_dat=ImageToGradient(Io)


I=rgb2gray(Io);% converting to grayscale
% I=imresize(I,[500 500]);
[imH,imW]=size(I);% obtaining image size

%% Set window 
% Window is a fundamental unit of examining input image
min_dim=min(imH,imW);% determining the smaller dimension
win_siz=ceil(0.025*min_dim);% window size
num_win=floor(imH/win_siz)*floor(imW / win_siz)*2;% number of windows

%% Generate sample points
% Note that quasi-random sampling method of Sobol sequence is used (can be
% substituted by regular grid or pure random methods)
% Note also that the range of sampling points is from win_siz/2 to 
% imH-win_siz/2 or imW-win_siz/2
sobol_seq=sobolset(2,'leap',1000);% initialization
samples=net(sobol_seq,num_win);% generating
x_sam=round(win_siz/4+samples(:,2)*(imW-win_siz/2));% scaling based on I
y_sam=round(win_siz/4+samples(:,1)*(imH-win_siz/2));% scaling based on I

%% Sample lineament
% Here, two different methods of Canny edge detection and Hough
% transformation were integratively applied
grad_dat=nan(num_win,4);
mm=0;
for ii=1:num_win
    xx=x_sam(ii);% center of sampling window x
    yy=y_sam(ii);% center of sampling window y
    
    % Sampling interval
    y_fm=round(max(1,yy-win_siz/2));
    y_to=round(min(imH,yy+win_siz/2));
    x_fm=round(max(1,xx-win_siz/2));
    x_to=round(min(imW,xx+win_siz/2));
    
    I_win=I(y_fm:y_to,x_fm:x_to);% image within the window

    % Perform edge detection on the window
    edges=edge(I_win,'canny');% Canny edge detection

    % Compute the Hough transform of the edge-detected image
    [H,theta,rho]=hough(edges,'rhoresolution',1);% <------RhoResolution?

    % Find peaks in the Hough transform
    peaks=houghpeaks(H,100,'threshold',ceil(0.25 * max(H(:)))); % only get one peak

    % Extract line segments from the peaks
    min_len=(x_to-x_fm)/2;% minimum lineament length
    lines=houghlines(edges,theta,rho,peaks,...
        'fillgap',2,'minlength',min_len); % win_siz/2

    % Accumulate sum of start and end points of all detected lines
    slp=[];
    for jj=1:length(lines)
        tmp=lines(jj).point1-lines(jj).point2;
        slp=[slp;tmp(2)/tmp(1)];% obtaining multiple slopes
    end
    
    % If at least one lineament is detected, obtain the representative
    if ~isempty(lines) && length(slp)>1
        mm=mm+length(lines);% <--- for debugging
        slp=median(slp);% selecting median slopes as a robust measure
        angles = atan(slp);

        dx = cos(angles);
        dy = sin(angles);

        grad_dat(ii,:)= [xx yy -dy dx];
    end
end
grad_dat=grad_dat(~isnan(grad_dat(:,4)),:);

%% draw the detected lines and its directions
figure('color','w','position',[200 200 800 700])
imagesc(imcomplement(Io))
hold on
quiver(grad_dat(:,1), grad_dat(:,2), grad_dat(:,4), -grad_dat(:,3), 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0);
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


