clc,clear,close all
% load inner part data
%load('20170604-3.mat');
 %load outer part data
% load('20170604-1-400.mat');
%load outerpart bio-data
% load('20170531-3-300.mat');

%load outersider part bio-data
load('20170531-3-700600.mat');
global k0; global LEDgap;global NA_led;global arraysize; global image_center
m1 = 128;n1 = 128;%the size of patch
% system parameter 
waveLength = 0.6297e-6;
k0 = 2*pi/waveLength; pix = 6.5e-6;    %pixel size of the CCD
NA = 0.1;
mag = 4.0;       %the magnification of objective
mag_image = 4;
pix = 6.5e-6;    %pixel size of the CCD
spsize = pix/mag; % sampling pixel size of the image plane with no magnification
psize = spsize/mag_image; % sampling pixel size of the object
m =m1*mag_image;n =n1*mag_image; % image size of the high resolution object

dkx = 2*pi/(psize*n);
dky = 2*pi/(psize*m);
cutoffFrequency = NA * k0;
kmax = pi/spsize;
[kxm, kym] = meshgrid(-kmax:kmax/((n1-1)/2):kmax,-kmax:kmax/((n1-1)/2):kmax);
CTF = double((kxm.^2+kym.^2)<cutoffFrequency^2); % coherent transfer function 

%aberration function
pupil=ones(m1,n1).*CTF;

%% defocus distance
z = 0e-6; % defocus distance
kzm = sqrt(k0^2-kxm.^2-kym.^2);
H = exp(1i.*z.*real(kzm)).*exp(-z.*imag(kzm));

arraysize = 17;
LEDgap = 4; % 4mm between adjacent LEDs
LEDheight = 113.5;% 108 mm between the LED matrix and the sample
theta=0; xshift=0;yshift=0;
% image_center = [0 0]*spsize*1e3;  % pixel numer of the patch center with respect to the raw image center

% for outside data
% image_center = 1.0.*[550 -750]*spsize*1e3+[0 0];  % true value around 2.4
%image_center = 1.0.*[750 -550]*spsize*1e3;
%for outsider biodata
image_center = 2.28.*[700 600]*spsize*1e3;    % 2.28 is around the true value

% for outside biodata
%image_center = 1.5.*[300 300]*spsize*1e3;

[kx,ky]=fpmangle(arraysize,xshift,yshift,theta,LEDheight,image_center);



NA_sys = NA + max(NA_led);

loop=9; %设定循环次数，次数小于10不用SA，大于10用SA
% [obj_pie,pupil_pie] = pie_solverexp(imSeqLowRes,m,n,arraysize,m1,n1,dkx,dky,kx,ky,CTF,H,loop);
[obj_pie,pupil_pie] = pie_solverexpv3(imSeqLowRes,m,n,arraysize,m1,n1,dkx,dky,kx,ky,CTF,H,loop);

figure(3)
subplot(221) ;imagesc(abs(obj_pie)); colormap(gray);
subplot(222) ;imagesc(angle(obj_pie)); colormap(gray);
subplot(223) ;imagesc(abs(pupil_pie)); colormap(gray);
subplot(224) ;imagesc(angle(pupil_pie)); colormap(gray);


% save('E:\我的快盘\matlab\fpm\simulation1\fpmsimu.mat','imSeqLowRes');







