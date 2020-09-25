clc,clear%,close all

% addpath('../.');
addpath(genpath('./function/'));  % Add funtion path with sub-folders
addpath(genpath('./data/'));  % Add funtion path with sub-folders

% load inner part data
% load('20170604-3.mat');
%load outer part data
load('20170604-1-400.mat');
%load outerpart bio-data
% load('20170531-3-300.mat');

%load outersider part bio-data
% load('20170531-3-700600.mat');
global k0; global LEDgap;global NA_led;global arraysize; global image_center
m1 = 128;n1 = 128;%the size of patch
% system parameter 
waveLength = 0.6297e-6;
k0 = 2*pi/waveLength; 
NA = 0.1;
mag = 4.0;       %the magnification of objective
mag_image = 4;
pix = 6.5e-6;    %pixel size of the CCD
spsize = pix/mag; % sampling pixel size of the image plane with no magnification
psize = spsize/mag_image; % sampling pixel size of the object
m = m1*mag_image;n = n1*mag_image; % image size of the high resolution object

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
% image_center = 2.4.*[550 -750]*spsize*1e3+[0 0];  % true value around 2.4
%image_center = 1.0.*[750 -550]*spsize*1e3;
%for outsider biodata
image_center = [700 600]*spsize*1e3;    % 2.28 is around the true value

% for outside biodata
image_center = [300 300]*spsize*1e3;    % true value 2.5

[kx,ky]=fpmangle(arraysize,xshift,yshift,theta,LEDheight,image_center);


NA_sys = NA + max(NA_led);
loop=10;
   
vlb=[-4;-4;-0*pi/180;113.5];
vub=[4;4;0*pi/180;113.5];
led_init=[0;0;0*pi/180;113.5];           
%[ledfit,fval,exitflag,output]=fmincon(@(zz)pie_solverexpv4(zz,imSeqLowRes,m,n,arraysize,m1,n1,dkx,dky,CTF,H,loop),led_init,[],[],[],[],vlb,vub);

  
%    opts.poscalibrate = 0;
opts.calbratetol = 1e-3;
poscost = @(zz) pie_solverexpv4(zz,imSeqLowRes,m,n,arraysize,m1,n1,dkx,dky,CTF,H,loop);
%   optsanneal = saoptimset('Display','off','TolFun',opts.calbratetol,'MaxIter',100);
optsanneal = saoptimset('Display','off','TolFun',opts.calbratetol,'MaxIter',100);
[ledfit,fval11, exitflag11, output] = simulannealbnd(poscost,led_init,vlb,vub,optsanneal);
%         cen_corr(:,i2)=cen1;

% display the output image  
loop=2;
[obj_pie,pupil_pie] = pie_solverexpv5(ledfit,imSeqLowRes,m,n,arraysize,m1,n1,dkx,dky,CTF,H,loop);

figure(3)
subplot(221) ;imagesc(abs(obj_pie)); colormap(gray);title('Amplitude');  
subplot(222) ;imagesc(angle(obj_pie)); colormap(gray);title('Phase');
subplot(223) ;imagesc(abs(pupil_pie)); colormap(gray);title('Pupil Amplitude');
subplot(224) ;imagesc(angle(pupil_pie)); colormap(gray);title('Pupil Phase');






