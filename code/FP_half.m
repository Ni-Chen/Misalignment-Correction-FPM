%{
----------------------------------------------------------------------------------------------------
Name: Fourier ptychography

Author:   Ni Chen (ni_chen@163.com)
Date:     Nov. 2016

Reference:
- G. Zheng, Fourier Ptychographic Imaging- A MATLAB® tutorial
 
Copyright (2016): Ni Chen (ni_chen@163.com)
This program is free software; you can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation; either version 2 of the License
, or (at your option) any later version.
----------------------------------------------------------------------------------------------------
%}

close all;
clear;
clc;

% FT = @(x) fftshift(fft2(ifftshift(x)));
FT = @(x) ifftshift(fft2(fftshift(x)));
iFT = @(x) ifftshift(ifft2(fftshift(x)));

%% Simulate the forward imaging process of Fourier ptychography
expName = 'sim_cameraman';

% file directory
addpath('./function/');
in_dir = './data/';

% out_dir = ['./data/', expName, '/'];
out_dir = './output/';

type = 'complex';
% type = 'amp';
% type = 'phase';
% Simulate the high resolution complex object
objAmp = double((imread([in_dir, 'lena.jpg'])));
phase = double((imread([in_dir, 'baboon.jpg'])));     % baboon, lena
% objAmp = Toxy(objAmp, 0, 1);
phase = Toxy(phase, 0, 0.5*pi);
object = objAmp.*exp(1i*phase);
% object = objAmp;
[Ny, Nx] = size(object);    % image size of the high resolution object

%% System parameters of the microscopy imaging system
N_led = 15;        % size of the LED array:15x15
LEDgap = 4;        % 4mm between adjacent LEDs
LEDheight = 90;    % 90mm between the LED matrix and the sample

pps = 2.75e-6;     % pixel size of the CCD
ppo_r = pps/6;     % final pixel size of the reconstructed object;
NA = 0.1;          % Numerical aperture
lambda = 0.532e-6;
z = 0e-3;

% Parameters of the hologram
pph = 8e-6;        % pixel size of the hologram (Fourier spectrum plane)

iterNum = 3; 
%%
k0 = 2*pi/lambda;

%% define propagation transfer function
[u, v] = meshgrid(((1:Ny)-Ny/2)*pph, ((1:Nx)-Nx/2)*pph);

% Fresnel, object defocus distance
H0 = exp(1i*2*pi/lambda*z)*exp(-1i*pi*lambda*z*(u.^2+v.^2));
% OR angular spectrum
% H0 = exp(1i*2*pi*sqrt((1/lambda^2-u.^2-v.^2).*double(sqrt(u.^2+v.^2)<1/lambda))*z);

% create the wave vectors  for the LED illumination
xlocation = zeros(1, N_led^2);
ylocation = zeros(1, N_led^2);

% Generate spatial positions of the LED elements
for i = 1:N_led  % From the top left to the bottom right
    xlocation(1, 1+N_led*(i-1):N_led+N_led*(i-1)) = (-(N_led-1)/2:1:(N_led-1)/2)*LEDgap; 
    ylocation(1, 1+N_led*(i-1):N_led+N_led*(i-1)) = ((N_led-1)/2-(i-1))*LEDgap; 
end

% Incident wave vectors for the 15x15 elements, assuming the object is located at (0,0)
kx_relative = -sin(atan(xlocation/LEDheight));
ky_relative = -sin(atan(ylocation/LEDheight));

% Magnification of the system
M = pps/ppo_r;

% image size of the CCD output
Nyl = round(Ny/M);
Nxl = round(Nx/M);

% imSeqLowRes = zeros(Nyl, Nxl, N_led^2); % output low-res image sequence
kx = k0*kx_relative;
ky = k0*ky_relative;

dkx = 2*pi/(ppo_r*Nx);
dky = 2*pi/(ppo_r*Ny);

cutoffFrequency = NA*k0;
kmax = pi/pps;
[kxm, kym] = meshgrid(-kmax:kmax/((Nxl-1)/2):kmax, -kmax:kmax/((Nyl-1)/2):kmax);
CTF = ((kxm.^2 + kym.^2) < cutoffFrequency^2);

sub_low_img = zeros(Nyl, Nxl, N_led^2);  % output low-res image sequence

%% Generating images according to the equation: G_out(kx, ky) = H_coh(kx,ky)G_obj(kx-kxn, ky-kyn)
for iImg = 1:N_led^2
        % Generate low resolution Fourier spectrum
        kxc = round((Nx+1)/2 + kx(1,iImg)/dkx);
        kyc = round((Ny+1)/2 + ky(1,iImg)/dky);

        kyl = round(kyc-(Nyl-1)/2);
        kyh = round(kyc+(Nyl-1)/2);
        kxl = round(kxc-(Nxl-1)/2);
        kxh = round(kxc+(Nxl-1)/2);
        
        objectFT = FT(object);
        
        imSeqLowFT = (Nyl/Ny)^2*objectFT(kyl:kyh, kxl:kxh).*CTF;  % (Ny_low/Ny)^2 is a scaling factor, to normalize the Fourier magnitude when changing the image size
        sub_low_img(:, :, iImg) = abs(iFT(imSeqLowFT));  % .^2 doesn't work?
        
        imwrite(uint8(sub_low_img(:, :, iImg)), [out_dir, num2str(iImg),'.tif']);
        disp([num2str(iImg/N_led^2*100),'% is finished.....']);
end

%% Recover the high resolution image
seq = gseq(N_led);  % define the order of recovery ,we start from the center(the 113th image) to the edge of the spectrum(the 225th image)

% Size of the hologram
Nyh = Ny;
Nxh = Nx;
hologram = FT(ones(Nyh, Nxh));  % initial guess of the spectrum


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for loop = 1:iterNum
    
    for iImg = 1:N_led^2
        iImgSeq = seq(iImg);
        
        % Update sub hologram with the low resolution spectrum
        kxc = round((Nxh+1)/2 + kx(1,iImgSeq)/dkx);
        kyc = round((Nyh+1)/2 + ky(1,iImgSeq)/dky);
        
        kyl = round(kyc-(Nyl-1)/2);
        kyh = round(kyc+(Nyl-1)/2);
        kxl = round(kxc-(Nxl-1)/2);
        kxh = round(kxc+(Nxl-1)/2);
        
        sub_holo = (Nyl/Ny)^2*hologram(kyl:kyh, kxl:kxh).*CTF;  % (Nyl/Ny)^2 is a scaling factor, to normalize the Fourier magnitude when changing the image size
        sub_img = iFT(sub_holo);        
        
%         sub_low = sub_low_img(:, :, iImgSeq);   
%         sub_img_prime = (Ny/Nyl)^2*sub_low.*exp(1i*angle(sub_img));
        
        sub_img_prime = (Ny/Nyl)^2*sub_low_img(:, :, iImgSeq).*exp(1i*angle(sub_img));
        sub_holo_prime = FT(sub_img_prime).*CTF;
        
        hologram(kyl:kyh, kxl:kxh) = (1-CTF).*hologram(kyl:kyh, kxl:kxh) + sub_holo_prime;        
        
        imshow(log(abs(hologram)), []);
        title('Fourier hologram');
        
        disp([num2str(iImg/N_led^2*100),'% is finished.....']);
    end
end

% temp = Toxy(log(abs(hologram)), 0, 255);
% figure
% imshow(log(abs(hologram)), []);

% title('hologram');
% imwrite((temp), [out_dir,'histology_spectrum.png']);


f = 200e-3;
objectRecover = dt_iLensFT(hologram, pph, pph, f, f, lambda);
% objectRecover = iFT(hologram);

reobj_amp = Toxy(abs(objectRecover), 0, 255); 
figure;
imshow(uint8(reobj_amp), []);
title('Reconstructed Obj');
imwrite(uint8(reobj_amp), [out_dir,'sim_fpm_', type, '_amp_whole.png']);
% print('-depsc', [out_dir, 'histology_amp.eps']);

reobj_phase = Toxy(angle(objectRecover),0, 255);
figure
imshow(uint8(reobj_phase), []);
title('Reconstructed Phase');
imwrite(uint8(reobj_phase), [out_dir,'sim_fpm_', type, '_phase_whole.png']);
% print('-depsc', [out_dir, 'histology_phase.eps']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hologram_half = FT(ones(Nyh, Nxh));  % initial guess of the spectrum
figure;
for loop = 1:iterNum
    
    for iImg = 1:N_led^2
        iImgSeq = seq(iImg);
        
        % Update sub hologram with the low resolution spectrum
        kxc = round((Nxh+1)/2 + kx(1,iImgSeq)/dkx);
        kyc = round((Nyh+1)/2 + ky(1,iImgSeq)/dky);
        
        kyl = round(kyc-(Nyl-1)/2);
        kyh = round(kyc+(Nyl-1)/2);
        kxl = round(kxc-(Nxl-1)/2);
        kxh = round(kxc+(Nxl-1)/2);
        
        sub_holo = (Nyl/Ny)^2*hologram_half(kyl:kyh, kxl:kxh).*CTF;  % (Nyl/Ny)^2 is a scaling factor, to normalize the Fourier magnitude when changing the image size
        sub_img = iFT(sub_holo);        
        
        if iImgSeq > N_led^2/2+8   % N_led^2/2+8
            sub_low = sub_low_img(:, :, N_led^2+1-iImgSeq);    % symetric spectrum with half images
            continue;    % use only half captured images
        else            
            sub_low = sub_low_img(:, :, iImgSeq);
        end     
        
        sub_img_prime = (Ny/Nyl)^2*sub_low.*exp(1i*angle(sub_img));
        sub_holo_prime = FT(sub_img_prime).*CTF;
        
        hologram_half(kyl:kyh, kxl:kxh) = (1-CTF).*hologram_half(kyl:kyh, kxl:kxh) + sub_holo_prime;        
        
        imshow(log(abs(hologram_half)), []);
        title('Fourier hologram');
        
        disp([num2str(iImg/N_led^2*100),'% is finished.....']);
    end
end

% temp = Toxy(log(abs(hologram_half)), 0, 255);
% figure
% imshow(log(abs(hologram_half)), []);

% title('hologram');
% imwrite((temp), [out_dir,'sim_fpm_complex_spectrum_half.png']);


f = 200e-3;
objectRecover_half = dt_iLensFT(hologram_half, pph, pph, f, f, lambda);
% objectRecover = iFT(hologram);

reobj_amp_half = Toxy(abs(objectRecover_half), 0, 255); 
figure;
imshow(uint8(reobj_amp_half), []);
title('Reconstructed Obj');
imwrite(uint8(reobj_amp_half), [out_dir,'sim_fpm_', type, '_amp_half.png']);
% print('-depsc', [out_dir, 'histology_amp_half.eps']);

reobj_phase_half = Toxy(angle(objectRecover_half)+pi,0, 255);
figure;
imshow(uint8(reobj_phase_half), []);
title('Reconstructed Phase');
imwrite(uint8(reobj_phase_half), [out_dir,'sim_fpm_', type, '_phase_half.png']);
% print('-depsc', [out_dir, 'histology_phase_half.eps']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compare%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


objAmp = Toxy(objAmp, 0, 255);
figure;
plot(objAmp(256,:), 'r','Linewidth',1.5);
hold on;
plot(reobj_amp(256,:), 'g-.','Linewidth',1.5);
hold on;
plot(reobj_amp_half(256,:), 'b--','Linewidth',1.5);
ylim([0 2^8-1]);
xlim([0 512]);
legend('Original', 'Whole images', 'Half number images');
set(gca,'fontsize',18);
print('-depsc', [out_dir, 'sim_fpm_', type, '_amp_diff.eps']);


reobj_phase_half = Toxy(angle(objectRecover_half), 0, 0.5*pi);
reobj_phase = Toxy(angle(objectRecover), 0, 0.5*pi);

phase = Toxy(phase, 0, 0.5*pi);
figure;
plot(phase(256,:), 'r','Linewidth',1.5);
hold on;
plot(reobj_phase(256,:), 'g-.','Linewidth',1.5);
hold on;
plot(reobj_phase_half(256,:), 'b--','Linewidth',1.5);
ylim([0 0.5*pi]);
xlim([0 512]);
legend('Original', 'Whole images', 'Half number images', 'SouthEast');
set(gca,'fontsize',18);
print('-depsc', [out_dir, 'sim_fpm_', type, '_phase_diff.eps']);

