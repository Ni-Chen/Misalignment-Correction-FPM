clc,clear,close all
%% simulate the forward imaging process of Fourier ptychography
% simulate the high resolution complex object
objectAmplitude = sqrt(double(imresize(rgb2gray(imread('4.2.03.tiff')),[511,511])));
% objectAmplitude = ones(256,256);

phase = double(imresize(imread('5.2.09.tiff'),[511,511]));
% phase = zeros(size(objectAmplitude));

phase = Toxy(phase, 0, 0.5*pi);

figure,imshow(phase,[]);title('phase');
object = objectAmplitude.*exp(1i.*phase);
figure,imshow(abs(object),[]);title('Input complex object');
%% create the wave vectors for the LED illumination
arraysize = 15; % size of LED array
xlocation = zeros(1,arraysize^2);
ylocation = zeros(1,arraysize^2);
LEDgap = 4; % 4mm between adjacent LEDs
LEDheight = 110; % 90 mm between the LED matrix and the sample
for i=1:arraysize % from top left to bottom right
    xlocation(1,1+arraysize*(i-1):15+arraysize*(i-1)) = (-(arraysize-1)/2:1:(arraysize-1)/2)*LEDgap;
    ylocation(1,1+arraysize*(i-1):15+arraysize*(i-1)) = ((arraysize-1)/2-(i-1))*LEDgap;
end
kx_relative = -sin(atan(xlocation/LEDheight)); % creat kx,ky wavevectors
ky_relative = -sin(atan(ylocation/LEDheight));
%% setup the parameters for the coherent imaging system
waveLength = 0.63e-6;
k0 = 2*pi/waveLength;
spsize = 6.5e-6/4; % sampling pixel size of the CCD             complex field
psize = spsize/4; % final pixel size of the reconstruction     即能恢复的物体复场的抽样间隔
NA = 0.1;
%% generate the low-pass filtered images
[m,n] = size(object); % image size of the high resolution object
m1 = floor(m/(spsize/psize));
n1 = floor(n/(spsize/psize)); % image size of the final output
imSeqLowRes = zeros(m1, n1, arraysize^2); %output low-res image sequence
kx = k0 * kx_relative;
ky = k0 * ky_relative;
cutoffFrequency = NA * k0;
kmax = pi/spsize;
[kxm, kym] = meshgrid(-kmax:kmax/((n1-1)/2):kmax,-kmax:kmax/((n1-1)/2):kmax);
CTF = double((kxm.^2+kym.^2)<cutoffFrequency^2); % coherent transfer function
dkx = 2*pi/(psize*n); % frequency sampling of reconstruction
dky = 2*pi/(psize*m);
objectFT = fftshift(fft2(object));
for tt=1:arraysize^2
    kxc = round((n+1)/2+kx(1,tt)/dkx);
    kyc = round((m+1)/2+ky(1,tt)/dky);
    kyl = round(kyc-(m1-1)/2);kyh = round(kyc+(m1-1)/2);
    kxl = round(kxc-(n1-1)/2);kxh = round(kxc+(n1-1)/2);
    imSeqLowFT = (m1/m)^2 * objectFT(kyl:kyh,kxl:kxh).*CTF;
    imSeqLowRes(:,:,tt) = Toxy((abs(ifft2(ifftshift(imSeqLowFT)))).^2,0,255);
end
figure;imshow(uint8(imSeqLowRes(:,:,113)),[]);

S = zeros(1,7); %记录SSIM值
for nn = 16:16:112;
I1 = uint8(imSeqLowRes(:,:,113-nn));
% figure,imshow(I1,[]);
I2 = uint8(imSeqLowRes(:,:,113+nn));
% figure,imshow(I2,[]);
S(nn/16) = ssim(I1,I2);
end

plot(S,'LineWidth',4);


% imwrite(I1,'amp_up_2.png','png');
% imwrite(I2,'amp_down_2.png','png');

% imwrite(I1,'pha_up_2.png','png');
% imwrite(I2,'pha_down_2.png','png');

% imwrite(I1,'com_up_2.png','png');
% imwrite(I2,'com_down_2.png','png');

