clc,clear
close all
addpath pic
%% simulate the forward imaging process of Fourier ptychography
% simulate the high resolution complex object
objectAmplitude = sqrt(double(imresize(rgb2gray(imread('4.2.03.tiff')),[511,511])));
% objectAmplitude = ones(511,511);

% phase = double(imresize(rgb2gray(imread('5.2.09.tiff')),[511,511]));
phase = zeros(size(objectAmplitude));

% phase = Toxy(phase, 0, 2*pi);
% figure,imshow(phase,[]);title('phase');
object = objectAmplitude.*exp(1i.*phase);
figure,imshow(abs(object),[]);title('Input complex object');
%% create the wave vectors for the LED illumination
arraysize = 15; % size of LED array
xlocation = zeros(1,arraysize^2);
ylocation = zeros(1,arraysize^2);
LEDgap = 4; % 4mm between adjacent LEDs
LEDheight = 110; % 110 mm between the LED matrix and the sample
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
m1 =floor( m/(spsize/psize) );
n1 =floor( n/(spsize/psize) ); % image size of the final output
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
    imSeqLowRes(:,:,tt) = abs(ifft2(ifftshift(imSeqLowFT)));
end
figure;imshow(uint8(imSeqLowRes(:,:,113)),[]);

%% SSIM self-learning based FPM
% s = zeros(arraysize);
% for nn=1:arraysize^2
%     s(nn) = SSIM1(imSeqLowRes(:,:,113),imSeqLowRes(:,:,nn));
% end
% ss = s>0.01; %该阈值可调整
% imagesc(ss);
% ss = ss';
%% recover the high resolution image
seq = gseq(arraysize); % define the order of recovery,we start from the center (the 113rd image) to the edge of the spectrum (the 225th image)
objectRecover = ones(m,n); % initial guess of the object
% objectRecover = imSeqLowRes(:,:,113);
objectRecoverFT = fftshift(fft2(objectRecover));
% objectRecoverFT = padarray(fftshift(fft2((imSeqLowRes(:,:,113)))),[(m-m1)/2, (m-m1)/2]);
loop = 15;
for tt=1:loop
    for i3=1:arraysize^2
        i2 = seq(i3);
%         if ss(i2) %是否使用self-learning
        if i2<=arraysize*(arraysize+1)/2   
        kxc = round((n+1)/2+kx(1,i2)/dkx);
        kyc = round((n+1)/2+ky(1,i2)/dky);
        kyl = round(kyc-(m1-1)/2);kyh = round(kyc+(m1-1)/2);
        kxl = round(kxc-(n1-1)/2);kxh = round(kxc+(n1-1)/2);
        lowResFT = (m1/m)^2 * objectRecoverFT(kyl:kyh,kxl:kxh).*CTF;
        im_lowRes = ifft2(ifftshift(lowResFT));
        im_lowRes = (m/m1)^2 * imSeqLowRes(:,:,i2).*exp(1i.*angle(im_lowRes));
        lowResFT = fftshift(fft2(im_lowRes)).*CTF;
        objectRecoverFT(kyl:kyh,kxl:kxh) = (1-CTF).*objectRecoverFT(kyl:kyh,kxl:kxh) + lowResFT;
        end
%         end
    end
end

%%频谱不补全
% objectRecover = ifft2(ifftshift(objectRecoverFT));
% objectRecoverFT3 = objectRecoverFT;

%%频谱补全 method1
% objectRecoverFT1 = rot90(objectRecoverFT,2);
% objectRecoverFT2 = objectRecoverFT;
% objectRecoverFT2(257:end,:) = 0;
% objectRecoverFT1(1:256,:) = 0;
% objectRecoverFT3 = conj(objectRecoverFT1) + objectRecoverFT2;
% objectRecover = ifft2(ifftshift(objectRecoverFT3));

%%频谱补全 method2
objectRecoverFT1 = rot90(objectRecoverFT,2);
objectRecoverFT2 = objectRecoverFT;
% objectRecoverFT2(257:end,:) = 0;
objectRecoverFT1(1:288,:) = 0;
objectRecoverFT3 = conj(objectRecoverFT1) + objectRecoverFT2;
objectRecover = ifft2(ifftshift(objectRecoverFT3));


imag_recover = uint8(Toxy((abs(objectRecover)).^2, 0, 255));
figure;imshow(imag_recover,[]);title('intensity recover');
% figure;imshow(Toxy(angle(objectRecover),0,2*pi),[]);title('phase recover');
figure;imshow(log(abs(objectRecoverFT3)+1),[]);



% imwrite(imag_recover, 'amp_full.png','png');
% imwrite(imag_recover, 'amp_half.png','png');
% imwrite(uint8(Toxy(log(abs(objectRecoverFT)+1), 0, 255)), 'amp_full_s.png','png');
% imwrite(uint8(Toxy(abs(objectRecover), 0, 255)), 'amp_half2full.png','png');

% imwrite(imag_recover, 'amp_half2full1.png','png');
% imwrite(uint8(Toxy(log(abs(objectRecoverFT3)+1), 0, 255)), 'amp_half2full1_s.png','png');

% imwrite(imag_recover, 'amp_half2full2.png','png');
% imwrite(uint8(Toxy(log(abs(objectRecoverFT3)+1), 0, 255)), 'amp_half2full2_s.png','png');


% imwrite(uint8(Toxy(angle(objectRecover), 0, 255)), 'pha_full.png','png');
% imwrite(uint8(Toxy(angle(objectRecover), 0, 255)), 'pha_half.png','png');


% imwrite(uint8(Toxy(angle(objectRecover), 0, 255)), 'com_pha_full.png','png');
% imwrite(uint8(Toxy(angle(objectRecover), 0, 255)), 'com_pha_half.png','png');

% imwrite(uint8(Toxy(abs(objectRecover), 0, 255)), 'com_amp_full.png','png');
% imwrite(uint8(Toxy(abs(objectRecover), 0, 255)), 'com_amp_half.png','png');

