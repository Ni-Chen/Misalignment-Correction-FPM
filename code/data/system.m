%% setup the parameters for the coherent imaging system
% m1 = 200;n1 = 200;            %the size of patch  ROOT
m1 = 128;n1 = 128;          %the size of patch USAF
waveLength = 0.6297e-6;
NA = 0.1;                     %objective lens NA
mag = 4.0;                    %objective lens magnification
LEDgap = 4;                 % 4mm between adjacent LEDs
% LEDheight = 113;            % 113 mm between the LED matrix and the sample   root
LEDheight = 113.5;          %high USAF
pix = 6.5e-6;                 %pixel size of the CCD
spsize = pix/mag;             % sampling pixel size of the CCD
mag_image = 4;                %the reconstruction image magnification
psize = spsize/mag_image;     % final pixel size of the reconstruction
m =m1*mag_image;n =n1*mag_image; % image size of the high resolution object

%% create the wave vectors for the LED illumination
arraysize = 17;%root / high USAF
xlocation = zeros(1,arraysize^2);
ylocation = zeros(1,arraysize^2);
for i=1:arraysize           % from top left to bottom right
    xlocation(1,1+arraysize*(i-1):arraysize+arraysize*(i-1)) = (-(arraysize-1)/2:1:(arraysize-1)/2)*LEDgap;
    ylocation(1,1+arraysize*(i-1):arraysize+arraysize*(i-1)) = ((arraysize-1)/2-(i-1))*LEDgap;
end
image_center = [0 0]*spsize*1e3;  % pixel numer of the patch center with respect to the raw image center
kx_relative = sin(atan((image_center(2)-xlocation)/LEDheight)); % creat kx,ky wavevectors
ky_relative = sin(atan((image_center(1)-ylocation)/LEDheight));
NA_led = sqrt(kx_relative.^2 + ky_relative.^2);
k0 = 2*pi/waveLength;
kx = k0 * kx_relative;
ky = k0 * ky_relative;
dkx = 2*pi/(psize*n);
dky = 2*pi/(psize*m);
cutoffFrequency = NA * k0;
kmax = pi/spsize;
[kxm, kym] = meshgrid(-kmax:kmax/((n1-1)/2):kmax,-kmax:kmax/((n1-1)/2):kmax);
CTF = double((kxm.^2+kym.^2)<cutoffFrequency^2);   % coherent transfer function
NA_syn = NA + max(NA_led);                         % synthesis NA

%% defocus distance
z = 0; % defocus distance
kzm = sqrt(k0^2-kxm.^2-kym.^2);
H = exp(1i.*z.*real(kzm)).*exp(-z.*imag(kzm));