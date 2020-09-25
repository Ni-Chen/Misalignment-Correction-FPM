function [kx,ky]=fpmangle(arraysize,xshift,yshift,theta,LEDheight,image_center)

global k0;global LEDgap;global NA_led;

% the shifts of LED in the fourier plane
%arraysize = 17;
xlocation1 = zeros(1,arraysize^2);
ylocation1 = zeros(1,arraysize^2);
% rlocation = zeros(1,arraysize^2);
%LEDgap = 4; % 4mm between adjacent LEDs
%LEDheight = 113.5;% 108 mm between the LED matrix and the sample

%LEDheight = 113.5;% 108 mm between the LED matrix and the sample
for i=1:arraysize % from top left to bottom right
    xlocation1(1,1+arraysize*(i-1):arraysize+arraysize*(i-1)) = (-(arraysize-1)/2:1:(arraysize-1)/2)*LEDgap;
    ylocation1(1,1+arraysize*(i-1):arraysize+arraysize*(i-1)) = ((arraysize-1)/2-(i-1))*LEDgap;
    
    
end

xlocation=cos(theta).*xlocation1+sin(theta).*ylocation1+xshift;
ylocation=-sin(theta).*xlocation1+cos(theta).*ylocation1+yshift;

kx_relative = sin(atan((image_center(2)-xlocation)/LEDheight)); % creat kx,ky wavevectors
ky_relative = sin(atan((image_center(1)-ylocation)/LEDheight));
NA_led = sqrt(kx_relative.^2 + ky_relative.^2);
kx = k0 * kx_relative;
ky = k0 * ky_relative;

end

