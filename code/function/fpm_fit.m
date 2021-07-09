function errorfpmfit = fpm_fit(zz, n, dkx, dky, seq, cen_corr, mode, arraysize1)

    global k0; global LEDgap; global NA_led; global arraysize; global image_center

    xshift = zz(1); yshift = zz(2); theta = zz(3); LEDheight = zz(4);

    % the shifts of LED in the fourier plane
    %arraysize = 17;
    xlocation1 = zeros(1, arraysize^2);
    ylocation1 = zeros(1, arraysize^2);
    rlocation = zeros(1, arraysize^2);
    %LEDgap = 4; % 4mm between adjacent LEDs
    %LEDheight = 113.5;% 108 mm between the LED matrix and the sample

    %LEDheight = 113.5;% 108 mm between the LED matrix and the sample
    for i = 1:arraysize % from top left to bottom right
        xlocation1(1, 1 + arraysize * (i - 1):arraysize + arraysize * (i - 1)) = (-(arraysize - 1) / 2:1:(arraysize - 1) / 2) * LEDgap;
        ylocation1(1, 1 + arraysize * (i - 1):arraysize + arraysize * (i - 1)) = ((arraysize - 1) / 2 - (i - 1)) * LEDgap;

    end

    xlocation = cos(theta) .* xlocation1 + sin(theta) .* ylocation1 + xshift;
    ylocation = -sin(theta) .* xlocation1 + cos(theta) .* ylocation1 + yshift;

    kx_relative = sin(atan((image_center(2) - xlocation) / LEDheight)); % creat kx,ky wavevectors
    ky_relative = sin(atan((image_center(1) - ylocation) / LEDheight));
    NA_led = sqrt(kx_relative.^2 + ky_relative.^2);
    kx = k0 * kx_relative;
    ky = k0 * ky_relative;

    if mode == 1

        errorfpmfit = 0;

        for i3 = 1:1:arraysize1^2
            i2 = seq(i3);
            kxc = ((n + 1) / 2 + kx(1, i2) / dkx); % delete the round
            kyc = ((n + 1) / 2 + ky(1, i2) / dky);
            errorfpmfit = errorfpmfit + (kxc - cen_corr(1, i2)).^2 + (kyc - cen_corr(2, i2)).^2;
        end

    else
        errorfpmfit = zeros(2, arraysize^2);

        for i3 = 1:1:arraysize^2
            i2 = seq(i3);
            kxc = ((n + 1) / 2 + kx(1, i2) / dkx); % delete the round
            kyc = ((n + 1) / 2 + ky(1, i2) / dky);
            % errorfpmfit=errorfpmfit+(kxc-cen_corr(1,i2)).^2+(kyc-cen_corr(2,i2)).^2;
            errorfpmfit(1, i2) = kxc;
            errorfpmfit(2, i2) = kyc;
        end

    end

end
