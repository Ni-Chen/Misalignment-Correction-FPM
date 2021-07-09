function errorpie1 = pie_solver_coarse(zz, imSeqLowRes, m, n, arraysize, m1, n1, dkx, dky, CTF, H, loop)

    global k0; global LEDgap; global NA_led; global arraysize; global image_center

    seq = gseq(9, 9, arraysize); % define the order of recovery
    objectRecoverFT = padarray(fftshift(fft2((imSeqLowRes(:, :, 145)))), [(m - m1) / 2, (m - m1) / 2]);

    xshift = zz(1); yshift = zz(2); theta = zz(3); LEDheight = zz(4);

    [kx, ky] = fpmangle(arraysize, xshift, yshift, theta, LEDheight, image_center);

    dkshift = [6; 6];
    pupil = 1;
    errorpie = zeros(1, loop);

    cen_wrong = zeros(2, arraysize^2);

    for tt = 1:1:loop
        tt
        errorpie(1, tt) = 0;

        for i3 = 1:1:arraysize^2
            i2 = seq(i3);
            % i2=i3;

            if tt < 20
                kxc = ((n + 1) / 2 + kx(1, i2) / dkx); % delete the round
                kyc = ((n + 1) / 2 + ky(1, i2) / dky); % delete the round
                cen_wrong(:, i2) = [kxc; kyc];
                cen_corr = cen_wrong;

                % with no sa correction
                kyl = round(kyc - (m1 - 1) / 2); kyh = round(kyc + (m1 - 1) / 2);
                kxl = round(kxc - (n1 - 1) / 2); kxh = round(kxc + (n1 - 1) / 2);
                lowResFT_1 = (m1 / m)^2 * objectRecoverFT(kyl:kyh, kxl:kxh) .* CTF .* pupil .* H;

            else

                % cen_corr(1,i)
                kxc = cen_corr(1, i2);
                kyc = cen_corr(2, i2);

                % kxc = round((n+1)/2+kx(1,i2)/dkx);
                % kyc = round((n+1)/2+ky(1,i2)/dky);

                % with sa procedure
                I_mea = imSeqLowRes(:, :, i2);

                cen = [kxc; kyc];
                cen11 = cen;

                opts.poscalibrate = 0;
                opts.calbratetol = 1e-3;
                poscost = @(ss) sum(sum((abs((m1 / m)^2 .* ifft2(ifftshift((m1 / m)^2 .* CTF .* H .* pupil .* objectftsamp(objectRecoverFT, ss, m1)))).^2 - I_mea).^2));
                optsanneal = saoptimset('Display', 'off', 'TolFun', opts.calbratetol);
                [cen1, fval11, exitflag11, output] = simulannealbnd(poscost, cen, cen11 - dkshift, cen11 + dkshift, optsanneal);
                cen_corr(:, i2) = cen1;

                kyl = round(cen1(2) - (m1 - 1) / 2); kyh = round(cen1(2) + (m1 - 1) / 2);
                kxl = round(cen1(1) - (n1 - 1) / 2); kxh = round(cen1(1) + (n1 - 1) / 2);
                lowResFT_1 = (m1 / m)^2 .* objectRecoverFT(kyl:kyh, kxl:kxh) .* CTF .* pupil .* H;

            end

            % disp([kxc;kyc]);

            % lowResFT_1 =objectRecoverFT(kyl:kyh,kxl:kxh).*CTF.*pupil;
            im_lowRes = ifft2(ifftshift(lowResFT_1));

            errorpie(1, tt) = errorpie(1, tt) + sum(sum(((m1 / m).^2 * abs(im_lowRes) - imSeqLowRes(:, :, i2)).^2));
            im_lowRes = (m / m1)^2 * sqrt(imSeqLowRes(:, :, i2)) .* exp(1i .* angle(im_lowRes)); % removed (m/m1)^2 *
            %  lowResFT_2 = fftshift(fft2(im_lowRes)).*CTF./pupil./H;
            lowResFT_2 = fftshift(fft2(im_lowRes));

            % sequential GS update
            %objectRecoverFT(kyl:kyh,kxl:kxh) = objectRecoverFT(kyl:kyh,kxl:kxh) + ...
            %    1.*conj(CTF.*pupil)./max(max(abs(pupil).^2)).*(lowResFT_2 - lowResFT_1);
            %pupil = pupil + 1.*conj(objectRecoverFT(kyl:kyh,kxl:kxh))./...
            %    max(max(abs(objectRecoverFT(kyl:kyh,kxl:kxh)).^2)).*(lowResFT_2 - lowResFT_1).*CTF;

            alpha = 1; beta = 1e0;
            % second order update
            objectRecoverFT(kyl:kyh, kxl:kxh) = objectRecoverFT(kyl:kyh, kxl:kxh) + ...
                1 / max(max(abs(pupil))) * abs(pupil) .* conj(pupil) .* (lowResFT_2 - lowResFT_1) ./ (abs(pupil).^2 + alpha);
            pupil = pupil + 1 / max(max(abs(objectRecoverFT(kyl:kyh, kxl:kxh)))) * abs(objectRecoverFT(kyl:kyh, kxl:kxh)) .* ...
                conj(objectRecoverFT(kyl:kyh, kxl:kxh)) .* (lowResFT_2 - lowResFT_1) ./ (abs(objectRecoverFT(kyl:kyh, kxl:kxh)).^2 + beta) .* CTF;

        end

        objectRecover = ifft2(ifftshift(objectRecoverFT));

    end

    objectRecover = ifft2(ifftshift(objectRecoverFT));

    errorpie1 = errorpie(1, tt) ./ 1e4;
    %  errorpie1=[ojbectRecover,pupil];

    figure(1)
    subplot(221); imagesc(abs(objectRecover)); colormap(gray); title('Amplitude');
    subplot(222); imagesc(angle(objectRecover)); colormap(gray); title('Phase');
    subplot(223); imagesc(abs(pupil)); colormap(gray); title('Pupil amplitude');
    subplot(224); imagesc(angle(pupil)); colormap(gray); title('Pupil phase');

    figure(2)

    scatter(cen_wrong(1, :), cen_wrong(2, :), [], 'k*'); hold on
    scatter(cen_corr(1, :), cen_corr(2, :), [], 'r');
    legend('Captured locations', 'Corrected locations');
    hold off

    disp(errorpie1);
    disp(zz);

end
