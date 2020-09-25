clc,clear
close all

addpath(genpath('.\function\'));  % Add funtion path with sub-folders
addpath(genpath('.\data\'));  % Add funtion path with sub-folders

% load('20170531-3.mat');        %root  ��Ҫ����system�еĲ���
load('20170604-3.mat');      %High USAF ��Ҫ����system�еĲ���

%%self-learning
S = zeros(1,4);
s = zeros(17);
for nn=1:17^2
    S = SSIM1(imSeqLowRes(:,:,145),imSeqLowRes(:,:,nn)); %MATLAB 2015�汾���Ͽ�ʹ���Դ�ssim��������
    s(nn) = S(1);
end
figure
imagesc(s);
ss = s>0.012;
ss = ss';
NN = sum(double(ss(:)));
figure
imagesc(ss);

system();
%% recover the high resolution image
seq = gseq(arraysize); % define the order of recovery
objectRecoverFT = padarray(fftshift(fft2((imSeqLowRes(:,:,145)))),[(m-m1)/2, (m-m1)/2]);



loop = 15;pupil = CTF;alpha = 1;beta = 1e3;
for tt=1:loop
    for i3=1:arraysize^2
        i2 = seq(i3);
        %                 if ss(i2) %�Ƿ���self-learning
        
        if i2<=arraysize*(arraysize+1)/2
            kxc = round((n+1)/2+kx(1,i2)/dkx);
            kyc = round((n+1)/2+ky(1,i2)/dky);
            kyl = round(kyc-(m1-1)/2);kyh = round(kyc+(m1-1)/2);
            kxl = round(kxc-(n1-1)/2);kxh = round(kxc+(n1-1)/2);
            lowResFT_1 = (m1/m)^2 *objectRecoverFT(kyl:kyh,kxl:kxh).*pupil;
            im_lowRes = ifft2(ifftshift(lowResFT_1));
            im_lowRes = (m/m1)^2 *sqrt(imSeqLowRes(:,:,i2)).*im_lowRes./(abs(im_lowRes)+eps);
            lowResFT_2 = fftshift(fft2(im_lowRes));
            objectRecoverFT(kyl:kyh,kxl:kxh) = objectRecoverFT(kyl:kyh,kxl:kxh) + ...
                1/max(max(abs(pupil)))*abs(pupil).*conj(pupil).*(lowResFT_2 - lowResFT_1)./(abs(pupil).^2+alpha);
            pupil = pupil + 1/max(max(abs(objectRecoverFT(kyl:kyh,kxl:kxh))))*abs(objectRecoverFT(kyl:kyh,kxl:kxh)).*...
                conj(objectRecoverFT(kyl:kyh,kxl:kxh)).*(lowResFT_2 - lowResFT_1)./(abs(objectRecoverFT(kyl:kyh,kxl:kxh)).^2+beta).*CTF;
        end
        %                 end
    end
end


%% root Ƶ�ײ�ȫ
% objectRecoverFT(453:end,:)=0;
% objectRecoverFTT = objectRecoverFT(2:800,2:800);
% objectRecoverFT1 = rot90(objectRecoverFTT,2);
% objectRecoverFT2 = objectRecoverFTT;
% objectRecoverFT1(100:451,:) = 0; %root
% objectRecoverFT3 = conj(objectRecoverFT1) + objectRecoverFT2;
% objectRecoverFT(2:800,2:800) = objectRecoverFT3;
% objectRecover = ifft2(ifftshift(objectRecoverFT));


%% high USAF 
objectRecoverFT(291:end,:)=0; %high USAF
objectRecoverFTT = objectRecoverFT(2:512,2:512);
objectRecoverFT1 = rot90(objectRecoverFTT,2);
objectRecoverFT2 = objectRecoverFTT;
% objectRecoverFT1(100:256,:) = 0; %Ƶ�ײ�ȫ method1
objectRecoverFT1(100:289,:) = 0;   %Ƶ�ײ�ȫ method2
objectRecoverFT3 = conj(objectRecoverFT1) + objectRecoverFT2;
objectRecoverFT(2:512,2:512) = objectRecoverFT3;
objectRecover = ifft2(ifftshift(objectRecoverFT));


%%Ƶ�ײ���ȫ
% objectRecover = ifft2(ifftshift(objectRecoverFT));

objectRecover_imag = (abs(objectRecover)).^2 ;
figure,imagesc(imSeqLowRes(:,:,145)); title('LR');axis image; colormap gray; axis off
figure,imagesc(objectRecover_imag); title('HR');axis image; colormap gray; axis off
figure,imagesc(angle(objectRecover)); title('phase');axis image; colormap gray; axis off
figure,imagesc(log(abs(objectRecoverFT)));title('spectrum'); axis image; colormap gray; axis off
