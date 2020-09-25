% Author      : Ni Chen
% Date        : 2012/05/17
% Description : Extending an image pixel value to a [min_val, max_val] scope.
% Status      : gold
%--------------------------------------------------------------------------
function outImage = Toxy(oriImage, min_val, max_val)
    [Ny, Nx, Chan] = size(oriImage);
    outImage = zeros(Ny, Nx, Chan);
    maxI = max(max(double(oriImage(:))));
    minI = min(min(double(oriImage(:))));
        
    for c = 1:Chan
        tempImg = oriImage(:,:,c);           
%         maxI = max(max(double(tempImg)));
%         minI = min(min(double(tempImg)));
        if maxI == minI
            if minI~=0
                outImage(:,:,c) = double(tempImg)*max_val/minI;
            end
        else
            coeA = (max_val-min_val)/double(maxI-minI);
            coeB = min_val-minI*coeA;
            outImage(:,:,c) = double(tempImg).*coeA+coeB;
        end
    end 
end