function [rotationAngle,X,Y,finalImage,errors,finalOriginalImage] = ...
            alignTwoImages(image1,image2,angleGuess,alignInfo,originalImage,noRotation)


    if nargin < 3 || isempty(angleGuess) == 1
        angleGuess = 0;
    else
        angleGuess = mod(angleGuess,360);
    end
    
    angleGuess = angleGuess*pi/180;
    
    if nargin < 6 || isempty(noRotation)
        noRotation = false;
    end
    
    if ~isfield(alignInfo,'spacing') || isempty(alignInfo.spacing)
        spacing = .5;
    else
        spacing = alignInfo.spacing;
    end
    N = 180/spacing;
    
    
    if nargin < 5
        originalImage = [];
        finalOriginalImage = [];
    end
    
    errors = zeros(2,1);
    
    s = size(image1);
    
    
    if ~noRotation
        
        
        rp = regionprops(image2 > 0,'Centroid');
        midVal = round(s./2);
        vals = fliplr(round(rp(1).Centroid - midVal));
        if min(abs(vals)) > 0
            image2 = circshift(image2,-vals);
            if ~isempty(originalImage)
                originalImage = circshift(originalImage,-vals);
            end
        end
        
        thetas = linspace(0, 180-spacing, N);
        
        %Find fft of the Radon transform
        
        F1 = abs(fft(radon(image1, thetas)));
        F2 = abs(fft(radon(image2, thetas)));
        
        %Find the index of the correlation peak
        correlation = sum(conj(fft2(F1)) .* fft2(F2));
        peaks = real(ifft(correlation));
        peakIndex = find(peaks==max(peaks));
        
        if length(peakIndex) > 1
            peakIndex = peakIndex(1);
        end
        
        %Find rotation angle via quadratic interpolation
        if (peakIndex~=1) && (peakIndex ~= N)
            p=polyfit(thetas((peakIndex-1):(peakIndex+1)),peaks((peakIndex-1):(peakIndex+1)),2);
            rotationAngle = -.5*p(2)/p(1);
            errors(1) = polyval(p,rotationAngle);
        else
            if peakIndex == 1
                p = polyfit([thetas(end)-180,thetas(1),thetas(2)],peaks([N,1,2]),2);
                rotationAngle = -.5*p(2)/p(1);
                errors(1) = polyval(p,rotationAngle);
                if rotationAngle < 0
                    rotationAngle = 180 + rotationAngle;
                end
            else
                p = polyfit([thetas(end-1),thetas(end),180+thetas(1)],peaks([N-1,N,1]),2);
                rotationAngle = -.5*p(2)/p(1);
                errors(1) = polyval(p,rotationAngle);
                if rotationAngle >= 180
                    rotationAngle = rotationAngle - 180;
                end
            end
        end
        
        
        %Check to see if rotation angle is in the correct direction
        rA = rotationAngle*pi/180;
        test = dot([cos(rA),sin(rA)],[cos(angleGuess),sin(angleGuess)]);
        if test < 0
            rotationAngle = mod(rotationAngle-180,360);
        end
        rotationAngle = mod(rotationAngle,360);
        toRotate = mod(-rotationAngle,360);
        
        %Rotate Image & Crop to original Size
        rotatedImage = imrotate(originalImage,toRotate,'crop');
        
    else
        
        rotationAngle = mod(angleGuess,360);
        toRotate = mod(-rotationAngle,360);
        rotatedImage = imrotate(image2,toRotate,'crop');
        
    end
    
    [X,Y] = translationalAlignment(image1,rotatedImage);

    finalImage = [];
    
    
    if ~isempty(originalImage)

        T = maketform('affine',[1 0 0 ;0 1 0;X Y 1]);
        
        finalOriginalImage = imtransform(rotatedImage,T,'XData',[1 s(2)],'YData',[1 s(1)]);

    end
    
    
    
    
    
   