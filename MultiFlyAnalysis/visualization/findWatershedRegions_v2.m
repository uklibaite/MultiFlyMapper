function [watershedRegions,segments,v,obj,pRest,vals,vx,vy] = ...
            findWatershedRegions_v2(all_z,xx,LL,vSmooth,medianLength,pThreshold,minRest,obj,fitOnly,numGMM)


    if nargin < 4 || isempty(vSmooth)
        vSmooth = .5;
    end

    if nargin < 5 || isempty(medianLength)
        medianLength = 2;
    end
    
    if nargin < 6 || isempty(pThreshold)
        pThreshold = [.67 .33];
    end
    
    if nargin < 7 || isempty(minRest)
        minRest = 5;
    end
    
    if nargin < 9 || isempty(fitOnly)
        fitOnly = false;
    end
    
    if nargin < 10 || isempty(numGMM)
        numGMM = 2;
    end
    
    restLength = 5;
    dt = .01;
    numToTest = 50000;
    minMixing = .1;
    N = length(all_z(:,1));
    
    fprintf(1,'Processing Data\n');
    smooth_z = all_z;
    if medianLength > 0
        smooth_z(:,1) = medfilt1(all_z(:,1),medianLength);
        smooth_z(:,2) = medfilt1(all_z(:,2),medianLength);
    end
    
    
    %rangeVal = xx(end) - xx(1);
        
    %vx = gaussianfilterdata_derivative_zeropad(smooth_z(:,1),vSmooth,dt);
    %vy = gaussianfilterdata_derivative_zeropad(smooth_z(:,2),vSmooth,dt);
    vx = [0;diff(smooth_z(:,1))]./dt;
    vy = [0;diff(smooth_z(:,2))]./dt;
    v = sqrt(vx.^2+vy.^2);

    
    fprintf(1,'Fitting Mixture Model\n');
    if nargin < 8 || isempty(obj)
        figure
        obj = gmixPlot(sampleFromMatrix(log(v(v>0))./log(10),numToTest),numGMM,[],200,[],true,[],[],3);
        drawnow;
    end
    %qq = find(obj.PComponents > minMixing);

    [~,maxIdx] = max(obj.mu);
    %minIdx = qq(minIdx);
    posts = posterior(obj,log(v)./log(10));
    posts(v==0,maxIdx) = 0;
    pRest = 1 - posts(:,maxIdx);
    
    
    if ~fitOnly
        
        fprintf(1,'Finding Watershed Values\n');
        
        vals = round((smooth_z + max(xx))*length(xx)/(2*max(xx)));
        vals(vals<1) = 1;
        vals(vals>length(xx)) = length(xx);
        
        
        watershedValues = zeros(N,1);
        parfor i=1:N
            watershedValues(i) = diag(LL(vals(i,2),vals(i,1)));
        end
        diffValues = abs([0;diff(watershedValues)]) == 0;
        
        fprintf(1,'Finding Segments\n');
        L = max(LL(:));
        if length(pThreshold) == 1
            
            CC = largeBWConnComp(pRest > pThreshold & diffValues,minRest);
            
        else
            
            minVal = min(pThreshold);
            maxVal = max(pThreshold);
            
            CC = largeBWConnComp(pRest > minVal & diffValues,minRest);
            maxInRanges = zeros(CC.NumObjects,1);
            parfor i=1:CC.NumObjects
                maxInRanges(i) = max(pRest(CC.PixelIdxList{i}));
            end
            
            CC.NumObjects = sum(maxInRanges >= maxVal);
            CC.PixelIdxList = CC.PixelIdxList(maxInRanges >= maxVal);
            
        end
        
        segmentAssignments = zeros(size(CC.PixelIdxList));
        watershedRegions = zeros(N,1);
        for i=1:CC.NumObjects
            segmentAssignments(i) = mode(watershedValues(CC.PixelIdxList{i}));
            watershedRegions(CC.PixelIdxList{i}) = segmentAssignments(i);
        end
        
        for i=1:L
            CC = largeBWConnComp(watershedValues == i,restLength);
            for j=1:length(CC.PixelIdxList)
                watershedRegions(CC.PixelIdxList{j}) = i;
            end
        end
        
                 segments = cell(L,1);
        %         for i=1:L
        %             CC = bwconncomp(watershedRegions == i);
        %             segments{i} = CC.PixelIdxList;
        %         end
        
    else
        
        watershedRegions = [];
        vals = [];
        segments = [];
        
    end