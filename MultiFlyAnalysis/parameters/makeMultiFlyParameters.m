function [parameters] = makeMultiFlyParameters(parameters)

if nargin < 1
        parameters = [];
    end

%**********Segmentation Parameters*************
    
    
    imageSize = 75;
    
    
    minSize = 300;
    
    %DO NOT CHANGE THIS!!!!!!!!!!!!!!!!!!!!!!!
    splitLength = 300;
   
    medianN = 100;
    
    %minimal pixel value in segmentation after median subtraction
    minImageValue = 20;
    
    %size of image erosion for fly isolation
    imageErodeSize = 1;
    
    %location of ffmpeg on the system architecture
    ffmpegPath = '/opt/local/bin/';
    
    %location of mencoder on the system architecture
    mencoderPath = '/opt/local/bin/';  
    
    %location of tcprobe on the system architecture
    tcPath = '/usr/bin/';  
    
    %number of processors to use in parallel code
    numProcessors = 12;
    
    openingFilter = ones(3,3);
    
    defaultBodyThreshold = 170;
    
%**********Alignment Parameters*************


    %angular spacing used for Radon alignment
    spacing = 2;
    
    %eroding parameter in flipping detector
    erodeSize = 4;
    
    %accuracy (in pixels) of translational alignment
    fractionalPixelAccuracy = .2;
    
    %how many pixels to exclude from the edge in alignment
    edgeThreshold = 5;
    
    %width of region for front-back determination
    yRange = 25;
    
    %dilating of the background images for asymmetry-finding mask
    asymDilate = 25;
    
    %true if only aligning male images
    maleOnly = false;      
    
    %minimum length of a isolated sequence
    minDataLength = 100;
    
    %baseline value for median image size
    rescaleBaseline = 1300;
    
    %length of median filter to smooth oscillations in alignment
    centroidMedianFilterLength = 15;
    
    %maximum allowed body threshold
    maximumBodyThreshold = 140;
    
    %minimum allowed body threshold
    minimumBodyThreshold = 110;
    
    
    if ~isfield(parameters,'centroidMedianFilterLength') || isempty(parameters.centroidMedianFilterLength)
        parameters.centroidMedianFilterLength = centroidMedianFilterLength;
    end
    
    if ~isfield(parameters,'tcPath') || isempty(parameters.tcPath)
        parameters.tcPath = tcPath;
    end
    
    if ~isfield(parameters,'minimumBodyThreshold') || isempty(parameters.minimumBodyThreshold)
        parameters.minimumBodyThreshold = minimumBodyThreshold;
    end
    
    if ~isfield(parameters,'maximumBodyThreshold') || isempty(parameters.maximumBodyThreshold)
        parameters.maximumBodyThreshold = maximumBodyThreshold;
    end
    
    if ~isfield(parameters,'imageSize') || isempty(parameters.imageSize)
        parameters.imageSize = imageSize;
    end
    
    if ~isfield(parameters,'rescaleBaseline') || isempty(parameters.rescaleBaseline)
        parameters.rescaleBaseline = rescaleBaseline;
    end
    
    if ~isfield(parameters,'erodeSize') || isempty(parameters.erodeSize)
        parameters.erodeSize = erodeSize;
    end
    
    if ~isfield(parameters,'minDataLength') || isempty(parameters.minDataLength)
        parameters.minDataLength = minDataLength;
    end
    
    if ~isfield(parameters,'minSize') || isempty(parameters.minSize)
        parameters.minSize = minSize;
    end
    
    if ~isfield(parameters,'imageErodeSize') || isempty(parameters.imageErodeSize)
        parameters.imageErodeSize = imageErodeSize;
    end
    
    if ~isfield(parameters,'maleOnly') || isempty(parameters.maleOnly)
        parameters.maleOnly = maleOnly;
    end
    
    if ~isfield(parameters,'splitLength') || isempty(parameters.splitLength)
        parameters.splitLength = splitLength;
    end
    
    if ~isfield(parameters,'minImageValue') || isempty(parameters.minImageValue)
        parameters.minImageValue = minImageValue;
    end
    
    if ~isfield(parameters,'numProcessors') || isempty(parameters.numProcessors)
        parameters.numProcessors = numProcessors;
    end
   
    
    if ~isfield(parameters,'openingFilter') || isempty(parameters.openingFilter)
        parameters.openingFilter = openingFilter;
    end
    
    if ~isfield(parameters,'defaultBodyThreshold') || isempty(parameters.defaultBodyThreshold)
        parameters.defaultBodyThreshold = defaultBodyThreshold;
    end

    
    if ~isfield(parameters,'backgroundImages') || isempty(parameters.backgroundImages) || length(parameters.backgroundImages) < 2
        parameters.backgroundImages = cell(2,1);
        parameters.backgroundImages{1} = imread('female_bg.tiff');
        parameters.backgroundImages{2} = imread('male_bg.tiff');
    end
        
    if ~isfield(parameters,'ffmpegPath') || isempty(parameters.ffmpegPath)
        parameters.ffmpegPath = ffmpegPath;
    end
        
    if ~isfield(parameters,'spacing') || isempty(parameters.spacing)
        parameters.spacing = spacing;
    end    
    
    if ~isfield(parameters,'fractionalPixelAccuracy') || isempty(parameters.fractionalPixelAccuracy)
        parameters.fractionalPixelAccuracy = fractionalPixelAccuracy;
    end
    
    if ~isfield(parameters,'edgeThreshold') || isempty(parameters.edgeThreshold)
        parameters.edgeThreshold = edgeThreshold;
    end
    
    if ~isfield(parameters,'asymDilate') || isempty(parameters.asymDilate)
        parameters.asymDilate = asymDilate;
    end    
    
    if ~isfield(parameters,'yRange') || isempty(parameters.yRange)
        parameters.yRange = yRange;
    end
       
    if ~isfield(parameters,'medianN') || isempty(parameters.medianN)
        parameters.medianN = medianN;
    end
    
    if ~isfield(parameters,'mencoderPath') || isempty(parameters.mencoderPath)
        parameters.mencoderPath = mencoderPath;
    end


end

