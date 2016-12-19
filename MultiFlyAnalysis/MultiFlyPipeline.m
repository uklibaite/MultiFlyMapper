%% Multi Fly Pipeline
% Ugne Klibaite
% 1/11/16

% Follow along in MultiFlyDoc text file for details and modifications
% Make sure current folder is set to /MultiFlyAnalysis


cd '/Users/uklibaite/Desktop/MultiFlyAnalysis'
filePath = '/Volumes/DataDisk/oregonR_JaneliaRecorded/oregonR_03_guid_9a3f41b0bd620188f071d16e37b6dc6.avi'

%% parameters
addpath('./parameters');
parameters = makeMultiFlyParameters([]);
        
%open new matlabpool
if matlabpool('size') > 0
    matlabpool close
end
matlabpool(parameters.numProcessors);

%% Preprocessing
addpath('./preprocessing');
addpath('./utilities');

%create savePath if needed
[dir,name,~] = fileparts(filePath);
parameters.savePath = [dir '/' name '/'];

[isNotSavePath,~] = unix(['ls ' parameters.savePath]);
if isNotSavePath ~= 0
    unix(['mkdir ' parameters.savePath]);
end

[s F] = unix(['ls ' parameters.savePath '*mat']);

[flyData,parameters] = runPreprocessing(filePath,parameters); 

savePath = parameters.savePath;
parameters.dataPath = [parameters.savePath 'prigData.mat'];
parameters.filePath = [parameters.savePath 'prig.avi'];
parameters.dataNum = 1;

% DB CHECK
save([savePath 'preCheck.mat'],'flyData','parameters');


%% Segmentation and Alignment
addpath('./segmentation');
addpath('./alignment');
startFrame = 1;
endFrame = 100000;
[flyData,runData,parameters] = multiFlySegmentAndAlign(startFrame,endFrame,parameters);

save([parameters.savePath 'flyInfo.mat'],'flyData','runData','parameters');

% ~~~~~~~~~~~~ End of Individual Movie Processing ~~~~~~~~~~~~~~~~
% Bookkeeping for file names
% flies from same movie have same filename somewhere
% flies from different populations have those specified
cd /Volumes/DataDisk/cantonS_JaneliaRecorded/cantons_court/
fBoth = '/Volumes/DataDisk/cantonS_JaneliaRecorded/cantons_court/';
[s t1] = unix(['find . -name "*_male.avi"']);
[s t2] = unix(['find . -name "*_female.avi"']);

cd /Volumes/DataDisk/cantonS_JaneliaRecorded/cantons_female_7day/
fFemale = '/Volumes/DataDisk/cantonS_JaneliaRecorded/cantons_female_7day/';
[s t3] = unix(['find . -name "*SINGLE.avi"']);

cd /Volumes/DataDisk/cantonS_JaneliaRecorded/cantons_male/
fMale = '/Volumes/DataDisk/cantonS_JaneliaRecorded/cantons_male/';
[s t4] = unix(['find . -name "*SINGLE.avi"']);

t1 = strsplit(t1,'\n'); t1 = t1(1:end-1);
t2 = strsplit(t2,'\n'); t2 = t2(1:end-1);
t3 = strsplit(t3,'\n'); t3 = t3(1:end-1);
t4 = strsplit(t4,'\n'); t4 = t4(1:end-1);

for i = 1:length(t1)
    T1{i} = [fBoth t1{i}(2:end)];
    T2{i} = [fBoth t2{i}(2:end)];
end
for i = 1:length(t3)
    T3{i} = [fFemale t3{i}(2:end)];
end
for i = 1:length(t4)
    T4{i} = [fMale t4{i}(2:end)];
end

files = [T1 T2 T3 T4];

% Also, from now on keep everything corresponding to the same experiment in
% a single folder
fName = '/Volumes/DataDisk/cantonS_JaneliaRecorded/';

%% Pixel Statistics

cd '/Users/uklibaite/Desktop/MultiFlyAnalysis'
addpath('./pixelStatistics');

num = 100000;
thetas = 0:2:178;
[meanRadon,stdRadon,vidObjs] = findImagePixelStatistics(files,num,thetas,1);

% determine # of pixels
[h x] = hist(stdRadon(:),50);
semilogy(x,h)
pixels = find(stdRadon>279);
    
save([fName '/multiFlySubsetStats.mat'],'meanRadon','stdRadon','vidObjs','pixels','thetas')




%% PCA

% do PCA
batchSize = 20000;
scale = 1;
for i = 1:length(files)
    vidObjs{i} = VideoReader(files{i});
end
[mu, vecs, vals] = onlineImagePCA_radon(vidObjs,batchSize,scale,pixels,thetas); 

save([fName '/VecsVals.mat'],'pixels','thetas','vecs','mu','vals');


% save projections
unix(['mkdir ' fName '/ProjectionDir/']);

numProjections = 50;
meanValues = mu;
for i = 1:length(files)
    projections = find_PCA_projections(vidObjs(i),vecs(:,1:numProjections),...
        meanValues,pixels,thetas,numProjections,scale,batchSize);
    save([fName '/ProjectionDir/projections_' num2str(i) '.mat'],'projections');
end
for i = 1:length(files)
projName{i} = [fName '/ProjectionDir/projections_' num2str(i) '.mat'];
end


%% Wavelet Decomposition


%% Sub-Sampled Embedding
projectionDirectory = [fName 'ProjectionDir/'];
 
addpath(genpath('./utilities/'));
addpath(genpath('./tSNE/'));

[trainingSetData,trainingSetAmps,projectionFiles] = runEmbeddingSubSampling(projectionDirectory,[]);

save([projectionDirectory 'training_set_embeddings.mat'],'trainingSetAmps','trainingSetData','projectionFiles');


%% t-SNE

[trainingEmbedding,betas,P,errors] = run_tSne(trainingSetData,parameters);
save([projectionDirectory 'training_set_embeddings2.mat'],'trainingSetAmps','trainingSetData','projectionFiles','parameters','trainingEmbedding','betas','P','errors');


L = length(files);
fprintf(1,'Finding t-SNE Embedding for each file\n');
embeddingValues = cell(L,1);
outputStats = cell(L,1);
for i=1:L
    
    fprintf(1,'\t Finding Embbeddings for File #%4i out of %4i\n',i,L);
    
    load(projName{i},'projections');
    projections = projections(:,1:50);
    
    [embeddingValues{i},outputStats{i}] = ...
        findEmbeddings(projections,trainingSetData,trainingEmbedding,[]);
    
    clear projections
    
end

% replace not in convex hull w/ z-guesses
embeddingValuesGuesses = cell(L,1);
for i = 1:L
    eV = embeddingValues{i};
    zGuesses = outputStats{i}.zGuesses;
    inCV = outputStats{i}.inConvHull;
    eV(~inCV,:) = zGuesses(~inCV,:);
    embeddingValuesGuesses{i} = eV;
end


save([projectionDirectory 'Embeddings.mat'],'embeddingValues','embeddingValuesGuesses','outputStats','projName')
%
maxVal = max(max(abs(combineCells(embeddingValuesGuesses))));
maxVal = round(maxVal * 1.1);
    
sigma = 1.2;
numPoints = 501;
rangeVals = [-maxVal maxVal];

[xx,density] = findPointDensity(combineCells(embeddingValuesGuesses),sigma,numPoints,rangeVals);

densities = zeros(numPoints,numPoints,L);
for i=1:L
    [~,densities(:,:,i)] = findPointDensity(embeddingValuesGuesses{i},sigma,numPoints,rangeVals);
end


%% Visualization
% different populations
D1 = densities(:,:,1:15);
D2 = densities(:,:,16:30);
D3 = densities(:,:,31:40);
D4 = densities(:,:,41:51);

D1All = mean(D1,3);
D2All = mean(D2,3);
D3All = mean(D3,3);
D4All = mean(D4,3);

subplot(2,2,1); imagesc(D1All); axis equal; axis off; title('Male Courting');
subplot(2,2,2); imagesc(D2All); axis equal; axis off; title('Female Courting');
subplot(2,2,3); imagesc(D4All); axis equal; axis off; title('Male Walking');
subplot(2,2,4); imagesc(D3All); axis equal; axis off; title('Female Walking');


% Watershed
LL = watershed(-density,8);
LL2 = LL; LL2(mean(densities,3) < 1e-6) = -1;
LL3 = zeros(size(LL));
for i = 1:501
    for j = 1:501
        if LL2(i,j)==0
            LL3(i,j)=1;
        end
    end
end


[IDbounds1 IDbounds2] = find(LL3==1);
scatter(IDbounds2,IDbounds1); axis equal off xy


%% Brady Movies
addpath('./Visualization');

vSmooth = .5;
medianLength = 1;
pThreshold = [];
minRest = [];
obj = [];
fitOnly = true;
numGMM = 2;

allZ = combineCells(embeddingValuesGuesses);
[watershedRegions,segments,v,obj,pRest,vals,vx,vy] = ...
    findWatershedRegions_v2(allZ,xx,LL,vSmooth,medianLength,pThreshold,minRest,obj,fitOnly,numGMM);


GMMobj = obj;
fitOnly = false;
watershedRegions = cell(1,51);
segments = cell(1,51);
v = cell(1,51);

for i = 1:length(embeddingValuesGuesses)
    
    zvals = embeddingValuesGuesses{i};
    [watershedRegions{i},segments{i},v{i}] = ...
        findWatershedRegions_v2(zvals,xx,LL,vSmooth,medianLength,pThreshold,minRest,GMMobj,fitOnly,numGMM);
end

maxNum = max(max(LL));
startFrames = ones(51,1);
minLength = 15;

[groups,~,~] = makeGroupsAndSegments(watershedRegions,maxNum,startFrames,minLength);







movies = files;

idx = makeManyBradyMovies_avi([projectionDirectory '/Bradies2/'],movies,groups,[],[],1,'true');






