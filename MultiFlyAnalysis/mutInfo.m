%%
load('/Users/uklibaite/Desktop/Emory_MacPro/CosyneFigs/EmbeddingsByCell.mat');
load('/Users/uklibaite/Desktop/Emory_MacPro/CosyneFigs/CantonSWRInfo.mat');
load('/Users/uklibaite/Desktop/Emory_MacPro/CosyneFigs/WatershedsByCell.mat');



allF = combineCells(eVfemale);
allM = combineCells(eVmale);


EE = [allF allM];

sigma = 1.2;
numPoints = 101;
rangeVals = [-56 56];
xx = linspace(rangeVals(1),rangeVals(2),numPoints);

[~, dM] = findPointDensity(allM,sigma,numPoints,rangeVals);
[~, dF] = findPointDensity(allF,sigma,numPoints,rangeVals);


yy = xx;
vv = xx;
ww = xx;
[XX,YY,VV,WW] = ndgrid(xx,yy,vv,ww);

dx = xx(2)-xx(1);
sigma = 1.2;
G = exp(-.5.*(XX.^2 + YY.^2 + VV.^2 + WW.^2)./sigma^2) ./ (4*(pi^2)*sigma^4);

Gp = G*dx^4;
Gp = Gp/sum(Gp(:));

Z = zeros(size(G));


for i = 1:length(EE)
    if mod(i,1000)==0
        i
    end
    currE = EE(i,:);
    [mV idx1] = min(abs(xx-currE(1)));
    [mV idx2] = min(abs(xx-currE(2)));
    [mV idx3] = min(abs(xx-currE(3)));
    [mV idx4] = min(abs(xx-currE(4)));
    Z(idx1, idx2, idx3, idx4) = Z(idx1, idx2, idx3, idx4)+1;
end

Z = Z./sum(Z(:));

% have histogram, need to convolve with gaussian to get density
clear EE mPlaceAngTheta m2fHeadTheta eVmale eVmaleSingle eVfemale eVfemaleSingle Rmm 

density = fftshift(real(ifftn(fftn(Gp).*fftn(Z))));
density(density<0) = 0;
density = real(density);

% check similarity to dM
dm = imrotate(flipud(squeeze(sum(sum(density,1),2))),-90);
df = imrotate(flipud(squeeze(sum(sum(density,3),4))),-90);
% 

subplot(1,2,1); imagesc(dm); axis equal off
subplot(1,2,2); imagesc(dM); axis equal off

% sum1 = sum(sum(dM*dF));
% sumD = sum(density(:));
% 
% I = sumD*log2(sumD/sum1);

dM = dM/sum(dM(:));
dF = dF/sum(dF(:));

MI = zeros(length(xx),length(yy),length(vv),length(ww));
for x = 1:length(xx)
    x
    for y = 1:length(yy)
        for v = 1:length(vv)
            for w = 1:length(ww)
                
                MI(x,y,v,w) = density(x,y,v,w)*log2(density(x,y,v,w)/(dF(x,y)*dM(v,w)));
                
            end
        end
    end
end
MI(isnan(MI)) = 0;
MI(isinf(MI)) = 0;


fSum = imrotate(flipud(squeeze(sum(sum(MI,3),4))),-90);
mSum = imrotate(flipud(squeeze(sum(sum(MI,1),2))),-90);


subplot(2,2,1); imagesc(dF); axis equal off
subplot(2,2,2); imagesc(abs(fSum)); axis equal off
subplot(2,2,3); imagesc(dM); axis equal off
subplot(2,2,4); imagesc(abs(mSum)); axis equal off



tZ = -3000:100:0;
MIValFshift30 = zeros(size(tZ));
MIValMshift30 = zeros(size(tZ));

for i = 1:length(tZ)
    
    tau = tZ(i);
    i

    for k = 1:length(eVfemale)
        orgF = eVfemale{k};
        orgM = eVmale{k};
       eFadj{k} = orgF(1-tau:end,:);
       eMadj{k} = orgM(1:end+tau,:);
    end
    
    eF = combineCells(eFadj');
    eM = combineCells(eMadj');
    
    
    EE = [eF eM];

    
    Z = zeros(size(G));
    
    fprintf('\n calculating Z \n');
    
    for j = 1:length(EE)
 
        currE = EE(j,:);
        [mV idx1] = min(abs(xx-currE(1)));
        [mV idx2] = min(abs(xx-currE(2)));
        [mV idx3] = min(abs(xx-currE(3)));
        [mV idx4] = min(abs(xx-currE(4)));
        Z(idx1, idx2, idx3, idx4) = Z(idx1, idx2, idx3, idx4)+1;
    end
    
    Z = Z./sum(Z(:));
    
    density = fftshift(real(ifftn(fftn(Gp).*fftn(Z))));
    density(density<0) = 0;
    density = real(density);
    
    fprintf('calculating MI \n');
    
    MI = zeros(size(density));
    
    for x = 1:length(xx)
        fprintf(['Calculating MI of ' num2str(x) '\n']);
        for y = 1:length(yy)
            for v = 1:length(vv)
                for w = 1:length(ww)
                    
                    MI(x,y,v,w) = density(x,y,v,w)*log2(density(x,y,v,w)/(dF(x,y)*dM(v,w)));
                    
                end
            end
        end
    end
    
    MI(isnan(MI)) = 0;
    MI(isinf(MI)) = 0;

    
    MIValFshift30(i) = sum(MI(:));

end



for i = 1:length(tZ)
    
    tau = tZ(i);
    i

    for k = 1:length(eVfemale)
        orgF = eVfemale{k};
        orgM = eVmale{k};
       eFadj{k} = orgF(1:end+tau,:);
       eMadj{k} = orgM(1-tau:end,:);
       
    end
    
    eF = combineCells(eFadj');
    eM = combineCells(eMadj');
    
    
    EE = [eF eM];

    
    Z = zeros(size(G));
    
    fprintf('\n calculating Z \n');
    
    for j = 1:length(EE)
 
        currE = EE(j,:);
        [mV idx1] = min(abs(xx-currE(1)));
        [mV idx2] = min(abs(xx-currE(2)));
        [mV idx3] = min(abs(xx-currE(3)));
        [mV idx4] = min(abs(xx-currE(4)));
        Z(idx1, idx2, idx3, idx4) = Z(idx1, idx2, idx3, idx4)+1;
    end
    
    Z = Z./sum(Z(:));
    
    density = fftshift(real(ifftn(fftn(Gp).*fftn(Z))));
    density(density<0) = 0;
    density = real(density);
    
    fprintf('calculating MI \n');
    
    MI = zeros(size(density));
    
    for x = 1:length(xx)
        fprintf(['Calculating MI of ' num2str(x) '\n']);
        for y = 1:length(yy)
            for v = 1:length(vv)
                for w = 1:length(ww)
                    
                    MI(x,y,v,w) = density(x,y,v,w)*log2(density(x,y,v,w)/(dF(x,y)*dM(v,w)));
                    
                end
            end
        end
    end
    
    MI(isnan(MI)) = 0;
    MI(isinf(MI)) = 0;

    
    MIValMshift30(i) = sum(MI(:));

end

plot(MIValFshift30,'b'); hold on
plot(MIValMshift30,'r'); 

save('mutualInfo30.mat','MIValFshift30','MIValMshift30')

%% shuffled data

tZ = -30000:1000:0;
MIValFshiftShuffled = zeros(size(tZ));
MIValMshiftShuffled = zeros(size(tZ));

for i = 1:length(tZ)
    
    tau = tZ(i);
    i

    for k = 1:length(eVfemale)
        orgF = eVfemale{k};
        orgM = eVmale{k};
       eFadj{k} = orgF(1-tau:end,:);
       eMadj{k} = orgM(1:end+tau,:);
    end
    
    eFtemp = combineCells(eFadj');
    eFidx = randperm(length(eFtemp));
    eF = zeros(length(eFidx),2);
    for t = 1:length(eFidx)
        eF(t,:) = eFtemp(eFidx(t),:);
    end
    
    eM = combineCells(eMadj');
    
    
    EE = [eF eM];

    
    Z = zeros(size(G));
    
    fprintf('\n calculating Z \n');
    
    for j = 1:length(EE)
 
        currE = EE(j,:);
        [mV idx1] = min(abs(xx-currE(1)));
        [mV idx2] = min(abs(xx-currE(2)));
        [mV idx3] = min(abs(xx-currE(3)));
        [mV idx4] = min(abs(xx-currE(4)));
        Z(idx1, idx2, idx3, idx4) = Z(idx1, idx2, idx3, idx4)+1;
    end
    
    Z = Z./sum(Z(:));
    
    density = fftshift(real(ifftn(fftn(Gp).*fftn(Z))));
    density(density<0) = 0;
    density = real(density);
    
    fprintf('calculating MI \n');
    
    MI = zeros(size(density));
    
    for x = 1:length(xx)
        fprintf(['Calculating MI of ' num2str(x) '\n']);
        for y = 1:length(yy)
            for v = 1:length(vv)
                for w = 1:length(ww)
                    
                    MI(x,y,v,w) = density(x,y,v,w)*log2(density(x,y,v,w)/(dF(x,y)*dM(v,w)));
                    
                end
            end
        end
    end
    
    MI(isnan(MI)) = 0;
    MI(isinf(MI)) = 0;

    
    MIValFshiftShuffled(i) = sum(MI(:));

end



for i = 1:length(tZ)
    
    tau = tZ(i);
    i

    for k = 1:length(eVfemale)
        orgF = eVfemale{k};
        orgM = eVmale{k};
       eFadj{k} = orgF(1:end+tau,:);
       eMadj{k} = orgM(1-tau:end,:);
       
    end
    
    eF = combineCells(eFadj');
    eMtemp = combineCells(eMadj');
    eMidx = randperm(length(eMtemp));
    eM = zeros(length(eMidx),2);
    for t = 1:length(eMidx)
        eM(t,:) = eMtemp(eMidx(t),:);
    end

    EE = [eF eM];

    
    Z = zeros(size(G));
    
    fprintf('\n calculating Z \n');
    
    for j = 1:length(EE)
 
        currE = EE(j,:);
        [mV idx1] = min(abs(xx-currE(1)));
        [mV idx2] = min(abs(xx-currE(2)));
        [mV idx3] = min(abs(xx-currE(3)));
        [mV idx4] = min(abs(xx-currE(4)));
        Z(idx1, idx2, idx3, idx4) = Z(idx1, idx2, idx3, idx4)+1;
    end
    
    Z = Z./sum(Z(:));
    
    density = fftshift(real(ifftn(fftn(Gp).*fftn(Z))));
    density(density<0) = 0;
    density = real(density);
    
    fprintf('calculating MI \n');
    
    MI = zeros(size(density));
    
    for x = 1:length(xx)
        fprintf(['Calculating MI of ' num2str(x) '\n']);
        for y = 1:length(yy)
            for v = 1:length(vv)
                for w = 1:length(ww)
                    
                    MI(x,y,v,w) = density(x,y,v,w)*log2(density(x,y,v,w)/(dF(x,y)*dM(v,w)));
                    
                end
            end
        end
    end
    
    MI(isnan(MI)) = 0;
    MI(isinf(MI)) = 0;

    
    MIValMshiftShuffled(i) = sum(MI(:));

end


plot(MIValFshift,'b'); hold on
plot(MIValMshift,'r'); 
plot(MIValFshiftShuffled,'c');
plot(MIValMshiftShuffled,'m');
title('MI with time shift');
xlabel('t-Tau')
ylabel('Mutual Info');




subplot(2,2,1); imagesc(dF); axis equal off
subplot(2,2,2); imagesc(abs(fSum)); axis equal off
subplot(2,2,3); imagesc(dM); axis equal off
subplot(2,2,4); imagesc(abs(mSum)); axis equal off



save('firstPassMIs.mat', 'MIValFshift','MIValFshiftShuffled','MIValMshift','MIValMshiftShuffled')
plot(MIValFshift,'c'); hold on
plot(MIValMshift,'m'); 
plot(MIValFshiftShuffled,'c');
plot(MIVal

%% Calculate Mutual Info

% I = sum (p(z1,z2,z3,z4) * log*(p(z1,z2,z3,z4)/pm(z1,z1)*pf(z3,z4))
% I(tau) = 
% introduce time lags I(tau), etc

% dm and df 
% 
numPoints = 201;
rangeVals = [-56 56];
xx = linspace(rangeVals(1),rangeVals(2),numPoints);

allF = combineCells(eVfemale);
allM = combineCells(eVmale);



[~, dM] = findPointDensity(allM,sigma,numPoints,rangeVals);
[~, dF] = findPointDensity(allF,sigma,numPoints,rangeVals);


tZ = -1000:1000;
mInfo = zeros(size(tZ));


dM = dM/sum(dM(:));
dF = dF/sum(dF(:));

MI = zeros(length(xx),length(yy),length(vv),length(ww));
for x = 1:length(xx)
    x
    for y = 1:length(yy)
        for v = 1:length(vv)
            for w = 1:length(ww)
                
                MI(x,y,v,w) = density(x,y,v,w)*log2(density(x,y,v,w)/(dF(x,y)*dM(v,w)));
                
            end
        end
    end
end


MI(isnan(MI)) = 0;
MI(isinf(MI)) = 0;


fSum = imrotate(squeeze(sum(sum(MI,3),4)),90);
mSum = squeeze(sum(sum(MI,1),2));


subplot(2,2,1); imagesc(dF); axis xy equal off
subplot(2,2,2); imagesc(fSum);
subplot(2,2,3); imagesc(dM);
subplot(2,2,4); imagesc(mSum);







%%


for i = 1:length(tZ)
    i

    tau = tZ(i);
    
    
    eM = zeros(size(allM,1)+length(tZ),2);
    eF = zeros(size(allF,1)+length(tZ),2);
    
    
    eM(1001+tau:size(allM)+1000+tau,:) = allM;
    eF(1001:size(allF)+1000,:) = allF;
    
    EE = [eF eM];

    
    Z = zeros(size(G));
    
    
    for j = 1:length(EE)
 
        currE = EE(j,:);
        [mV idx1] = min(abs(xx-currE(1)));
        [mV idx2] = min(abs(xx-currE(2)));
        [mV idx3] = min(abs(xx-currE(3)));
        [mV idx4] = min(abs(xx-currE(4)));
        Z(idx1, idx2, idx3, idx4) = Z(idx1, idx2, idx3, idx4)+1;
    end
    
    Z = Z./sum(Z(:));
    
    density = fftshift(real(ifftn(fftn(Gp).*fftn(Z))));
    density(density<0) = 0;
    
    
    
    
    sum1 = sum(sum(dM*dF));
    sumD = sum(density(:));
    
    mInfo(i) = sumD*log2(sumD/sum1);
end










%% Using just watersheds

WR1 = [F M];

WRZ = zeros(100,100);
for i = 1:length(WR1)
    
    WRZ(WR1(i,1)+1,WR1(i,2)+1) = WRZ(WR1(i,1)+1,WR1(i,2)+1)+1;
end

[x y idx] = find(WRZ>20);









