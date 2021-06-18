
clear all;
% close all;

saveDir = fullfile('F:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise_final');
if ~exist(saveDir)
    mkdir(saveDir);
end

% randomizations = 20;
% neuronCount = 87;
randomizations = 10;
neuronCount = 87;


dataPaths{1} = fullfile('EP004','AW117','20200221-1');
dataPaths{2} = fullfile('EP004','AW117','20200221-2');
dataPaths{3} = fullfile('EP004','AW118','20200221-1');
dataPaths{4} = fullfile('EP004','AW118','20200221-2');
dataPaths{5} = fullfile('EP004','AW121','20200226-1');
dataPaths{6} = fullfile('EP004','AW121','20200226-2');
dataPaths{7} = fullfile('EP004','AW124','20200303-1');
dataPaths{8} = fullfile('EP004','AW124','20200303-2');
dataPaths{9} = fullfile('EP010','AW157','20201212-1');
dataPaths{10} = fullfile('EP010','AW157','20201212-2');
dataPaths{11} = fullfile('EP010','AW158','20201212-1');
dataPaths{12} = fullfile('EP010','AW158','20201212-2');

% unitDPandU = [5 41; 5 133; 5 5; 8 30; 3 54; 2 15; 10 22;...
%     12 1; 12 106; 12 35; 12 41; 12 59; 12 100];



%% Load response data of units of interest
disp('Creating masterLog...');

% unitsOfInterest = size(unitDPandU,1);
% for uoiU = 1:unitsOfInterest
%     dp = unitDPandU(uoiU,1);
%     neuron = unitDPandU(uoiU,2);
%     
%     dataFile = fullfile('F:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_final','AVMultiContrastDriftingGratingsWhiteNoiseData_selectivity.mat');
%     load(dataFile);
%     
%     neuronNumber = unitData(neuron).neuronNumber;
%     
%     masterLog(uoiU).trialResponse = unitData(neuron).trialResponse;
% end


unitsOfInterest = 0;
for dp = 1:length(dataPaths)
    dataFile = fullfile('F:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_final','AVMultiContrastDriftingGratingsWhiteNoiseData_selectivity.mat');
    load(dataFile);

    for u = 1:length(unitData)
        neuronNumber = unitData(u).neuronNumber;
        if (ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
                (ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) &&...
                ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits))) &&...
                (ismember(neuronNumber,responsiveUnits.orientationSelectiveUnitsAND) ||...
                ismember(neuronNumber,responsiveUnits.directionSelectiveUnitsAND)) &&...
                ismember(neuronNumber,[responsiveUnits.singleUnits responsiveUnits.multiUnits])

            unitsOfInterest = unitsOfInterest+1;
            masterLog(unitsOfInterest).trialResponse = unitData(u).trialResponse;
        end
    end
end

stimPath = dir(fullfile('F:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
load(fullfile(stimPath(end).folder,stimPath(end).name));
orientations = stimInfo.orientations;
contrasts = stimInfo.contrasts;
repeats = stimInfo.repeats;
quantScalar = 1000/(analysisParams.quantWindow(2)-analysisParams.quantWindow(1));
iMat = logical(eye(length(orientations)));
orientMat = circshift(iMat,3,2) | circshift(iMat,-3,2);
direcMat = circshift(iMat,6,2);

%% Run MLE decoder

tic
accuracyTrajectory = zeros(length(contrasts),2,neuronCount);
accuracySTD = zeros(length(contrasts),2,neuronCount);
choiceMap = zeros(length(contrasts),2,length(orientations),length(orientations));
orientationAcc = zeros(2,length(contrasts),neuronCount,randomizations);
directionAcc = zeros(2,length(contrasts),neuronCount,randomizations);
overallAcc = zeros(2,length(contrasts),neuronCount,randomizations);
for c = 1:length(contrasts)
    toc
    trajMean = zeros(2,neuronCount);
    trajSTD = zeros(2,neuronCount);
    
    
    for nc = 1:neuronCount
        repeatAccuracy = zeros(2,randomizations);
        tempChoiceMap = zeros(2,randomizations,length(orientations),length(orientations));
        
        for rando = 1:randomizations
            [c nc rando]
            toc
            randoUnits = randsample(1:unitsOfInterest,nc,'true');
            
            vInd = (c-1)*length(orientations);
            avInd = vInd + length(orientations)*length(contrasts);
            
            accuracyMat = zeros(2,length(orientations),length(orientations));
            for iDir = 1:length(orientations)
                for jDir = 1:length(orientations)
                    %                     [c nc rando iDir jDir]
                    
                    accuracy = zeros(2,repeats);
                    for trial = 1:repeats
                        neuronStats = zeros(nc,2,2);
                        probeData = [];
                        
                        iTestTrials = trial;
                        iTrainingTrials = setdiff(1:repeats,iTestTrials);
                        jTrainingTrials = 1:repeats;
                        
                        for neuron = 1:nc
                            chosenUnit = randoUnits(neuron);
                            
                            probeTrialV = masterLog(chosenUnit).trialResponse(vInd+iDir,trial)/quantScalar;
                            probeTrialAV = masterLog(chosenUnit).trialResponse(avInd+iDir,trial)/quantScalar;
                            probeData = [probeData [probeTrialV; probeTrialAV]];
                            
                            vData = masterLog(chosenUnit).trialResponse(vInd+iDir,iTrainingTrials)/quantScalar;
                            avData = masterLog(chosenUnit).trialResponse(avInd+iDir,iTrainingTrials)/quantScalar;
                            
                            neuronStats(neuron,1,1) = mean(vData);
                            neuronStats(neuron,1,2) = mean(avData);
                            
                            vData = masterLog(chosenUnit).trialResponse(vInd+jDir,jTrainingTrials)/quantScalar;
                            avData = masterLog(chosenUnit).trialResponse(avInd+jDir,jTrainingTrials)/quantScalar;
                            
                            neuronStats(neuron,2,1) = mean(vData);
                            neuronStats(neuron,2,2) = mean(avData);
                        end
                        
                        distData = neuronStats(:,:,1);
                        likelihoods = zeros(1,size(distData,2));
                        probe = probeData(1,:);
                        for d = 1:length(likelihoods)
                            prob = 0;
                            for nn = 1:size(distData,1)
                                post = pdf('Poisson',probe(nn),distData(nn,d));
                                if post>0
                                    prob = log(post) + prob;
                                end
                                [nn prob];
                            end
                            likelihoods(d) = prob;
                        end
                        [~, vEstimate] = max(likelihoods);
                        accuracy(1,trial) = vEstimate==1; %1 is iDir answer
                        
                        
                        distData = neuronStats(:,:,2);
                        likelihoods = zeros(1,size(distData,2));
                        probe = probeData(2,:);
                        for d = 1:length(likelihoods)
                            prob = 0;
                            for nn = 1:size(distData,1)
                                post = pdf('Poisson',probe(nn),distData(nn,d));
                                if post>0
                                    prob = log(post) + prob;
                                end
                                [nn prob];
                            end
                            likelihoods(d) = prob;
                        end
                        [~, avEstimate] = max(likelihoods);
                        accuracy(2,trial) = avEstimate==1; %1 is iDir answer
                        
                    end
                    
                    accuracyMat(1,iDir,jDir) = mean(accuracy(1,:));
                    accuracyMat(2,iDir,jDir) = mean(accuracy(2,:));
                end
            end
            squeezedV = squeeze(accuracyMat(1,:,:));
            squeezedAV = squeeze(accuracyMat(2,:,:));
            
            tempChoiceMap(1,rando,:,:) = squeezedV;
            tempChoiceMap(2,rando,:,:) = squeezedAV;
            
            repeatAccuracy(1,rando) = mean(squeezedV(~iMat));
            repeatAccuracy(2,rando) = mean(squeezedAV(~iMat));
            
            orientationAcc(1,c,nc,rando) = mean(squeezedV(orientMat));
            orientationAcc(2,c,nc,rando) = mean(squeezedAV(orientMat));
            
            directionAcc(1,c,nc,rando) = mean(squeezedV(direcMat));
            directionAcc(2,c,nc,rando) = mean(squeezedAV(direcMat));
            
            overallAcc(1,c,nc,rando) = mean(squeezedV(~iMat));
            overallAcc(2,c,nc,rando) = mean(squeezedAV(~iMat));
        end
        
        trajMean(:,nc) = [nanmean(repeatAccuracy(1,:)) nanmean(repeatAccuracy(2,:))]';
        trajSTD(:,nc) = [std(repeatAccuracy(1,:))/sqrt(randomizations) std(repeatAccuracy(2,:))/sqrt(randomizations)]';
        
    end
    
    accuracyTrajectory(c,:,:) = trajMean;
    accuracySTD(c,:,:) = trajSTD;
    
    choiceMap(c,1,:,:) = squeeze(mean(tempChoiceMap(1,:,:,:),2));
    choiceMap(c,2,:,:) = squeeze(mean(tempChoiceMap(2,:,:,:),2));
end
toc


blueMap = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1)];
redMap = [ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
map = [blueMap; redMap];
maxVal = 0;
for c = 1:length(contrasts)
    dif = abs(choiceMap(c,2,:,:) - choiceMap(c,1,:,:));
    maxVal = max([maxVal max(dif(:))]);
end
f1 = figure;
set(f1,'Position',[300 150 1000 600]);
for c = 1:length(contrasts)
    vv = squeeze(choiceMap(c,1,:,:));
    av = squeeze(choiceMap(c,2,:,:));
    
    subplot(4,length(contrasts),c);
    imagesc(vv,[0 1]);
    cm = colormap(gca,'parula');
    caxis([0 1]);
    xticks(1:3:12);
    yticks(1:3:12);
    xticklabels(orientations(xticks));
    yticklabels(orientations(yticks));
    title(['C = ' num2str(contrasts(c)) ', vis']);
    
    subplot(4,length(contrasts),c+length(contrasts));
    imagesc(av,[0 1]);
    cm = colormap(gca,'parula');
    caxis([0 1]);
    xticks(1:3:12);
    yticks(1:3:12);
    xticklabels(orientations(xticks));
    yticklabels(orientations(yticks));
    title(['C = ' num2str(contrasts(c)) ', audiovis']);
    
    subplot(4,length(contrasts),c+2*length(contrasts));
    imagesc(av-vv,[-maxVal maxVal]);
    colormap(gca,map);
    caxis([-maxVal maxVal]);
    xticks(1:3:12);
    yticks(1:3:12);
    xticklabels(orientations(xticks));
    yticklabels(orientations(yticks));
    title(['C = ' num2str(contrasts(c)) ', difference']);
    
    subplot(4,length(contrasts),c+3*length(contrasts));hold on;
    iMat = logical(eye(length(orientations)));
    scatter(vv(~iMat),av(~iMat));
    line([0 1],[0 1],'Color',[0 0 0]);
    xlim([0 1]);ylim([0 1]);
    delt = mean(av(~iMat)-vv(~iMat));
    [h p] = ttest(vv(~iMat),av(~iMat));
    xlabel('SVM accuracy (visual)');
    ylabel('SVM accuracy (audiovisual)');
    title(['\Delta = ' num2str(delt) ', p = ' num2str(p)]);
end





maxAcc = squeeze(accuracyTrajectory(:,:,end))';
maxStd = squeeze(accuracySTD(:,:,end))';
lineV = maxAcc(1,:);
lineAV = maxAcc(2,:);
areaV = maxStd(1,:);
areaAV = maxStd(2,:);

figure;hold on;
hA = area(contrasts,[lineV-areaV; 2*areaV]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[lineAV-areaAV; 2*areaAV]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];

plot(contrasts,lineV,'Color',[0 0 0]);
plot(contrasts,lineAV,'Color',[0 0 1]);



orientTrajV = squeeze(mean(orientationAcc(1,:,neuronCount,:),4));
orientSTDv = squeeze(std(orientationAcc(1,:,neuronCount,:),[],4))./sqrt(randomizations);
orientTrajAV = squeeze(mean(orientationAcc(2,:,neuronCount,:),4));
orientSTDav = squeeze(std(orientationAcc(2,:,neuronCount,:),[],4))./sqrt(randomizations);

figure;hold on;
hA = area(contrasts,[orientTrajV-orientSTDv; 2*orientSTDv]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[orientTrajAV-orientSTDav; 2*orientSTDav]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];

plot(contrasts,orientTrajV,'Color',[0 0 0]);
plot(contrasts,orientTrajAV,'Color',[0 0 1]);


direcTrajV = squeeze(mean(directionAcc(1,:,neuronCount,:),4));
direcSTDv = squeeze(std(directionAcc(1,:,neuronCount,:),[],4))./sqrt(randomizations);
direcTrajAV = squeeze(mean(directionAcc(2,:,neuronCount,:),4));
direcSTDav = squeeze(std(directionAcc(2,:,neuronCount,:),[],4))./sqrt(randomizations);

figure;hold on;
hA = area(contrasts,[direcTrajV-direcSTDv; 2*direcSTDv]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[direcTrajAV-direcSTDav; 2*direcSTDav]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];

plot(contrasts,direcTrajV,'Color',[0 0 0]);
plot(contrasts,direcTrajAV,'Color',[0 0 1]);