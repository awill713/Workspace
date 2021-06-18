
clear all;
% close all;

saveDir = fullfile('F:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise_final');
if ~exist(saveDir)
    mkdir(saveDir);
end

randomizations = 10;
neuronCount = 87;
% randomizations = 5;
% neuronCount = 12;


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
% unitLog = [];
% for uoiU = 1:unitsOfInterest
%     dp = unitDPandU(uoiU,1);
%     u = unitDPandU(uoiU,2);
%
%     dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_final','AVMultiContrastDriftingGratingsWhiteNoiseData_selectivity.mat');
%     load(dataFile);
%
%     neuronNumber = unitData(u).neuronNumber;
%
%     masterLog(uoiU) = unitData(u);
%     unitLog = [unitLog; dp u neuronNumber];
% end


unitsOfInterest = 0;
% unitLog = [];
% masterLog = zeros(120,10,0);
for dp = 1:length(dataPaths)
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_final','AVMultiContrastDriftingGratingsWhiteNoiseData_selectivity.mat');
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
            %             masterLog = cat(3,masterLog,unitData(u).trialResponse);
            %             unitLog = [unitLog; dp u neuronNumber];
        end
    end
end

stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
load(fullfile(stimPath(end).folder,stimPath(end).name));
orientations = stimInfo.orientations;
contrasts = stimInfo.contrasts;
repeats = stimInfo.repeats;
quantScalar = 1000/(analysisParams.quantWindow(2)-analysisParams.quantWindow(1));


%% Run MLE decoder
tic
accuracyTrajectory = zeros(length(contrasts),2,neuronCount);
accuracySTD = zeros(length(contrasts),2,neuronCount);
choiceMap = zeros(length(contrasts),2,length(orientations),length(orientations));
orientationAcc = zeros(2,length(contrasts),neuronCount,randomizations);
directionAcc = zeros(2,length(contrasts),neuronCount,randomizations);
overallAcc = zeros(2,length(contrasts),neuronCount,randomizations);
for c = 1:length(contrasts)
    trajMean = zeros(2,neuronCount);
    trajSTD = zeros(2,neuronCount);
    
    for nc = 1:neuronCount
        repeatAccuracy = zeros(2,randomizations);
        
        for rando = 1:randomizations
            %             disp(['Contrast = ' num2str(c) ', ' num2str(nc) '/' num2str(neuronCount) ' units, repeat ' num2str(rando)]);
            [c nc rando]
            randoUnits = randsample(1:unitsOfInterest,nc,'true');
            
            %             for rr = 1:length(randoUnits)
            %                 neuronsIncluded(c,nc,rando,rr,:) = unitLog(randoUnits(rr),:);
            %             end
            
            correctTrialsV = 0;
            correctTrialsAV = 0;
            totalTrials = 0;
            for o = 1:length(orientations)
                for trial = 1:repeats
                    %                     disp(['Contrast = ' num2str(c) ', ' num2str(nc) '/' num2str(neuronCount) ' units, repeat ' num2str(rando), ', o = ' num2str(o) ', trial ' num2str(trial)]);
                    
                    %                     if exist('neuronStats','var')
                    %                         clear neuronStats
                    %                     end
                    
                    %                     neuronStats = zeros(nc,length(orientations),2);
                    neuronStats = zeros(nc,2,2);
                    probeData = [];
                    
                    for neuron = 1:nc
                        chosenUnit = randoUnits(neuron);
                        
                        orthInd = [mod(o+3-1,length(orientations)) mod(o-3-1,length(orientations))]+1;
                        
                        orientsOfInterest = [o orthInd(1)];
                        for oo = 1:2 %orientOfInterest and ortho
                            
                            if oo ==1 %orientOfInterest
                                orient = orientsOfInterest(oo);
                                vInd = (c-1)*length(orientations) + orient;
                                avInd = (c-1)*length(orientations) + orient +length(orientations)*length(contrasts);
                                
                                trialsIncluded  = setdiff(1:repeats,trial);
                                probeTrialV = masterLog(chosenUnit).trialResponse(vInd,trial)/quantScalar;
                                probeTrialAV = masterLog(chosenUnit).trialResponse(avInd,trial)/quantScalar;
                                probeData = [probeData [probeTrialV; probeTrialAV]];
                                
                                vData = masterLog(chosenUnit).trialResponse(vInd,trialsIncluded)/quantScalar;
                                avData = masterLog(chosenUnit).trialResponse(avInd,trialsIncluded)/quantScalar;
                                
                                neuronStats(neuron,oo,1) = median(vData);
                                neuronStats(neuron,oo,2) = median(avData);
                            else %ortho
                                orient = orthInd;
                                vInd = (c-1)*length(orientations) + orient;
                                avInd = (c-1)*length(orientations) + orient +length(orientations)*length(contrasts);
                                
                                trialsIncluded  = 1:repeats;
                                
                                vData = reshape(masterLog(chosenUnit).trialResponse(vInd,trialsIncluded)/quantScalar,[1 2*repeats]);
                                avData = reshape(masterLog(chosenUnit).trialResponse(avInd,trialsIncluded)/quantScalar,[1 2*repeats]);
                                
                                neuronStats(neuron,oo,1) = median(vData);
                                neuronStats(neuron,oo,2) = median(avData);
                            end
                        end
                    end
                    
                    %                     vEstimate = maximumLikelihoodFunction(neuronStats(:,:,1),probeData(1,:));
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
                    if vEstimate==1
                        correctTrialsV = correctTrialsV+1;
                    end
                    
                    %                     avEstimate = maximumLikelihoodFunction(neuronStats(:,:,2),probeData(2,:));
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
                    if avEstimate==1
                        correctTrialsAV = correctTrialsAV+1;
                    end
                    totalTrials = totalTrials+1;
                end
            end
            repeatAccuracy(1,rando) = correctTrialsV / totalTrials;
            repeatAccuracy(2,rando) = correctTrialsAV / totalTrials;
        end
        
        trajMean(:,nc) = [mean(repeatAccuracy(1,:)) mean(repeatAccuracy(2,:))]';
        trajSTD(:,nc) = [std(repeatAccuracy(1,:))/sqrt(randomizations) std(repeatAccuracy(2,:))/sqrt(randomizations)]';
        
    end
    accuracyTrajectory(c,:,:) = trajMean;
    accuracySTD(c,:,:) = trajSTD;
    
    
end
toc

for c = 1:length(contrasts)
    lineV = squeeze(accuracyTrajectory(c,1,:))';
    lineAV = squeeze(accuracyTrajectory(c,2,:))';
    areaV = squeeze(accuracySTD(c,1,:))';
    areaAV = squeeze(accuracySTD(c,2,:))';
    
    figure;hold on;
    hA = area(1:neuronCount,[lineV-areaV; 2*areaV]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    hA = area(1:neuronCount,[lineAV-areaAV; 2*areaAV]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 1];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    
    plot(1:neuronCount,lineV,'Color',[0 0 0]);
    plot(1:neuronCount,lineAV,'Color',[0 0 1]);
    
    title(['Contrast ' num2str(contrasts(c))]);
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
