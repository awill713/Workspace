
clear all;
% close all;

saveDir = fullfile('E:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise');
if ~exist(saveDir)
    mkdir(saveDir);
end
selectivityDir = fullfile('E:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise','D prime - individual neuron');


dataPaths{1} = fullfile('EP004','AW117','20200221-1');
dataPaths{2} = fullfile('EP004','AW117','20200221-2');
dataPaths{3} = fullfile('EP004','AW118','20200221-1');
dataPaths{4} = fullfile('EP004','AW118','20200221-2');
dataPaths{5} = fullfile('EP004','AW121','20200226-1');
dataPaths{6} = fullfile('EP004','AW121','20200226-2');
dataPaths{7} = fullfile('EP004','AW124','20200303-1');
dataPaths{8} = fullfile('EP004','AW124','20200303-2');

neuronCount = 20;
%which neurons to include
onlySingleUnits = 0;
allNeurons = 0;
lightResponsive = 1;
orientationSelective = 1;
soundResponsive = 1;

%Sort units based on OSI and DSI
disp('Sorting units by OSI and DSI...');
selectivityFile = fullfile(selectivityDir,'SelectivityIndices.mat');
load(selectivityFile);
for c = 1:size(dsiMat,2)
    sortedDSI{c} = sortrows(dsiMat{c},1,'descend');
    sortedOSI{c} = sortrows(osiMat{c},1,'descend');
end

%Create masterLog
disp('Creating masterLog...');
for dp = 1:length(dataPaths)
    %     [randomize dp]
    %Load data
    dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Recat','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('E:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    temp = struct2cell(unitData);
    masterLog{dp} = squeeze(temp(6,1,:))';
end


orientations = stimInfo.orientations;
contrasts = stimInfo.contrasts;
repeats = stimInfo.repeats;
accuracyTrajectory = zeros(length(contrasts),2,neuronCount);
for c = 1:length(contrasts)
    tempTrajectory = zeros(2,neuronCount);
    for nc = 1:neuronCount
        disp(['Contrast = ' num2str(c) ', ' num2str(nc) '/' num2str(neuronCount) ' units']);
        
        correctTrialsV = 0;
        correctTrialsAV = 0;
        totalTrials = 0;
        for o = 1:length(orientations)
            for trial = 1:repeats
                disp(['Contrast = ' num2str(c) ', ' num2str(nc) '/' num2str(neuronCount) ' units, o = ' num2str(o) ', trial ' num2str(trial)]);
        
                if exist('neuronStats','var')
                    clear neuronStats
                end
                probeData = [];
                
                for neuron = 1:nc
                    neuronDP = sortedDSI{c}(neuron,5);
                    neuronU = sortedDSI{c}(neuron,7);
                    
                    for orient = 1:length(orientations)
                        vInd = (c-1)*length(orientations) + orient;
                        avInd = (c-1)*length(orientations) + orient +length(orientations)*length(contrasts);
                        
                        trialsIncluded = 1:repeats;
                        if orient==o
                            trialsIncluded  = setdiff(trialsIncluded,trial);
                            
                            probeTrialV = masterLog{neuronDP}{neuronU}(vInd,trial);
                            probeTrialAV = masterLog{neuronDP}{neuronU}(avInd,trial);
                            probeData = [probeData [probeTrialV; probeTrialAV]];
                        end
                        
                        vData = masterLog{neuronDP}{neuronU}(vInd,trialsIncluded);
                        avData = masterLog{neuronDP}{neuronU}(avInd,trialsIncluded);
                        
                        neuronStats(neuron,orient,1) = fitdist(vData','Normal');
                        neuronStats(neuron,orient,2) = fitdist(avData','Normal');
                    end
                end
                
                vEstimate = maximumLikelihoodFunction(neuronStats(:,:,1),probeData(1,:));
                if vEstimate==o
                    correctTrialsV = correctTrialsV+1;
                end
                
                avEstimate = maximumLikelihoodFunction(neuronStats(:,:,2),probeData(2,:));
                if avEstimate==o
                    correctTrialsAV = correctTrialsAV+1;
                end
                totalTrials = totalTrials+1;
            end
        end
        tempTrajectory(1,nc) = correctTrialsV / totalTrials;
        tempTrajectory(2,nc) = correctTrialsAV / totalTrials;
        
    end
    accuracyTrajectory(c,:,:) = tempTrajectory;
    
   
end

for c = 1:length(contrasts)
    figure;
    hold on;
    plot(squeeze(accuracyTrajectory(c,1,:)),'Color',[0 0 0]);
    plot(squeeze(accuracyTrajectory(c,2,:)),'Color',[0 0 1]);
    xlabel('Number of units included');
    ylabel('MLE model accuracy');
    title(['Contrast ' num2str(contrasts(c))]);
end

meanMap = squeeze(mean(choiceMap{1,5},3));
figure;imagesc(meanMap); 
statsMat = [];
for c = 1:5
    clear tempV tempAV tempDiff
    for s=1:size(choiceMap{1,c},3)
        tempV(s) = mean(diag(choiceMap{1,c}(:,:,s)));
        tempAV(s) = mean(diag(choiceMap{2,c}(:,:,s)));
        tempDiff(s) = tempAV(s) - tempV(s);
    end
    vAcc(c) = mean(tempV);
    vAccStd(c) = std(tempV)/sqrt(length(tempV));
    avAcc(c) = mean(tempAV);
    avAccStd(c) = std(tempAV)/sqrt(length(tempAV));
    diffAcc(c) = mean(tempDiff);
    diffAccStd(c) = std(diffAcc)/sqrt(length(diffAcc));
    diffAccStd(c) = std(tempDiff)/sqrt(length(tempDiff));
    
    statsMat = [statsMat; [tempV' tempAV']];
    
%     [c vAcc avAcc]
end
figure;hold on;
hA = area(contrasts,[vAcc-vAccStd; 2*vAccStd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[avAcc-avAccStd; 2*avAccStd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];

plot(contrasts,vAcc,'Color',[0 0 0]);
plot(contrasts,avAcc,'Color',[0 0 1]);

figure;hold on;
hA = area(contrasts,[diffAcc-diffAccStd; 2*diffAccStd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,diffAcc,'Color',[0 0 0]);

[p,~,~] = anova2(statsMat,8,'off')

