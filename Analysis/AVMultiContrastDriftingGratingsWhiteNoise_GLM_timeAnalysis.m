
clear;

experiment = 'EP010';
mouseID = 'AW159';
session = 'Session2';
date = '20201213-2';
stimPath = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'StimInfo',[date '_' mouseID '_AVmultiContrastDriftingGratingsWhiteNoise_stimInfo']);

analysisWindow = [-200 1200]; %ms relative to stimulus onset
frBinWidth = 10; %ms


dataFolder = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'*AV_driftingGratingsMultiContrast_whiteNoise*'));

movementFile = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'Video data','movementData.mat');
unitDataFile = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'AVMultiContrastDriftingGratingsWhiteNoise_Video_final','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');

newDir = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'AVMultiContrastDriftingGratingsWhiteNoise_GLM_TimeCourse');
if ~exist(newDir)
    mkdir(newDir);
end


load(stimPath);
uniqueEvents = size(stimInfo.index,1);
indices = stimInfo.index(:,1);
index = stimInfo.index;
repeats = stimInfo.repeats;
orientations = stimInfo.orientations;
contrasts = stimInfo.contrasts;
stimDuration = stimInfo.stimDuration;

load(unitDataFile);

load(movementFile);
frameRate = movementData.frameRate;

totalUnits = length(dataFiles);

binCount = analysisWindow(2)-analysisWindow(1)-frBinWidth+1; %sliding window
binEdges = analysisWindow(1)+frBinWidth:1:analysisWindow(2);

for n = 1:totalUnits
    neuronNumber = unitData(n).neuronNumber;
    if ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits) ||...
            ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) ||...
            ismember(neuronNumber,responsiveUnitsGLM.movementResponsiveUnits) ||...
            ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits)
        n
        
        nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
        
        eventTimes = nData.Events;
        spikeTimes = nData.SpikeData(1,:);
        
        
        glmCoefficients = zeros(7,binCount); %intercept + 3 unique predictors + their pairwise interactions
        for b=1:binCount
%             [n b]
            glmLight = zeros(uniqueEvents*repeats,1);
            glmSound = zeros(uniqueEvents*repeats,1);
            glmMove = zeros(uniqueEvents*repeats,1);
            glmFR = zeros(uniqueEvents*repeats,1);
            
            evCount = 0;
            for u = 1:uniqueEvents
                lightType = index(u,3);
                soundType = boolean(index(u,4));
                
                eventID = indices(u);
                eventsOfInterest = find(stimInfo.order==eventID);
                for ev = 1:length(eventsOfInterest)
                    evCount = evCount+1;
                    
                    time = eventTimes(eventsOfInterest(ev));
                    alignedSpikes = (spikeTimes-time)/1000; %temporal resolution of cheetah is in microseconds, converting to milliseconds
                    timeMax = binEdges(b);
                    timeMin = timeMax-frBinWidth;
                    
                    trialSpikes = length(find(alignedSpikes>timeMin & alignedSpikes<timeMax));
                    
                    evFrame = movementData.eventFrames(eventsOfInterest(ev));
                    binEndFrame = evFrame + round(binEdges(b)*frameRate/1000);
                    binStartFrame = evFrame + round((binEdges(b)-frBinWidth)*frameRate/1000);
                    movement = mean(movementData.frameMovement(binStartFrame:binEndFrame));
                    
                    glmLight(evCount) = lightType;
                    glmSound(evCount) = soundType;
                    glmMove(evCount) = movement;
                    glmFR(evCount) = trialSpikes;
                end
            end
            
            one = min(glmMove);
            two = prctile(glmMove,25);
            three = prctile(glmMove,50);
            four = prctile(glmMove,75);
            five = max(glmMove);
            glmMoveDisc = discretize(glmMove,[one two three four five]);
            
            glmOutput = fitglm([glmLight glmSound glmMoveDisc],glmFR,'interactions','distr','poisson','link','log');
            coeff = glmOutput.Coefficients.Estimate;
            glmCoefficients(:,b) = coeff;
            
        end
        unitGLMTimeCourse(n).coefficients = glmCoefficients;
    end
end

[moveQuants moveBins] = discretize(movementData.frameMovement,100);
movementStats.mode = moveBins(mode(moveQuants));
movementStats.minimum = min(movementData.frameMovement);
movementStats.lowerQuart = prctile(movementData.frameMovement,25);
movementStats.average = prctile(movementData.frameMovement,50);
movementStats.upperQuart = prctile(movementData.frameMovement,75);
movementStats.maximum = max(movementData.frameMovement);

% 
% lightTrace = 1000/frBinWidth * exp(glmCoefficients(1,:) + glmCoefficients(2,:)*1 + glmCoefficients(4,:)*1 + glmCoefficients(7,:)*1*1);
% lightSoundTrace = 1000/frBinWidth * exp(glmCoefficients(1,:) + glmCoefficients(2,:)*1 + glmCoefficients(3,:)*boolean(1) + glmCoefficients(4,:)*1 + glmCoefficients(5,:) + glmCoefficients(6,:)*boolean(1)*1 + glmCoefficients(7,:)*boolean(1)*1);
% figure;hold on;
% plot(binEdges,lightTrace);
% plot(binEdges,lightSoundTrace);
% 
% lightTrace = 1000/frBinWidth * exp(glmCoefficients(1,:) + glmCoefficients(2,:)*1 + glmCoefficients(4,:)*1 + glmCoefficients(7,:)*1*1);
% lightMoveTrace = 1000/frBinWidth * exp(glmCoefficients(1,:) + glmCoefficients(2,:)*1 + glmCoefficients(4,:)*3 + glmCoefficients(7,:)*1*3);
% figure;hold on;
% plot(binEdges,lightTrace);
% plot(binEdges,lightMoveTrace);
% 
% lightTrace = 1000/frBinWidth * exp(glmCoefficients(1,:) + glmCoefficients(2,:)*0 + glmCoefficients(4,:)*1 + glmCoefficients(7,:)*0*1);
% lightHigh = 1000/frBinWidth * exp(glmCoefficients(1,:) + glmCoefficients(2,:)*0.25 + glmCoefficients(4,:)*1 + glmCoefficients(7,:)*0.25*1);
% figure;hold on;
% plot(binEdges,lightTrace);
% plot(binEdges,lightHigh);

glmAnalysisParams.binCount = binCount;
glmAnalysisParams.binEdges = binEdges;
glmAnalysisParams.frBinWidth = frBinWidth;
glmAnalysisParams.analysiWindow = analysisWindow;

save(fullfile(newDir,'glmTimeCourseData_10msBin.mat'),'unitGLMTimeCourse','movementStats','glmAnalysisParams');

lightTrain = zeros(0,binCount);
lightSoundTrain = zeros(0,binCount);
lightMoveTrain = zeros(0,binCount);
lightSoundMoveTrain = zeros(0,binCount);
for n = 1:totalUnits
    neuronNumber = unitData(n).neuronNumber;
    if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
            (ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits) &&...
            ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits))
        
        coeff = unitGLMTimeCourse(n).coefficients;
        
        lTrain = 1000/frBinWidth * exp(coeff(1,:) + coeff(2,:)*1 + coeff(3,:)*boolean(0) + coeff(4,:)*1 + coeff(5,:)*1*boolean(0) + coeff(6,:)*1*1 + coeff(7,:)*boolean(0)*1);
        lsTrain = 1000/frBinWidth * exp(coeff(1,:) + coeff(2,:)*1 + coeff(3,:)*boolean(1) + coeff(4,:)*1 + coeff(5,:)*1*boolean(1) + coeff(6,:)*1*1 + coeff(7,:)*boolean(1)*1);
        lmTrain = 1000/frBinWidth * exp(coeff(1,:) + coeff(2,:)*1 + coeff(3,:)*boolean(0) + coeff(4,:)*3 + coeff(5,:)*1*boolean(0) + coeff(6,:)*1*3 + coeff(7,:)*boolean(0)*3);
        lsmTrain = 1000/frBinWidth * exp(coeff(1,:) + coeff(2,:)*1 + coeff(3,:)*boolean(1) + coeff(4,:)*3 + coeff(5,:)*1*boolean(1) + coeff(6,:)*1*3 + coeff(7,:)*boolean(1)*3);
        
        lightTrain = [lightTrain; lTrain];
        lightSoundTrain = [lightSoundTrain; lsTrain];
        lightMoveTrain = [lightMoveTrain; lmTrain];
        lightSoundMoveTrain = [lightSoundMoveTrain; lsmTrain];
    end
end

figure;hold on;
plot(binEdges,mean(lightTrain));
plot(binEdges,mean(lightSoundTrain));
legend({'Light','Light + sound'});

figure;hold on;
plot(binEdges,mean(lightSoundTrain));
plot(binEdges,mean(lightSoundMoveTrain));
legend({'Light+sound','Light+sound+move'});

figure;hold on;
plot(binEdges,mean(lightTrain));
plot(binEdges,mean(lightMoveTrain));
plot(binEdges,mean(lightSoundMoveTrain));
legend({'Light','Light + move','Light+sound+move'});