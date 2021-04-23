
clear all;
% close all;

dataPaths{1} = fullfile('EP004','AW117','20200221-2');
dataPaths{2} = fullfile('EP004','AW118','20200221-1');
dataPaths{3} = fullfile('EP004','AW121','20200226-1');
dataPaths{4} = fullfile('EP004','AW124','20200303-2');
dataPaths{5} = fullfile('EP010','AW157','20201212-2');
dataPaths{6} = fullfile('EP010','AW158','20201212-2');

unitDPandU = [3 41; 3 133; 3 5; 4 30; 2 54; 1 15; 5 22;...
    6 1; 6 106; 6 35; 6 41; 6 59; 6 100];


%% Load response data of units of interest
disp('Creating masterLog...');
unitsOfInterest = size(unitDPandU,1);
for uoiU = 1:unitsOfInterest;
    dp = unitDPandU(uoiU,1);
    u = unitDPandU(uoiU,2);
    
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_final','AVMultiContrastDriftingGratingsWhiteNoiseData_selectivity.mat');
    load(dataFile);
    
    masterLog(uoiU) = unitData(u);
end

stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{1},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
load(fullfile(stimPath(end).folder,stimPath(end).name));
orientations = stimInfo.orientations;
contrasts = stimInfo.contrasts;
repeats = stimInfo.repeats;
quantScalar = 1000/(analysisParams.quantWindow(2)-analysisParams.quantWindow(1));


%% Run population decoder
choiceMap = cell(2,length(contrasts));

for c = 1:length(contrasts)
    tempChoiceMap = zeros(2,length(orientations),length(orientations));
    
    for o = 1:length(orientations)
        
        mleClassification = zeros(2,length(orientations));
        randomClassification = zeros(2,length(orientations));
        
        for trial = 1:repeats
            [c o trial]
            
            if exist('neuronStats','var')
                clear neuronStats randStats
            end
            probeData = [];
            
            for neuron = 1:unitsOfInterest
                
                for orient = 1:length(orientations)
                    vInd = (c-1)*length(orientations) + orient;
                    avInd = (c-1)*length(orientations) + orient +length(orientations)*length(contrasts);
                    
                    trialsIncluded = 1:repeats;
                    if orient==o
                        trialsIncluded  = setdiff(trialsIncluded,trial);
                        
                        probeTrialV = masterLog(neuron).trialResponse(vInd,trial)/quantScalar;
                        probeTrialAV = masterLog(neuron).trialResponse(avInd,trial)/quantScalar;
                        probeData = [probeData [probeTrialV; probeTrialAV]];
                    end
                    
                    vData = masterLog(neuron).trialResponse(vInd,trialsIncluded)/quantScalar;
                    avData = masterLog(neuron).trialResponse(avInd,trialsIncluded)/quantScalar;
                    
                    neuronStats(neuron,orient,1) = fitdist(vData','Poisson');
                    neuronStats(neuron,orient,2) = fitdist(avData','Poisson');
                end
            end
            
            vEstimate = maximumLikelihoodFunction(neuronStats(:,:,1),probeData(1,:));
            mleClassification(1,vEstimate) = mleClassification(1,vEstimate) +1;
            
            avEstimate = maximumLikelihoodFunction(neuronStats(:,:,2),probeData(2,:));
            mleClassification(2,avEstimate) = mleClassification(2,avEstimate) +1;
            
        end
        
        tempChoiceMap(1,o,:) = mleClassification(1,:) / sum(mleClassification(1,:));
        tempChoiceMap(2,o,:) = mleClassification(2,:) / sum(mleClassification(2,:));
    end
    choiceMap{1,c} = cat(3,choiceMap{1,c},squeeze(tempChoiceMap(1,:,:)));
    size(choiceMap{1,c});
    choiceMap{2,c} = cat(3,choiceMap{2,c},squeeze(tempChoiceMap(2,:,:)));
end


%% Plot data

meanMap = squeeze(choiceMap{1,5});
figure;imagesc(meanMap); 
for c = 1:5
    clear tempV tempAV tempDiff

        tempV = mean(diag(choiceMap{1,c}));
        tempAV = mean(diag(choiceMap{2,c}));
        tempDiff = tempAV - tempV;

    vAcc(c) = (tempV);
    avAcc(c) = (tempAV);
    diffAcc(c) = (tempDiff);

end
figure;hold on;
plot(contrasts,vAcc,'Color',[0 0 0]);
plot(contrasts,avAcc,'Color',[0 0 1]);

figure;hold on
plot(contrasts,diffAcc,'Color',[0 0 0]);

