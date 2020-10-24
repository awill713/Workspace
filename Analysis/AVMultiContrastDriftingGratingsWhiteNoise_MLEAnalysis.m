
clear all;
% close all;

saveDir = fullfile('D:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise');
if ~exist(saveDir)
    mkdir(saveDir);
end

dataPaths{1} = fullfile('EP004','AW117','20200221-1');
dataPaths{2} = fullfile('EP004','AW117','20200221-2');
dataPaths{3} = fullfile('EP004','AW118','20200221-1');
dataPaths{4} = fullfile('EP004','AW118','20200221-2');
dataPaths{5} = fullfile('EP004','AW121','20200226-1');
dataPaths{6} = fullfile('EP004','AW121','20200226-2');
dataPaths{7} = fullfile('EP004','AW124','20200303-1');
dataPaths{8} = fullfile('EP004','AW124','20200303-2');

%which neurons to include
onlySingleUnits = 0;
allNeurons = 0;
lightResponsive = 1;
orientationSelective = 1;
soundResponsive = 1;


for dp = 1:length(dataPaths)
    %     [randomize dp]
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Recat','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    repeats = stimInfo.repeats;
    
    if ~exist('choiceMap','var')
        choiceMap = cell(2,length(contrasts)); %v and av, by contrasts
        randomMap = cell(2,length(contrasts)); %v and av, by contrasts, trials shuffled
    end
    
    %     analysisWindow = analysisParams.analysisWindow;
    %     quantWindow = analysisParams.quantWindow;
    %     quantWindow = [0 250];
    %     frBinWidth = analysisParams.frBinWidth;
    %     quantFirstBin = quantWindow(1)-analysisWindow(1)-frBinWidth+1;
    %     quantLastBin = quantWindow(2)-analysisWindow(1)-frBinWidth+1;
    
    neuronCount = 0;
    neuronNumberLog = [];
    for u  = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        if (allNeurons || (lightResponsive && ismember(neuronNumber,responsiveUnits.lightResponsiveUnits))...
                || (orientationSelective && ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits))...
                || (soundResponsive && ismember(neuronNumber,responsiveUnits.soundResponsiveUnits)))...
                && (~onlySingleUnits || unitData(u).type==1) && unitData(u).type~=0
            
            neuronCount = neuronCount+1;
            neuronNumberLog = [neuronNumberLog; neuronCount neuronNumber u];
        end
    end
    
    if neuronCount>0
        
        for c = 1:length(contrasts)
            tempChoiceMap = zeros(2,length(orientations),length(orientations));
            
            for o = 1:length(orientations)
                [dp c o]
                mleClassification = zeros(2,length(orientations));
                randomClassification = zeros(2,length(orientations));
                
                
                
                for trial = 1:repeats
                
                    if exist('neuronStats','var')
                        clear neuronStats randStats
                    end
                    probeData = [];
                    
                    for neuron = 1:neuronCount
                        neuronU = neuronNumberLog(neuron,3);
                        
                        for orient = 1:length(orientations)
                            vInd = (c-1)*length(orientations) + orient;
                            avInd = (c-1)*length(orientations) + orient +length(orientations)*length(contrasts);
                            
                            trialsIncluded = 1:repeats;
                            if orient==o
                                trialsIncluded  = setdiff(trialsIncluded,trial);
                                
                                probeTrialV = unitData(neuronU).trialResponse(vInd,trial);
                                probeTrialAV = unitData(neuronU).trialResponse(avInd,trial);
                                probeData = [probeData [probeTrialV; probeTrialAV]];
                            end
                            
                            vData = unitData(neuronU).trialResponse(vInd,trialsIncluded);
                            avData = unitData(neuronU).trialResponse(avInd,trialsIncluded);
                            
                            neuronStats(neuron,orient,1) = fitdist(vData','Normal');
                            neuronStats(neuron,orient,2) = fitdist(avData','Normal');
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
            size(choiceMap{1,c})
            choiceMap{2,c} = cat(3,choiceMap{2,c},squeeze(tempChoiceMap(2,:,:)));
        end
    end
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

     
