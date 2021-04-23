
% clear all;
% close all;
clearvars -except dsiData osiData osiSorted

saveDir = fullfile('D:\Electrophysiology','EP010','MetaData - AVMultiContrastDriftingGratingsWhiteNoise','MLE model','Individual neurons - recat');
if ~exist(saveDir)
    mkdir(saveDir);
end

% dataPaths{1} = fullfile('EP004','AW117','20200221-1');
% dataPaths{1} = fullfile('EP004','AW117','20200221-2');
% dataPaths{1} = fullfile('EP004','AW118','20200221-1');
% dataPaths{4} = fullfile('EP004','AW118','20200221-2');
% dataPaths{1} = fullfile('EP004','AW121','20200226-1');
% dataPaths{6} = fullfile('EP004','AW121','20200226-2');
% dataPaths{1} = fullfile('EP004','AW124','20200303-1');
% dataPaths{1} = fullfile('EP004','AW124','20200303-2');
dataPaths{1} = fullfile('EP010','AW157','20201212-1');
% dataPaths{1} = fullfile('EP010','AW157','20201212-2');
% dataPaths{1} = fullfile('EP010','AW158','20201212-1');
% dataPaths{1} = fullfile('EP010','AW158','20201212-2');

u = 6;

% dataPaths{13} = fullfile('EP010','AW159','20201213-1');
% dataPaths{14} = fullfile('EP010','AW159','20201213-2');
% dataPaths{15} = fullfile('EP010','AW162','20210102-1');
% dataPaths{16} = fullfile('EP010','AW162','20210102-2');
% dataPaths{17} = fullfile('EP010','AW163','20210102-1');
% dataPaths{18} = fullfile('EP010','AW163','20210102-2');
% dataPaths{19} = fullfile('EP010','AW164','20210105');
% dataPaths{20} = fullfile('EP010','AW165','20210106-1');
% dataPaths{21} = fullfile('EP010','AW165','20210106-2');

%which neurons to include
onlySingleUnits = 0;
allNeurons = 0;
lightResponsive = 1;
orientationSelective = 1;
soundResponsive = 1;


dp = 1;

%Load data
dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_final','AVMultiContrastDriftingGratingsWhiteNoiseData_selectivity.mat');
load(dataFile);
stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
load(fullfile(stimPath(end).folder,stimPath(end).name));

orientations = stimInfo.orientations;
contrasts = stimInfo.contrasts;
repeats = stimInfo.repeats;
quantScalar = 1000/(analysisParams.quantWindow(2)-analysisParams.quantWindow(1));

if ~exist('choiceMap','var')
    choiceMap = cell(2,length(contrasts)); %v and av, by contrasts
    accuracy = cell(2,length(contrasts));
    orientationMap = cell(1,length(contrasts));
    directionMap = cell(1,length(contrasts));
end

figure;
for c = 1:5
    vIndBase = (c-1)*length(orientations) + 1;
    avIndBase = (c-1)*length(orientations) + 1 + length(orientations)*length(contrasts);
    
    vResp = unitData(u).meanResponse(vIndBase:vIndBase+length(orientations)-1,4);
    avResp = unitData(u).meanResponse(avIndBase:avIndBase+length(orientations)-1,4);
    %                 if c==1
    %                     prefOrientV = randi(length(orientations));
    %                     prefOrientAV = randi(length(orientations));
    %                 else
    [~, prefOrientV] = max(vResp);
    [~, prefOrientAV] = max(avResp);
    %                 end
    
    orthIndV = [mod(prefOrientV+3-1,length(orientations)) mod(prefOrientV-3-1,length(orientations))]+1;
    orthIndAV = [mod(prefOrientAV+3-1,length(orientations)) mod(prefOrientAV-3-1,length(orientations))]+1;
    
    vvs = unitData(u).trialResponse(vIndBase + prefOrientV - 1,:);
    vvO = unitData(u).trialResponse(vIndBase + [orthIndV] - 1,:);
    avs = unitData(u).trialResponse(avIndBase + prefOrientAV - 1,:);
    avO = unitData(u).trialResponse(avIndBase + [orthIndAV] - 1,:);
    
    subplot(2,5,c);
    hold on;
    histogram(vvs,0:5:120);
    histogram(vvO,0:5:120);
    subplot(2,5,c+5);
    hold on;
    histogram(avs,0:5:120);
    histogram(avO,0:5:120);
end



neuronNumber = unitData(u).neuronNumber;
%         if ismember(neuronNumber, responsiveUnitsGLM.lightSoundInteractUnits)...
%                 && ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits)
if (ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
        (ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) &&...
        ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits))) &&...
        (ismember(neuronNumber,responsiveUnits.directionSelectiveUnits) ||...
        ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits)) &&...
        ismember(neuronNumber,[responsiveUnits.singleUnits responsiveUnits.multiUnits])
    
    
    
    for c = 1:length(contrasts)
        
        %%Overall accuracy
        tempChoiceMap = zeros(2,length(orientations),length(orientations));
        for testDir = 1:length(orientations)
            [dp u c testDir]
            
            mleClassification = zeros(2,length(orientations));
            
            
            for trial = 1:repeats
                if exist('neuronStats','var')
                    clear neuronStats
                end
                
                
                for orient = 1:length(orientations)
                    
                    vInd = (c-1)*length(orientations) + orient;
                    avInd = (c-1)*length(orientations) + orient + length(orientations)*length(contrasts);
                    
                    trialsIncluded = 1:repeats;
                    if orient==testDir
                        trialsIncluded  = setdiff(trialsIncluded,trial);
                        
                        probeTrialV = unitData(u).trialResponse(vInd,trial)/quantScalar;
                        probeTrialAV = unitData(u).trialResponse(avInd,trial)/quantScalar;
                        probeData = [probeTrialV; probeTrialAV];
                    end
                    
                    vData = unitData(u).trialResponse(vInd,trialsIncluded)/quantScalar;
                    avData = unitData(u).trialResponse(avInd,trialsIncluded)/quantScalar;
                    
                    neuronStats(orient,1) = fitdist(vData','Poisson');
                    neuronStats(orient,2) = fitdist(avData','Poisson');
                end
                
                vEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,1),probeData(1));
                mleClassification(1,vEstimate) = mleClassification(1,vEstimate)+1;
                
                avEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,2),probeData(2));
                mleClassification(2,avEstimate) = mleClassification(2,avEstimate) +1;
            end
            
            tempChoiceMap(1,testDir,:) = mleClassification(1,:) / sum(mleClassification(1,:));
            tempChoiceMap(2,testDir,:) = mleClassification(2,:) / sum(mleClassification(2,:));
        end
        choiceMap{1,c} = cat(3,choiceMap{1,c},squeeze(tempChoiceMap(1,:,:)));
        choiceMap{2,c} = cat(3,choiceMap{2,c},squeeze(tempChoiceMap(2,:,:)));
        
        accuracy{1,c} = [accuracy{1,c}; mean(diag(squeeze(tempChoiceMap(1,:,:)))) dp u neuronNumber size(accuracy{1,c},1)+1];
        accuracy{2,c} = [accuracy{2,c}; mean(diag(squeeze(tempChoiceMap(2,:,:)))) dp u neuronNumber size(accuracy{2,c},1)+1];
        
        %%Orientation accuracy
        vIndBase = (c-1)*length(orientations) + 1;
        avIndBase = (c-1)*length(orientations) + 1 + length(orientations)*length(contrasts);
        
        vResp = unitData(u).meanResponse(vIndBase:vIndBase+length(orientations)-1,4);
        avResp = unitData(u).meanResponse(avIndBase:avIndBase+length(orientations)-1,4);
        %                 if c==1
        %                     prefOrientV = randi(length(orientations));
        %                     prefOrientAV = randi(length(orientations));
        %                 else
        [~, prefOrientV] = max(vResp);
        [~, prefOrientAV] = max(avResp);
        %                 end
        
        orthIndV = [mod(prefOrientV+3-1,length(orientations)) mod(prefOrientV-3-1,length(orientations))]+1;
        orthIndAV = [mod(prefOrientAV+3-1,length(orientations)) mod(prefOrientAV-3-1,length(orientations))]+1;
        oppIndV = mod(prefOrientV+6-1,length(orientations))+1;
        oppIndAV = mod(prefOrientAV+6-1,length(orientations))+1;
        
        orientClassification = zeros(2,2);
        for trial = 1:repeats
            if exist('neuronStats','var')
                clear neuronStats
            end
            vIndBase = (c-1)*length(orientations);
            avIndBase = (c-1)*length(orientations) + length(orientations)*length(contrasts);
            
            %             prefSamplesV = unitData(u).trialResponse(vIndBase + [prefOrientV oppIndV],:);
            %             prefSamplesAV = unitData(u).trialResponse(avIndBase + [prefOrientAV oppIndAV],:);
            %             probeTrialV = prefSamplesV(trial)/quantScalar;
            %             probeTrialAV = prefSamplesAV(trial)/quantScalar;
            
            probeTrialV = unitData(u).trialResponse(vIndBase + prefOrientV,trial)/quantScalar;
            probeTrialAV = unitData(u).trialResponse(avIndBase +prefOrientAV,trial)/quantScalar;
            probeData = [probeTrialV; probeTrialAV];
            
            %             prefDataV = reshape(prefSamplesV(setdiff(1:(2*repeats),trial))/quantScalar,[1 2*repeats-1]);
            %             prefDataAV = reshape(prefSamplesAV(setdiff(1:(2*repeats),trial))/quantScalar,[ 1 2*repeats-1]);
            
            prefDataV = unitData(u).trialResponse(vIndBase+prefOrientV,setdiff(1:repeats,trial))/quantScalar;
            prefDataAV = unitData(u).trialResponse(avIndBase+prefOrientAV,setdiff(1:repeats,trial))/quantScalar;
            neuronStats(1,1) = fitdist(prefDataV','Poisson');
            neuronStats(1,2) = fitdist(prefDataAV','Poisson');
            
            orthDataV = reshape(unitData(u).trialResponse(vIndBase+orthIndV,:)/quantScalar,[1 2*repeats]);
            orthDataAV = reshape(unitData(u).trialResponse(avIndBase+orthIndAV,:)/quantScalar,[1 2*repeats]);
            neuronStats(2,1) = fitdist(orthDataV','Poisson');
            neuronStats(2,2) = fitdist(orthDataAV','Poisson');
            
            vEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,1),probeData(1));
            orientClassification(1,vEstimate) = orientClassification(1,vEstimate)+1;
            
            avEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,2),probeData(2));
            orientClassification(2,avEstimate) = orientClassification(2,avEstimate)+1;
        end
        if ismember(neuronNumber, responsiveUnits.orientationSelectiveUnits)
            orientationMap{1,c} = [orientationMap{1,c}; orientClassification(:,1)'./repeats dp u neuronNumber];
        end
        
        %%Direction accuracy
        vIndBase = (c-1)*length(orientations) + 1;
        avIndBase = (c-1)*length(orientations) + 1 + length(orientations)*length(contrasts);
        
        vResp = unitData(u).meanResponse(vIndBase:vIndBase+length(orientations)-1,4);
        avResp = unitData(u).meanResponse(avIndBase:avIndBase+length(orientations)-1,4);
        %                 if c==1
        %                     prefOrientV = randi(length(orientations));
        %                     prefOrientAV = randi(length(orientations));
        %                 else
        [~, prefOrientV] = max(vResp);
        [~, prefOrientAV] = max(avResp);
        %                 end
        
        oppIndV = mod(prefOrientV+6-1,length(orientations))+1;
        oppIndAV = mod(prefOrientAV+6-1,length(orientations))+1;
        
        directClassification = zeros(2,2);
        for trial = 1:repeats
            if exist('neuronStats','var')
                clear neuronStats
            end
            vIndBase = (c-1)*length(orientations);
            avIndBase = (c-1)*length(orientations) + length(orientations)*length(contrasts);
            
            probeTrialV = unitData(u).trialResponse(vIndBase + prefOrientV,trial)/quantScalar;
            probeTrialAV = unitData(u).trialResponse(avIndBase +prefOrientAV,trial)/quantScalar;
            probeData = [probeTrialV; probeTrialAV];
            
            prefDataV = unitData(u).trialResponse(vIndBase+prefOrientV,setdiff(1:repeats,trial))/quantScalar;
            prefDataAV = unitData(u).trialResponse(avIndBase+prefOrientAV,setdiff(1:repeats,trial))/quantScalar;
            neuronStats(1,1) = fitdist(prefDataV','Poisson');
            neuronStats(1,2) = fitdist(prefDataAV','Poisson');
            
            oppDataV = unitData(u).trialResponse(vIndBase+oppIndV,:)/quantScalar;
            oppDataAV = unitData(u).trialResponse(avIndBase+oppIndAV,:)/quantScalar;
            neuronStats(2,1) = fitdist(oppDataV','Poisson');
            neuronStats(2,2) = fitdist(oppDataAV','Poisson');
            
            vEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,1),probeData(1));
            directClassification(1,vEstimate) = directClassification(1,vEstimate)+1;
            
            avEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,2),probeData(2));
            directClassification(2,avEstimate) = directClassification(2,avEstimate)+1;
        end
        if ismember(neuronNumber, responsiveUnits.directionSelectiveUnits)
            directionMap{1,c} = [directionMap{1,c}; directClassification(:,1)'./repeats dp u neuronNumber];
        end
        
    end
    
end

if ismember(neuronNumber, responsiveUnits.orientationSelectiveUnits)
    vvv = zeros(1,5);
    avav = zeros(1,5);
    for c = 1:5
        vvv(c) = orientationMap{1,c}(1,1);
        avav(c) = orientationMap{1,c}(1,2);
    end
    figure; hold on;
    plot(contrasts,vvv,'Color',[0 0 0],'LineWidth',2);
    plot(contrasts,avav,'Color',[0 0 1],'LineWidth',2);
    xticks(contrasts);
    xlabel('Contrast');
    ylabel('Accuracy');
end
